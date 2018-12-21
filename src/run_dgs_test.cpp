/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setfill, std::setw
#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include <functional>
#include <fcntl.h>
#include <tuple>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/mpi.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/distributed_db.hpp>

#include <dynamic_graph_store/graphstore_utilities.hpp>
#include <dynamic_graph_store/baseline/baseline.hpp>
#include <dynamic_graph_store/baseline/baseline_map.hpp>
#include <dynamic_graph_store/degawarerhh/degawarerhh.hpp>
#include <dynamic_graph_store/dist_dynamic_graphstore.hpp>
#include <dynamic_graph_store/edge_list_buffer.hpp>

/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = unsigned char;
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              havoqgt::distributed_db::segment_manager_type>;

using baselinemap_type      = graphstore::graphstore_baseline_map<vertex_id_type,
                                                                  vertex_property_type,
                                                                  edge_property_type,
                                                                  havoqgt::distributed_db::segment_manager_type>;

enum : size_t { middle_high_degree_threshold = 8 }; // must be more than 1
using degawarerhh_type  = graphstore::degawarerhh<vertex_id_type,
                                                  vertex_property_type,
                                                  edge_property_type,
                                                  havoqgt::distributed_db::segment_manager_type,
                                                  middle_high_degree_threshold>;


/// --- Choose graph store type --- ///
/// baseline_type, baselinemap_type or degawarerhh_type
using gstore_type = degawarerhh_type;

using dist_gstore_type = graphstore::dist_dynamic_graphstore<gstore_type>;
using visitor_type = dg_visitor<dist_gstore_type>;
//using dg_visitor_queue_type = havoqgt::visitor_queue<visitor_type,
//                                                     havoqgt::detail::visitor_priority_queue,
//                                                     dist_gstore_type, char>;

using edge_list_buffer_type = graphstore::edge_list_buffer<uint64_t, bool>;

void usage() {
  if (havoqgt::mpi_comm_rank() == 0) {
    std::cerr << "Usage: \n"
              << " -o <string>   - base filename to create segmentfiles\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_options(int argc, char **argv, std::string &segment_file_prefix) {
  if (havoqgt::mpi_comm_rank() == 0) {
    std::cout << "MPI initialized with " << havoqgt::mpi_comm_size() << " ranks." << std::endl;
  }

  char c;
  while ((c = getopt(argc, argv, "o:h")) != -1) {
    switch (c) {
      case 'o': {
        segment_file_prefix = optarg;
        break;
      }
      case 'h':usage();
        break;
    }
  }
}

void generate_edgelist(edge_list_buffer_type &edgelist) {
  uint64_t num_edges_per_rank = (1ULL << 20) * 16ULL / havoqgt::mpi_comm_size();
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(havoqgt::mpi_comm_rank()) * 3ULL,
                                    20, num_edges_per_rank,
                                    0.57, 0.19, 0.19, 0.05, true, true);
  edgelist.clear();
  for (const auto &edge : rmat) {
    edgelist.emplace_back(edge.first, edge.second, false);
  }
}

using eddge_holder_type = typename edge_list_buffer_type::value_type;
using mpi_buf_vec_type = std::vector<std::vector<eddge_holder_type>>;
void global_sort(edge_list_buffer_type &edgelist) {
  int mpi_rank = havoqgt::mpi_comm_rank();
  int mpi_size = havoqgt::mpi_comm_size();

  mpi_buf_vec_type send_buf_vec(mpi_size);
  mpi_buf_vec_type recev_buf_vec(mpi_size);

  for (const auto &edge : edgelist) {
    const int target_rank = std::get<0>(edge) % mpi_size;
    send_buf_vec[target_rank].emplace_back(edge);
  }
  edgelist.clear();
  MPI_Barrier(MPI_COMM_WORLD);

  havoqgt::mpi_all_to_all(send_buf_vec, recev_buf_vec, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  for (auto &vec_per_rank : recev_buf_vec) {
    for (auto &edge : vec_per_rank) {
      edgelist.emplace_back(edge);
    }
  }

  graphstore::sort_edge_list_buffer(edgelist);
}

void remove_duplicated_edges(edge_list_buffer_type &edgelist) {
  auto pre_edge(edgelist[0]);
  size_t tail = 1;
  for (size_t i = 1; i < edgelist.size(); ++i) {
    if (pre_edge == edgelist[i]) continue;
    edgelist[tail++] = edgelist[i];
    pre_edge = edgelist[i];
  }
  std::cout << "#duplcates " << edgelist.size() - tail << std::endl;
  edgelist.resize(tail);
}

void dump_all_edges(dist_gstore_type &dist_graph, edge_list_buffer_type &buffer) {
  buffer.clear();
  for (auto vrtx = dist_graph.vertices_begin(), vrtx_end = dist_graph.vertices_end();
       vrtx != vrtx_end;
       ++vrtx) {
    for (auto edge = dist_graph.adjacent_edge_begin(vrtx.source_vertex()),
             edge_end = dist_graph.adjacent_edge_end(vrtx.source_vertex());
         edge != edge_end;
         ++edge) {
      buffer.emplace_back(vrtx.source_vertex().id(), edge.target_vertex().id(), false);
    }
  }
}

int main(int argc, char **argv) {

  havoqgt::init(&argc, &argv);
  {
    int mpi_rank = havoqgt::mpi_comm_rank();

    MPI_Barrier(MPI_COMM_WORLD);

    /// --- init segment file --- ///
    std::string segment_file_prefix = "/dev/shm/gstore";
    parse_options(argc, argv, segment_file_prefix);
    size_t graph_capacity_gb_per_rank = 2;
    havoqgt::distributed_db ddb(havoqgt::db_create(), segment_file_prefix.c_str(), graph_capacity_gb_per_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    gstore_type gstore(ddb.get_segment_manager());
    dist_gstore_type dist_graph(&gstore);
    visitor_type::set_graph_ref(&dist_graph);
    char dummy;
    auto dg_visitor_queue = havoqgt::create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(&dist_graph, dummy);
    MPI_Barrier(MPI_COMM_WORLD);

    edge_list_buffer_type input_edgelist;
    generate_edgelist(input_edgelist);
    MPI_Barrier(MPI_COMM_WORLD);

    /// --- dynamic graph construction --- ///
    if (mpi_rank == 0) std::cout << "--- Start dynamic graph construction ---" << std::endl;
    dg_visitor_queue.dynamic_graph_construction(input_edgelist.begin(), input_edgelist.end());
    MPI_Barrier(MPI_COMM_WORLD);

    /// --- Check result --- ///
    if (mpi_rank == 0) std::cout << "--- start dumping edges ---" << std::endl;
    edge_list_buffer_type dumped_edgelist;
    dump_all_edges(dist_graph, dumped_edgelist);
    graphstore::sort_edge_list_buffer(dumped_edgelist);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) std::cout << "--- start global sort ---" << std::endl;
    global_sort(input_edgelist);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) std::cout << "--- start removing duplicated edges ---" << std::endl;
    remove_duplicated_edges(input_edgelist);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) std::cout << "--- start checking result ---" << std::endl;
    if (dumped_edgelist.size() != input_edgelist.size()) {
      std::cout << dumped_edgelist.size() << " != " << input_edgelist.size() << std::endl;
      exit(1);
    }

    for (size_t i = 0; i < dumped_edgelist.size(); ++i) {
      if (dumped_edgelist[i] != input_edgelist[i]) {
        std::cout << "different edge" << std::endl;
        std::cout << std::get<0>(dumped_edgelist[i]) << " " << std::get<1>(dumped_edgelist[i]) << "\n"
                  << std::get<0>(input_edgelist[i]) << " " << std::get<1>(input_edgelist[i]) << std::endl;
        exit(1);
      }
    }
  }

  std::cout << "Finished the test successfully" << std::endl;

  return 0;
}
