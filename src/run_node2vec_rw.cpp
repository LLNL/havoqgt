/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <mpi.h>
#include <node2vec_rw/node2vec_rw.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

struct node2vec_rw_option {
  std::size_t walk_length{80}; // Length of walk per source
  std::size_t num_walkers_per_vertex{10}; // Number of walks per source
  double p{1.0}; // Return hyperparameter
  double q{1.0}; // Inout hyperparameter
  bool small_edge_weight_variance = false; // This option is not used for now (will be implemented soon)
  unsigned int rnd_seed{113}; // This option is not used for now (will be implemented soon)
  bool verbose{true}; // This option is not used for now (will be implemented soon)
};

// Graph data class
// The actual implemented is include/havoqgt/delegate_partitioned_graph.hpp
using graph_type = node2vec_rw::node2vec_rw<>::graph_type;

// Vertex data is a 64 bits of class which has some informaion for distributed memory
// One can get a vertex ID with uint64_t format using locator_to_label() function in graph_type
// The actual implemented is include/havoqgt/impl/vertex_locator.hpp
using vertex_type = node2vec_rw::node2vec_rw<>::vertex_type;

// In HavoqGT, edge data are stored separately from the graph data
// The actual implementaion is include/havoqgt/impl/edge_data.hpp
using edge_data_type = node2vec_rw::node2vec_rw<>::edge_weight_data_type;

void parse_cmd_line(int argc, char **argv,
                    std::string *graph_filename,
                    std::string *backup_filename,
                    node2vec_rw_option *option,
                    std::string *walk_history_out_file_name) {
  int c;
  while ((c = getopt(argc, argv, "g:b:l:w:p:q:s:r:vo:")) != -1) {
    switch (c) {
      case 'g':*graph_filename = optarg;
        break;

      case 'b':*backup_filename = optarg;
        break;

      case 'l':option->walk_length = std::stoull(optarg);
        break;

      case 'w': option->num_walkers_per_vertex = std::stoll(optarg);
        break;

      case 'p':option->p = std::stod(optarg);
        break;

      case 'q':option->q = std::stod(optarg);
        break;

      case 's':option->small_edge_weight_variance = std::stoul(optarg);
        break;

      case 'r':option->rnd_seed = std::stoul(optarg);
        break;

      case 'v': option->verbose = true;
        break;

      case 'o':*walk_history_out_file_name = optarg;
        break;

    }
  }

  if (graph_filename->empty() && backup_filename->empty()) {
    std::cerr << "Graph file name and backup file name are empty" << std::endl;
    std::abort();
  }
}

void dump_walk_history(const graph_type &graph,
                       const std::vector<std::vector<vertex_type>> &walk_history,
                       const std::string &walk_history_out_file_name,
                       const bool truncate_file) {
  int mpi_rank = -1;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  std::string fout_name(walk_history_out_file_name + "-" + std::to_string(mpi_rank));
  std::ofstream ofs(fout_name,
                    (truncate_file) ? std::ofstream::trunc : std::ofstream::app);
  for (auto &per_wk_history : walk_history) {
    for (auto &v : per_wk_history) {
      ofs << graph.locator_to_label(v) << " ";
    }
    ofs << "\n";
  }
  ofs.close();
}

int main(int argc, char **argv) {
  havoqgt::init(&argc, &argv);
  {
    std::string graph_filename;
    std::string backup_filename;
    node2vec_rw_option option;
    std::string walk_history_out_file_name;
    parse_cmd_line(argc, argv, &graph_filename, &backup_filename, &option, &walk_history_out_file_name);

    // Copy files from backup_filename to graph_filename
    if (!backup_filename.empty()) {
      havoqgt::distributed_db::transfer(backup_filename.c_str(), graph_filename.c_str());
    }

    // Load graph data
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_filename.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;

    // Load edge data
    // If it does not exist, set 1.0 to all edges
    edge_data_type *edge_data = ddb.get_segment_manager()->find<edge_data_type>("graph_edge_data_obj").first;
    if (!edge_data) {
      edge_data = new edge_data_type(*graph);
      edge_data->reset(1.0);
    }

    node2vec_rw::node2vec_rw<> n2v_rw(*graph,
                                      *edge_data,
                                      option.small_edge_weight_variance,
                                      option.walk_length,
                                      option.p,
                                      option.q,
                                      MPI_COMM_WORLD,
                                      option.verbose);

    // Start mini-batched random walk
    int mini_batch_size = 1024; // Run 1024 walkers from each process
    bool first_walk_history_dump = true;
    for (int n = 0; n < option.num_walkers_per_vertex; ++n) {
      auto itr = graph->vertices_begin();
      while (true) {
        std::vector<vertex_type> start_vertices;
        for (int b = 0; b < mini_batch_size && itr != graph->vertices_end(); ++b, ++itr) {
          start_vertices.push_back(*itr);
        }

        std::vector<std::vector<vertex_type>> walk_history = n2v_rw.run_walker(start_vertices);
        if (!walk_history_out_file_name.empty()) {
          dump_walk_history(*graph, walk_history, walk_history_out_file_name, first_walk_history_dump);
          first_walk_history_dump = false;
        }

        // Global termination check
        const char local_finished = (itr == graph->vertices_end());
        if (havoqgt::mpi_all_reduce(local_finished, std::logical_and<char>{}, MPI_COMM_WORLD)) {
          break;
        }
      }
    }
  }

  return 0;
}