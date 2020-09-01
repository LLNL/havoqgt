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

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

typedef havoqgt::delegate_partitioned_graph<havoqgt::distributed_db::allocator<>> graph_type;
typedef uint8_t edge_data_type;
typedef havoqgt::distributed_db::allocator<edge_data_type> edge_data_allocator_type;
typedef graph_type::vertex_data<uint16_t, std::allocator<uint16_t>> bfs_level_data_type;
typedef graph_type::vertex_data<graph_type::vertex_locator, std::allocator<graph_type::vertex_locator>>
    bfs_parent_data_type;

static constexpr const char *k_data_store_path = "/tmp/havoqgt_test_bfs";
static constexpr const char *k_graph_key_name = "graph_obj";

/// \brief Ingest edge lists
void ingest_edges(const std::vector<std::string> &edge_list_file_names) {
  havoqgt::parallel_edge_list_reader<edge_data_type> pelr(edge_list_file_names, true);
  bool has_edge_data = pelr.has_edge_data();
  if (has_edge_data) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  havoqgt::distributed_db ddb(havoqgt::db_create(), k_data_store_path);
  graph_type::edge_data<edge_data_type, edge_data_allocator_type> dummy_edge_data(ddb.get_allocator());
  graph_type *graph = ddb.get_manager()->construct<graph_type>(k_graph_key_name)
      (ddb.get_allocator(), MPI_COMM_WORLD, pelr, pelr.max_vertex_id(), 1024, 1, 8192, dummy_edge_data);
}

/// \brief Run BFS from vertex 0
void run_bfs(const std::size_t max_vertex_id,
             const std::function<void(const graph_type &,
                                      const std::size_t,
                                      const bfs_level_data_type &,
                                      const bfs_parent_data_type &)> &validate_bfs_result) {
  havoqgt::distributed_db ddb(havoqgt::db_open_read_only(), k_data_store_path);
  auto graph = ddb.get_manager()->find<graph_type>(k_graph_key_name).first;
  bfs_level_data_type bfs_level_data(*graph);
  bfs_level_data.reset(std::numeric_limits<uint16_t>::max());
  bfs_parent_data_type bfs_parent_data(*graph);
  havoqgt::comm_world().barrier();

  auto source = graph->label_to_locator(0);
  havoqgt::breadth_first_search(graph, bfs_level_data, bfs_parent_data, source);

  validate_bfs_result(*graph, max_vertex_id, bfs_level_data, bfs_parent_data);
}

/// \brief Validate BFS results on a chain graph
void validate_chain_bfs_result(const graph_type &graph,
                               const std::size_t max_vertex_id,
                               const bfs_level_data_type &level_data,
                               const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (size_t i = 0; i <= max_vertex_id; ++i) {
    graph_type::vertex_locator vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      size_t level = level_data[vertex];
      if (i != level) {
        std::cout << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        graph_type::vertex_locator parent = graph.label_to_locator(i - 1);
        if (parent_data[vertex] != parent) {
          std::cout << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
    havoqgt::comm_world().barrier();
  }
}

/// \brief Validate BFS results on a star graph
void validate_star_bfs_result(const graph_type &graph,
                              const std::size_t max_vertex_id,
                              const bfs_level_data_type &level_data,
                              const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (size_t i = 0; i <= max_vertex_id; ++i) {
    graph_type::vertex_locator vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      size_t level = level_data[vertex];
      if ((i == 0 && level != 0) || (i > 0 && level != 1)) {
        std::cout << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        graph_type::vertex_locator parent = graph.label_to_locator(0);
        if (parent_data[vertex] != parent) {
          std::cout << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
    havoqgt::comm_world().barrier();
  }
}

/// \brief Validate BFS results on a binary tree graph
void validate_binary_tree_bfs_result(const graph_type &graph,
                                     const std::size_t max_vertex_id,
                                     const bfs_level_data_type &level_data,
                                     const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (size_t i = 0; i <= max_vertex_id; ++i) {
    graph_type::vertex_locator vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      size_t level = level_data[vertex];
      if (level != (size_t)std::log2(i + 1)) {
        std::cout << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        graph_type::vertex_locator parent = graph.label_to_locator((i - 1) / 2);
        if (parent_data[vertex] != parent) {
          std::cout << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
    havoqgt::comm_world().barrier();
  }
}

int main(int argc, char **argv) {

  havoqgt::init(&argc, &argv);
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    havoqgt::cout_rank0() << "MPI initialized with " << mpi_size << " ranks." << std::endl;

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_10");
      ingest_edges(f);
      run_bfs(9, validate_chain_bfs_result);
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_1024");
      ingest_edges(f);
      run_bfs(1023, validate_chain_bfs_result);
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_10");
      ingest_edges(f);
      run_bfs(9, validate_star_bfs_result);
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_1024");
      ingest_edges(f);
      run_bfs(1023, validate_star_bfs_result);
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_7");
      ingest_edges(f);
      run_bfs(6, validate_binary_tree_bfs_result);
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_1023");
      ingest_edges(f);
      run_bfs(1022, validate_binary_tree_bfs_result);
    }

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded all BFS tests" << std::endl;
  } // End of MPI

  return 0;
}