// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <bfs.hpp>
#include <havoqgt/mpi.hpp>
#include <util.hpp>
#include <ingest_edges.hpp>

using namespace havoqgt::test;

/// \brief Validate BFS results on a chain graph
void validate_chain_bfs_result(const graph_type &          graph,
                               const uint64_t              max_vertex_id,
                               const bfs_level_data_type & level_data,
                               const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    const auto vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      const auto level = level_data[vertex];
      if (i != level) {
        std::cerr << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        const auto parent = graph.label_to_locator(i - 1);
        if (parent_data[vertex] != parent) {
          std::cerr << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
  }
}

/// \brief Validate BFS results on a star graph
void validate_star_bfs_result(const graph_type &          graph,
                              const uint64_t              max_vertex_id,
                              const bfs_level_data_type & level_data,
                              const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    const auto vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      const auto level = level_data[vertex];
      if ((i == 0 && level != 0) || (i > 0 && level != 1)) {
        std::cerr << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        const auto parent = graph.label_to_locator(0);
        if (parent_data[vertex] != parent) {
          std::cerr << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
  }
}

/// \brief Validate BFS results on a binary tree graph
void validate_binary_tree_bfs_result(const graph_type &          graph,
                                     const uint64_t              max_vertex_id,
                                     const bfs_level_data_type & level_data,
                                     const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    const auto vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      const auto level = level_data[vertex];
      if (level != (uint64_t)std::log2(i + 1)) {
        std::cerr << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }

      if (i > 0) {
        const auto parent = graph.label_to_locator((i - 1) / 2);
        if (parent_data[vertex] != parent) {
          std::cerr << "Unexpected BFS parent" << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  havoqgt::init(&argc, &argv);
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    havoqgt::cout_rank0() << "MPI initialized with " << mpi_size << " ranks."
                          << std::endl;

    havoqgt::cout_rank0() << "Data store path:\t"
                          << gen_test_dir_path(k_test_name) << std::endl;

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_10");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(9, validate_chain_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_1024-0");
      f.emplace_back("./datasets/chain_1024-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(1023, validate_chain_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_10");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(9, validate_star_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_1024-0");
      f.emplace_back("./datasets/star_1024-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(1023, validate_star_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_7");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(6, validate_binary_tree_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_1023-0");
      f.emplace_back("./datasets/binary_tree_1023-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_bfs(1022, validate_binary_tree_bfs_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded all BFS tests" << std::endl;
  }  // End of MPI

  return 0;
}