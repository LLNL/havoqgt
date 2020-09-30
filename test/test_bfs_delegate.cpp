// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <bfs.hpp>
#include <havoqgt/mpi.hpp>
#include <ingest_edges.hpp>
#include <util.hpp>

using namespace havoqgt::test;

/// \brief Validates BFS results
void validate_bfs_result(const graph_type &          graph,
                         const std::size_t           max_vertex_id,
                         const bfs_level_data_type & level_data,
                         const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  // Load the correct BFS levels
  std::unordered_map<uint64_t, uint64_t> level_table;
  {
    std::ifstream ifs("./datasets/bfs_level_rmat_s17");
    if (!ifs) {
      std::cerr << "Cannot open" << std::endl;
      CHK_MPI(MPI_Abort(MPI_COMM_WORLD, -1));
    }

    uint64_t vid;
    uint64_t level;
    while (ifs >> vid >> level) {
      level_table[vid] = level;
    }
  }

  // Validate BFS level
  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    graph_type::vertex_locator vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      size_t level = level_data[vertex];
      if (level_table[i] != level) {
        std::cerr << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
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

    std::vector<std::string> f;
    f.emplace_back("./datasets/edge_list_rmat_s17-0");
    f.emplace_back("./datasets/edge_list_rmat_s17-1");
    f.emplace_back("./datasets/edge_list_rmat_s17-2");
    f.emplace_back("./datasets/edge_list_rmat_s17-3");
    f.emplace_back("./datasets/edge_list_rmat_s17-4");
    f.emplace_back("./datasets/edge_list_rmat_s17-5");
    f.emplace_back("./datasets/edge_list_rmat_s17-6");
    f.emplace_back("./datasets/edge_list_rmat_s17-7");
    static constexpr uint64_t delegate_degree_threshold = 1024 * 10;
    ingest_unweighted_edges(f, gen_test_dir_path(k_test_name),
                            delegate_degree_threshold);

    static constexpr uint64_t max_vertex_id = 131070;
    run_bfs(max_vertex_id, validate_bfs_result);

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded the BFS with delegate test"
                          << std::endl;

    havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
  }  // End of MPI

  return 0;
}