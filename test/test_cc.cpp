// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <cc.hpp>
#include <havoqgt/mpi.hpp>
#include <ingest_edges.hpp>
#include <util.hpp>

using namespace havoqgt::test;

/// \brief Validate BFS results on a chain graph
void validate_cc_result(const graph_type &graph, const uint64_t max_vertex_id,
                        const cc_data_type &cc_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  const auto ans = cc_data[graph.label_to_locator(0)];
  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    const auto vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      if (cc_data[vertex] != ans) {
        std::cerr << "Unexpected CC value" << std::endl;
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

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_10");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(9, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/chain_1024-0");
      f.emplace_back("./datasets/chain_1024-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(1023, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_10");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(9, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/star_1024-0");
      f.emplace_back("./datasets/star_1024-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(1023, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_7");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(6, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    {
      std::vector<std::string> f;
      f.emplace_back("./datasets/binary_tree_1023-0");
      f.emplace_back("./datasets/binary_tree_1023-1");
      ingest_unweighted_edges(f, gen_test_dir_path(k_test_name));
      run_cc(1022, validate_cc_result);
      havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
    }

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded all CC tests" << std::endl;
  }  // End of MPI

  return 0;
}