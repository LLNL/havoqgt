// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <cc.hpp>
#include <havoqgt/mpi.hpp>
#include <ingest_edges.hpp>
#include <util.hpp>

using namespace havoqgt::test;

/// \brief Validates BFS results
void validate_cc_result(const graph_type &graph, const uint32_t max_vertex_id,
                        const cc_data_type &cc_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  // Read the correct CC ID
  std::vector<uint32_t> ans_cc_id(max_vertex_id + 1,
                                  std::numeric_limits<uint32_t>::max());

  {
    std::ifstream ifs("./datasets/cc_rmat_s17");
    if (!ifs) {
      std::cerr << "Cannot open" << std::endl;
      CHK_MPI(MPI_Abort(MPI_COMM_WORLD, -1));
    }

    uint32_t vid;
    uint32_t cc_id;
    while (ifs >> vid >> cc_id) {
      ans_cc_id[vid] = cc_id;
    }
  }

  // Check CC ID
  for (uint32_t v = 0; v <= max_vertex_id; ++v) {
    graph_type::vertex_locator vertex = graph.label_to_locator(v);
    if (vertex.owner() == mpi_rank) {
      const auto cc_id = graph.locator_to_label(cc_data[vertex]);
      // std::cout << cc_id << " == " << ans_cc_id[v] << std::endl;
      if (cc_id != ans_cc_id[v]) {
        std::cerr << "Unexpected CC ID" << std::endl;
        std::cerr << "Vertex " << v << " is in CC " << cc_id << std::endl;
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

    static constexpr uint32_t max_vertex_id = 131070;
    run_cc(max_vertex_id, validate_cc_result);

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded the CC with delegate test" << std::endl;

    havoqgt::distributed_db::remove(gen_test_dir_path(k_test_name));
  }  // End of MPI

  return 0;
}