// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/kronecker_edge_generator.hpp>
#include <havoqgt/triangle_count_global.hpp>

#include <assert.h>
#include <unistd.h>
#include <algorithm>
#include <deque>
#include <fstream>  // std::ifstream
#include <functional>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <iostream>
#include <utility>

// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;

typedef uint64_t gt_tc_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr << "Usage: file1 file2 \n"
              << "file1          - Edge list file for first graph\n"
              << "file2          - Edge list file for second graph\n\n";
  }
}

int main(int argc, char** argv) {
  using allocator_type = std::allocator<void>;
  using graph_type     = havoqgt::delegate_partitioned_graph<allocator_type>;

  init(&argc, &argv);

  int mpi_rank = comm_world().rank();
  int mpi_size = comm_world().size();

  if (argc != 3) {
    usage();
    exit(-1);
  }

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
  }
  comm_world().barrier();

  uint64_t    delegate_threshold = 1024 * 1024 * 1024;
  std::string input_filename1    = argv[1];
  std::string input_filename2    = argv[2];
  uint64_t    partition_passes   = 1;

  uint64_t chunk_size = 8 * 1024;
  bool     scramble   = true;

  // Setup Kronecker generator
  kronecker_edge_generator<gt_tc_type> kron(input_filename1, input_filename2,
                                            scramble, false);

  if (mpi_rank == 0) {
    std::cout << "Generating  graph." << std::endl;
  }
  graph_type graph(allocator_type(), MPI_COMM_WORLD, kron, kron.max_vertex_id(),
                   delegate_threshold, partition_passes, chunk_size);

  graph.print_graph_statistics();

  if (comm_world().rank() == 0) {
    std::cout << "\n\n\n";
    std::cout << "Benchmark results for:" << std::endl
              << "A = " << input_filename1 << std::endl
              << "B = " << input_filename2 << std::endl;
  }
  comm_world().barrier();

  uint64_t global_tc = triangle_count_global(graph);
  uint64_t global_gt_tc(0);
  for (auto edgetuple : kron) {
    global_gt_tc += std::get<2>(edgetuple);
  }
  global_gt_tc /= 2;  // Each edge is counted twice from kron stream
  global_gt_tc /= 3;  // Each triangle is counted 3 times, once per edge

  global_gt_tc = comm_world().all_reduce(global_gt_tc, MPI_SUM);
  if (comm_world().rank() == 0) {
    if (global_tc == global_gt_tc) {
      std::cout << "!!PASSED!!" << std::endl;
      // std::cout << "Global Triangle Count = " << global_gt_tc / 3 <<
      // std::endl;
    } else {
      std::cout << "?? FAILED ??" << std::endl;
      std::cout << global_gt_tc << " != " << global_tc << std::endl;
    }
  }

  comm_world().barrier();
  return 0;
}
