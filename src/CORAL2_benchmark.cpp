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
 * License. http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
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
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by the United States Government or Lawrence Livermore National Security, LLC.
 * The views and opinions of authors expressed herein do not necessarily state
 * or reflect those of the United States Government or Lawrence Livermore
 * National Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

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

  uint64_t chunk_size = 64 * 1024;
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
