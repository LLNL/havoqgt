// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/rmat_edge_generator.hpp>

#include <assert.h>
#include <algorithm>
#include <deque>
#include <fstream>  // std::ifstream
#include <functional>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <iostream>
#include <utility>
#include <ygm/mailbox_p2p_nrroute.hpp>
#include <ygm/mpi.hpp>

// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;
using havoqgt::detail::hash32;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr << "Usage: -s <int> -d <int> -o <string>\n"
              << " -s <int>    - RMAT graph Scale (default 17)\n"
              << " -o <string> - output graph base filename\n"
              << " -h          - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, uint64_t& scale,
                    uint64_t& delegate_threshold) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  scale              = 17;
  delegate_threshold = 1024;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:d:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's':
        scale = atoll(optarg);
        break;
      case 'd':
        delegate_threshold = atoll(optarg);
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {
  init(&argc, &argv);
  int mpi_rank = comm_world().rank();
  int mpi_size = comm_world().size();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
  }
  comm_world().barrier();

  uint64_t num_vertices = 1;
  uint64_t vert_scale;
  uint64_t delegate_threshold;

  parse_cmd_line(argc, argv, vert_scale, delegate_threshold);

  num_vertices <<= vert_scale;

  if (mpi_rank == 0) {
    std::cout << "Building Graph500" << std::endl
              << "Building graph Scale: " << vert_scale << std::endl;
  }

  // Generate RMAT graph
  uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                    vert_scale, num_edges_per_rank, 0.57, 0.19,
                                    0.19, 0.05, true, false);

  //
  //
  uint64_t global_desired_num_edges = num_vertices * 16;
  std::map<uint64_t, uint64_t> vertex_degree;
  std::map<uint64_t, uint64_t> vertex_high_degree;
  uint64_t local_recv(0);
  uint64_t local_gen(0);
  auto     recvr = [&vertex_degree, &vertex_high_degree, &local_recv,
                &delegate_threshold](auto* mail, bool bcast, uint64_t v) {
    if (bcast) {
      std::cout << whoami() << " - recv_bcast" << std::endl;
      vertex_high_degree[v] = 0;
    } else {
      ++local_recv;
      if (vertex_degree[v]++ > delegate_threshold) {
        if (vertex_high_degree.count(v) == 0) {
          vertex_high_degree[v] = 0;
          std::cout << whoami() << " - send_bcast" << std::endl;
          mail->send_bcast(v);
        }
      }
    }

  };

  uint64_t local_gen_count = global_desired_num_edges / comm_world().size();
  if (global_desired_num_edges % comm_world().size()) ++local_gen_count;

  {
    ygm::mailbox_p2p_nrroute<uint64_t, decltype(recvr)> mailbox(recvr,
                                                                1024 * 1024);

    auto rmat_itr = rmat.begin();
    for (uint64_t i = 0; i < local_gen_count; ++i) {
      local_gen += 2;
      uint64_t ei = rmat_itr->first;
      if (vertex_high_degree.count(ei) > 0) {
        vertex_high_degree[ei]++;
      } else {
        mailbox.send(ei % comm_world().size(), ei);
      }
      uint64_t ej = rmat_itr->second;
      if (vertex_high_degree.count(ej) > 0) {
        vertex_high_degree[ej]++;
      } else {
        mailbox.send(ej % comm_world().size(), ej);
      }
      ++rmat_itr;
    }
  }

  uint64_t global_generated = comm_world().all_reduce(local_gen, MPI_SUM);
  uint64_t global_recv      = comm_world().all_reduce(local_recv, MPI_SUM);
  size_t   global_min_dsize =
      comm_world().all_reduce(vertex_high_degree.size(), MPI_MIN);
  if (global_min_dsize != vertex_high_degree.size()) {
    std::cout << whoami()
              << " - delegate counts don't match:   global_min_dsize = "
              << global_min_dsize
              << "vertex_high_degree.size() = " << vertex_high_degree.size()
              << std::endl;
    abort();
  } else {
    if (mpi_rank == 0) {
      std::cout << "Delegate Counts match!" << std::endl;
    }
  }

  for (const auto& kv : vertex_high_degree) {
    auto gkv = comm_world().all_reduce(kv.second, MPI_SUM);
    if (kv.first % comm_world().size() == comm_world().rank()) {
      vertex_degree[kv.second] += gkv;
    }
  }

  uint64_t local_sum(0), local_xor(0);
  for (const auto& kv : vertex_degree) {
    local_sum += kv.second;
    local_xor ^= kv.first ^ kv.second;
  }

  uint64_t global_sum = comm_world().all_reduce(local_sum, MPI_SUM);
  uint64_t global_xor = comm_world().all_reduce(local_xor, MPI_LXOR);

  if (mpi_rank == 0) {
    std::cout << "global_generated = " << global_generated << std::endl;
    std::cout << "global_recv = " << global_recv << std::endl;
    std::cout << "global_sum = " << global_sum << std::endl;
    std::cout << "global_xor = " << global_xor << std::endl;
  }

  comm_world().barrier();
  return 0;
}
