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
                    std::string& output_filename) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_output_filename = false;
  scale                      = 17;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:o:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's':
        scale = atoll(optarg);
        break;
      case 'o':
        found_output_filename = true;
        output_filename       = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_output_filename) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  init(&argc, &argv);
  {
    std::string output_filename;
    int         mpi_rank = comm_world().rank();
    int         mpi_size = comm_world().size();

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks."
                << std::endl;
    }
    comm_world().barrier();

    uint64_t num_vertices = 1;
    uint64_t vert_scale;

    parse_cmd_line(argc, argv, vert_scale, output_filename);

    num_vertices <<= vert_scale;
    if (mpi_rank == 0) {
      std::cout << "Building Graph500" << std::endl
                << "Building graph Scale: " << vert_scale << std::endl;
    }

    // Generate RMAT graph
    uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
    havoqgt::rmat_edge_generator rmat(
        uint64_t(5489) + uint64_t(mpi_rank) * 3ULL, vert_scale,
        num_edges_per_rank, 0.57, 0.19, 0.19, 0.05, true, false);

    //
    //
    uint64_t global_desired_num_edges = num_vertices * 16;
    using edge_type                   = std::pair<uint64_t, uint64_t>;
    std::deque<edge_type> local_edges;
    uint64_t              local_edges_generated(0);
    uint64_t              local_edges_recv(0);
    auto recvr = [&local_edges, &local_edges_recv](auto* mail, bool bcast,
                                                   const edge_type& edge) {
      local_edges.push_back(edge);
      ++local_edges_recv;
    };

    auto     rmat_itr = rmat.begin();
    uint64_t global_num_edges(0);
    do {
      {
        ygm::mailbox_p2p_nrroute<edge_type, decltype(recvr)> mailbox(
            recvr, 1024 * 1024);
        uint64_t local_gen_count =
            (global_desired_num_edges - global_num_edges) / comm_world().size();
        if ((global_desired_num_edges - global_num_edges) % comm_world().size())
          ++local_gen_count;
        for (uint64_t i = 0; i < local_gen_count; ++i) {
          ++local_edges_generated;
          edge_type edge = *rmat_itr;
          ++rmat_itr;
          if (edge.first != edge.second) {
            if (edge.first > edge.second) std::swap(edge.first, edge.second);
            size_t owner = (hash32(edge.first) ^ hash32(edge.second)) %
                           comm_world().size();
            mailbox.send(owner, edge);
          }
        }
      }
      std::sort(local_edges.begin(), local_edges.end());
      local_edges.erase(std::unique(local_edges.begin(), local_edges.end()),
                        local_edges.end());
      global_num_edges =
          comm_world().all_reduce((uint64_t)local_edges.size(), MPI_SUM);

    } while (global_num_edges < global_desired_num_edges);

    uint64_t global_edges_generated =
        comm_world().all_reduce(local_edges_generated, MPI_SUM);
    uint64_t global_edges_recv =
        comm_world().all_reduce(local_edges_recv, MPI_SUM);

    if (mpi_rank == 0) {
      std::cout << "Total Edges gen = " << global_num_edges << std::endl;
      std::cout << "Edges should be = " << num_vertices * 16 << std::endl;
      std::cout << "global_edges_generated = " << global_edges_generated
                << std::endl;
      std::cout << "% duplicate = "
                << double(global_edges_generated) / double(global_num_edges)
                << std::endl;
      // std::cout << "comm_nl().size() = " << comm_nl().size() << std::endl;
    }

    //
    // Write out edges
    {
      std::vector<char> vec(16 * 1024 * 1024);
      std::stringstream fname;
      fname << output_filename << "_" << mpi_rank << "_of_" << mpi_size;
      std::ofstream ofs(fname.str().c_str());
      ofs.rdbuf()->pubsetbuf(&vec.front(), vec.size());
      // for(size_t i=0; i< comm_nl().size(); ++i) {
      //    if(i == comm_nl().rank()) {
      for (auto& edge : local_edges) {
        ofs << edge.first << " " << edge.second << "\n"
            << edge.second << " " << edge.first << "\n";
      }
      //    }

      // comm_world().barrier();
      //  }
    }

    comm_world().barrier();
  }  // END Main MPI
  ;
  return 0;
}
