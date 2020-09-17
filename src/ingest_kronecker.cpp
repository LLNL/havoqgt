// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/kronecker_edge_generator.hpp>

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

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

typedef uint64_t edge_data_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr
        << "Usage: -o <string> -d <int> file1 file2 \n"
        << " -o <string>   - output graph base filename (required)\n"
        << " -b <string>   - backup graph base filename \n"
        << " -d <int>      - delegate threshold (Default is 1048576)\n"
        << " -h            - print help and exit\n"
        << " -p <int>      - number of Low & High partition passes (Default is "
           "1)\n"
        << " -f <float>    - Gigabytes reserved per rank (Default is 0.25)\n"
        << " -c <int>      - Edge partitioning chunk size (Defulat is 8192)\n"
        << " -u <bool>     - Treat edgelist as undirected (Default is 0)\n"
        << "file1          - Edge list file for first graph\n"
        << "file2          - Edge list file for second graph\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& output_filename,
                    std::string& backup_filename, uint64_t& delegate_threshold,
                    std::string& input_filename1, std::string& input_filename2,
                    double& gbyte_per_rank, uint64_t& partition_passes,
                    uint64_t& chunk_size, bool& undirected) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_output_filename = false;
  delegate_threshold         = 1048576;
  gbyte_per_rank             = 0.25;
  partition_passes           = 1;
  chunk_size                 = 8 * 1024;
  undirected                 = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "o:d:p:f:c:b:u:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 'd':
        delegate_threshold = atoll(optarg);
        break;
      case 'o':
        found_output_filename = true;
        output_filename       = optarg;
        break;
      case 'b':
        backup_filename = optarg;
        break;
      case 'p':
        partition_passes = atoll(optarg);
        break;
      case 'f':
        gbyte_per_rank = atof(optarg);
        break;
      case 'c':
        chunk_size = atoll(optarg);
        break;
      case 'u':
        undirected = atoi(optarg);
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

  if (argc - optind != 2) {
    usage();
    exit(-1);
  }

  input_filename1 = argv[optind++];
  input_filename2 = argv[optind];
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  init(&argc, &argv);
  {
    std::string output_filename;
    std::string backup_filename;

    {  // Build Distributed_DB
      int mpi_rank = comm_world().rank();
      int mpi_size = comm_world().size();

      if (mpi_rank == 0) {
        std::cout << "MPI initialized with " << mpi_size << " ranks."
                  << std::endl;
      }
      comm_world().barrier();

      uint64_t    delegate_threshold;
      std::string input_filename1, input_filename2;
      uint64_t    partition_passes;
      double      gbyte_per_rank;
      uint64_t    chunk_size;
      bool        undirected;
      bool        scramble = true;

      parse_cmd_line(argc, argv, output_filename, backup_filename,
                     delegate_threshold, input_filename1, input_filename2,
                     gbyte_per_rank, partition_passes, chunk_size, undirected);

      if (mpi_rank == 0) {
        std::cout << "Ingesting graphs" << std::endl;
      }

      distributed_db ddb(db_create(), output_filename.c_str());

      typedef distributed_db::allocator<edge_data_type> edge_data_allocator;
      graph_type::edge_data<edge_data_type, edge_data_allocator> edge_data(ddb.get_allocator());

      // Setup Kronecker generator
      kronecker_edge_generator<edge_data_type> kron(
          input_filename1, input_filename2, scramble, undirected);
      bool has_edge_data = kron.has_edge_data();

      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }
      graph_type* graph = ddb.get_manager()->construct<graph_type>("graph_obj")(
          ddb.get_allocator(), MPI_COMM_WORLD, kron, kron.max_vertex_id(),
          delegate_threshold, partition_passes, chunk_size, edge_data);

      if (has_edge_data) {
        auto edge_data_ptr = ddb.get_manager()->
            construct<graph_type::edge_data<edge_data_type, edge_data_allocator>>("graph_edge_data_obj")(edge_data);
      }

      comm_world().barrier();
      if (mpi_rank == 0) {
        std::cout << "Graph Ready, Calculating Stats. " << std::endl;
      }

      // TODO: implement get_size() and get_free_memory() in Metall
      // for (int i = 0; i < mpi_size; i++) {
      //   if (i == mpi_rank) {
      //     double percent = double(segment_manager->get_free_memory()) /
      //                      double(segment_manager->get_size());
      //     std::cout << "[" << mpi_rank << "] "
      //               << segment_manager->get_free_memory() << "/"
      //               << segment_manager->get_size() << " = " << percent
      //               << std::endl;
      //   }
      //   comm_world().barrier();
      // }

      //    graph->print_graph_statistics();

      //
      // Calculate max degree
      //    uint64_t max_degree(0);
      //    for (auto citr = graph->controller_begin(); citr !=
      //    graph->controller_end(); ++citr) {
      //      max_degree = std::max(max_degree, graph->degree(*citr));
      //    }

      //    uint64_t global_max_degree = havoqgt_all_reduce(max_degree,
      //    std::greater<uint64_t>(), MPI_COMM_WORLD);

      comm_world().barrier();

      if (mpi_rank == 0) {
        //      std::cout << "Max Degree = " << global_max_degree << std::endl;
      }

      comm_world().barrier();
    }  // Complete build distributed_db
    if (backup_filename.size() > 0) {
      distributed_db::transfer(output_filename.c_str(),
                               backup_filename.c_str());
    }
    comm_world().barrier();
    if (comm_nl().rank() == 0) {
      sync();
    }
    comm_world().barrier();
  }  // END Main MPI
  ;
  return 0;
}
