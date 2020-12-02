// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include <havoqgt/cache_utilities.hpp>

#include <assert.h>
#include <unistd.h>
#include <algorithm>
#include <deque>
#include <fstream>  // std::ifstream
#include <functional>
#include <iostream>
#include <utility>

// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

typedef double                                    edge_data_type;
typedef distributed_db::allocator<edge_data_type> edge_data_allocator_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr
        << "Usage: -o <string> -d <int> [file ...]\n"
        << " -o <string>   - output graph base filename (required)\n"
        << " -b <string>   - backup graph base filename \n"
        << " -d <int>      - delegate threshold (Default is 1048576)\n"
        << " -h            - print help and exit\n"
        << " -p <int>      - number of Low & High partition passes (Default is "
           "1)\n"
        << " -f <float>    - Gigabytes reserved per rank (Default is 0.25)\n"
        << " -c <int>      - Edge partitioning chunk size (Defulat is 8192)\n"
        << " -u <bool>     - Treat edgelist as undirected (Default is 0)\n"
        << "[file ...] - list of edge list files to ingest\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& output_filename,
                    std::string& backup_filename, uint64_t& delegate_threshold,
                    std::vector<std::string>& input_filenames,
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
  input_filenames.clear();
  gbyte_per_rank   = 0.25;
  partition_passes = 1;
  chunk_size       = 8 * 1024;
  undirected       = false;

  int  c;
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

  for (int index = optind; index < argc; index++) {
    /// std::cout << "Input file = " << argv[index] << std::endl;
    input_filenames.push_back(argv[index]);
  }
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

      uint64_t                 delegate_threshold;
      std::vector<std::string> input_filenames;
      uint64_t                 partition_passes;
      double                   gbyte_per_rank;
      uint64_t                 chunk_size;
      bool                     undirected;

      parse_cmd_line(argc, argv, output_filename, backup_filename,
                     delegate_threshold, input_filenames, gbyte_per_rank,
                     partition_passes, chunk_size, undirected);

      if (mpi_rank == 0) {
        std::cout << "Ingesting graph from " << input_filenames.size()
                  << " files." << std::endl;
      }

      distributed_db ddb(db_create(), output_filename.c_str());

      auto edge_data_ptr =
          ddb.get_manager()
              ->construct<graph_type::edge_data<edge_data_type,
                                                edge_data_allocator_type>>(
                  "graph_edge_data_obj")(ddb.get_allocator());

      // Setup edge list reader
      havoqgt::parallel_edge_list_reader<edge_data_type> pelr(input_filenames,
                                                              undirected);
      bool has_edge_data = pelr.has_edge_data();

      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }
      graph_type* graph = ddb.get_manager()->construct<graph_type>("graph_obj")(
          ddb.get_allocator(), MPI_COMM_WORLD, pelr, pelr.max_vertex_id(),
          delegate_threshold, partition_passes, chunk_size, *edge_data_ptr);

      if (!has_edge_data) {
        ddb.get_manager()
            ->destroy<graph_type::edge_data<edge_data_type,
                                            edge_data_allocator_type>>(
                "graph_edge_data_obj");
      }

      comm_world().barrier();
      if (mpi_rank == 0) {
        std::cout << "Graph Ready, Calculating Stats. " << std::endl;
      }

      // TODO: implement get_size() and get_free_memory() in Metall
      // for (int i = 0; i < mpi_size; i++) {
      //  if (i == mpi_rank) {
      //    double percent = double(ddb.get_manager()->get_free_memory()) /
      //    double(ddb.get_manager()->get_size());
      //    std::cout << "[" << mpi_rank << "] " <<
      //    ddb.get_manager()->get_free_memory()
      //              << "/" << ddb.get_manager()->get_size() << " = " <<
      //              percent << std::endl;
      //  }
      //  comm_world().barrier();
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
  } //END Main MPI
  return 0;
}
