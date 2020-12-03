// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/hdf5_edge_list_reader.hpp>
#include <havoqgt/distributed_db.hpp>

#include <iostream>
#include <assert.h>
#include <functional>
#include <unistd.h>
#include <memory>

// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

typedef double edge_data_type;
typedef std::allocator<edge_data_type> edge_data_allocator_type;

void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -o <string> -d <int> [file ...]\n"
              << " -o <string>   - output graph base filename (required)\n"
              << " -s <string>   - key name of source vertex list (required)\n"
              << " -t <string>   - key name of destination vertex list (required)\n"
              << " -b <string>   - backup graph base filename \n"
              << " -d <int>      - delegate threshold (Default is 1048576)\n"
              << " -h            - print help and exit\n"
              << " -p <int>      - number of Low & High partition passes (Default is 1)\n"
              << " -f <float>    - Gigabytes reserved per rank (Default is 0.25)\n"
              << " -c <int>      - Edge partitioning chunk size (Defulat is 8192)\n"
              << " -u <bool>     - Treat edgelist as undirected (Default is 0)\n"
              << "[file ...] - list of edge list files to ingest\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& output_filename, std::string& backup_filename,
                    uint64_t& delegate_threshold, std::vector< std::string >& input_filenames,
                    double& gbyte_per_rank, uint64_t& partition_passes, uint64_t& chunk_size,
                    bool& undirected, std::string &src_key, std::string &dst_key) {
  if(comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_output_filename = false;
  delegate_threshold = 1048576;
  input_filenames.clear();
  gbyte_per_rank = 0.25;
  partition_passes = 1;
  chunk_size = 8*1024;
  undirected = false;

  int c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "o:d:p:f:c:b:u:hs:t: ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 'd':
        delegate_threshold = atoll(optarg);
        break;
      case 'o':
        found_output_filename = true;
        output_filename = optarg;
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
      case 's':
        src_key = optarg;
        break;
      case 't':
        dst_key = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_output_filename || src_key.empty() || dst_key.empty()) {
    usage();
    exit(-1);
  }

  for (int index = optind; index < argc; index++) {
    ///std::cout << "Input file = " << argv[index] << std::endl;
    input_filenames.push_back(argv[index]);
  }
}

auto read_edge(const std::vector<std::string> &input_filenames,
               const std::string &src_key,
               const std::string &dst_key,
               const bool undirected) {
  const int mpi_rank = havoqgt::comm_world().rank();
  const int mpi_size = havoqgt::comm_world().size();

  std::vector<std::tuple<uint64_t, uint64_t, edge_data_type>> edge_pair_list;
  for (std::size_t i = 0; i < input_filenames.size(); ++i) {
    if (i % mpi_size == mpi_rank) {
      havoqgt::hdf5_edge_list_reader reader;
      if (!reader.read(input_filenames[i], src_key, dst_key)) {
        std::cerr << "Failed to read edge: " << input_filenames[i] << std::endl;
      }
      const auto &read_edge_lists = reader.edges();
      assert(read_edge_lists.first.size() == read_edge_lists.second.size());
      for (std::size_t e = 0; e < read_edge_lists.first.size(); ++e) {
        edge_pair_list.emplace_back(read_edge_lists.first[e], read_edge_lists.first[e], edge_data_type());
        if (undirected)
          edge_pair_list.emplace_back(read_edge_lists.first[e], read_edge_lists.first[e], edge_data_type());
      }
    }
  }

  return edge_pair_list;
}

uint64_t find_global_max_vertex_id(const std::vector<std::tuple<uint64_t, uint64_t, edge_data_type>> &edge_list) {
  uint64_t local_max_vertex = 0;
  for (const auto &edge : edge_list) {
    local_max_vertex = std::max(std::get<0>(edge), local_max_vertex);
    local_max_vertex = std::max(std::get<1>(edge), local_max_vertex);
  }

  const auto global_max_vertex = havoqgt::mpi_all_reduce(local_max_vertex, std::greater<uint64_t>(), MPI_COMM_WORLD);

  return global_max_vertex;
}

int main(int argc, char** argv) {

  init(&argc, &argv);
  {
    std::string                output_filename;
    std::string                backup_filename;

    { // Build Distributed_DB
      int mpi_rank = comm_world().rank();
      int mpi_size = comm_world().size();

      if (mpi_rank == 0) {
        std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      }
      comm_world().barrier();

      uint64_t                   delegate_threshold;
      std::vector< std::string > input_filenames;
      uint64_t                   partition_passes;
      double                     gbyte_per_rank;
      uint64_t                   chunk_size;
      bool                       undirected;
      std::string src_key;
      std::string dst_key;

      parse_cmd_line(argc, argv, output_filename, backup_filename, delegate_threshold, input_filenames, gbyte_per_rank, partition_passes, chunk_size, undirected, src_key, dst_key);

      if (mpi_rank == 0) {
        std::cout << "Ingesting graph from " << input_filenames.size() << " files." << std::endl;
      }

      distributed_db ddb(db_create(), output_filename.c_str());
      graph_type::edge_data<edge_data_type, edge_data_allocator_type> dummy_edge_data;

      //Setup edge list reader
      const auto edge_list = read_edge(input_filenames, src_key, dst_key, undirected);
      const auto max_vertex_id = find_global_max_vertex_id(edge_list);
      const bool has_edge_data = false;

      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }
      graph_type *graph = ddb.get_manager()->construct<graph_type>
          ("graph_obj")
          (ddb.get_allocator(), MPI_COMM_WORLD,edge_list, max_vertex_id, delegate_threshold, partition_passes, chunk_size, dummy_edge_data);

      comm_world().barrier();
      if (mpi_rank == 0) {
        std::cout << "Graph Ready, Calculating Stats. " << std::endl;
      }

      comm_world().barrier();
    } // Complete build distributed_db
    if(backup_filename.size() > 0) {
      distributed_db::transfer(output_filename.c_str(), backup_filename.c_str());
    }
    comm_world().barrier();
  } //END Main MPI

  return 0;
}
