// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/upper_triangle_edge_generator.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <iostream>
#include <assert.h>
#include <deque>
#include <utility>
#include <algorithm>
#include <functional>
#include <fstream>      // std::ifstream


// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;

void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -s <int> -d <int> -o <string>\n"
         << " -s <int>    - RMAT graph Scale (default 17)\n"
         << " -d <int>    - delegate threshold (Default is 1048576)\n"
         << " -o <string> - output graph base filename\n"
         << " -b <string>   - backup graph base filename \n"
         << " -p <int>    - number of Low & High partition passes (Default is 1)\n"
         << " -f <float>  - Gigabytes reserved per rank (Default is 0.25)\n"
         << " -c <int>      - Edge partitioning chunk size (Defulat is 8192)\n"
         << " -h          - print help and exit\n\n";
         
  }
}

void parse_cmd_line(int argc, char** argv, uint64_t& scale, uint64_t& delegate_threshold, 
                    std::string& output_filename, std::string& backup_filename, double& gbyte_per_rank, 
                    uint64_t& partition_passes, uint64_t& chunk_size) {
  if(comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }
  
  bool found_output_filename = false;
  scale = 17;
  delegate_threshold = 1048576;
  gbyte_per_rank = 0.25;
  partition_passes = 1;
  chunk_size = 8*1024;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:d:o:b:p:f:c:h ")) != -1) {
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
      default:
         std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
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
  typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

  int mpi_rank(0), mpi_size(0);

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

    uint64_t      num_vertices = 1;
    uint64_t      vert_scale;
    uint64_t      hub_threshold;
    uint64_t      partition_passes;
    double        gbyte_per_rank;
    uint64_t      chunk_size;
        
    parse_cmd_line(argc, argv, vert_scale, hub_threshold, output_filename, backup_filename, 
                   gbyte_per_rank, partition_passes, chunk_size);

    num_vertices <<= vert_scale;
    if (mpi_rank == 0) {
      std::cout << "Building Graph500"<< std::endl
        << "Building graph Scale: " << vert_scale << std::endl
        << "Hub threshold = " << hub_threshold << std::endl
        << "File name = " << output_filename << std::endl
        << "Reserved Gigabytes per Rank = " << gbyte_per_rank << std::endl
        << "High/Low partition passes = " << partition_passes << std::endl; 
    }

    distributed_db ddb(db_create(), output_filename.c_str(), gbyte_per_rank);

    //Generate RMAT graph
    uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
    havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                      vert_scale, num_edges_per_rank,
                                      0.57, 0.19, 0.19, 0.05, true, true);


    if (mpi_rank == 0) {
      std::cout << "Generating new graph." << std::endl;
    }
    auto graph = ddb.get_manager()->construct<graph_type>
        ("graph_obj")
        (ddb.get_allocator(), MPI_COMM_WORLD, rmat, rmat.max_vertex_id(), hub_threshold, partition_passes, chunk_size);


    comm_world().barrier();
    if (mpi_rank == 0) {
      std::cout << "Graph Ready, Calculating Stats. " << std::endl;
    }

    // TODO: implement get_size() and get_free_memory() in Metall
    // for (int i = 0; i < mpi_size; i++) {
    //  if (i == mpi_rank) {
    //    double percent = double(ddb.get_manager()->get_free_memory()) /
    //    double(ddb.get_manager()->get_size());
    //    std::cout << "[" << mpi_rank << "] " << ddb.get_manager()->get_free_memory()
    //              << "/" << ddb.get_manager()->get_size() << " = " << percent << std::endl;
    //  }
    //  comm_world().barrier();
    // }

    graph->print_graph_statistics();


    //
    // Calculate max degree
    uint64_t max_degree(0);
    for (auto citr = graph->controller_begin(); citr != graph->controller_end(); ++citr) {
      max_degree = std::max(max_degree, graph->degree(*citr));
    }

    uint64_t global_max_degree = mpi_all_reduce(max_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);

    comm_world().barrier();

    if (mpi_rank == 0) {
      std::cout << "Max Degree = " << global_max_degree << std::endl;
    }

    comm_world().barrier();
    } // Complete build distributed_db
    if(backup_filename.size() > 0) {
      distributed_db::transfer(output_filename.c_str(), backup_filename.c_str());
    }
    comm_world().barrier();
    if(comm_nl().rank() == 0) {
      sync();
    }
    comm_world().barrier();
  } //END Main MPI
  ;
  return 0;
}
