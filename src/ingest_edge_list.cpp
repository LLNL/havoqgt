/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <iostream>
#include <assert.h>
#include <deque>
#include <utility>
#include <algorithm>
#include <functional>
#include <fstream>      // std::ifstream
#include <unistd.h>


// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

typedef uint8_t edge_data_type;
  
void usage()  {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -o <string> -d <int> [file ...]\n"
         << " -o <string>   - output graph base filename (required)\n"
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
                    bool& undirected) {
  if(havoqgt_env()->world_comm().rank() == 0) {
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

   for (int index = optind; index < argc; index++) {
     ///std::cout << "Input file = " << argv[index] << std::endl;
     input_filenames.push_back(argv[index]);
   }
}


int main(int argc, char** argv) {

  int mpi_rank(0), mpi_size(0);

  havoqgt_init(&argc, &argv);
  {
    std::string                output_filename;
    std::string                backup_filename;
 
    { // Build Distributed_DB
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();
    havoqgt::get_environment();
    
    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      havoqgt::get_environment().print();
    }
    havoqgt_env()->world_comm().barrier();

    uint64_t                   delegate_threshold;
    std::vector< std::string > input_filenames;
    uint64_t                   partition_passes;
    double                     gbyte_per_rank;
    uint64_t                   chunk_size;
    bool                       undirected;
   
    parse_cmd_line(argc, argv, output_filename, backup_filename, delegate_threshold, input_filenames, gbyte_per_rank, partition_passes, chunk_size, undirected);

    if (mpi_rank == 0) {
      std::cout << "Ingesting graph from " << input_filenames.size() << " files." << std::endl;
    }

    havoqgt::distributed_db ddb(havoqgt::db_create(), output_filename.c_str(), gbyte_per_rank);

    segment_manager_t* segment_manager = ddb.get_segment_manager();
    bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

    graph_type::edge_data<edge_data_type, bip::allocator<edge_data_type, segment_manager_t>> edge_data(alloc_inst); 

    //Setup edge list reader
    havoqgt::parallel_edge_list_reader<edge_data_type> pelr(input_filenames, undirected);
    bool has_edge_data = pelr.has_edge_data();
    has_edge_data = true; // TODO: force to create edge data 

    if (mpi_rank == 0) {
      std::cout << "Generating new graph." << std::endl;
    }
    graph_type *graph = segment_manager->construct<graph_type>
        ("graph_obj")
        (alloc_inst, MPI_COMM_WORLD, pelr, pelr.max_vertex_id(), delegate_threshold, partition_passes, chunk_size,
         edge_data); 
    
    if (has_edge_data) {
      if (mpi_rank == 0) {  
        std::cout << "Graph with edge data." << std::endl;
      }
      graph_type::edge_data<edge_data_type, bip::allocator<edge_data_type, segment_manager_t>>* edge_data_ptr
      = segment_manager->construct<graph_type::edge_data<edge_data_type, bip::allocator<edge_data_type, segment_manager_t>>>
          ("graph_edge_data_obj")
          (edge_data);
    }

    havoqgt_env()->world_comm().barrier();
    if (mpi_rank == 0) {
      std::cout << "Graph Ready, Calculating Stats. " << std::endl;
    }

    for (int i = 0; i < mpi_size; i++) {
      if (i == mpi_rank) {
        double percent = double(segment_manager->get_free_memory()) /
        double(segment_manager->get_size());
        std::cout << "[" << mpi_rank << "] " << segment_manager->get_free_memory()
                  << "/" << segment_manager->get_size() << " = " << percent << std::endl;
      }
      havoqgt_env()->world_comm().barrier();
    }

//    graph->print_graph_statistics();

    //
    // Calculate max degree
//    uint64_t max_degree(0);
//    for (auto citr = graph->controller_begin(); citr != graph->controller_end(); ++citr) {
//      max_degree = std::max(max_degree, graph->degree(*citr));
//    }

//    uint64_t global_max_degree = havoqgt::mpi::mpi_all_reduce(max_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);

    havoqgt_env()->world_comm().barrier();

    if (mpi_rank == 0) {
//      std::cout << "Max Degree = " << global_max_degree << std::endl;
    }

    havoqgt_env()->world_comm().barrier();
    } // Complete build distributed_db
    if(backup_filename.size() > 0) {
      distributed_db::transfer(output_filename.c_str(), backup_filename.c_str());
    }
    havoqgt_env()->world_comm().barrier();
    if(havoqgt_env()->node_local_comm().rank() == 0) {
      sync();
    }
    havoqgt_env()->world_comm().barrier();
  } //END Main MPI
  havoqgt_finalize();
  return 0;
}
