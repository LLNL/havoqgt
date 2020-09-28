// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT


#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/kth_core.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <havoqgt/distributed_db.hpp>
#include <assert.h>

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>

#include <boost/interprocess/managed_heap_memory.hpp>

using namespace havoqgt;

void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
         << " -i <string>   - input graph base filename (required)\n"
         << " -b <string>   - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename, std::string& backup_filename) {
  if(comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }
  
  bool found_input_filename = false;
  
  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:b:h ")) != -1) {
     switch (c) {
       case 'h':  
         prn_help = true;
         break;
      case 'i':
         found_input_filename = true;
         input_filename = optarg;
         break;
      case 'b':
         backup_filename = optarg;
         break;
      default:
         std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
         prn_help = true;
         break;
     }
   } 
   if (prn_help || !found_input_filename) {
     usage();
     exit(-1);
   }
}

int main(int argc, char** argv) {
  typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

  int mpi_rank(0), mpi_size(0);

  havoqgt::init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  
  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
     //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  std::string graph_input;
  std::string backup_filename;
  
  parse_cmd_line(argc, argv, graph_input, backup_filename);


  MPI_Barrier(MPI_COMM_WORLD);
  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  distributed_db ddb(db_open_read_only(), graph_input.c_str());

  auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }
  //graph->print_graph_statistics();
  MPI_Barrier(MPI_COMM_WORLD);
 
  //graph_type::vertex_data<kth_core_data, std::allocator<kth_core_data> >  vert_kcore_data(*graph);

  // uint64_t count_alive = 0;
  // uint64_t core = 0;
  // do {
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   double time_start = MPI_Wtime();
  //   count_alive = kth_core(*graph, vert_kcore_data, core++);
  //   MPI_Barrier(MPI_COMM_WORLD);
  //   double time_end = MPI_Wtime();
  //
  //   if(mpi_rank == 0){
  //     std::cout << "Core " << core-1 << " size = " << count_alive << ", time = " << time_end - time_start << std::endl;
  //   }
  // } while(count_alive);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    kth_core(*graph);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    if(mpi_rank == 0){
      std::cout << "Total Time = " << time_end - time_start << std::endl;
    }
  
  }  // END Main MPI
  ;

  return 0;
}
