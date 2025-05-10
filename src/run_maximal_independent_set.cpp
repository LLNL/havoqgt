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
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
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
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>

//#include <havoqgt/single_source_shortest_path.hpp>
#include <independent_set/initialize.hpp>
//#include <independent_set/maximal_independent_set.cpp>
#include <independent_set/maximal_independent_set_asynchronous.cpp>
#include <independent_set/utility.hpp>

#include <assert.h>

#include <algorithm>
#include <deque>
#include <functional>
#include <string>
#include <utility>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

using namespace havoqgt;

//typedef double edge_data_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -i <string>   - input graph base filename (required)\n"
              << " -b <string>   - backup graph base filename.  If set, "
                 "\"input\" graph will be deleted if it exists\n"
              << " -s <int>      - Source vertex of BFS (Default is 0)\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename,
                    std::string& backup_filename, uint64_t& source_vertex) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;
  source_vertex = 0;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:b:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's':
        source_vertex = atoll(optarg);
        break;
      case 'i':
        found_input_filename = true;
        input_filename = optarg;
        break;
      case 'b':
        backup_filename = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
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
  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
  typedef havoqgt::delegate_partitioned_graph
    <typename segment_manager_t::template allocator<void>::type> graph_type;

  typedef typename graph_type::vertex_locator vloc_type;
  typedef typename graph_type::vertex_iterator vitr_type;

  int mpi_rank(0), mpi_size(0);

  havoqgt::init(&argc, &argv);
  {
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    
    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks."
                << std::endl;
      // print_system_info(false);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_filename;
    uint64_t    source_vertex = 0;

    parse_cmd_line(argc, argv, graph_input, backup_filename, source_vertex);

    MPI_Barrier(MPI_COMM_WORLD);
    if (backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

    graph_type* graph =
        ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    //auto edge_data_qry =
    //    ddb.get_segment_manager()
    //        ->find<graph_type::edge_data<
    //            edge_data_type,
    //            bip::allocator<edge_data_type, segment_manager_t>>>(
    //            "graph_edge_data_obj");
    //if (edge_data_qry.second == false) {
    //  if (mpi_rank == 0) {
    //    std::cout << "ERROR, edge weights not found" << std::endl;
    //  }
    //  abort();
    //}
    //auto edge_data_ptr = edge_data_qry.first;

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph loaded ready" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // maximal independent set
    {
    double time_start = MPI_Wtime();	    
    double time_end = MPI_Wtime();
    double elapsed_time = time_end - time_start;

    std::vector<double> step_time = {0, 0, 0, 0};

    using VertexID = uint64_t;
    using VertexPriority = uint64_t;
    using VertexIDPriorityMap = std::unordered_map<VertexID, VertexPriority>;

    using VertexPriorityCollection = graph_type::vertex_data<VertexPriority, 
      std::allocator<VertexPriority> >;

    using VertexIDPriorityMapCollection = 
      graph_type::vertex_data<VertexIDPriorityMap, 
      std::allocator<VertexIDPriorityMap> >;

    // application 

    time_start = MPI_Wtime();
  
    VertexPriorityCollection vertex_priority_list(*graph);
    VertexIDPriorityMapCollection higher_priority_neighbor_map(*graph);
    VertexIDPriorityMapCollection lower_priority_neighbor_map(*graph);

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // graph statistic
    
    double time_start_0 = MPI_Wtime();

    maximal_independent_set::graph_statistic<graph_type, VertexID,
      VertexPriority>(graph);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_end_0 = MPI_Wtime();
    step_time[0] = time_end_0 - time_start_0;
    if (mpi_rank == 0) {
      std::cout << "Elapsed time 0 (graph statistic) : " << step_time[0] << 
        std::endl;
    }

    // initialize vertex priorities 
     
    double time_start_1 = MPI_Wtime();

    maximal_independent_set::initialize_priority<graph_type, VertexID, 
      VertexPriority, VertexPriorityCollection, 
      VertexIDPriorityMapCollection> 
      (graph, vertex_priority_list, higher_priority_neighbor_map, 
      lower_priority_neighbor_map);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_end_1 = MPI_Wtime();
    step_time[1] = time_end_1 - time_start_1;
    if (mpi_rank == 0) {
      std::cout << "Elapsed time 1 (priority initialization) : " << 
        step_time[1] << 
        std::endl;
    }

    // compute maximal independent set
    
    double time_start_2 = MPI_Wtime(); 

    //maximal_independent_set::maximal_independent_set<graph_type,
    maximal_independent_set::maximal_independent_set_asynchronous<graph_type, 
      VertexID, VertexPriority, VertexPriorityCollection, 
      VertexIDPriorityMapCollection>    
      (graph, vertex_priority_list, higher_priority_neighbor_map, 
      lower_priority_neighbor_map);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_end_2 = MPI_Wtime();
    step_time[2] = time_end_2 - time_start_2;
    if (mpi_rank == 0) {
      std::cout << "Elapsed time 2 (maximal independent set computation) : " << 
        step_time[2] << 
        std::endl;
    }

    // total time

    MPI_Barrier(MPI_COMM_WORLD);
    time_end = MPI_Wtime();
    elapsed_time = time_end - time_start;
    step_time[3] = elapsed_time;
    if (mpi_rank == 0) {
      std::cout << "Total elapsed time: " << elapsed_time << std::endl;
    }

    // end application

    // results
    
    std::vector<size_t> local_vertex_status_count = {0, 0, 0};   
   
    for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
      ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_priority_list[vertex] == maximal_independent_set::in) {
        local_vertex_status_count[0]++;	      
      } else if (vertex_priority_list[vertex] == maximal_independent_set::out) {
        local_vertex_status_count[1]++;  	      
      } else {
	local_vertex_status_count[2]++;      
      }	      
    } // for		   

    for(vitr_type vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_priority_list[vertex] == maximal_independent_set::in) {
        local_vertex_status_count[0]++;	      
      } else if (vertex_priority_list[vertex] == maximal_independent_set::out) {
        local_vertex_status_count[1]++;  	      
      } else {
	local_vertex_status_count[2]++;      
      }  	    
    } // for		   

    MPI_Barrier(MPI_COMM_WORLD); 
    
    havoqgt::mpi_all_reduce_inplace(local_vertex_status_count, 
      std::plus<size_t>(), MPI_COMM_WORLD);    

    if (mpi_rank == 0) {
      std::cout << 
	"#IN vertices (maximal independent set): " << 
	local_vertex_status_count[0] <<
        ", #OUT vertices: " << local_vertex_status_count[1] <<
	", #UNDECIDED vertices (should be 0): " << 
	local_vertex_status_count[2] << std::endl;
    }	    
    
    }; // end of maximal independent set	     

  }; // end of main MPI  

  return 0;
}
