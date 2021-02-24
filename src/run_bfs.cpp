// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT


#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/distributed_db.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <assert.h>

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>

using namespace havoqgt;

void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
         << " -i <string>   - input graph base filename (required)\n"
         << " -b <string>   - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
         << " -s <int>      - Source vertex of BFS (Default is 0)\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename, std::string& backup_filename, uint64_t& source_vertex) {
  if(comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
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
  }
  MPI_Barrier(MPI_COMM_WORLD);


  std::string graph_input;
  std::string backup_filename;
  uint64_t source_vertex = 0;

  parse_cmd_line(argc, argv, graph_input, backup_filename, source_vertex);


  MPI_Barrier(MPI_COMM_WORLD);
  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  distributed_db ddb(db_open_read_only(), graph_input);

  auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }
  //graph->print_graph_statistics();
  MPI_Barrier(MPI_COMM_WORLD);


  // BFS Experiments
  {

    graph_type::vertex_data<uint16_t, std::allocator<uint16_t> >                      bfs_level_data(*graph);
    graph_type::vertex_data<graph_type::vertex_locator, std::allocator<graph_type::vertex_locator> >  bfs_parent_data(*graph);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "BFS data allocated.  Starting BFS from vertex " << source_vertex << std::endl;
    }

    //  Run BFS experiments
    double time(0);
    int count(0);
    uint64_t isource = source_vertex;
      graph_type::vertex_locator source = graph->label_to_locator(isource);
      uint64_t global_degree(0);
      do {
        uint64_t local_degree = 0;
        source = graph->label_to_locator(isource);
        if (source.is_delegate()) {
          break;
        }
        if (uint32_t(mpi_rank) == source.owner()) {
          local_degree = graph->degree(source);
        }
        global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(),
            MPI_COMM_WORLD);
        if(global_degree == 0) ++isource;
      } while (global_degree == 0);
      if (uint32_t(mpi_rank) == source.owner()) {
        if(isource != source_vertex) {
          std::cout << "Vertex " << source_vertex << " has a degree of 0.   New source vertex = " << isource << std::endl;
        } else {
          std::cout << "Starting vertex = " << isource << std::endl;
        }
        std::cout << "delegate? = " << source.is_delegate() << std::endl;
        std::cout << "local_id = " << source.local_id() << std::endl;
        std::cout << "degree = " << graph->degree(source) << std::endl;
      }

      bfs_level_data.reset(std::numeric_limits<uint16_t>::max());

      MPI_Barrier(MPI_COMM_WORLD);
      double time_start = MPI_Wtime();
      havoqgt::breadth_first_search(graph, bfs_level_data, bfs_parent_data,
          source);
      MPI_Barrier(MPI_COMM_WORLD);
      double time_end = MPI_Wtime();

      uint64_t visited_total(0);
      for (uint64_t level = 0; level < std::numeric_limits<uint16_t>::max(); ++level) {
        uint64_t local_count(0);
        graph_type::vertex_iterator vitr;
        for (vitr = graph->vertices_begin();
             vitr != graph->vertices_end();
             ++vitr) {
          if (bfs_level_data[*vitr] == level) {
            ++local_count;
          }
        }

        // Count the controllers!
        graph_type::controller_iterator citr;
        for (citr = graph->controller_begin();
             citr != graph->controller_end();
             ++citr) {
          if (bfs_level_data[*citr] == level) {
            ++local_count;
          }
        }

        uint64_t global_count = mpi_all_reduce(local_count,
          std::plus<uint64_t>(), MPI_COMM_WORLD);
        visited_total += global_count;
        if (mpi_rank == 0 && global_count > 0) {
          std::cout << "Level " << level << ": " << global_count << std::endl;
        }
        if (global_count == 0) {
          break;
        }
      }  // end for level

      if (mpi_rank == 0) {
        if (visited_total > 1) {
          std::cout
            << "Visited total = " << visited_total << std::endl
            << "BFS Time = " << time_end - time_start << std::endl;
          time += time_end - time_start;
          ++count;
        }
      }
    if (mpi_rank == 0) {
      std::cout << "Count BFS = " << count << std::endl;
      std::cout << "AVERAGE BFS = " << time / double(count) << std::endl;
    }
  }  // End BFS Test
  }  // END Main MPI
  ;

  return 0;
}
