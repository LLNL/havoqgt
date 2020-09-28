// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

#include <havoqgt/single_source_shortest_path.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <assert.h>
#include <havoqgt/distributed_db.hpp>

#include <algorithm>
#include <deque>
#include <functional>
#include <string>
#include <utility>

using namespace havoqgt;

typedef double edge_data_type;

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
  typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

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

    distributed_db ddb(havoqgt::db_open_read_only(), graph_input.c_str());

    auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    typedef distributed_db::allocator<edge_data_type> edge_data_allocator;
    auto edge_data_qry = ddb.get_manager()->find<graph_type::edge_data<edge_data_type, edge_data_allocator>>("graph_edge_data_obj");
    if (edge_data_qry.second == false) {
      if (mpi_rank == 0) {
        std::cout << "ERROR, edge weights not found" << std::endl;
      }
      abort();
    }
    auto edge_data_ptr = edge_data_qry.first;

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // SSSP Experiment
    {
      graph_type::vertex_data<edge_data_type, std::allocator<edge_data_type>>
          sssp_path_data(*graph);

      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0) {
        std::cout << "SSSP data allocated.  Starting SSSP from vertex "
                  << source_vertex << std::endl;
      }

      //  Run SSSP experiments
      double                     time(0);
      int                        count(0);
      uint64_t                   isource = source_vertex;
      graph_type::vertex_locator source = graph->label_to_locator(isource);
      uint64_t                   global_degree(0);
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
        if (global_degree == 0) ++isource;
      } while (global_degree == 0);
      if (uint32_t(mpi_rank) == source.owner()) {
        if (isource != source_vertex) {
          std::cout << "Vertex " << source_vertex
                    << " has a degree of 0.   New source vertex = " << isource
                    << std::endl;
        } else {
          std::cout << "Starting vertex = " << isource << std::endl;
        }
        // std::cout << "delegate? = " << source.is_delegate() << std::endl;
        // std::cout << "local_id = " << source.local_id() << std::endl;
        std::cout << "degree = " << graph->degree(source) << std::endl;
      }

      sssp_path_data.reset(std::numeric_limits<edge_data_type>::max());

      MPI_Barrier(MPI_COMM_WORLD);
      double time_start = MPI_Wtime();
      havoqgt::single_source_shortest_path(*graph, sssp_path_data,
                                           *edge_data_ptr, source);
      MPI_Barrier(MPI_COMM_WORLD);
      double time_end = MPI_Wtime();

      uint64_t local_visited_count = 0;
      double   local_max_distance = 0.0f;
      // process low degree verts
      for (auto vitr = graph->vertices_begin(); vitr != graph->vertices_end();
           ++vitr) {
        if (sssp_path_data[*vitr] !=
            std::numeric_limits<edge_data_type>::max()) {
          ++local_visited_count;
          local_max_distance =
              std::max(local_max_distance, sssp_path_data[*vitr]);
        }
      }
      // process high-degree controllers
      for (auto citr = graph->controller_begin();
           citr != graph->controller_end(); ++citr) {
        if (sssp_path_data[*citr] !=
            std::numeric_limits<edge_data_type>::max()) {
          ++local_visited_count;
          local_max_distance =
              std::max(local_max_distance, sssp_path_data[*citr]);
        }
      }

      uint64_t global_count_visited = mpi_all_reduce(
          local_visited_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t global_max_distance = mpi_all_reduce(
          local_max_distance, std::greater<double>(), MPI_COMM_WORLD);
      if (mpi_rank == 0) {
        std::cout << "Total vertices reached/visited = " << global_count_visited
                  << std::endl;
        std::cout << "Max distance = " << global_max_distance << std::endl;
        std::cout << "Elapsed time = " << time_end - time_start << std::endl;
      }
    }  // End BFS Test
  }    // END Main MPI
  ;

  return 0;
}
