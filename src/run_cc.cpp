// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/connected_components.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <assert.h>
#include <havoqgt/distributed_db.hpp>

#include <algorithm>
#include <deque>
#include <functional>
#include <string>
#include <utility>

#include <boost/interprocess/managed_heap_memory.hpp>

using namespace havoqgt;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr
        << "Usage: -i <string> -s <int>\n"
        << " -i <string>   - input graph base filename (required)\n"
        << " -b <string>   - backup graph base filename.  If set, "
           "\"input\" graph will be deleted if it exists\n"
        << " -w <double>   - Filter edge weights greater than weight value \n"
        << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename,
                    std::string& backup_filename, bool& has_edge_filter,
                    double& edge_filter, std::string& output_base_filename,
                    size_t& print_threshold) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;
  has_edge_filter           = false;
  print_threshold           = 0;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:b:w:o:p:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 'i':
        found_input_filename = true;
        input_filename       = optarg;
        break;
      case 'b':
        backup_filename = optarg;
        break;
      case 'w':
        has_edge_filter = true;
        edge_filter     = atof(optarg);
        break;
      case 'o':
        output_base_filename = optarg;
        break;
      case 'p':
        print_threshold = atoi(optarg);
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
  typedef double                                    edge_data_type;
  typedef distributed_db::allocator<edge_data_type> edge_data_allocator_type;

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
    std::string output_base_filename;
    bool        has_edge_filter;
    double      edge_filter;
    size_t      print_threshold;

    parse_cmd_line(argc, argv, graph_input, backup_filename, has_edge_filter,
                   edge_filter, output_base_filename, print_threshold);

    MPI_Barrier(MPI_COMM_WORLD);
    if (backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    distributed_db ddb(db_open_read_only(), graph_input.c_str());

    auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    // graph->print_graph_statistics();
    MPI_Barrier(MPI_COMM_WORLD);

    graph_type::vertex_data<graph_type::vertex_locator,
                            std::allocator<graph_type::vertex_locator>>
        cc_data(*graph);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    if (!has_edge_filter) {
      connected_components(graph, cc_data);
    } else {
      auto finder = ddb.get_manager()
                        ->find<graph_type::edge_data<edge_data_type,
                                                     edge_data_allocator_type>>(
                            "graph_edge_data_obj");
      if (finder.second == false) {
        throw std::runtime_error("Edge weights not found");
      }
      auto edge_weights = finder.first;
      auto v_predicate  = [](const auto& vi) { return true; };
      auto e_predicate  = [edge_weights, edge_filter](auto e) {
        if ((*edge_weights)[e] > edge_filter) {
          return false;
        } else {
          return true;
        }
      };
      connected_components(graph, cc_data, v_predicate, e_predicate);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();

    if (print_threshold > 0) {
      std::map<uint64_t, uint64_t> cc_count;
      for (auto vitr = graph->vertices_begin(); vitr != graph->vertices_end();
           ++vitr) {
        if (graph->degree(*vitr) > 0) {
          cc_count[graph->locator_to_label(cc_data[*vitr])]++;
        }
      }

      for (auto citr = graph->controller_begin();
           citr != graph->controller_end(); ++citr) {
        if (graph->degree(*citr) > 0) {
          cc_count[graph->locator_to_label(cc_data[*citr])]++;
        }
      }

      auto     cc_count_itr = cc_count.begin();
      uint64_t largest_cc   = 0;
      uint64_t num_ccs      = 0;
      size_t   ccs_printed(0);
      while (!detail::global_iterator_range_empty(cc_count_itr, cc_count.end(),
                                                  MPI_COMM_WORLD)) {
        uint64_t local_next_cc = (cc_count_itr != cc_count.end())
                                     ? cc_count_itr->first
                                     : std::numeric_limits<uint64_t>::max();
        uint64_t global_next_cc = mpi_all_reduce(
            local_next_cc, std::less<uint64_t>(), MPI_COMM_WORLD);
        assert(global_next_cc != std::numeric_limits<uint64_t>::max());
        uint64_t local_count =
            (local_next_cc == global_next_cc) ? (cc_count_itr++)->second : 0;
        uint64_t global_cc_count =
            mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
        if (mpi_rank == 0) {
          std::cout << "CC " << global_next_cc << ", size = " << global_cc_count
                    << std::endl;
        }
        largest_cc = std::max(global_cc_count, largest_cc);
        num_ccs++;
        if (++ccs_printed >= print_threshold) {
          break;
        }
      }
      if (mpi_rank == 0) {
        std::cout << "Num CCs = " << num_ccs
                  << ", largest CC (approx) = " << largest_cc
                  << ", Traversal Time = " << time_end - time_start
                  << std::endl;
      }
    }
    if (output_base_filename.size() > 0) {
      std::stringstream output_filename;
      output_filename << output_base_filename << "_" << mpi_rank;
      std::ofstream ofs(output_filename.str().c_str());

      for (auto vitr = graph->vertices_begin(); vitr != graph->vertices_end();
           ++vitr) {
        if (graph->degree(*vitr) > 0) {
          ofs << graph->locator_to_label(*vitr) << " "
              << graph->locator_to_label(cc_data[*vitr]) << "\n";
        }
      }

      for (auto citr = graph->controller_begin();
           citr != graph->controller_end(); ++citr) {
        if (graph->degree(*citr) > 0) {
          ofs << graph->locator_to_label(*citr) << " "
              << graph->locator_to_label(cc_data[*citr]) << "\n";
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }  // END Main MPI
  return 0;
}
