/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */
#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setfill, std::setw
#include <random>
#include <chrono>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/graphstore/csr/csr_graph.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#include "dynamicgraphstore_bench.hpp"
#include "page_rank_bench.hpp"

using index_type      = uint64_t;
using vertex_id_type     = uint64_t;

using vertex_prop_type = graphstore::utility::packed_pair<double, double>;
using edge_prop_type = unsigned char;
using graphstore_type = csr_graph::prop_csr_graph_container<index_type, vertex_prop_type, vertex_id_type, edge_prop_type>;

template <typename graphstore_type, typename vertex_id_type>
static void run_page_rank_sync (graphstore_type& graphstore, const std::size_t num_vertices, const double damping_factor, const int num_loops)
{

  /// ---- init ---- ///
  auto tic_init = graphstore::utility::duration_time();
  graphstore.init_vertex_property(vertex_prop_type(0.0, 0.0));
  std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

  for (int k = 0; k < num_loops; ++k) {
    auto tic_step = graphstore::utility::duration_time();
    for (auto vrt_itr = graphstore.vertices_begin(), vrt_end = graphstore.vertices_end();
         vrt_itr != vrt_end;
         ++vrt_itr) {
      const double pr = graphstore.vertex_property_data(vrt_itr).first;
      const size_t degree = graphstore.degree(vrt_itr);
      for (auto adj_itr = graphstore.adjacent_edge_begin(vrt_itr), adj_end = graphstore.adjacent_edge_end(vrt_itr);
           adj_itr != adj_end;
           ++adj_itr) {
        graphstore.vertex_property_data(*adj_itr).second += (pr / degree) * damping_factor;
      }
    }
    for (auto vrt_itr = graphstore.vertices_begin(), vrt_end = graphstore.vertices_end();
         vrt_itr != vrt_end;
         ++vrt_itr) {
      graphstore.vertex_property_data(vrt_itr).first = graphstore.vertex_property_data(vrt_itr).second + (1.0 - damping_factor) / num_vertices;
      graphstore.vertex_property_data(vrt_itr).second = 0.0;
    }
    std::cout << "\n[" << k << "] : progress (sec.):\t"  << graphstore::utility::duration_time_sec(tic_step) << std::endl;
  }
}

/// Avoid linker errors with template function
template void run_page_rank_sync<graphstore_type, vertex_id_type>(graphstore_type&, const std::size_t, const double, const int);


std::vector<std::string> fname_edge_list_;
vertex_id_type max_vertex_id_ = 0;
size_t num_edges_ = 0;
int num_loop_ = 10;
double damping_factor_ = 1.0;

void parse_options(int argc, char **argv)
{

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  char c;

  while ((c = getopt (argc, argv, "g:S:o:E:v:m:l:")) != -1) {
    switch (c) {
      case 'v':
        max_vertex_id_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'm':
        num_edges_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'd':
        damping_factor_ = boost::lexical_cast<double>(optarg);
        break;

      case 'l':
        num_loop_ = boost::lexical_cast<int>(optarg);
        break;

      case 'E':
      {
        std::string fname(optarg);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while (std::getline(fin, line)) {
          fname_edge_list_.push_back(line);
        }
        break;
      }
    }
  }

  if (fname_edge_list_.empty()) HAVOQGT_ERROR_MSG("Edge lists are not given");
  for (auto itr : fname_edge_list_) {
    std::cout << "Load edge list from " << itr << std::endl;
  }

}


int main(int argc, char* argv[])
{

  parse_options(argc, argv);

  graphstore_type* graph;
  std::cout << "\n--- Initializing edgelist ---" << std::endl;
  graphstore::utility::print_time();
  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    havoqgt::get_environment();
    if (mpi_rank == 0) {
      havoqgt::parallel_edge_list_reader edge_list(fname_edge_list_);
      std::cout << "\n--- Constructing csr graph ---" << std::endl;
      graphstore::utility::print_time();
      graph = new graphstore_type(max_vertex_id_ + 1, num_edges_);
      graph->construct(edge_list);
    }
  }

  run_page_rank<graphstore_type, vertex_id_type>(*graph, max_vertex_id_ + 1, num_edges_, damping_factor_, num_loop_);

  delete graph;
}
