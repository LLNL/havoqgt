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
#include "bfs_bench.hpp"

using index_type      = uint64_t;
using vertex_type     = uint64_t;
#if BFS_USE_BITMAP
using graphstore_type = csr_graph::csr_graph_container<index_type, vertex_type>;
#else
using vertex_prop_type = bool;
using edge_prop_type = unsigned char;
#if 1
using graphstore_type = csr_graph::prop_csr_graph_container<index_type, vertex_prop_type, vertex_type, edge_prop_type>;
#else
using graphstore_type = csr_graph::hyper_prop_csr_graph_container<index_type, vertex_prop_type, vertex_type, edge_prop_type>;
#endif
#endif

template <typename graphstore_type, typename vertex_type>
void run_bfs_sync (graphstore_type& graph,
                   trv_inf<vertex_type>& inf,
                   std::queue<vertex_type>& frontier_queue,
                   std::queue<vertex_type>& next_queue,
                   vertex_type& root_vrtx)
{

  /// ---- init inf ---- ///
  auto tic_init = graphstore::utility::duration_time();
#if BFS_USE_BITMAP
  inf.init(true);
  uint64_t* const visited = inf.visited;
  visited[bitmap_global_pos(root_vrtx)] |= 0x1 << bitmap_local_pos(root_vrtx);
#else
  inf.init(false);
  bool init(false);
  graph.init_vertex_property(init);
  graph.vertex_property(root_vrtx) = true;
#endif
  /// inf.tree[root_vrtx] = root_vrtx;
  std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

  /// --- BFS main loop -------- ///
  size_t level = 0;
  while (true) {
    std::cout << "Lv. " << level << ", size of frontier queue =\t" << frontier_queue.size() << std::endl;

    /// loop for current frontier
    while (!frontier_queue.empty()) {
      vertex_type src = frontier_queue.front();
      // std::cout << "src " << src << std::endl;
      frontier_queue.pop();
      ++inf.count_visited_vertices;


      /// push adjacent vertices to the next queue
      for (auto edge = graph.adjacencylist(src), end = graph.adjacencylist_end(src); edge != end; ++edge) {
        const vertex_type dst = *edge;
#if BFS_USE_BITMAP
          const bool is_visited = visited[bitmap_global_pos(dst)] & (0x1 << bitmap_local_pos(dst));
#else
        bool& is_visited = graph.vertex_property(dst);
#endif
        if (!is_visited) {
          next_queue.push(dst);
#if BFS_USE_BITMAP
          visited[bitmap_global_pos(dst)] |= 0x1 << bitmap_local_pos(dst);
#else
          graphstore.vertex_property_data(dst) = true;
#endif
        }
        ++(inf.count_visited_edges);
      }
    }  /// end of loop for a frontier

    if (next_queue.empty()) break; /// termination condition
    frontier_queue.swap(next_queue);
    ++level;
  } /// end of BFS loop
}

/// Avoid linker errors with template function
template void run_bfs_sync<graphstore_type, vertex_type>(graphstore_type&,
                                                         trv_inf<vertex_type>&,
                                                         std::queue<vertex_type>&,
                                                         std::queue<vertex_type>&,
                                                         vertex_type&);


std::string fname_graph_;
std::string fname_segmentfile_;
size_t segment_size_log2_ = 30;
std::vector<std::string> fname_edge_list_;
vertex_type max_vertex_id_ = 0;
size_t num_edges_ = 0;
std::vector<vertex_type> source_list_;

void parse_options(int argc, char **argv)
{

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  char c;

  while ((c = getopt (argc, argv, "g:S:o:E:v:m:r:")) != -1) {
    switch (c) {
      case 'g':
        fname_graph_ = optarg;
        break;
      case 'o':
        fname_segmentfile_ = optarg;
        break;
      case 'S':
        segment_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;
      case 'v':
        max_vertex_id_ = boost::lexical_cast<size_t>(optarg);
        break;
      case 'm':
        num_edges_ = boost::lexical_cast<size_t>(optarg);
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

      case 'r':
      {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_list_.push_back(boost::lexical_cast<size_t>(buf));
        break;
      }
    }
  }

  if (!fname_edge_list_.empty())
    std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
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

  /// ---------- Graph Traversal --------------- ///
  std::cout << "\n<Run BFS>" << std::endl;
  if (source_list_.empty())
    generate_bfs_sources(4, max_vertex_id_, source_list_);
  run_bfs(*graph, max_vertex_id_, max_vertex_id_, source_list_);

  delete graph;
}
