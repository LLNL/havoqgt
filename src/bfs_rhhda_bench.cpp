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
#include <utility>
#include <algorithm>
#include <functional>
#include <fcntl.h>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include <havoqgt/graphstore/graphstore_rhhda.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#include "dynamicgraphstore_bench.hpp"
#include "bfs_bench.hpp"

#define DEBUG_MODE 0
#if DEBUG_MODE
 std::ofstream ofs_edges;
#endif

enum : size_t {
  midle_high_degree_threshold = 2 // must be more or equal than 1
};

/// --- typenames --- ///
using vertex_type           = uint64_t;
using vertex_meta_data_type = bool;
using edge_weight_type      = unsigned char;
using graphstore_type       = graphstore::graphstore_rhhda<vertex_type,
                                                           vertex_meta_data_type,
                                                           edge_weight_type,
                                                           segment_manager_type,
                                                           midle_high_degree_threshold>;

/// --- global variables --- ///
vertex_type max_vertex_id_ = 0;
size_t num_edges_ = 0;


template <typename Edges>
void constract_graph(mapped_file_type& mapped_file,
                     segment_manager_type *const segment_manager,
                     graphstore_type& graph_store,
                     Edges& edges, const size_t chunk_size)
{
  std::cout << "-- Disp status of before generation --" << std::endl;
  std::cout << "segment size (GB); " << get_segment_size(segment_manager) << std::endl;
  print_system_mem_usages();


  size_t count_inserted = 0;
  size_t loop_cnt = 0;
  double construction_time = 0;

  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  request_vector_type<vertex_type> update_request_vec = request_vector_type<vertex_type>();
  update_request_vec.reserve(chunk_size);

  auto global_start = graphstore::utility::duration_time();
  while (edges_itr != edges_itr_end) {
    std::cout << "[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);

    unsigned char dummy_weight = 0;
    auto local_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      max_vertex_id_ = std::max(max_vertex_id_, edge.first);
      max_vertex_id_ = std::max(max_vertex_id_, edge.second);
      count_inserted += graph_store.insert_edge(edge.first, edge.second, dummy_weight);
    }
    std::cout << "shrink to fit low table" << std::endl;
    graph_store.shrink_to_fit_low_table();
    double t = graphstore::utility::duration_time_sec(local_start);
    construction_time += t;
    std::cout << "progress (sec.): " << t << std::endl;

    ++loop_cnt;
  }
  std::cout << "sync mmap" << std::endl;
  flush_mmmap(mapped_file);
  std::cout << "sync di-mmap" << std::endl;
  sync_dimmap();
  const double whole_construction_time = graphstore::utility::duration_time_sec(global_start);

  num_edges_ = count_inserted;

  std::cout << "\n-- All edge updations done --" << std::endl;
  std::cout << "inserted edges : " << count_inserted << std::endl;
  std::cout << "construction time (insertion only) : " << construction_time << std::endl;
  std::cout << "whole construction time : " << whole_construction_time << std::endl;
  std::cout << "segment size (GB); " << get_segment_size(segment_manager) << std::endl;
  print_system_mem_usages();

}


template <typename graphstore_type, typename vertex_type>
  static void run_bfs_sync (
      graphstore_type& graph,
      trv_inf<vertex_type>& inf,
      std::queue<vertex_type>& frontier_queue,
      std::queue<vertex_type>& next_queue,
      vertex_type& start_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
#if BFS_USE_BITMAP
    inf.init(true);
    inf.is_visited[start_vrtx] = true;
#else
    inf.init(false);
    for (auto itr = graph.begin_low_edges(); !itr.is_end(); ++itr) {
      itr->value.first = false;
    }
    for (auto itr = graph.begin_mid_high_edges(); !itr.is_end(); ++itr) {
      itr->value.first = false;
    }
    graph.vertex_meta_data(start_vrtx) = true;
#endif
    std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

    /// --- BFS main loop -------- ///
    size_t level = 0;
    while (true) {
      std::cout << "Lv. " << level << ", size of frontier queue =\t" << frontier_queue.size() << std::endl;

      /// loop for current frontier
      while (!frontier_queue.empty()) {
        vertex_type src = frontier_queue.front();
        frontier_queue.pop();
        ++inf.count_visited_vertices;

        size_t count_visited_edges = 0;
        /// push adjacent vertices to the next queue
        for (auto edge = graph.find_low_edge(src); !edge.is_end(); ++edge) {
          const vertex_type dst = edge->second;
#if BFS_USE_BITMAP
          bool& is_visited = inf.is_visited[dst];
#else
          bool& is_visited = graph.vertex_meta_data(dst);
#endif
          if (!is_visited) {
            next_queue.push(dst);
            is_visited = true;
          }
          ++count_visited_edges;
        }
        inf.count_visited_edges += count_visited_edges;
        if (count_visited_edges > 0) continue;

        for (auto edge = graph.find_mid_high_edge(src); !edge.is_end(); ++edge) {
          const vertex_type dst = edge->key;
#if BFS_USE_BITMAP
          bool& is_visited = inf.is_visited[dst];
#else
          bool& is_visited = graph.vertex_meta_data(dst);
#endif
          if (!is_visited) {
            next_queue.push(dst);
            is_visited = true;
          }
          ++count_visited_edges;
        }
        inf.count_visited_edges += count_visited_edges;

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



/// --- option variables --- ///
std::string fname_segmentfile_;
std::vector<std::string> fname_edge_list_;
size_t segment_size_log2_ = 30;
std::vector<vertex_type> source_list_;

void parse_options(int argc, char **argv)
{

  char c;

  while ((c = getopt (argc, argv, "S:o:E:r:")) != -1) {
    switch (c) {
      case 'S':
        segment_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'o':
        fname_segmentfile_ = optarg;
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

}


int main(int argc, char** argv) {

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  parse_options(argc, argv);
  std::cout << "Midle-high degree threshold = " << midle_high_degree_threshold << std::endl;
  if (!fname_edge_list_.empty())
    std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
  for (auto itr : fname_edge_list_) {
    std::cout << "Load edge list from " << itr << std::endl;
  }


  /// --- create a segument file --- ///
  boost::interprocess::file_mapping::remove(fname_segmentfile_.c_str());
  std::cout << "\n<<Construct segment>>" << std::endl;
  std::cout << "Create and map a segument file" << std::endl;
  uint64_t graph_capacity = std::pow(2, segment_size_log2_);
  mapped_file_type mapped_file = mapped_file_type(
                                   boost::interprocess::create_only,
                                   fname_segmentfile_.c_str(),
                                   graph_capacity);
  std::cout << "Call posix_fallocate\n";
  fallocate(fname_segmentfile_.c_str(), graph_capacity, mapped_file);

  /// --- Get a segument manager --- ///
  std::cout << "\n<Get a segment manager>" << std::endl;
  segment_manager_type* segment_manager = mapped_file.get_segment_manager();
  print_system_mem_usages();


  /// --- Allocate graphstore_rhhda --- ///
  std::cout << "\n<Allocate graphstore_rhhda>" << std::endl;
  graphstore_type graph_store(segment_manager);

  /// --- Graph Construction --- ////
  std::cout << "\n<Construct graph>" << std::endl;
  havoqgt::havoqgt_init(&argc, &argv);
  {
    havoqgt::get_environment();

    havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_);
    constract_graph(
          mapped_file,
          segment_manager,
          graph_store,
          edgelist,
          static_cast<uint64_t>(std::pow(10, 6)));
  }


  /// ---------- Graph Traversal --------------- ///
  std::cout << "\n<Run BFS>" << std::endl;
  if (source_list_.empty())
    generate_bfs_sources(4, max_vertex_id_, source_list_);
  run_bfs(graph_store, max_vertex_id_, num_edges_, source_list_);

  return 0;
}
