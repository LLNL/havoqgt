#ifndef BFS_HPP
#define BFS_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <utility>
#include <vector>
#include <stack>
#include <queue>
#include <limits>
#include <fstream>
#include <ctime>

#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/graph_traversal/graph_trav_info.hpp>

template <typename graph_type, typename vertex_type, int type>
class bfs_core {

 public:

  bfs_core()=delete;

  static void run_bfs_sync (
      graph_type& graph,
      graph_trv_info::trv_inf<vertex_type>& inf,
      std::queue<vertex_type>& frontier_queue,
      std::queue<vertex_type>& next_queue,
      vertex_type& start_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
    inf.init(true);
    inf.is_visited[start_vrtx] = true;
    inf.tree[start_vrtx] = start_vrtx;
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
        for (auto edge = graph.begin_adjlist(src), end = graph.end_adjlist(src); edge != end; ++edge) {
          const vertex_type dst = *edge;
          if (!inf.is_visited[dst]) {
            next_queue.push(dst);
            inf.tree[dst] = src;
            inf.is_visited[dst] = true;
          }
          ++(inf.count_visited_edges);
        }
      }  /// end of loop for a frontier

      if (next_queue.empty()) break; /// termination condition
      frontier_queue.swap(next_queue);
      ++level;
    } /// end of BFS loop
  }
};

template <typename graph_type, typename vertex_type>
class bfs_core <graph_type, vertex_type, 1> {

 public:

  bfs_core()=delete;

  static void run_bfs_sync (
      graph_type& graph,
      graph_trv_info::trv_inf<vertex_type>& inf,
      std::queue<vertex_type>& frontier_queue,
      std::queue<vertex_type>& next_queue,
      vertex_type& start_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
    inf.init(false);
    for (auto itr = graph.begin_low_edges(); !itr.is_end(); ++itr) {
      itr->value.first = false;
    }
    for (auto itr = graph.begin_mid_high_edges(); !itr.is_end(); ++itr) {
      itr->value.first = false;
    }
    graph.vertex_meta_data(start_vrtx) = true;
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
          bool& is_visited = graph.vertex_meta_data(dst);
          if (!is_visited) {
            next_queue.push(dst);
            /// inf.tree[dst] = src;
            is_visited = true;
          }
          ++count_visited_edges;
        }
        inf.count_visited_edges += count_visited_edges;
        if (count_visited_edges > 0) continue;

        for (auto edge = graph.find_mid_high_edge(src); !edge.is_end(); ++edge) {
          const vertex_type dst = edge->key;
          bool& is_visited = graph.vertex_meta_data(dst);
          if (!is_visited) {
            next_queue.push(dst);
            /// inf.tree[dst] = src;
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
};


template <typename graph_type, typename vertex_type>
class bfs_core <graph_type, vertex_type, 2> {

 public:

  bfs_core()=delete;

  static void run_bfs_sync (
      graph_type& graph,
      graph_trv_info::trv_inf<vertex_type>& inf,
      std::queue<vertex_type>& frontier_queue,
      std::queue<vertex_type>& next_queue,
      vertex_type& start_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
    inf.init(false);
    graph.find(start_vrtx)->second.first = true;
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

        /// push adjacent vertices to the next queue
        auto adjlist_vec = graph.find(src)->second.second;
        for (const auto edge : adjlist_vec) {
          const vertex_type dst = edge.first;
          bool& is_visited = graph.find(dst)->second.first;
          if (!is_visited) {
            next_queue.push(dst);
            /// inf.tree[dst] = src;
            is_visited = true;
          }
          ++(inf.count_visited_edges);
        }
      }  /// end of loop for a frontier

      if (next_queue.empty()) break; /// termination condition
      frontier_queue.swap(next_queue);
      ++level;
    } /// end of BFS loop
  }
};


template <typename graph_type, typename vertex_type, int type>
void bfs_sync(graph_type& graph, vertex_type start_vrtx, size_t max_vertex_id, size_t num_edges)
{

  /// ---- Initialization ----- ////
  auto tic_init = graphstore::utility::duration_time();
  graph_trv_info::trv_inf<vertex_type> inf(max_vertex_id, num_edges);
  std::queue<vertex_type> frontier_queue;
  std::queue<vertex_type> next_queue;
  frontier_queue.push(start_vrtx);
  std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

  /// --- BFS main loop -------- ///
  auto tic = graphstore::utility::duration_time();
  bfs_core<graph_type, vertex_type, type>::run_bfs_sync(graph, inf, frontier_queue, next_queue, start_vrtx);
  double duration_sec = graphstore::utility::duration_time_sec(tic);

  std::cout << "-----------" << std::endl;
  std::cout << "#visited vertices:\t" << inf.count_visited_vertices << std::endl;
  std::cout << "#visited edges:\t"    << inf.count_visited_edges    << std::endl;
  std::cout << "BFS time (sec.):\t"  << duration_sec << std::endl;
  std::cout << "MegaTEPS:\t" << num_edges / duration_sec / (1ULL << 20) << std::endl;
}

void print_statics(std::vector<double>& mteps)
{
  const double ave_mteps = std::accumulate(mteps.begin(), mteps.end(), 0.0) / mteps.size();
  std::sort(mteps.begin(), mteps.end());
  const double med_mteps = (mteps[mteps.size()/2] + mteps[(mteps.size()- 1) / 2]) / 2.0;

  std::cout << "average MTEPS =\t" << ave_mteps << std::endl;
  std::cout << "median MTEPS =\t" << med_mteps << std::endl;
}

#endif // BFS_HPP

