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

template <typename graph_type, typename vertex_type, bool is_rhhda>
class bfs_core {

 public:

  bfs_core()=delete;

  static void run_bfs (graph_type& graph, graph_trv_info::trv_inf<vertex_type>& inf, std::queue<vertex_type>& frontier_queue, std::queue<vertex_type>& next_queue)
  {
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
class bfs_core <graph_type, vertex_type, true> {

 public:

  bfs_core()=delete;

  static void run_bfs (graph_type& graph, graph_trv_info::trv_inf<vertex_type>& inf, std::queue<vertex_type>& frontier_queue, std::queue<vertex_type>& next_queue)
  {
    /// --- BFS main loop -------- ///
    size_t level = 0;
    while (true) {
      std::cout << "Lv. " << level << ", size of frontier queue =\t" << frontier_queue.size() << std::endl;

      /// loop for current frontier
      while (!frontier_queue.empty()) {
        vertex_type src = frontier_queue.front();
  //      std::cout << "src " << src << std::endl;
        frontier_queue.pop();
        ++inf.count_visited_vertices;

        /// push adjacent vertices to the next queue
        for (auto edge = graph.find_low_edge(src); !edge.is_end(); ++edge) {
          const vertex_type dst = edge->second;
          if (!inf.is_visited[dst]) {
            next_queue.push(dst);
            inf.tree[dst] = src;
            inf.is_visited[dst] = true;
          }
          ++(inf.count_visited_edges);
        }

        for (auto edge = graph.find_mid_high_edge(src); !edge.is_end(); ++edge) {
          const vertex_type dst = edge->key;
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


template <typename graph_type, typename vertex_type, bool is_rhhda>
void bfs_sync(graph_type& graph, vertex_type start_vrtx, size_t max_vertex_id, size_t num_edges)
{

  /// ---- Initialization ----- ////
  graph_trv_info::trv_inf<vertex_type> inf(max_vertex_id, num_edges);
  std::queue<vertex_type> frontier_queue;
  std::queue<vertex_type> next_queue;

  /// ---- Set source vertex ---- ///
  frontier_queue.push(start_vrtx);

  inf.is_visited[start_vrtx] = true;
  inf.tree[start_vrtx] = start_vrtx;

  /// --- BFS main loop -------- ///
  auto tic = graphstore::utility::duration_time();
  bfs_core<graph_type, vertex_type, is_rhhda>::run_bfs(graph, inf, frontier_queue, next_queue);
  auto duration_misec = graphstore::utility::duration_time_sec(tic);

  std::cout << "-----------" << std::endl;
  std::cout << "#visited vertices:\t" << inf.count_visited_vertices << std::endl;
  std::cout << "#visited edges:\t"    << inf.count_visited_edges    << std::endl;
  std::cout << "BFS time (misec):\t"  << duration_misec << std::endl;
  std::cout << "MegaTEPS:\t" << num_edges / (duration_misec / 1000000.0) / (1ULL << 20) << std::endl;
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

