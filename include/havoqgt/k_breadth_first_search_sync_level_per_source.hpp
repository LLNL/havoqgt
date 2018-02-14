/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */
//
// Created by Iwabuchi, Keita on 12/14/17.
//

#ifndef HAVOQGT_K_BREADTH_FIRST_SEARCH_SYNC_LEVEL_PER_SOURCE_HPP
#define HAVOQGT_K_BREADTH_FIRST_SEARCH_SYNC_LEVEL_PER_SOURCE_HPP

#include <iostream>
#include <iomanip>
#include <array>
#include <limits>
#include <unordered_map>

#include <boost/container/deque.hpp>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/detail/bitmap.hpp>


namespace havoqgt
{
using level_t = uint16_t;

template <typename _graph_t, int _k_num_sources>
struct kbfs_type
{
  using graph_t = _graph_t;
  static constexpr int k_num_sources = _k_num_sources;
  using k_level = std::array<level_t, k_num_sources>;
  using k_bitmap = typename havoqgt::detail::static_bitmap<k_num_sources>;
  static constexpr level_t unvisited_level = std::numeric_limits<level_t>::max();


  struct vertex_data
  {
    using level_array_t = typename graph_t::template vertex_data<k_level, std::allocator<k_level>>;
    using visited_sources_bitmap_t = typename graph_t::template vertex_data<k_bitmap, std::allocator<k_bitmap>>;
    using visited_flag_t = typename graph_t::template vertex_data<bool, std::allocator<bool>>;

    explicit vertex_data(const graph_t &graph)
      : level(graph),
        visited_sources_bitmap(graph),
        new_visited_sources_bitmap(graph),
        visited_flag(graph),
        new_visited_flag(graph) { }

    void init()
    {
      k_level initial_value;
      initial_value.fill(unvisited_level);
      level.reset(initial_value);

      visited_sources_bitmap.reset(k_bitmap());
      new_visited_sources_bitmap.reset(k_bitmap());

      visited_flag.reset(false);
      new_visited_flag.reset(false);
    }

    level_array_t level;
    visited_sources_bitmap_t visited_sources_bitmap;
    visited_sources_bitmap_t new_visited_sources_bitmap;
    visited_flag_t visited_flag; // whether visited at the previous level
    visited_flag_t new_visited_flag;
  };
};

// Should add to graph object?
level_t g_current_level;

template <typename Graph, int num_k_sources>
class kbfs_visitor
{
 private:
  using kbfs_t = kbfs_type<Graph, num_k_sources>;
  static constexpr int index_level = 0;
  static constexpr int index_bitmap = 1;
  static constexpr int index_next_bitmap = 2;
  static constexpr int index_visit_flag = 3;
  static constexpr int index_next_visit_flag = 4;

 public:
  typedef typename Graph::vertex_locator vertex_locator;

  kbfs_visitor()
    : vertex(),
      visit_bitmap() { }

  explicit kbfs_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      visit_bitmap() { }

#pragma GCC diagnostic pop

  kbfs_visitor(vertex_locator _vertex, typename kbfs_t::k_bitmap _visit_bitmap)
    : vertex(_vertex),
      visit_bitmap(_visit_bitmap) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //

    if (!std::get<index_visit_flag>(alg_data)[vertex]) return false; // TODO: erase?

    const auto &bitmap = std::get<index_bitmap>(alg_data)[vertex];
    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(kbfs_visitor(eitr.target(), bitmap));
    }

    return true; // trigger bcast from the masters (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const
  {
    // -------------------------------------------------- //
    // This function applies sent data to the vertex (apply step)
    // -------------------------------------------------- //

    bool updated(false);
    for (size_t k = 0; k < kbfs_t::k_num_sources; ++k) { // TODO: change to actual num bit
      if (std::get<index_level>(alg_data)[vertex][k] != kbfs_t::unvisited_level)
        continue; // Already visited by source k

      if (!visit_bitmap.get(k)) continue; // No visit request from source k

      std::get<index_level>(alg_data)[vertex][k] = g_current_level;
      std::get<index_next_bitmap>(alg_data)[vertex].set(k);
      std::get<index_next_visit_flag>(alg_data)[vertex] = true; // TODO: doesn't need?
      updated = true;
    }

    // return true to call pre_visit() at the master of delegated vertices required by visitor_queue.hpp 274
    // (label:FLOW2)
    return updated && vertex.is_delegate();
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    // return true, to bcast bitmap, to visit non-master delegates, sent by label:FLOW2
    if (!vertex.get_bcast()) {
      // const int mpi_rank = havoqgt_env()->world_comm().rank();
      // std::cout << "trigger bcast " << mpi_rank << " : " << g.master(vertex) << " : " << vertex.owner() << " : " << vertex.local_id() << " : " << visit_bitmap.get(0) << std::endl;
      return true;
    }

    // handle bcast triggered by the line above
    // visitors issued by label:FLOW1 paths this condition
    if (pre_visit(alg_data)) {
      // const int mpi_rank = havoqgt_env()->world_comm().rank();
      // std::cout << "recived bcast " << mpi_rank << " : " << g.master(vertex) << " : " << vertex.owner() << " : " << vertex.local_id() << " : " << visit_bitmap.get(0) << std::endl;
      return false;
    }
    // FIXME:
    // if I'm not the owner of a vertex, I recieve the vertex with bitmap == 1 ??

    // for the case, triggered by label:FLOW1
    // scatter current status to neighbors (this line is supposed to be executed on only non-master delegates)
    // const int mpi_rank = havoqgt_env()->world_comm().rank();
    // std::cout << "queue visitor " << mpi_rank << " : " << g.master(vertex) << " : " << vertex.owner() << " : " << vertex.local_id() << " : " << visit_bitmap.get(0) << std::endl;
    const auto &bitmap = std::get<index_bitmap>(alg_data)[vertex];
    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(kbfs_visitor(eitr.target(), bitmap));
    }

    return false;
  }

  friend inline bool operator>(const kbfs_visitor &v1, const kbfs_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const kbfs_visitor &v1, const kbfs_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  typename kbfs_t::k_bitmap visit_bitmap;
} __attribute__ ((packed));

namespace detail
{
template <typename TGraph, int k_num_sources>
void set_bfs_sources(const std::vector<typename TGraph::vertex_locator> source_list,
                     typename kbfs_type<TGraph, k_num_sources>::vertex_data &vertex_data)
{
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  const double time_start = MPI_Wtime();
  for (size_t k = 0; k < source_list.size(); ++k) {
    if (source_list[k].owner() == mpi_rank || source_list[k].is_delegate()) {
      vertex_data.level[source_list[k]][k] = 0; // source vertices are visited at level 0
      typename kbfs_type<TGraph, k_num_sources>::k_bitmap visited_bitmap;
      visited_bitmap.set(k);
      vertex_data.visited_sources_bitmap[source_list[k]] = visited_bitmap;
      vertex_data.visited_flag[source_list[k]] = true;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  const double time_end = MPI_Wtime();
}

template <typename TGraph, int k_num_sources, typename iterator_t>
void update_local_vertex_data(iterator_t vitr, iterator_t end,
                              typename kbfs_type<TGraph, k_num_sources>::vertex_data &vertex_data)
{
  for (; vitr != end; ++vitr) {
    vertex_data.visited_sources_bitmap[*vitr] = vertex_data.new_visited_sources_bitmap[*vitr];
    vertex_data.visited_flag[*vitr] = vertex_data.new_visited_flag[*vitr];
    vertex_data.new_visited_flag[*vitr] = false;
  }
}

template <typename TGraph, int k_num_sources, typename iterator_t>
std::vector<size_t>
count_visited_vertices(iterator_t vitr, iterator_t end,
                       const typename kbfs_type<TGraph, k_num_sources>::vertex_data &vertex_data,
                       const level_t level)
{
  std::vector<size_t> num_vertices(k_num_sources, 0);
  for (; vitr != end; ++vitr) {
    for (size_t k = 0; k < k_num_sources; ++k) {
      if (vertex_data.level[*vitr][k] == level) ++num_vertices[k];
    }
  }
  return num_vertices;
}
} // namespace detail


template <typename TGraph, int k_num_sources>
void k_breadth_first_search_level_per_source(TGraph *g,
                                             typename kbfs_type<TGraph, k_num_sources>::vertex_data &vertex_data,
                                             std::vector<typename TGraph::vertex_locator> source_list)
{
  typedef kbfs_visitor<TGraph, k_num_sources> visitor_type;

  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  auto alg_data = std::forward_as_tuple(vertex_data.level,
                                        vertex_data.visited_sources_bitmap, vertex_data.new_visited_sources_bitmap,
                                        vertex_data.visited_flag, vertex_data.new_visited_flag);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // std::cout << std::setprecision(2);

  // Set BFS sources (level 0)
  detail::set_bfs_sources<TGraph, k_num_sources>(source_list, vertex_data);
  std::vector<size_t> num_total_visited_vertices(source_list.size(), 1);
  g_current_level = 1;

  if (mpi_rank == 0) {
    std::cout << "[Level]\t" << "[Traversal time]\t" <<"[Local Update]\t" << "[Visited vertices]" << std::endl;
  }

  while (g_current_level < std::numeric_limits<level_t>::max()) { // level
    if (mpi_rank == 0) std::cout << g_current_level << "\t";
    // ------------------------------ Traversal ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      vq.init_visitor_traversal_new(); // init_visit -> queue
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) std::cout << time_end - time_start << "\t";
    }

    // ------------------------------ Update BFS information ------------------------------ //
    std::vector<size_t> num_local_visited_vertices(source_list.size(), 0);
    {
      const double time_start = MPI_Wtime();
      detail::update_local_vertex_data<TGraph, k_num_sources>(g->vertices_begin(), g->vertices_end(), vertex_data);
      detail::update_local_vertex_data<TGraph, k_num_sources>(g->delegate_vertices_begin(), g->delegate_vertices_end(), vertex_data);
      std::vector<size_t> cnt1 = detail::count_visited_vertices<TGraph, k_num_sources>(g->vertices_begin(), g->vertices_end(),
                                                                                       vertex_data, g_current_level);
      std::vector<size_t> cnt2 = detail::count_visited_vertices<TGraph, k_num_sources>(g->controller_begin(), g->controller_end(),
                                                                                       vertex_data, g_current_level);
      for (size_t i = 0; i < num_local_visited_vertices.size(); ++i)
        num_local_visited_vertices[i] = cnt1[i] + cnt2[i];
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) std::cout << time_end - -time_start << "\t";
    }

    std::vector<size_t> num_global_visited_vertices;
    mpi_all_reduce(num_local_visited_vertices, num_global_visited_vertices, std::plus<size_t>(), MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      for (auto n : num_global_visited_vertices)
        std::cout << n << "\t";
      std::cout << std::endl;
    }

    // Terminate condition
    if (std::accumulate(num_global_visited_vertices.begin(), num_global_visited_vertices.end(), 0ULL) == 0) break;

    for (size_t i = 0; i < num_total_visited_vertices.size(); ++i)
      num_total_visited_vertices[i] += num_global_visited_vertices[i];

    ++g_current_level;
  }

  if (mpi_rank == 0) {
    std::cout << "# total visited vertices from each sources: " << num_total_visited_vertices[0] << std::endl;

    /// Simple validation code
    for (size_t i = 0; i < num_total_visited_vertices.size() - 1; ++i) {
      if (num_total_visited_vertices[i] != num_total_visited_vertices[i+1]) {
        std::cerr << "# total visited vertices do not much" << std::endl;
        std::abort();
      }
    }
  }
  // std::cout << std::setprecision(6);
}
} //end namespace havoqgt

#endif //HAVOQGT_K_BREADTH_FIRST_SEARCH_SYNC_LEVEL_PER_SOURCE_HPP
