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

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/detail/bitmap.hpp>


namespace havoqgt
{

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class k_breadth_first_search_vertex_data
{
 public:
  using k_level_t = std::array<level_t, k_num_sources>;
  using k_bitmap_t = typename havoqgt::detail::static_bitmap<k_num_sources>;

 private:
  using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;
  using vertex_data_visited_level_t = typename graph_t::template vertex_data<k_level_t, std::allocator<k_level_t>>;
  using vertex_data_visited_sources_bitmap_t = typename graph_t::template vertex_data<k_bitmap_t, std::allocator<k_bitmap_t>>;
  using queue_status_t = typename graph_t::template vertex_data<uint8_t, std::allocator<uint8_t>>;

  static constexpr level_t unvisited_level = std::numeric_limits<level_t>::max();
  static constexpr uint8_t k_not_in_queue_bit = 0x0;
  static constexpr uint8_t k_in_frontier_bit = 0x1;
  static constexpr uint8_t k_in_next_queue_bit = 0x2;

 public:
  explicit k_breadth_first_search_vertex_data(const graph_t &graph)
    : m_graph(graph),
      m_level(graph),
      m_visited_bitmap(graph),
      m_tmp_visited_bitmap(graph),
      m_queue_status(graph)
  {
#ifdef DEBUG
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
#endif
  }

  void reset()
  {
    k_level_t initial_value;
    initial_value.fill(unvisited_level);
    m_level.reset(initial_value);

    m_visited_bitmap.reset(k_bitmap_t());
    m_tmp_visited_bitmap.reset(k_bitmap_t());

    m_queue_status.reset(k_not_in_queue_bit);
  }

  bool in_frontier(const typename graph_t::vertex_locator vertex)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return (m_queue_status[vertex] & k_in_frontier_bit);
  }

  k_bitmap_t& visited_bitmap(const typename graph_t::vertex_locator vertex)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return m_visited_bitmap[vertex];
  }

  bool visited_by(const typename graph_t::vertex_locator vertex, const size_t k)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return (m_level[vertex][k] != unvisited_level);
  }

  void visit(const typename graph_t::vertex_locator vertex, const size_t k, const level_t level)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    m_level[vertex][k] = level;
    m_tmp_visited_bitmap[vertex].set(k);
    m_queue_status[vertex] |= k_in_next_queue_bit;
  }

  void set_source(const typename graph_t::vertex_locator vertex, const size_t k)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    m_level[vertex][k] = 0;
    m_visited_bitmap[vertex].set(k);
    m_tmp_visited_bitmap[vertex].set(k);
    m_queue_status[vertex] = k_in_frontier_bit;
  }

  void prepare_for_next_iteration()
  {
    prepare_for_next_iteration_helper(m_graph.vertices_begin(), m_graph.vertices_end());
    prepare_for_next_iteration_helper(m_graph.delegate_vertices_begin(), m_graph.delegate_vertices_end());
  }

  std::vector<size_t> count_visited_vertices(const level_t level, const size_t num_souces)
  {
    std::vector<size_t> cnt_vrtx = count_visited_vertices_helper(level, num_souces,
                                                                 m_graph.vertices_begin(), m_graph.vertices_end());
    std::vector<size_t> cnt_ctrl = count_visited_vertices_helper(level, num_souces,
                                                                 m_graph.controller_begin(), m_graph.controller_end());
    std::transform(cnt_vrtx.begin(), cnt_vrtx.end(), cnt_ctrl.begin(), cnt_vrtx.begin(), std::plus<size_t>());

    return cnt_vrtx;
  }

  k_level_t& level(const typename graph_t::vertex_locator vertex)
  {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return m_level[vertex];
  }

 private:
  template <typename iterator_t>
  std::vector<size_t> count_visited_vertices_helper(const level_t level, const size_t num_souces,
                                                    iterator_t vitr, iterator_t end)
  {
    std::vector<size_t> num_vertices(num_souces, 0);
    for (; vitr != end; ++vitr) {
      for (size_t k = 0; k < num_souces; ++k) {
        if (m_level[*vitr][k] == level) ++num_vertices[k];
      }
    }
    return num_vertices;
  }

  template <typename iterator_t>
  void prepare_for_next_iteration_helper(iterator_t vitr, iterator_t end)
  {
    for (; vitr != end; ++vitr) {
      m_visited_bitmap[*vitr] = m_tmp_visited_bitmap[*vitr];
      m_queue_status[*vitr] = (m_queue_status[*vitr] & k_in_next_queue_bit) ? k_in_frontier_bit : k_not_in_queue_bit;
    }
  }

  const graph_t &m_graph;
  vertex_data_visited_level_t m_level;
  vertex_data_visited_sources_bitmap_t m_visited_bitmap;
  vertex_data_visited_sources_bitmap_t m_tmp_visited_bitmap;
  queue_status_t m_queue_status;
#ifdef DEBUG
  int m_mpi_rank;
#endif
};

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
constexpr level_t  k_breadth_first_search_vertex_data<segment_manager_t, level_t, k_num_sources>::unvisited_level;

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
constexpr uint8_t k_breadth_first_search_vertex_data<segment_manager_t, level_t, k_num_sources>::k_not_in_queue_bit;

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
constexpr uint8_t k_breadth_first_search_vertex_data<segment_manager_t, level_t, k_num_sources>::k_in_frontier_bit;

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
constexpr uint8_t k_breadth_first_search_vertex_data<segment_manager_t, level_t, k_num_sources>::k_in_next_queue_bit;


template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class k_breadth_first_search
{
 private:
  using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;
  using self_type = k_breadth_first_search<segment_manager_t, level_t, k_num_sources>;

 public:
  using vertex_data_t = k_breadth_first_search_vertex_data<segment_manager_t, level_t, k_num_sources>;
  class kbfs_visitor;

  k_breadth_first_search<segment_manager_t, level_t, k_num_sources>(graph_t& graph)
  : m_graph(graph),
    m_vertex_data(graph)
  { }

  size_t run(const std::vector<typename graph_t::vertex_locator>& source_list)
  {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    std::cout << std::setprecision(2);

    m_vertex_data.reset();
    set_source(source_list);

    level_t current_level(1);
    auto alg_data = std::forward_as_tuple(m_vertex_data, current_level);
    auto vq = create_visitor_queue<kbfs_visitor, havoqgt::detail::visitor_priority_queue>(&m_graph, alg_data);

    if (mpi_rank == 0) {
      std::cout << "[Level]\t" << "[Traversal time]\t" << "[Local Update]\t" << "[Visited vertices]" << std::endl;
    }

    std::vector<size_t> num_total_visited_vertices(source_list.size(), 1); // 1 is for source vertices
    while (current_level < std::numeric_limits<level_t>::max()) { // level
      if (mpi_rank == 0) std::cout << current_level << "\t";
      // ------------------------------ Traversal ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        vq.init_visitor_traversal_new();
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << time_end - time_start << "\t";
      }

      // ------------------------------ Update BFS information ------------------------------ //
      std::vector<size_t> num_visited_vertices;
      {
        const double time_start = MPI_Wtime();
        m_vertex_data.prepare_for_next_iteration();
        num_visited_vertices = m_vertex_data.count_visited_vertices(current_level, source_list.size());
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << time_end - -time_start << "\t";
      }

      // ------------------------------ Terminate condition ------------------------------ //
      {
        mpi_all_reduce_inplace(num_visited_vertices, std::plus<size_t>(), MPI_COMM_WORLD);
        if (mpi_rank == 0) {
          for (auto n : num_visited_vertices)
            std::cout << n << "\t";
          std::cout << std::endl;
        }

        // Terminate condition
        if (std::accumulate(num_visited_vertices.begin(), num_visited_vertices.end(), 0ULL) == 0) break;
      }

      std::transform(num_total_visited_vertices.begin(), num_total_visited_vertices.end(), num_visited_vertices.begin(),
                     num_total_visited_vertices.begin(), std::plus<size_t>());
      ++current_level;
    }

    if (mpi_rank == 0) {
      std::cout << "# total visited vertices from each sources: " << num_total_visited_vertices[0] << std::endl;

      /// Simple validation code
      std::sort(num_total_visited_vertices.begin(), num_total_visited_vertices.end());
      if (*(num_total_visited_vertices.begin()) != *(num_total_visited_vertices.end()-1)) {
        std::cerr << "# total visited vertices do not much: "
                  << *(num_total_visited_vertices.begin()) << " " << *(num_total_visited_vertices.end()-1) << std::endl;
        std::abort();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << std::fixed;
    return num_total_visited_vertices[0];
  }

  vertex_data_t& vertex_data()
  {
    return m_vertex_data;
  }

 private:
  void set_source(const std::vector<typename graph_t::vertex_locator>& source_list)
  {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    for (size_t k = 0; k < source_list.size(); ++k) {
      if (source_list[k].owner() == static_cast<uint32_t>(mpi_rank) || source_list[k].is_delegate()) {
        m_vertex_data.set_source(source_list[k], k);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

 private:
  graph_t& m_graph;
  vertex_data_t m_vertex_data;
};

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class k_breadth_first_search<segment_manager_t, level_t, k_num_sources>::kbfs_visitor
{
 private:
  enum index
  {
    vertex_data = 0,
    current_level = 1
  };

 public:
  typedef typename graph_t::vertex_locator vertex_locator;

  kbfs_visitor()
    : vertex(),
      visit_bitmap() { }

  explicit kbfs_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      visit_bitmap() { }

#pragma GCC diagnostic pop

  kbfs_visitor(vertex_locator _vertex, typename vertex_data_t::k_bitmap_t _visit_bitmap)
    : vertex(_vertex),
      visit_bitmap(_visit_bitmap) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    auto &vertex_data = std::get<index::vertex_data>(alg_data);

    if (!vertex_data.in_frontier(vertex)) return false;

    const auto &bitmap = vertex_data.visited_bitmap(vertex);
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
    auto &vertex_data = std::get<index::vertex_data>(alg_data);

    bool updated(false);
    for (size_t k = 0; k < k_num_sources; ++k) { // TODO: change to actual num bit
      if (vertex_data.visited_by(vertex, k)) continue; // Already visited by source k

      if (!visit_bitmap.get(k)) continue; // No visit request from source k

      vertex_data.visit(vertex, k, std::get<index::current_level>(alg_data));
      updated = true;
    }

    // return true to call pre_visit() at the master of delegated vertices required by visitor_queue.hpp 274
    // (label:FLOW2)
    return updated && vertex.is_delegate();
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
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
    const auto &bitmap = std::get<index::vertex_data>(alg_data).visited_bitmap(vertex);
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
  typename vertex_data_t::k_bitmap_t visit_bitmap;
} __attribute__ ((packed));

} //end namespace havoqgt

#endif //HAVOQGT_K_BREADTH_FIRST_SEARCH_SYNC_LEVEL_PER_SOURCE_HPP
