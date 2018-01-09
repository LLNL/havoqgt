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
#include <boost/container/deque.hpp>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/k_visit_bitmap.hpp>

#ifdef NUM_SOURCES
constexpr int k_num_sources = NUM_SOURCES;
#else
constexpr int k_num_sources = 8;
#endif

using level_t = uint16_t;
using visit_bitmap_t = typename havoqgt::k_visit_bitmap<k_num_sources>;

struct kbfs_result_t {
  level_t level;
  size_t num_visited_vertices;
};
// Should add to graph object?
level_t g_current_level;

namespace havoqgt {

template<typename Graph, typename VisitBitmap>
class kbfs_visitor {
 private:
  static constexpr int index_level = 0;
  static constexpr int index_bitmap = 1;
  static constexpr int index_next_bitmap = 2;
  static constexpr int index_visit_flag = 3;
  static constexpr int index_next_visit_flag = 4;

 public:
  typedef typename Graph::vertex_locator                 vertex_locator;

  kbfs_visitor()
    : vertex(),
      // m_level(0),
      m_visit_bitmap() { }

  explicit kbfs_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      // m_level(0),
      m_visit_bitmap()
  { }

  kbfs_visitor(vertex_locator _vertex, uint16_t _source_no)
    : vertex(_vertex),
      // m_level(0),
      m_visit_bitmap()
  {
    m_visit_bitmap.set(_source_no);
  }

#pragma GCC diagnostic pop
  kbfs_visitor(vertex_locator _vertex, visit_bitmap_t _visit_bitmap)
  // kbfs_visitor(vertex_locator _vertex, uint16_t _level, visit_bitmap_t _visit_bitmap)
    : vertex(_vertex),
      // m_level(_level),
      m_visit_bitmap(_visit_bitmap) { }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //

    if (!std::get<index_visit_flag>(alg_data)[vertex]) return false;

    // const auto g_current_level = std::get<index_level>(alg_data)[vertex];
    const auto& current_bitmap = std::get<index_bitmap>(alg_data)[vertex];
    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      kbfs_visitor new_visitor(neighbor, current_bitmap);
      // kbfs_visitor new_visitor(neighbor, g_current_level + 1, current_bitmap);
      vis_queue->queue_visitor(new_visitor);
    }

    return true;
  }

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    // -------------------------------------------------- //
    // This function apply sent data for this vertex (apply step)
    // -------------------------------------------------- //

    visit_bitmap_t& next_visit_bitmap = std::get<index_next_bitmap>(alg_data)[vertex];

    for (size_t i = 0; i < visit_bitmap_t::num_bit; ++i) {
      bool next_visited_bit = next_visit_bitmap.get(i);
      if (!next_visited_bit && m_visit_bitmap.get(i)) {
        next_visit_bitmap.set(i);
        std::get<index_next_visit_flag>(alg_data)[vertex] = true;
        std::get<index_level>(alg_data)[vertex][i] = g_current_level + 1;
      }
    }
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    return false;
  }

  // uint16_t level() const {  return m_level; }
  typename VisitBitmap::value_type visit_bitmap() const {  return m_visit_bitmap; }

  friend inline bool operator>(const kbfs_visitor& v1, const kbfs_visitor& v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const kbfs_visitor &v1, const kbfs_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  // uint16_t       m_level;
  visit_bitmap_t  m_visit_bitmap;
} __attribute__ ((packed));


template <typename TGraph, typename LevelData, typename VisitBitmapData, typename VisitFlagData>
kbfs_result_t k_breadth_first_search_level_per_source(TGraph *g,
                                                 LevelData &level_data,
                                                 VisitBitmapData &visit_bitmap,
                                                 VisitBitmapData &next_visit_bitmap,
                                                 VisitFlagData &visit_flag,
                                                 VisitFlagData &next_visit_flag,
                                                 std::vector<typename TGraph::vertex_locator> source_list)
{
  typedef  kbfs_visitor<TGraph, VisitBitmapData>    visitor_type;
  auto alg_data = std::forward_as_tuple(level_data, visit_bitmap, next_visit_bitmap, visit_flag, next_visit_flag);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  g_current_level = 0;
  size_t num_total_visited_vertices(0);

  // Set (pre_visit) BFS sources
  int mpi_rank = 0;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  for (int i = 0; i < source_list.size(); ++i) {
    if (source_list[i].owner() == mpi_rank) {
      vq.queue_visitor(visitor_type(source_list[i], i)); // -> call pre_visit
      visit_bitmap[source_list[i]] = next_visit_bitmap[source_list[i]];
      visit_flag[source_list[i]] = true;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  ++g_current_level;
  while (g_current_level < std::numeric_limits<level_t>::max()) { // level
    // ------------------------------ Traversal ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      vq.init_visitor_traversal_new(); // init_visit -> queue
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) {
        std::cout << "Traversal time: " << time_end - time_start << std::endl;
      }
    }

    // ------------------------------ Local update ------------------------------ //
    size_t num_visited_vertices = 0;
    {
      const double time_start = MPI_Wtime();
      for (auto vitr = g->vertices_begin(), end = g->vertices_end(); vitr != end; ++vitr) {
        visit_bitmap[*vitr] = next_visit_bitmap[*vitr];
        num_visited_vertices += next_visit_flag[*vitr];
        visit_flag[*vitr] = next_visit_flag[*vitr];
        next_visit_flag[*vitr] = false;
      }

      for (auto citr = g->controller_begin(), end = g->controller_end(); citr != end; ++citr) {
        visit_bitmap[*citr] = next_visit_bitmap[*citr];
        num_visited_vertices += next_visit_flag[*citr];
        visit_flag[*citr] = next_visit_flag[*citr];
        next_visit_flag[*citr] = false;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) {
        std::cout << "Local Update time: " << time_end - -time_start << std::endl;
      }
    }

    num_visited_vertices = mpi_all_reduce(num_visited_vertices, std::plus<size_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // Just in case

    if (num_visited_vertices == 0) break;
    num_total_visited_vertices += num_visited_vertices;

    if (mpi_rank == 0) {
      std::cout << "Level: " << g_current_level << std::endl;
      std::cout << "# visited vertices (not unique): " << num_visited_vertices << std::endl;
    }
    ++g_current_level;
  }

  if (mpi_rank == 0) {
    std::cout << "# total visited vertices (not unique): " << num_total_visited_vertices << std::endl;
  }

  return kbfs_result_t{g_current_level, num_total_visited_vertices};
}
} //end namespace havoqgt

#endif //HAVOQGT_K_BREADTH_FIRST_SEARCH_SYNC_LEVEL_PER_SOURCE_HPP
