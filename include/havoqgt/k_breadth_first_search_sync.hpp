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


#ifndef HAVOQGT_MPI_K_BREADTH_FIRST_SEARCH_HPP
#define HAVOQGT_MPI_K_BREADTH_FIRST_SEARCH_HPP

#include <boost/container/deque.hpp>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/detail/bitmap.hpp>

#ifdef NUM_SOURCES
constexpr int k_num_sources = NUM_SOURCES;
#else
constexpr int k_num_sources = 64;
#endif

using visit_bitmap_t = typename havoqgt::detail::static_bitmap<k_num_sources>;

template <size_t num_sources>
bool is_contain(const visit_bitmap_t &lhs, const visit_bitmap_t &rhs)
{
  for (size_t i = 0; i < visit_bitmap_t::size; ++i) {
    if ((lhs.bitmap[i] | rhs.bitmap[i]) > lhs.bitmap[i]) return false;
  }
  return true;
}

namespace havoqgt {

template<typename Graph, typename VisitBitmap>
class bfs_visitor {
 private:
  static constexpr int index_max_level = 0;
  static constexpr int index_next_max_level = 1;
  static constexpr int index_bitmap = 2;
  static constexpr int index_next_bitmap = 3;
  static constexpr int index_visit_flag = 4;
  static constexpr int index_next_visit_flag = 5;

 public:
  typedef typename Graph::vertex_locator                 vertex_locator;

  bfs_visitor()
    : vertex(),
      m_level(0),
      m_visit_bitmap() { }

  explicit bfs_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      m_level(0),
      m_visit_bitmap()
  { }

  bfs_visitor(vertex_locator _vertex, uint16_t _source_no)
    : vertex(_vertex),
      m_level(0),
      m_visit_bitmap()
  {
    m_visit_bitmap.set(_source_no);
  }

#pragma GCC diagnostic pop
  bfs_visitor(vertex_locator _vertex, uint16_t _level, visit_bitmap_t _visit_bitmap)
    : vertex(_vertex),
      m_level(_level),
      m_visit_bitmap(_visit_bitmap) { }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //

    if (!std::get<index_visit_flag>(alg_data)[vertex]) return false;

    const auto current_level = std::get<index_max_level>(alg_data)[vertex];
    const auto& current_bitmap = std::get<index_bitmap>(alg_data)[vertex];
    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      bfs_visitor new_visitor(neighbor, current_level + 1, current_bitmap);
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

    const bool already_visited = is_contain(next_visit_bitmap, m_visit_bitmap);
    if (!already_visited) {
      next_visit_bitmap |= m_visit_bitmap;
      std::get<index_next_visit_flag>(alg_data)[vertex] = true;
      auto& next_level = std::get<index_next_max_level>(alg_data)[vertex];
      if (next_level < level()) next_level = level();
    }
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    return false;
  }

  uint16_t level() const {  return m_level; }
  typename VisitBitmap::value_type visit_bitmap() const {  return m_visit_bitmap; }

  friend inline bool operator>(const bfs_visitor& v1, const bfs_visitor& v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const bfs_visitor &v1, const bfs_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  uint16_t       m_level;
  visit_bitmap_t  m_visit_bitmap;
} __attribute__ ((packed));


template <typename TGraph, typename LevelData, typename VisitBitmapData, typename VisitFlagData>
uint16_t k_breadth_first_search(TGraph* g,
                                LevelData& level_data,
                                LevelData& next_level_data,
                                VisitBitmapData& visit_bitmap,
                                VisitBitmapData& next_visit_bitmap,
                                VisitFlagData& visit_flag,
                                VisitFlagData& next_visit_flag,
                                std::vector<typename TGraph::vertex_locator> source_list) {
  typedef  bfs_visitor<TGraph, VisitBitmapData>    visitor_type;
  auto alg_data = std::forward_as_tuple(level_data, next_level_data,
                                        visit_bitmap, next_visit_bitmap,
                                        visit_flag, next_visit_flag);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // Set (pre_visit) BFS sources
  int mpi_rank = 0;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  for (uint16_t i = 0; i < source_list.size(); ++i) {
    if (source_list[i].owner() == mpi_rank) {
      vq.queue_visitor(visitor_type(source_list[i], i)); // -> call pre_visit
      level_data[source_list[i]] = next_level_data[source_list[i]];
      visit_bitmap[source_list[i]] = next_visit_bitmap[source_list[i]];
      visit_flag[source_list[i]] = true;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  uint16_t level = 0;
  while (level < std::numeric_limits<uint16_t>::max()) { // level
    if (mpi_rank == 0) std::cout << "==================== " << level + 1 << " ====================" << std::endl;

    vq.init_visitor_traversal_new(); // init_visit -> queue

    bool visited_next_vertex = false;
    for (auto vitr = g->vertices_begin(), end = g->vertices_end(); vitr != end; ++vitr) {
      level_data[*vitr] = next_level_data[*vitr];
      visit_bitmap[*vitr] = next_visit_bitmap[*vitr];
      visited_next_vertex |= next_visit_flag[*vitr];
      visit_flag[*vitr] = next_visit_flag[*vitr];
      next_visit_flag[*vitr] = false;
    }

    // Count the controllers!
    for (auto citr = g->controller_begin(), end = g->controller_end(); citr != end; ++citr) {
      level_data[*citr] = next_level_data[*citr];
      visit_bitmap[*citr] = next_visit_bitmap[*citr];
      visited_next_vertex |= next_visit_flag[*citr];
      visit_flag[*citr] = next_visit_flag[*citr];
      next_visit_flag[*citr] = false;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const char wk = static_cast<char>(visited_next_vertex);
    const char global = mpi_all_reduce(wk, std::logical_or<char>(), MPI_COMM_WORLD);
    if (!global) break;
    MPI_Barrier(MPI_COMM_WORLD); // Just in case
    ++level;
  }

  return level;
}



} //end namespace havoqgt


#endif //HAVOQGT_MPI_T_K_BREADTH_FIRST_SEARCH_HPP
