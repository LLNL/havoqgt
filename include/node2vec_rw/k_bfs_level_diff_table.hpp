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

#ifndef HAVOQGT_K_BFS_LEVEL_DIFF_TABLE_HPP
#define HAVOQGT_K_BFS_LEVEL_DIFF_TABLE_HPP

#include <vector>
#include <memory>
#include <scoped_allocator>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <node2vec_rw/hash.hpp>
#include <node2vec_rw/k_bfs.hpp>

namespace node2vec_rw {

template <typename graph_allocator_type>
class k_bfs_level_diff_table {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;

  struct vertex_locator_hash {
    std::size_t operator()(const vertex_type &v) const noexcept {
      return hash<decltype(v.local_id())>{}(v.local_id());
    }
  };
  using per_vertex_level_diff_type = std::unordered_map<vertex_type,
                                                        std::vector<int8_t>,
                                                        vertex_locator_hash>;
  using internal_table_type = typename graph_type::template vertex_data<per_vertex_level_diff_type,
                                                                        std::allocator<per_vertex_level_diff_type>>;

 public:

  k_bfs_level_diff_table(const graph_type &graph, const std::size_t k_size)
      : m_graph(graph),
        m_table(graph) {}

  int8_t get(const vertex_type &source, const vertex_type &target, const std::size_t k) const {
    assert(m_table[source].count(target));
    return m_table[source].at(target).at(k);
  }

  void set(const vertex_type &source, const vertex_type &target, const std::size_t k, const int8_t value) {
    m_table[source][target].resize(k + 1);
    assert(value == -1 || value == 0 || value == 1);
    m_table[source][target][k] = value;
  }

  std::size_t k_size(const vertex_type &source, const vertex_type &target) const {
    assert(m_table[source].count(target));
    return m_table[source].at(target).size();
  }

  bool has_entry(const vertex_type &vertex) const {
    return m_table[vertex].size() > 0;
  }

 private:
  const graph_type &m_graph;
  internal_table_type m_table;
};

/// \NOTE this algorithm assumes an undirected graph
template <typename graph_allocator_type, typename level_type>
class level_diff_compute_visitor {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;
  using k_bfs_level_table_type = k_bfs_level_table<graph_type, level_type>;
  using k_bfs_level_diff_table_type = k_bfs_level_diff_table<graph_allocator_type>;

  enum algo_data {
    level_table = 0,
    level_diff_table = 1
  };

 public:

  level_diff_compute_visitor() = default;

  explicit level_diff_compute_visitor(vertex_type _vertex)
      : vertex(_vertex),
        source_vertex(),
        level(),
        k() {}

  level_diff_compute_visitor(vertex_type _vertex, vertex_type _source_vertex, level_type _level, uint16_t _k)
      : vertex(_vertex), source_vertex(_source_vertex), level(_level), k(_k) {}

  template <typename VisitorQueueHandle, typename algorithm_data_type>
  bool init_visit(graph_type &g, VisitorQueueHandle vis_queue, algorithm_data_type &algorithm_data) const {
    return visit(g, vis_queue, algorithm_data);
  }

  template <typename algorithm_data_type>
  bool pre_visit(algorithm_data_type &algorithm_data) const {
    const auto &level_table = std::get<algo_data::level_table>(algorithm_data);
    const int8_t diff = level - level_table.get_level(vertex, k);
    auto &level_diff_table = std::get<algo_data::level_diff_table>(algorithm_data);
    level_diff_table.set(vertex, source_vertex, k, diff);

    return false;
  }

  template <typename VisitorQueueHandle, typename algorithm_data_type>
  bool visit(graph_type &g, VisitorQueueHandle vis_queue, algorithm_data_type &algorithm_data) const {
    const auto &level_table = std::get<algo_data::level_table>(algorithm_data);

    if (!level_table.visited(vertex, 0)) return false;

    for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      for (uint16_t k = 0; k < level_table.k_size(); ++k) {
        assert(level_table.visited(vertex, k));
        auto neighbor = eitr.target();
        level_diff_compute_visitor new_visitor(neighbor, vertex, level_table.get_level(vertex, k), k);
        vis_queue->queue_visitor(new_visitor);
      }
    }
    
    return true;
  }

  friend inline bool operator>(const level_diff_compute_visitor &v1, const level_diff_compute_visitor &v2) {
    return false;
  }

  friend inline bool operator<(const level_diff_compute_visitor &v1, const level_diff_compute_visitor &v2) {
    return false;
  }

  vertex_type vertex;
  vertex_type source_vertex;
  level_type level;
  uint16_t k;
};

template <typename graph_allocator_type, typename level_type>
auto construct_k_bfs_level_diff_table(havoqgt::delegate_partitioned_graph<graph_allocator_type> &graph,
                                      const k_bfs_level_table<graph_allocator_type, level_type> &level_table) {
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  auto level_diff_table = std::make_unique<k_bfs_level_diff_table<graph_allocator_type>>
      (const_cast<graph_type &>(graph), level_table.k_size());

  auto alg_data = std::forward_as_tuple(level_table, *level_diff_table);

  using vistor_type = level_diff_compute_visitor<graph_allocator_type, level_type>;
  auto vq = havoqgt::create_visitor_queue<vistor_type, havoqgt::detail::visitor_priority_queue>(&graph, alg_data);
  vq.init_visitor_traversal();

  return level_diff_table;
}
}

#endif //HAVOQGT_K_BFS_LEVEL_DIFF_TABLE_HPP
