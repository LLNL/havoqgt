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

#ifndef HAVOQGT_NODE2VEC_RW_VISITOR_HPP
#define HAVOQGT_NODE2VEC_RW_VISITOR_HPP

#include <vector>
#include <memory>
#include <scoped_allocator>

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/mpi.hpp>
#include <node2vec_rw/random.hpp>
#include <node2vec_rw/neighbor_selector.hpp>
#include <node2vec_rw/hash.hpp>
#include <node2vec_rw/global_shuffle.hpp>

namespace node2vec_rw {

template <typename graph_allocator_type, typename edge_weight_data_type>
class rw_visitor_data {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;
  using neighbor_selector_type = neighbor_selector_with_rejection_sampling<graph_allocator_type,
                                                                           edge_weight_data_type>;
 public:
  using walk_id_type = uint64_t;
  using walk_length_type = uint16_t;
  using walk_log_type = std::tuple<walk_id_type, walk_length_type, vertex_type>;
  using walk_log_store_type = std::vector<walk_log_type>;

  enum walk_log {
    walk_id = 0,
    step_no = 1,
    visited_vertex = 2
  };

  rw_visitor_data(graph_type &graph,
                  const edge_weight_data_type &edge_weight_data,
                  const bool small_edge_weight_variance,
                  const k_bfs_level_diff_table<graph_allocator_type> &level_diff_table,
                  const walk_length_type walk_length,
                  const double p,
                  const double q,
                  const MPI_Comm mpi_comm)
      : m_graph(graph),
        m_level_diff_table(level_diff_table),
        m_walk_length(walk_length),
        m_walker_id(0),
        m_rnd_generator(std::random_device{}()),
        m_neighbor_selector(graph, edge_weight_data, level_diff_table, p, q, small_edge_weight_variance),
        m_walk_log_store(),
        m_mpi_comm(mpi_comm) {
    clear();
  }

  bool visited(const vertex_type &vertex) const {
    return m_level_diff_table.has_entry(vertex);
  }

  walk_length_type walk_length() const {
    return m_walk_length;
  }

  walk_id_type gen_walker_id() {
    const walk_id_type max = (1ULL << 40ULL) - 1;
    if ((m_walker_id & max) == max) {
      std::cerr << "Too many walkers are generated" << std::endl;
      ::MPI_Abort(m_mpi_comm, -1);
    }
    return m_walker_id++;
  }

  void store(const walk_id_type walker_id, const walk_length_type walk_step, const vertex_type &vertex) {
    m_walk_log_store.push_back(std::make_tuple(walker_id, walk_step, vertex));
  }

  vertex_type select_neighbor(const vertex_type &prev_vertex,
                              const vertex_type &current_vertex,
                              const bool first_step) {
    return m_neighbor_selector.select(prev_vertex, current_vertex, first_step, &m_rnd_generator);
  }

  void clear() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(m_mpi_comm, &mpi_rank));
    m_walker_id = ((walk_id_type)mpi_rank << 40ULL);
    m_walk_log_store.clear();
  }

  std::vector<std::vector<vertex_type>> get_walk_history() {
    return convert_to_walk_history(shuffle_walk_log());
  }

 private:
  walk_log_store_type shuffle_walk_log() const {
    walk_log_store_type shuffled_walk_log;
    int mpi_size = 0;
    int mpi_rank = 0;
    CHK_MPI(MPI_Comm_size(m_mpi_comm, &mpi_size));
    CHK_MPI(MPI_Comm_rank(m_mpi_comm, &mpi_rank));

    auto walk_log_partitioner = [mpi_size, this](const walk_log_type &data) -> int {
      const uint64_t walk_id = std::get<walk_log::walk_id>(data);
      const int rank = (int)(walk_id >> 40ULL);
      if (rank >= mpi_size) {
        std::cerr << "Invalid walk ID" << std::endl;
        MPI_Abort(m_mpi_comm, -1);
      }
      return rank;
    };

    node2vec_rw::mpi_large_all_to_all(m_walk_log_store.begin(), m_walk_log_store.end(),
                                      walk_log_partitioner,
                                      m_mpi_comm,
                                      std::back_inserter(shuffled_walk_log));

    std::sort(shuffled_walk_log.begin(), shuffled_walk_log.end(),
              [](const walk_log_type &lhd, const walk_log_type &rhd) {
                if (std::get<walk_log::walk_id>(lhd) != std::get<walk_log::walk_id>(rhd))
                  return std::get<walk_log::walk_id>(lhd) < std::get<walk_log::walk_id>(rhd);
                return std::get<walk_log::step_no>(lhd) < std::get<walk_log::step_no>(rhd);
              });

    return shuffled_walk_log;
  }

  std::vector<std::vector<vertex_type>> convert_to_walk_history(walk_log_store_type &&walk_log_store) {
    assert(walk_log_store.size() % m_walk_length == 0);
    const std::size_t num_walks = walk_log_store.size() / m_walk_length;
    std::vector<std::vector<vertex_type>> walk_history(num_walks);
    for (std::size_t w = 0; w < num_walks; ++w) {
      walk_history[w].reserve(m_walk_length);
      for (std::size_t l = 0; l < m_walk_length; ++l) {
        walk_history[w].push_back(std::get<walk_log::visited_vertex>(walk_log_store[w * m_walk_length + l]));
      }
    }
    return walk_history;
  }

  const graph_type &m_graph;
  const k_bfs_level_diff_table<graph_allocator_type> &m_level_diff_table;
  const walk_length_type m_walk_length;
  walk_id_type m_walker_id;
  rand_1024 m_rnd_generator;
  neighbor_selector_type m_neighbor_selector;
  walk_log_store_type m_walk_log_store;
  MPI_Comm m_mpi_comm;
};

template <typename graph_allocator_type, typename edge_weight_type>
class node2vec_rw_visitor {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;
  using algorithm_data_type = rw_visitor_data<graph_allocator_type, edge_weight_type>;

 public:
  using walk_id_type = typename algorithm_data_type::walk_id_type;
  using walk_length_type = typename algorithm_data_type::walk_length_type;

  node2vec_rw_visitor() = default;

  explicit node2vec_rw_visitor(const vertex_type _vertex)
      : vertex(_vertex),
        previous_vertex(),
        wk_id(),
        walk_step(0) {}

  node2vec_rw_visitor(const vertex_type _vertex,
                      const vertex_type _previous_vertex,
                      const walk_id_type _wk_id,
                      const walk_length_type _walk_length)
      : vertex(_vertex),
        previous_vertex(_previous_vertex),
        wk_id(_wk_id),
        walk_step(_walk_length) {}

  template <typename VisitorQueueHandle>
  bool init_visit(graph_type &g, VisitorQueueHandle vis_queue, algorithm_data_type &algorithm_data) const {
    if (!algorithm_data.visited(vertex)) return false;

    const auto walker_id = algorithm_data.gen_walker_id();

    algorithm_data.store(walker_id, 0, vertex);

    const auto next_vertex = algorithm_data.select_neighbor(vertex, vertex, true);

    node2vec_rw_visitor new_visitor(next_vertex, vertex, walker_id, 2);
    vis_queue->queue_visitor(new_visitor);

    return true;
  }

  bool pre_visit(algorithm_data_type &algorithm_data) const {
    assert(walk_step <= algorithm_data.walk_length());

    algorithm_data.store(wk_id, walk_step - 1, vertex);

    if (walk_step == algorithm_data.walk_length()) // This vertex is the last visit
      return false;

    return true;
  }

  template <typename VisitorQueueHandle>
  bool visit(graph_type &g, VisitorQueueHandle vis_queue, algorithm_data_type &algorithm_data) const {
    assert(walk_step <= algorithm_data.walk_length());

    const auto next_vertex = algorithm_data.select_neighbor(previous_vertex, vertex, false);

    node2vec_rw_visitor new_visitor(next_vertex, vertex, wk_id, walk_step + 1);
    vis_queue->queue_visitor(new_visitor);

    return true;
  }

  friend inline bool operator>(const node2vec_rw_visitor &v1, const node2vec_rw_visitor &v2) {
    return false;
  }

  friend inline bool operator<(const node2vec_rw_visitor &v1, const node2vec_rw_visitor &v2) {
    return false;
  }

  vertex_type vertex;
  vertex_type previous_vertex;
  walk_id_type wk_id;
  walk_length_type walk_step;
};
}

#endif //HAVOQGT_NODE2VEC_RW_VISITOR_HPP
