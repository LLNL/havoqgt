/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */
#ifndef HAVOQGT_NODE2VEC_RW_HPP
#define HAVOQGT_NODE2VEC_RW_HPP

#include <vector>
#include <algorithm>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/mpi.hpp>
#include <node2vec_rw/k_bfs.hpp>
#include <node2vec_rw/k_bfs_level_diff_table.hpp>
#include <node2vec_rw/neighbor_selector.hpp>
#include <node2vec_rw/node2vec_rw_visitor.hpp>

namespace node2vec_rw {
namespace detail {
using default_graph_allocator_type
= typename havoqgt::distributed_db::segment_manager_type::template allocator<void>::type;
}

template <typename _edge_weight_type = double,
          typename _graph_allocator_type = detail::default_graph_allocator_type,
          typename _edge_data_allocator_type = std::allocator<_edge_weight_type>>
class node2vec_rw {
 public:
  using edge_weight_type = _edge_weight_type;
  using graph_allocator_type = _graph_allocator_type;
  using edge_data_allocator_type = _edge_data_allocator_type;
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;
  using edge_weight_data_type = typename graph_type::template edge_data<edge_weight_type, edge_data_allocator_type>;

 private:
  using bfs_level_type = uint16_t;
  using k_bfs_level_diff_table_type = k_bfs_level_diff_table<graph_allocator_type>;
  using neighbor_selector_type = neighbor_selector_with_rejection_sampling<graph_allocator_type, edge_weight_data_type>;
  using rw_visitor_type = node2vec_rw_visitor<graph_allocator_type, edge_weight_data_type>;
  using rw_visitor_data_type = rw_visitor_data<graph_allocator_type, edge_weight_data_type>;
  using rw_visitor_queue_type = havoqgt::visitor_queue<rw_visitor_type,
                                                       havoqgt::detail::visitor_priority_queue,
                                                       graph_type,
                                                       rw_visitor_data_type &>;

 public:
  node2vec_rw(graph_type &graph,
              const edge_weight_data_type &edge_weight_data,
              const bool small_edge_weight_variance,
              const uint16_t walk_length,
              const double p,
              const double q,
              const MPI_Comm mpi_comm,
              const bool verbose = false)
      : m_graph(graph),
        m_edge_weight_data(edge_weight_data),
        m_ptr_k_bfs_level_diff_table(nullptr),
        m_rw_visitor_data(nullptr),
        m_rw_visitor_queue(nullptr),
        m_mpi_comm(mpi_comm),
        m_mpi_rank(0),
        m_mpi_size(0),
        m_verbose(verbose) {
    CHK_MPI(MPI_Comm_rank(m_mpi_comm, &m_mpi_rank));
    CHK_MPI(MPI_Comm_size(m_mpi_comm, &m_mpi_size));

    make_kbfs_diff_table();
    setup_rw_visitor(small_edge_weight_variance, walk_length, p, q);
  }

  std::vector<std::vector<vertex_type>> run_walker(const std::vector<vertex_type> &start_vertices) {
    m_rw_visitor_data->clear();
    m_rw_visitor_queue->init_visitor_traversal(start_vertices);
    return m_rw_visitor_data->get_walk_history();
  }

 private:
  void make_kbfs_diff_table() {
    k_bfs<graph_allocator_type, bfs_level_type> k_bfs(m_graph);

    std::mt19937 rand(static_cast<uint32_t>(m_graph.max_global_vertex_id() + m_mpi_size));
    const int num_sources = std::max(static_cast<int>(std::log10(m_graph.max_global_vertex_id())), 4);
    std::vector<vertex_type> source_candidate_list;
    for (int i = 0; i < num_sources; ++i) {
      source_candidate_list.emplace_back(m_graph.label_to_locator(rand() % m_graph.max_global_vertex_id()));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const double k_bfs_time_start = MPI_Wtime();
    auto k_bfs_level_table = k_bfs.run(source_candidate_list);
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_verbose && m_mpi_rank == 0) {
      std::cout << "k-BFS took:\t" << MPI_Wtime() - k_bfs_time_start << std::endl;
    }

    const double k_bfs_level_diff_compute_time_start = MPI_Wtime();
    m_ptr_k_bfs_level_diff_table = construct_k_bfs_level_diff_table(m_graph, k_bfs_level_table);
    MPI_Barrier(MPI_COMM_WORLD);
    if (m_verbose && m_mpi_rank == 0) {
      std::cout << "k-BFS level diff computing took:\t" << MPI_Wtime() - k_bfs_level_diff_compute_time_start
                << std::endl;
    }
  }

  void setup_rw_visitor(const bool small_edge_weight_variance, const uint16_t walk_length, const double p, const double q) {

    m_rw_visitor_data = std::make_unique<rw_visitor_data_type>(m_graph,
                                                               m_edge_weight_data,
                                                               small_edge_weight_variance,
                                                               *m_ptr_k_bfs_level_diff_table,
                                                               walk_length,
                                                               p,
                                                               q,
                                                               m_mpi_comm);
    m_rw_visitor_queue = std::make_unique<rw_visitor_queue_type>(&m_graph, std::ref(*m_rw_visitor_data));
  }

  graph_type &m_graph;
  const edge_weight_data_type &m_edge_weight_data;
  std::unique_ptr<k_bfs_level_diff_table<graph_allocator_type>> m_ptr_k_bfs_level_diff_table;
  std::unique_ptr<rw_visitor_data_type> m_rw_visitor_data;
  std::unique_ptr<rw_visitor_queue_type> m_rw_visitor_queue;
  MPI_Comm m_mpi_comm;
  int m_mpi_rank;
  int m_mpi_size;
  bool m_verbose;
};

} // namespace node2vec_rw

#endif //HAVOQGT_NODE2VEC_RW_HPP
