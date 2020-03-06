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


#ifndef HAVOQGT_K_BFS_HPP
#define HAVOQGT_K_BFS_HPP

#include <iostream>
#include <cassert>
#include <memory>

#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/mpi.hpp>

namespace node2vec_rw {

template <typename graph_allocator_type, typename level_type>
class k_bfs_level_table {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using single_bfs_level_table = typename graph_type::template vertex_data<level_type, std::allocator<level_type>>;
  using table_type = std::vector<single_bfs_level_table>;
  using vertex_type = typename graph_type::vertex_locator;

 public:
  k_bfs_level_table(const graph_type &graph, const std::size_t k_size)
      : m_graph(graph),
        m_table(k_size, single_bfs_level_table(m_graph)) {
    for (auto &elem : m_table) {
      elem.reset(unvisited_level());
    }
  }

  level_type get_level(const vertex_type &source, const std::size_t k) const {
    assert(m_table.size() > k);
    return m_table[k][source];
  }

  void set_level(const vertex_type &source, const std::size_t k, const level_type level) {
    assert(level != unvisited_level());
    m_table[k][source] = level;
  }

  std::size_t k_size() const {
    return m_table.size();
  }

  bool visited(const vertex_type &vertex, const std::size_t k) const {
    assert(m_table.size() > k);
    return m_table[k][vertex] != unvisited_level();
  }

  single_bfs_level_table& bfs_level_data(const std::size_t k) {
    assert(m_table.size() > k);
    return m_table[k];
  }

  const single_bfs_level_table& bfs_level_data(const std::size_t k) const {
    assert(m_table.size() > k);
    return m_table[k];
  }

  void statistic() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    for (std::size_t k = 0; k < m_table.size(); ++k) {
      auto &bfs_level_data = m_table[k];

      std::size_t visited_total(0);
      for (level_type level = 0; level < std::numeric_limits<level_type>::max(); ++level) {
        std::size_t local_count(0);
        for (auto vitr = m_graph.vertices_begin(); vitr != m_graph.vertices_end(); ++vitr) {
          if (bfs_level_data[*vitr] == level) {
            ++local_count;
          }
        }
        for (auto vitr = m_graph.controller_begin(); vitr != m_graph.controller_end(); ++vitr) {
          if (bfs_level_data[*vitr] == level) {
            ++local_count;
          }
        }

        const std::size_t global_count = havoqgt::mpi_all_reduce(local_count, std::plus<std::size_t>(), MPI_COMM_WORLD);
        visited_total += global_count;
        if (mpi_rank == 0 && global_count > 0) {
          std::cout << "Level\t" << level << "\t" << global_count << std::endl;
        }
        if (global_count == 0) {
          break;
        }
      }  // end for level
      if (mpi_rank == 0) {
        std::cout << "Total\t" << visited_total << std::endl;
      }
    }
  }

 private:
  static constexpr level_type unvisited_level() {
    return std::numeric_limits<level_type>::max();
  }

  const graph_type &m_graph;
  table_type m_table;
};

template <typename graph_allocator_type, typename level_type>
class k_bfs {
 private:
  using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
  using vertex_type = typename graph_type::vertex_locator;
  using bfs_level_table_type = typename graph_type::template vertex_data<level_type, std::allocator<level_type>>;
  using k_bfs_level_table_type = k_bfs_level_table<graph_allocator_type, level_type>;

 public:
  explicit k_bfs(const graph_type &graph)
      : m_graph(graph) {}

  k_bfs_level_table_type run(const std::vector<vertex_type> &source_candidate_list) const {
    return run_kbfs(source_candidate_list);
  }

 private:
  k_bfs_level_table_type run_kbfs(const std::vector<vertex_type> &source_candidate_list) const {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    auto source_list = select_non_zero_degree_source(source_candidate_list);

    // -- Init k-BFS level table -- //
    k_bfs_level_table_type k_bfs_level_table(m_graph, source_list.size());
    MPI_Barrier(MPI_COMM_WORLD);

    // -- k-BFS core -- //
    for (int k = 0; k < source_list.size(); ++k) {
      if (mpi_rank == 0) {
        std::cout << "\nk-BFS " << k << std::endl;
      }
      run_single_bfs(source_list[k], &k_bfs_level_table.bfs_level_data(k));
      MPI_Barrier(MPI_COMM_WORLD);
    }

    return k_bfs_level_table;
  }

  template <typename bfs_level_data_type>
  void run_single_bfs(const vertex_type &source, bfs_level_data_type* bfs_level_data) const {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    typename graph_type::template vertex_data<vertex_type, std::allocator<vertex_type>> bfs_parent_data(m_graph);

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    havoqgt::breadth_first_search(&m_graph, *bfs_level_data, bfs_parent_data, source);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "BFS took: " << time_end - time_start << std::endl;
    }
    show_level_info(*bfs_level_data);
  }

  auto select_non_zero_degree_source(const std::vector<vertex_type> &source_candidate_list) const {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    std::vector<vertex_type> source_locator_list;

    for (auto candidate_itr = source_candidate_list.begin(), end = source_candidate_list.end();
         candidate_itr != end; ++candidate_itr) {
      vertex_type source = *candidate_itr;
      while (true) {
        if (source.is_delegate()) {
          break;
        }

        const std::size_t local_degree = (uint32_t(mpi_rank) == source.owner()) ? m_graph.degree(source) : 0;
        const std::size_t
            global_degree = havoqgt::mpi_all_reduce(local_degree, std::greater<std::size_t>(), MPI_COMM_WORLD);

        const bool already_exist = (std::find(source_candidate_list.begin(), candidate_itr, source) != candidate_itr);

        if (global_degree > 0 && !already_exist) break;

        source = m_graph.label_to_locator(m_graph.locator_to_label(source) + 1);
      }

      if (uint32_t(mpi_rank) == source.owner()) {
        std::cout << std::endl;
        if (source != *candidate_itr) {
          std::cout << "Vertex " << candidate_itr->local_id() << " has a degree of 0." << std::endl;
        }
        std::cout << "Source vertex\t" << m_graph.locator_to_label(source) << std::endl;
        std::cout << "delegate?\t" << source.is_delegate() << std::endl;
        std::cout << "local_id\t" << source.local_id() << std::endl;
        std::cout << "degree\t" << m_graph.degree(source) << std::endl;
      }

      source_locator_list.emplace_back(std::move(source));
    }

    return source_locator_list;
  }

  void show_level_info(const bfs_level_table_type &bfs_level_data) const {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    uint64_t visited_total(0);
    for (uint64_t level = 0; level < std::numeric_limits<uint16_t>::max(); ++level) {
      uint64_t local_count(0);
      for (auto vitr = m_graph.vertices_begin(); vitr != m_graph.vertices_end(); ++vitr) {
        if (bfs_level_data[*vitr] == level) {
          ++local_count;
        }
      }

      for (auto vitr = m_graph.controller_begin(); vitr != m_graph.controller_end(); ++vitr) {
        if (bfs_level_data[*vitr] == level) {
          ++local_count;
        }
      }

      uint64_t global_count = havoqgt::mpi_all_reduce(local_count,
                                                      std::plus<uint64_t>(), MPI_COMM_WORLD);
      visited_total += global_count;
      if (mpi_rank == 0 && global_count > 0) {
        std::cout << "Level\t" << level << ":\t" << global_count << std::endl;
      }
      if (global_count == 0) {
        break;
      }
    }  // end for level

    if (mpi_rank == 0) {
      std::cout << "Visited total\t" << visited_total << std::endl;
    }
  }

  const graph_type &m_graph;
};

}
#endif //HAVOQGT_K_BFS_HPP
