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
// Created by Iwabuchi, Keita on 2/4/18.
//

#ifndef HAVOQGT_EXACT_ECCENTRICITY_HPP
#define HAVOQGT_EXACT_ECCENTRICITY_HPP

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <functional>
#include <unordered_map>
#include <cassert>
#include <limits>

#include <havoqgt/environment.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>
#include "k_breadth_first_search_sync_level_per_source.hpp"


namespace havoqgt
{
template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity_vertex_data
{
 private:
  using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;
  using ecc_t = typename graph_t::template vertex_data<level_t, std::allocator<level_t>>;
  using just_solved_flag_t = typename graph_t::template vertex_data<bool, std::allocator<bool>>;

 public:
  explicit exact_eccentricity_vertex_data(const graph_t &graph)
    : m_graph(graph),
      m_lower(graph),
      m_upper(graph),
      m_just_solved(graph) { }

  void init()
  {
    m_lower.reset(std::numeric_limits<level_t>::min());
    m_upper.reset(std::numeric_limits<level_t>::max());
    reset();
  }

  void reset()
  {
    m_just_solved.reset(false);
  }

  level_t &lower(const typename graph_t::vertex_locator &vertex)
  {
    return m_lower[vertex];
  }

  level_t &upper(const typename graph_t::vertex_locator &vertex)
  {
    return m_upper[vertex];
  }

  void set_just_solved(const typename graph_t::vertex_locator &vertex)
  {
    m_just_solved[vertex] = true;
  }

  bool get_just_solved(const typename graph_t::vertex_locator &vertex)
  {
    return m_just_solved[vertex];
  }

 private:
  const graph_t &m_graph;
  ecc_t m_lower;
  ecc_t m_upper;
  just_solved_flag_t m_just_solved;
};

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity
{
 private:
  using self_type = exact_eccentricity<segment_manager_t, level_t, k_num_sources>;
  using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;
  using vertex_locator_t = typename graph_t::vertex_locator;
  using source_selection_strategy_t = std::function<uint64_t(const vertex_locator_t)>;
  using source_selection_strategy_list_t = std::vector<source_selection_strategy_t>;

 public:
  using kbfs_t = k_breadth_first_search<segment_manager_t, level_t, k_num_sources>;
  using kbfs_vertex_data_t = typename kbfs_t::vertex_data_t;
  using ecc_vertex_data_t = exact_eccentricity_vertex_data<segment_manager_t, level_t, k_num_sources>;

  class pruning_visitor;


  explicit exact_eccentricity<segment_manager_t, level_t, k_num_sources>(graph_t &graph)
    : m_graph(graph),
      m_kbfs(graph),
      m_ecc_vertex_data(graph),
      m_source_info_holder(),
      m_source_selection_strategy_list()
  {
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return degree_score(vertex);
    });
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return shell_score(vertex);
    });
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return diff_score(vertex);
    });
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return level2_score(vertex);
    });
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return min_lower(vertex);
    });
    m_source_selection_strategy_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
      return max_upper(vertex);
    });
  }


  void run()
  {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    m_ecc_vertex_data.init();
    m_source_info_holder.init(m_source_selection_strategy_list.size());

    while (true) {
      adaptively_select_source();
      if (mpi_rank == 0) {
        std::cout << "#souces: " << m_source_info_holder.source_list.size() << std::endl;

        std::cout << "ID: ";
        for (auto v : m_source_info_holder.source_list) std::cout << m_graph.locator_to_label(v) << " ";
        std::cout << std::endl;

        std::cout << "Strategy: ";
        for (auto s : m_source_info_holder.strategy_list) std::cout << s << " ";
        std::cout << std::endl;
      }

      m_kbfs.run(m_source_info_holder.source_list);

      compute_ecc_k_source();

      size_t num_solved;
      size_t num_unsolved;
      {
        const double time_start = MPI_Wtime();
        m_ecc_vertex_data.reset();
        std::tie(num_solved, num_unsolved) = bound_ecc();
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Bound: " << time_end - time_start << std::endl;
          std::cout << "#num_solved: " << num_solved << std::endl;
          std::cout << "#num_unsolved: " << num_unsolved << std::endl;
        }
      }
      if (num_unsolved == 0) break;

      size_t num_pruned;
      {
        const double time_start = MPI_Wtime();
        num_pruned = prun_single_degree_vertices();
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Pruning: " << time_end - time_start << std::endl;
          std::cout << "#pruned: " << num_pruned << std::endl;
        }
      }

      if (num_unsolved - num_pruned == 0) break;
    }
  }


 private:

  class source_info_holder_t {

    public:
    void init(const size_t num_strategy)
    {
      source_list.resize(k_num_sources); // dummy locators
      num_solved_list.resize(k_num_sources, 0);

      for (size_t i = 0; i < k_num_sources; ++i)
        strategy_list.push_back(i % num_strategy);
    }

    bool add_source(vertex_locator_t source, int strategy)
    {
      auto ret = std::find(source_list.begin(), source_list.end(), source);
      if (ret != source_list.end()) return false;

      source_list.emplace_back(source);
      num_solved_list.emplace_back(0);
      strategy_list.emplace_back(strategy);

      return true;
    }

    size_t num_source() const
    {
#ifdef DEBUG
      assert(source_list.size() == num_solved_list.size() && source_list.size() == strategy_list.size());
#endif
      return source_list.size();
    }

    std::vector<vertex_locator_t> source_list;
    std::vector<size_t> num_solved_list;
    std::vector<int> strategy_list;
  };

  // -------------------------------------------------------------------------------------------------------------- //
  // select source
  // -------------------------------------------------------------------------------------------------------------- //
  uint64_t hash_vertex_id(const uint64_t vid)
  {
    return detail::hash_nbits(vid, 64);
  }

  template <typename iterator_t>
  void select_source_in_local(const size_t max_num_sources, const bool consider_all_vertices,
                              const std::vector<std::function<uint64_t(const vertex_locator_t)>> &score_calculater_list,
                              std::vector<vertex_locator_t> &candidate_list,
                              iterator_t vitr, iterator_t end)
  {
    auto is_candidate = [this] (const vertex_locator_t& locator) -> bool {
      return (m_kbfs.vertex_data().visited_by(locator, 0)) && (m_ecc_vertex_data.lower(locator) == m_ecc_vertex_data.upper(locator));
    };
    auto is_prior = [this, &score_calculater_list](const vertex_locator_t& lhd, const vertex_locator_t& rhd) -> bool {
      for (auto& score_calculator : score_calculater_list) {
        const uint64_t score_lhd = score_calculator(lhd);
        const uint64_t score_rhd = score_calculator(rhd);
        if (score_lhd != score_rhd) return (score_lhd > score_rhd);
      }
      // --- Final tie braker--- //
      const uint64_t vid_lhd = m_graph.locator_to_label(lhd);
      const uint64_t vid_rhd = m_graph.locator_to_label(rhd);
      const uint64_t h_lhd = hash_vertex_id(vid_lhd);
      const uint64_t h_rhd = hash_vertex_id(vid_rhd);
      return (h_lhd > h_rhd);
    };

    for (; vitr != end; ++vitr) {
      auto locator = *vitr;
      if (!consider_all_vertices && !is_candidate(locator)) continue;
      if (candidate_list.size() < max_num_sources) {
        candidate_list.emplace_back(locator);
      } else {
        auto min_itr = std::min_element(candidate_list.begin(), candidate_list.end(),
                                        [&is_prior](const vertex_locator_t& lhd,
                                                    const vertex_locator_t& rhd) -> bool {
                                          return !is_prior(lhd, rhd);
                                        });
        if (is_prior(locator, *min_itr)) {
          *min_itr = locator; // kick out lowest priority element
        }
      }
    }
  }

  std::vector<vertex_locator_t>
  select_source_in_global(const size_t max_num_sources,
                          const std::vector<std::function<uint64_t(const vertex_locator_t)>> &score_calculater_list,
                          const std::vector<vertex_locator_t> &candidate_list)
  {
    // -------------------- Construct and exchange score lists in global -------------------- //
    std::vector<std::vector<uint64_t>> global_score_matrix(score_calculater_list.size());
    for (size_t i = 0; i < score_calculater_list.size(); ++i) {
      std::vector<uint64_t> local_score_list;
      for (auto locator : candidate_list) {
        local_score_list.emplace_back(score_calculater_list[i](locator));
      }
      mpi_all_gather(local_score_list, global_score_matrix[i], MPI_COMM_WORLD);
    }

    std::vector<uint64_t> local_vid_list;
    for (auto locator : candidate_list) {
      local_vid_list.emplace_back(m_graph.locator_to_label(locator));
    }
    std::vector<uint64_t> global_vid_list;
    havoqgt::mpi_all_gather(local_vid_list, global_vid_list, MPI_COMM_WORLD);

    // -------------------- Select high priority sources based on the exchanged scores -------------------- //
    std::vector<std::unordered_map<uint64_t, uint64_t>> global_score_table_list(global_score_matrix.size());
    for (size_t i = 0; i < global_score_matrix.size(); ++i) {
      auto& global_score_table = global_score_table_list[i];
      auto& global_score_list = global_score_matrix[i];
      for (size_t j = 0; j < global_vid_list.size(); ++j) {
        global_score_table[global_vid_list[j]] = global_score_list[j];
      }
    }

    const size_t final_size = std::min(static_cast<size_t>(max_num_sources), global_vid_list.size());
    std::partial_sort(global_vid_list.begin(), global_vid_list.begin() + final_size, global_vid_list.end(),
                      [this, &global_score_table_list](const uint64_t vid1, const uint64_t vid2) -> bool {
                        for (auto& score_table : global_score_table_list) {
                          const uint64_t score1 = score_table[vid1];
                          const uint64_t score2 = score_table[vid2];
                          if (score1 != score2) return (score1 > score2);
                        }
                        const uint64_t h1 = hash_vertex_id(vid1);
                        const uint64_t h2 = hash_vertex_id(vid2);
                        return (h1 > h2);
                      });
    global_vid_list.resize(final_size);

    std::vector<vertex_locator_t> selected_source_list;
    for (uint64_t vid : global_vid_list) {
      selected_source_list.emplace_back(m_graph.label_to_locator(vid));
    }

    return selected_source_list;
  }

  std::vector<vertex_locator_t>
  select_source(const size_t max_num_sources, const bool consider_all_vertices,
                std::vector<std::function<uint64_t(const vertex_locator_t)>> score_calculater_list)
  {

    std::vector<vertex_locator_t> local_candidate_list;
    select_source_in_local(max_num_sources, consider_all_vertices,
                           score_calculater_list, local_candidate_list,
                           m_graph.vertices_begin(), m_graph.vertices_end());
    select_source_in_local(max_num_sources, consider_all_vertices,
                           score_calculater_list, local_candidate_list,
                           m_graph.controller_begin(), m_graph.controller_end());

    return select_source_in_global(max_num_sources, score_calculater_list, local_candidate_list);
  }

  uint64_t degree_score(const vertex_locator_t vertex)
  {
    return m_graph.degree(vertex);
  }

  uint64_t shell_score(const vertex_locator_t vertex)
  {
    uint64_t total(0);
    for (size_t k = 0; k < m_source_info_holder.num_source(); ++k) {
#ifdef DEBUG
      if (!m_kbfs.vertex_data().visited_by(vertex, k)) {
        std::abort();
      }
#endif
      total += m_kbfs.vertex_data().level(vertex)[k];
    }
    return total;
  }

  uint64_t diff_score(const vertex_locator_t vertex)
  {
    return m_ecc_vertex_data.upper(vertex) - m_ecc_vertex_data.lower(vertex);
  }

  uint64_t level2_score(const vertex_locator_t vertex)
  {
    uint64_t total(0);
    for (size_t k = 0; k < m_source_info_holder.num_source(); ++k) {
      total += (m_kbfs.vertex_data().level(vertex)[k] == 2);
    }
    return total;
  }

  uint64_t min_lower(const vertex_locator_t vertex)
  {
#ifdef DEBUG
    if(!m_kbfs.vertex_data().visited_by(vertex, 0)) {
      std::abort();
    }
#endif
    return std::numeric_limits<level_t>::max() - m_ecc_vertex_data.lower(vertex);
  }

  uint64_t max_upper(const vertex_locator_t vertex)
  {
#ifdef DEBUG
    if(!m_kbfs.vertex_data().visited_by(vertex, 0)) {
      std::abort();
    }
#endif
    return m_ecc_vertex_data.upper(vertex);
  }

  void adaptively_select_source()
  {
    std::vector<size_t> num_solved_by_strategy(m_source_selection_strategy_list.size(), 0);
    for (size_t i = 0; i < m_source_info_holder.num_source(); ++i) {
      num_solved_by_strategy[m_source_info_holder.strategy_list[i]] = m_source_info_holder.num_solved_list[i] + 1;
    }
    const size_t num_total_solved = std::accumulate(num_solved_by_strategy.begin(), num_solved_by_strategy.end(), 0ULL);

    source_info_holder_t new_source_info_holder;
    for (int strategy_id = 0; strategy_id < m_source_selection_strategy_list.size(); ++strategy_id) {
      const size_t num_additiional_souces = (static_cast<double>(num_solved_by_strategy[strategy_id]) / num_total_solved) * k_num_sources;
      if (num_additiional_souces == 0) continue;

      const size_t new_num_sources = new_source_info_holder.num_source() + num_additiional_souces;
      const bool consider_all_vertices = false;
      auto &strategy = m_source_selection_strategy_list[strategy_id];
      auto source_candidate_list = select_source(new_num_sources, consider_all_vertices, {strategy});

      // ----- Merge source ----- //
      for (auto candidate : source_candidate_list) {
        new_source_info_holder.add_source(candidate, strategy_id);
        if (new_source_info_holder.num_source() == new_num_sources) break;
      }
    }

    m_source_info_holder = std::move(new_source_info_holder);
    assert(m_source_info_holder.num_source() > 0);
  }

  // -------------------------------------------------------------------------------------------------------------- //
  // compute_ecc_k_source
  // -------------------------------------------------------------------------------------------------------------- //
  template <typename iterator_t>
  void compute_ecc_k_source_helper(std::vector<level_t> &k_source_ecc, iterator_t vitr, iterator_t end)
  {
    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue;
      for (size_t k = 0; k < k_source_ecc.size(); ++k) {
        k_source_ecc[k] = std::max(m_kbfs.vertex_data().level(*vitr)[k], k_source_ecc[k]);
      }
    }
  }

  void compute_ecc_k_source()
  {
    std::vector<level_t> k_source_ecc(m_source_info_holder.num_source(), std::numeric_limits<level_t>::min());
    compute_ecc_k_source_helper(k_source_ecc, m_graph.vertices_begin(), m_graph.vertices_end());
    compute_ecc_k_source_helper(k_source_ecc, m_graph.delegate_vertices_begin(), m_graph.delegate_vertices_end());

    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    mpi_all_reduce_inplace(k_source_ecc, std::greater<level_t>(), MPI_COMM_WORLD);
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      auto& source = m_source_info_holder.source_list[k];
      if (source.owner() == static_cast<uint32_t>(mpi_rank) || source.is_delegate()) {
        m_ecc_vertex_data.lower(source) = m_ecc_vertex_data.upper(source) = k_source_ecc[k];
      }
    }

    if (mpi_rank == 0) {
      std::cout << "ECC for k sources: ";
      for (level_t ecc : k_source_ecc) std::cout << ecc << " ";
      std::cout << std::endl;
    }
  }

  // -------------------------------------------------------------------------------------------------------------- //
  // bound_ecc
  // -------------------------------------------------------------------------------------------------------------- //
  template <typename iterator_t>
  std::pair<size_t, size_t> bound_ecc_helper(iterator_t vitr, iterator_t end)
  {
    size_t num_solved(0);
    size_t num_unsolved(0);

    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue; // skip unvisited vertices

      level_t &current_lower = m_ecc_vertex_data.lower(*vitr);
      level_t &current_upper = m_ecc_vertex_data.upper(*vitr);
      if (current_lower == current_upper) continue; // Exact ecc has been already found

      for (size_t k = 0; k < m_source_info_holder.num_source(); ++k) {
        const level_t level = m_kbfs.vertex_data().level(*vitr)[k];
        const level_t k_ecc = m_ecc_vertex_data.lower(m_source_info_holder.source_list[k]); // TODO: this is not efficient

        current_lower = std::max(current_lower, std::max(level, static_cast<uint16_t>(k_ecc - level)));
        current_upper = std::min(current_upper, static_cast<uint16_t>(k_ecc + level));

        if (current_lower == current_upper) {
          ++m_source_info_holder.num_solved_list[k];
          m_ecc_vertex_data.set_just_solved(*vitr);
          break;
        }
      }
      num_solved += m_ecc_vertex_data.get_just_solved(*vitr);
      num_unsolved += !m_ecc_vertex_data.get_just_solved(*vitr);
    }

    return std::make_pair(num_solved, num_unsolved);
  }

  std::pair<size_t, size_t> bound_ecc()
  {
    const auto ret_vrtx = bound_ecc_helper(m_graph.vertices_begin(), m_graph.vertices_end());
    const auto ret_ctrl = bound_ecc_helper(m_graph.controller_begin(), m_graph.controller_end());

    size_t num_solved = ret_vrtx.first + ret_ctrl.first;
    size_t num_unsolved = ret_vrtx.second + ret_ctrl.second;

    num_solved = mpi_all_reduce(num_solved, std::plus<size_t>(), MPI_COMM_WORLD);
    num_unsolved = mpi_all_reduce(num_unsolved, std::plus<size_t>(), MPI_COMM_WORLD);

    return std::make_pair(num_solved, num_unsolved);
  }

  // -------------------------------------------------------------------------------------------------------------- //
  // prun_single_degree_vertices
  // -------------------------------------------------------------------------------------------------------------- //
  size_t prun_single_degree_vertices()
  {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    size_t local_num_pruned = 0;
    auto alg_data = std::forward_as_tuple(m_ecc_vertex_data, local_num_pruned);
    auto vq = create_visitor_queue<pruning_visitor, havoqgt::detail::visitor_priority_queue>(&m_graph, alg_data);
    vq.init_visitor_traversal_new();
    MPI_Barrier(MPI_COMM_WORLD);

    const size_t global_num_pruned = mpi_all_reduce(local_num_pruned, std::plus<level_t>(), MPI_COMM_WORLD);

    return global_num_pruned;
  }

  graph_t& m_graph;
  kbfs_t m_kbfs;
  ecc_vertex_data_t m_ecc_vertex_data;
  source_info_holder_t m_source_info_holder;
  source_selection_strategy_list_t m_source_selection_strategy_list;
};

template <typename segment_manager_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity<segment_manager_t, level_t, k_num_sources>::pruning_visitor
{
 private:
  enum index {
    ecc_data = 0,
    count_num_pruned = 1
  };

 public:
  pruning_visitor()
    : vertex(),
      ecc() { }

  explicit pruning_visitor(vertex_locator_t _vertex)
    : vertex(_vertex),
      ecc() { }

#pragma GCC diagnostic pop

  pruning_visitor(vertex_locator_t _vertex, level_t _cee)
    : vertex(_vertex),
      ecc(_cee) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (!std::get<index::ecc_data>(alg_data).get_just_solved(vertex)) return false;

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      if (!eitr.target().is_delegate()) // Don't send visitor to delegates because they are obviously not degree 1 vertices
        vis_queue->queue_visitor(pruning_visitor(eitr.target(), std::get<index::ecc_data>(alg_data).lower(vertex)));
    }

    return true; // trigger bcast from masters of delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const
  {
    return true;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    if (vertex.is_delegate()) { // for case, label:FLOW1
      init_visit(g, vis_queue, alg_data);
    } else {
      // ----- Update ecc of 1 degree vertices ----- //
      if (g.degree(vertex) == 1 && std::get<index::ecc_data>(alg_data).lower(vertex) != std::get<index::ecc_data>(alg_data).upper(vertex)) {
        std::get<index::ecc_data>(alg_data).lower(vertex) = std::get<index::ecc_data>(alg_data).upper(vertex) = ecc + 1;
        ++(std::get<index::count_num_pruned>(alg_data));
      }
    }

    return false;
  }

  friend inline bool operator>(const pruning_visitor &v1, const pruning_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const pruning_visitor &v1, const pruning_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator_t vertex;
  level_t ecc;
} __attribute__ ((packed));


//namespace detail
//{
#if 0
// -------------------------------------------------------------------------------------------------------------- //
// propagate ecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename Graph>
class propagate_ecc_visitor
{
 private:
  enum index {
    ecc_data = 0,
    count_num_solved = 1
  };

 public:
  typedef typename Graph::vertex_locator vertex_locator;

  propagate_ecc_visitor()
    : vertex(),
      ecc() { }

  explicit propagate_ecc_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      ecc() { }

#pragma GCC diagnostic pop

  propagate_ecc_visitor(vertex_locator _vertex, level_t _cee)
    : vertex(_vertex),
      ecc(_cee) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (!std::get<index::ecc_data>(alg_data).just_solved_flag[vertex]) return false;

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(propagate_ecc_visitor(eitr.target(), std::get<index::ecc_data>(alg_data).lower[vertex]));
    }

    return true; // trigger bcast from masters of delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const
  {
    return true;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    if (vertex.get_bcast()) { // for case, label:FLOW1, where delegates recieved request from their master node
      init_visit(g, vis_queue, alg_data);
    } else {
      // ----- Update ecc of 1 degree vertices ----- //
      auto& lower = std::get<index::ecc_data>(alg_data).lower[vertex];
      auto& upper = std::get<index::ecc_data>(alg_data).upper[vertex];
      if (lower != upper) {
        lower = std::max(lower, static_cast<level_t>(std::max(1, ecc - 1)));
        upper = std::min(upper, static_cast<level_t>(ecc + 1));

        std::get<index::count_num_solved>(alg_data) += (lower == upper);
      }
    }

    return false;
  }

  friend inline bool operator>(const propagate_ecc_visitor &v1, const propagate_ecc_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const propagate_ecc_visitor &v1, const propagate_ecc_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  level_t ecc;
} __attribute__ ((packed));


template <typename graph_t, uint32_t k_num_sources>
size_t propagate_ecc(graph_t *const graph,
                     typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  typedef propagate_ecc_visitor<graph_t> visitor_type;
  size_t local_num_solved = 0;
  auto alg_data = std::forward_as_tuple(ecc_vertex_data, local_num_solved);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(graph, alg_data);
  vq.init_visitor_traversal_new();
  MPI_Barrier(MPI_COMM_WORLD);

  const size_t global_num_solved = mpi_all_reduce(local_num_solved, std::plus<level_t>(), MPI_COMM_WORLD);
  return global_num_solved;
}

#endif

// -------------------------------------------------------------------------------------------------------------- //
// compute_eecc
// -------------------------------------------------------------------------------------------------------------- //
//std::pair<size_t, size_t> compute_eecc()
//{
//  int mpi_rank(0), mpi_size(0);
//  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
//  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
//
//  // ------------------------------ Compute ECC for k sources ------------------------------ //
//  std::vector<level_t> k_source_ecc(source_locator_list.size(), 0);
//  {
//    const double time_start = MPI_Wtime();
//    detail::compute_ecc_k_source<graph_t, k_num_sources>(kbfs_vertex_data,
//                                                         graph->vertices_begin(), graph->vertices_end(), k_source_ecc);
//    detail::compute_ecc_k_source<graph_t, k_num_sources>(kbfs_vertex_data,
//                                                         graph->controller_begin(), graph->controller_end(), k_source_ecc);
//    mpi_all_reduce_inplace(k_source_ecc, std::greater<level_t>(), MPI_COMM_WORLD);
//    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
//      auto& locator = source_locator_list[k];
//      if (locator.owner() == static_cast<uint32_t>(mpi_rank) || locator.is_delegate()) {
//        ecc_vertex_data.lower[locator] = ecc_vertex_data.upper[locator] = k_source_ecc[k];
//      }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    const double time_end = MPI_Wtime();
//    if (mpi_rank == 0) {
//      std::cout << "Compute ECC for k sources: " << time_end - time_start << std::endl;
//      std::cout << "ECC: ";
//      for (level_t ecc : k_source_ecc) std::cout << ecc << " ";
//      std::cout << std::endl;
//    }
//  }
//
//  std::vector<size_t> num_solved(k_source_ecc.size(), 0);
//  // ------------------------------ Bound ECC ------------------------------ //
//  {
//    const double time_start = MPI_Wtime();
//    ecc_vertex_data.just_solved_flag.reset(false);
//    {
//      std::vector<size_t> ret = detail::bound_ecc<graph_t, k_num_sources>(graph, kbfs_vertex_data, k_source_ecc, ecc_vertex_data,
//                                                           graph->vertices_begin(), graph->vertices_end());
//      std::transform(ret.begin(), ret.end(), num_solved.begin(), num_solved.begin(), std::plus<size_t>);
//    }
//
//    {
//      auto ret = detail::bound_ecc<graph_t, k_num_sources>(graph, kbfs_vertex_data, k_source_ecc, ecc_vertex_data,
//                                                           graph->controller_begin(), graph->controller_end());
//      std::transform(ret.begin(), ret.end(), num_solved.begin(), num_solved.begin(), std::plus<size_t>);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    const double time_end = MPI_Wtime();
//
//    mpi_all_reduce_inplace(num_solved, std::plus<size_t>, MPI_COMM_WORLD);
//    if (mpi_rank == 0) {
//      std::cout << "Update ECC bounds: " << time_end - time_start << std::endl;
//      std::cout << "# bounded vertices: " << std::accumulate(num_solved.begin(), num_solved.end(), 0) << std::endl;
//    }
//  }
//
//  // ------------------------------ Propagate just solved ecc ------------------------------ //
//  {
//    const double time_start = MPI_Wtime();
//    const size_t global_num_solved = propagate_ecc<graph_t, k_num_sources>(graph, ecc_vertex_data);
//    MPI_Barrier(MPI_COMM_WORLD);
//    const double time_end = MPI_Wtime();
//    if (mpi_rank == 0) std::cout << "Propagation: " << time_end - time_start << std::endl;
//    if (mpi_rank == 0) std::cout << "# solved: " << global_num_solved << std::endl;
//  }
//
//  // ------------------------------ Pruning ------------------------------ //
//  {
//    const double time_start = MPI_Wtime();
//    const size_t global_num_pruned = prun_single_degree_vertices<graph_t, k_num_sources>(graph, ecc_vertex_data);
//    MPI_Barrier(MPI_COMM_WORLD);
//    const double time_end = MPI_Wtime();
//    if (mpi_rank == 0) std::cout << "Pruning: " << time_end - time_start << std::endl;
//    if (mpi_rank == 0) std::cout << "# pruned: " << global_num_pruned << std::endl;
//    global_num_solved_vertices += global_num_pruned;
//  }
//
//  return std::make_pair(global_num_solved_vertices, global_num_unsolved_vertices);
//}

#if 0
// -------------------------------------------------------------------------------------------------------------- //
// compute_distance_score
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, uint32_t k_num_sources, typename iterator_t>
void compute_distance_score_histgram_helper(const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                            const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                            const size_t max_distance,
                                            std::vector<size_t> &histgram, iterator_t vitr, iterator_t end)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    const uint64_t x = ecc_vertex_data.upper[*vitr] - ecc_vertex_data.lower[*vitr];
    if (x > max_distance) ++histgram[max_distance + 1];
    else ++histgram[x];
  }
}

template <typename graph_t, uint32_t k_num_sources>
std::vector<size_t> compute_distance_score_histgram(const graph_t *const graph,
                                                    const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                    const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                                    const size_t max_distance)
{
  std::vector<size_t> histgram(max_distance + 2, 0);
  compute_distance_score_histgram_helper<graph_t, k_num_sources>(kbfs_vertex_data, ecc_vertex_data, max_distance,
                                                                 histgram, graph->vertices_begin(), graph->vertices_end());
  compute_distance_score_histgram_helper<graph_t, k_num_sources>(kbfs_vertex_data, ecc_vertex_data, max_distance,
                                                                 histgram, graph->controller_begin(), graph->controller_end());

  mpi_all_reduce_inplace(histgram, std::plus<size_t>(), MPI_COMM_WORLD);

  return histgram;
}

// -------------------------------------------------------------------------------------------------------------- //
// collect information about unsolved vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, uint32_t k_num_sources, typename iterator_t>
void dump_unsolved_vertices_info_helper(const graph_t *const graph,
                                        const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                        const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                        iterator_t vitr, iterator_t end, std::ofstream &ofs)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    if (ecc_vertex_data.upper[*vitr] == ecc_vertex_data.lower[*vitr]) continue;
    ofs << graph->locator_to_label(*vitr) << " " << graph->degree(*vitr) << " " << kbfs_vertex_data.level[*vitr]
        << " " << ecc_vertex_data.upper[*vitr] << " " << ecc_vertex_data.lower[*vitr] << std::endl;
  }
}

template <typename graph_t, uint32_t k_num_sources>
void dump_unsolved_vertices_info(const graph_t *const graph,
                                 const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                 const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                 const std::string &output_prefix)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::string output_path = output_prefix + "_" + std::to_string(mpi_rank);
  std::ofstream ofs(output_path);
  dump_unsolved_vertices_info_helper(graph, kbfs_vertex_data, ecc_vertex_data,
                                     graph->vertices_begin(), graph->vertices_end(), ofs);
  dump_unsolved_vertices_info_helper(graph, kbfs_vertex_data, ecc_vertex_data,
                                     graph->controller_begin(), graph->controller_end(), ofs);
}

template <typename graph_t, uint32_t k_num_sources, typename iterator_t>
void collect_unsolved_vertices_statistics_helper(const graph_t *const graph,
                                                 const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                 const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                                 const size_t num_sources,
                                                 std::map<size_t, size_t> &degree_count,
                                                 std::map<size_t, size_t> &level_count,
                                                 std::map<size_t, size_t> &lower_count,
                                                 std::map<size_t, size_t> &upper_count,
                                                 iterator_t vitr, iterator_t end)
{
  auto count_up = [](const size_t key, std::map<size_t, size_t>& table){
    if (table.count(key) == 0) table[key] = 0;
    ++table[key];
  };

  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    if (ecc_vertex_data.upper[*vitr] == ecc_vertex_data.lower[*vitr]) continue;
    {
      size_t degree = graph->degree(*vitr);
      if (degree > 10) {
        degree = std::pow(10, static_cast<uint64_t>(std::log10(degree)));
      }
      count_up(degree, degree_count);
    }
    {
      level_t average_level = 0;
      for (size_t k = 0; k < num_sources; ++k) {
        average_level += kbfs_vertex_data.level[*vitr][k];
      }
      average_level /= num_sources;
      count_up(average_level, level_count);
    }
    {
      count_up(ecc_vertex_data.lower[*vitr], lower_count);
      count_up(ecc_vertex_data.upper[*vitr], upper_count);
    }
  }
}

template <typename graph_t, uint32_t k_num_sources>
void collect_unsolved_vertices_statistics(const graph_t *const graph,
                                          const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                          const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                          const size_t num_sources)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::map<size_t, size_t> degree_count;
  std::map<size_t, size_t> level_count;
  std::map<size_t, size_t> lower_count;
  std::map<size_t, size_t> upper_count;

  collect_unsolved_vertices_statistics_helper<graph_t, k_num_sources>(graph,
                                                                      kbfs_vertex_data, ecc_vertex_data,
                                                                      num_sources,
                                                                      degree_count, level_count, lower_count, upper_count,
                                                                      graph->vertices_begin(), graph->vertices_end());
  collect_unsolved_vertices_statistics_helper<graph_t, k_num_sources>(graph,
                                                                      kbfs_vertex_data, ecc_vertex_data,
                                                                      num_sources,
                                                                      degree_count, level_count, lower_count, upper_count,
                                                                      graph->controller_begin(), graph->controller_end());

  auto all_gather = [](std::map<size_t, size_t>& table){
    std::vector<size_t> first;
    std::vector<size_t> second;
    for (auto elem : table) {
      first.emplace_back(elem.first);
      second.emplace_back(elem.second);
    }
    std::vector<size_t> gl_first;
    mpi_all_gather(first, gl_first, MPI_COMM_WORLD);
    std::vector<size_t> gl_second;
    mpi_all_gather(second, gl_second, MPI_COMM_WORLD);

    table.clear();
    for (size_t i = 0; i < gl_first.size(); ++i) {
      if (table.count(gl_first[i]) == 0) table[gl_first[i]] = 0;
      table[gl_first[i]] += gl_second[i];
    }
  };

  all_gather(degree_count);
  all_gather(level_count);
  all_gather(lower_count);
  all_gather(upper_count);

  auto print = [](std::map<size_t, size_t>& table) {
    for (auto elem : table) {
      std::cout << elem.first << " " << elem.second << ", ";
    }
    std::cout << std::endl;
  };

  if (mpi_rank == 0) {
    std::cout << "degree: "; print(degree_count);
    std::cout << "level: "; print(level_count);
    std::cout << "lower: "; print(lower_count);
    std::cout << "upper: "; print(upper_count);
  }
}
#endif

} // namespace havoqgt
#endif //HAVOQGT_EXACT_ECCENTRICITY_HPP
