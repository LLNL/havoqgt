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

#include <havoqgt/environment.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>
#include "k_breadth_first_search_sync_level_per_source.hpp"


namespace havoqgt
{

namespace eecc_source_select_mode_tag
{
struct random
{
  explicit random(const size_t _max_sources)
    : max_sources(_max_sources) { }

  const size_t max_sources;
};


struct max_lower
{
  explicit max_lower(const size_t _max_sources)
    : max_sources(_max_sources) { }

  const size_t max_sources;
};


struct high_degree
{
  explicit high_degree(const size_t _max_sources)
    : max_sources(_max_sources) { }

  const size_t max_sources;
};


struct level2
{
  explicit level2(const size_t _max_sources)
    : max_sources(_max_sources) { }

  const size_t max_sources;
};


struct far_and_hdeg
{
  far_and_hdeg(const size_t _max_num_far_sources, const size_t _max_num_hdeg_sources)
    : max_num_far_sources(_max_num_far_sources),
      max_num_hdeg_sources(_max_num_hdeg_sources) { }

  const size_t max_num_far_sources;
  const size_t max_num_hdeg_sources;
};
};

template <typename graph_t, int k_num_sources>
struct eecc_type
{
  struct vertex_data
  {
    using vertex_data_ecc_t = typename graph_t::template vertex_data<level_t, std::allocator<level_t>>;

    explicit vertex_data(const graph_t &graph)
      : lower(graph),
        upper(graph) { }

    void init()
    {
      lower.reset(std::numeric_limits<level_t>::min());
      upper.reset(std::numeric_limits<level_t>::max());
    }

    vertex_data_ecc_t lower;
    vertex_data_ecc_t upper;
  };
};

namespace detail
{
// -------------------------------------------------------------------------------------------------------------- //
// hash vertex id
// -------------------------------------------------------------------------------------------------------------- //
uint64_t hash_vertex_id(const uint64_t vid)
{
  return shifted_n_hash16(vid, 64);
}

// -------------------------------------------------------------------------------------------------------------- //
// compute_ecc_k_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void compute_ecc_k_source(const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                          iterator_t vitr, iterator_t end,
                          std::vector<level_t> &k_source_ecc)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      k_source_ecc[k] = std::max(kbfs_vertex_data.level[*vitr][k], k_source_ecc[k]);
    }
  }
}

// -------------------------------------------------------------------------------------------------------------- //
// select_source_in_local
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void select_source_in_local(const graph_t *const graph,
                            const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                            const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                            const size_t max_num_sources,
                            const bool consider_all_vertices,
                            const std::vector<std::function<uint64_t(const typename graph_t::vertex_locator)>> &score_calculater_list,
                            std::vector<typename graph_t::vertex_locator> &candidate_list,
                            iterator_t vitr, iterator_t end)
{
  auto is_candidate = [&graph, &kbfs_vertex_data, &ecc_vertex_data] (const typename graph_t::vertex_locator locator) -> bool {
    return ((kbfs_vertex_data.level[locator][0] != kbfs_type<graph_t, k_num_sources>::unvisited_level) &&
            (ecc_vertex_data.lower[locator] != ecc_vertex_data.upper[locator]));
  };
  auto is_prior = [graph, &score_calculater_list](const typename graph_t::vertex_locator locator1, const typename graph_t::vertex_locator locator2) -> bool {
    for (auto& score_calculator : score_calculater_list) {
      const uint64_t score1 = score_calculator(locator1);
      const uint64_t score2 = score_calculator(locator2);
      if (score1 != score2) return (score1 > score2);
    }
    // --- Final tie braker--- //
    const uint64_t vid1 = graph->locator_to_label(locator1);
    const uint64_t vid2 = graph->locator_to_label(locator2);
    const uint64_t h1 = hash_vertex_id(vid1);
    const uint64_t h2 = hash_vertex_id(vid2);
    return (h1 > h2);
  };

  for (; vitr != end; ++vitr) {
    auto locator = *vitr;
    if (!consider_all_vertices && !is_candidate(locator)) continue;
    if (candidate_list.size() < max_num_sources) {
      candidate_list.emplace_back(locator);
    } else {
      auto min_itr = std::min_element(candidate_list.begin(), candidate_list.end(),
                                      [&is_prior](const typename graph_t::vertex_locator x1,
                                                  const typename graph_t::vertex_locator x2) -> bool {
        return !is_prior(x1, x2);
      });
      if (is_prior(locator, *min_itr)) {
        *min_itr = locator; // kick out lowest priority element
      }
    }
  }
}

template <typename graph_t, int k_num_sources>
std::vector<typename graph_t::vertex_locator>
select_source_in_global(const graph_t *const graph,
                        const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                        const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                        const size_t max_num_sources,
                        const std::vector<std::function<uint64_t(const typename graph_t::vertex_locator)>> &score_calculater_list,
                        const std::vector<typename graph_t::vertex_locator> &candidate_list)
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
    local_vid_list.emplace_back(graph->locator_to_label(locator));
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
                    [&global_score_table_list](const uint64_t vid1, const uint64_t vid2) -> bool {
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

  std::vector<typename graph_t::vertex_locator> selected_source_list;
  for (uint64_t vid : global_vid_list) {
    selected_source_list.emplace_back(graph->label_to_locator(vid));
  }

  return selected_source_list;
}

template <typename graph_t, int k_num_sources>
std::vector<typename graph_t::vertex_locator>
select_source(const graph_t *const graph,
              const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
              const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
              const size_t max_num_sources, const bool consider_all_vertices,
              std::vector<std::function<uint64_t(const typename graph_t::vertex_locator)>> score_calculater_list)
{

  std::vector<typename graph_t::vertex_locator> local_candidate_list;
  select_source_in_local<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                 max_num_sources, consider_all_vertices,
                                                 score_calculater_list, local_candidate_list,
                                                 graph->vertices_begin(), graph->vertices_end());
  select_source_in_local<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                 max_num_sources, consider_all_vertices,
                                                 score_calculater_list, local_candidate_list,
                                                 graph->controller_begin(), graph->controller_end());

  return select_source_in_global<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                         max_num_sources, score_calculater_list, local_candidate_list);
}

// -------------------------------------------------------------------------------------------------------------- //
// bound_ecc
// -------------------------------------------------------------------------------------------------------------- //
struct pair_hash
{
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &x) const
  {
    const auto h1 = std::hash<T1>()(x.first);
    const auto h2 = std::hash<T2>()(x.second);
    return h1 ^ h2;
  }
};

template <typename graph_t, int k_num_sources, typename iterator_t>
std::pair<size_t, size_t> bound_ecc(const graph_t *const graph,
                                    const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                    const std::vector<level_t> &k_source_ecc,
                                    typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                    iterator_t vitr, iterator_t end)
{
  size_t num_bounded_vertices(0);
  size_t num_non_bounded_vertices(0);

  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level)
      continue; // skip unvisited vertices

    level_t &current_lower = ecc_vertex_data.lower[*vitr];
    level_t &current_upper = ecc_vertex_data.upper[*vitr];
    if (current_lower == current_upper) continue; // Exact ecc has been already found

    bool bounded = false;
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      const level_t level = kbfs_vertex_data.level[*vitr][k];

      current_lower = std::max(current_lower, std::max(level, static_cast<uint16_t>(k_source_ecc[k] - level)));
      current_upper = std::min(current_upper, static_cast<uint16_t>(k_source_ecc[k] + level));

      if (current_lower == current_upper) {
        bounded = true;
        break;
      }
    }
    num_bounded_vertices += bounded;
    num_non_bounded_vertices += !bounded;
  }

  return std::make_pair(num_bounded_vertices, num_non_bounded_vertices);
}

} // namespace detail


// -------------------------------------------------------------------------------------------------------------- //
// compute_eecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
size_t compute_eecc(const graph_t *const graph,
                    const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                    const std::vector<typename graph_t::vertex_locator> &source_locator_list,
                    typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  // ------------------------------ Compute ECC for k sources ------------------------------ //
  std::vector<level_t> k_source_ecc(source_locator_list.size(), 0);
  {
    const double time_start = MPI_Wtime();
    detail::compute_ecc_k_source<graph_t, k_num_sources>(kbfs_vertex_data,
                                                         graph->vertices_begin(), graph->vertices_end(), k_source_ecc);
    detail::compute_ecc_k_source<graph_t, k_num_sources>(kbfs_vertex_data,
                                                         graph->controller_begin(), graph->controller_end(), k_source_ecc);
    mpi_all_reduce_inplace(k_source_ecc, std::greater<level_t>(), MPI_COMM_WORLD);
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      if (source_locator_list[k].owner() == static_cast<uint32_t>(mpi_rank) || source_locator_list[k].is_delegate()) {
        ecc_vertex_data.lower[source_locator_list[k]] = ecc_vertex_data.upper[source_locator_list[k]] = k_source_ecc[k];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    const double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Compute ECC for k sources: " << time_end - time_start << std::endl;
      std::cout << "ECC: ";
      for (level_t ecc : k_source_ecc) std::cout << ecc << " ";
      std::cout << std::endl;
    }
  }

  // ------------------------------ Bound ECC ------------------------------ //
  size_t num_bounded_vertices(0);
  size_t num_non_bounded_vertices(0);
  {
    const double time_start = MPI_Wtime();
    {
      auto ret = detail::bound_ecc<graph_t, k_num_sources>(graph, kbfs_vertex_data, k_source_ecc,
                                                           ecc_vertex_data, graph->vertices_begin(), graph->vertices_end());
      num_bounded_vertices += ret.first;
      num_non_bounded_vertices += ret.second;
    }

    {
      auto ret = detail::bound_ecc<graph_t, k_num_sources>(graph, kbfs_vertex_data, k_source_ecc,
                                                           ecc_vertex_data, graph->controller_begin(), graph->controller_end());
      num_bounded_vertices += ret.first;
      num_non_bounded_vertices += ret.second;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const size_t global_num_bounded_vertices = mpi_all_reduce(num_bounded_vertices,
                                                              std::plus<size_t>(),
                                                              MPI_COMM_WORLD);
    const size_t global_num_non_bounded_vertices = mpi_all_reduce(num_non_bounded_vertices,
                                                                  std::plus<size_t>(),
                                                                  MPI_COMM_WORLD);
    const double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Update ECC bounds: " << time_end - time_start << std::endl;
      std::cout << "# bounded vertices: " << global_num_bounded_vertices << std::endl;
      std::cout << "# non-bounded vertices: " << global_num_non_bounded_vertices << std::endl;
    }
  }

  const size_t global_num_non_bounded_vertices = mpi_all_reduce(num_non_bounded_vertices, std::logical_or<size_t>(),
                                                                MPI_COMM_WORLD);
  return global_num_non_bounded_vertices;
}

// -------------------------------------------------------------------------------------------------------------- //
// select_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
std::vector<typename graph_t::vertex_locator>
select_initial_source(const graph_t *const graph,
                      const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                      const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                      const size_t max_num_sources)
{
  auto degree_score = [graph](const typename graph_t::vertex_locator locator) -> uint64_t {
    return graph->degree(locator);
  };
  const bool consider_all_vertices = true;
  return detail::select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                       max_num_sources, consider_all_vertices, {degree_score});
}

template <typename graph_t, int k_num_sources>
std::vector<typename graph_t::vertex_locator>
select_source(const graph_t *const graph,
              const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
              const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
              const size_t max_num_sources, const size_t previous_num_sources)
{
  auto degree_score = [graph](const typename graph_t::vertex_locator locator) -> uint64_t {
    return graph->degree(locator);
  };

  auto shell_score = [&kbfs_vertex_data, &previous_num_sources](const typename graph_t::vertex_locator locator) -> uint64_t {
//    level_t max(0);
//    for (int k = 0; k < previous_num_sources; ++k) {
//      max = std::max(kbfs_vertex_data.level[locator][k], max);
//    }
//    return max;

    uint64_t total(0);
    for (size_t k = 0; k < previous_num_sources; ++k) {
#ifdef DEBUG
      if (kbfs_vertex_data.level[locator][k] == kbfs_type<graph_t, k_num_sources>::unvisited_level) {
        std::abort();
      }
#endif
      total += kbfs_vertex_data.level[locator][k];
    }
    return total;
  };

  auto diff_score = [&ecc_vertex_data](const typename graph_t::vertex_locator locator) -> uint64_t {
    return ecc_vertex_data.upper[locator] - ecc_vertex_data.lower[locator];
  };

  auto level2_score = [&kbfs_vertex_data, &previous_num_sources](const typename graph_t::vertex_locator locator) -> uint64_t {
    uint64_t total(0);
    for (size_t k = 0; k < previous_num_sources; ++k) {
      total += (kbfs_vertex_data.level[locator][k] == 2);
    }
    return total;
  };

  const bool consider_all_vertices = false;
  auto source_list1 = detail::select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                                    max_num_sources / 2, consider_all_vertices,
                                                                    {degree_score, shell_score, diff_score, level2_score});
  auto source_list2 = detail::select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                                    max_num_sources / 4 * 3, consider_all_vertices,
                                                                    {shell_score, diff_score, level2_score, degree_score});
  auto source_list3 = detail::select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                                    max_num_sources, consider_all_vertices,
                                                                    {diff_score, level2_score, degree_score, shell_score});

  for (auto locator : source_list2) {
    if (std::find(source_list1.begin(), source_list1.end(), locator) == source_list1.end()) {
      source_list1.push_back(locator);
      if (source_list1.size() == max_num_sources / 4 * 3) break;
    }
  }

  for (auto locator : source_list3) {
    if (std::find(source_list1.begin(), source_list1.end(), locator) == source_list1.end()) {
      source_list1.push_back(locator);
      if (source_list1.size() == max_num_sources) break;
    }
  }

  return source_list1;
}

// -------------------------------------------------------------------------------------------------------------- //
// pluning
// -------------------------------------------------------------------------------------------------------------- //
template <typename Graph, int num_k_sources>
class pluning_visitor
{
 private:
  using kbfs_t = kbfs_type<Graph, num_k_sources>;

 public:
  typedef typename Graph::vertex_locator vertex_locator;

  pluning_visitor()
    : vertex(),
      ecc() { }

  explicit pluning_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      ecc() { }

#pragma GCC diagnostic pop

  pluning_visitor(vertex_locator _vertex, level_t _cee)
    : vertex(_vertex),
      ecc(_cee) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (std::get<0>(alg_data).level[vertex][0] == kbfs_t::unvisited_level) return false; // skip unvisited vertices
    if (std::get<1>(alg_data).lower[vertex] != std::get<1>(alg_data).upper[vertex]) return false;

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      if (!eitr.target().is_delegate()) // Don't send visitor to delegates because they are obviously not degree 1 vertices
        vis_queue->queue_visitor(pluning_visitor(eitr.target(), std::get<1>(alg_data).lower[vertex]));
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
    if (vertex.is_delegate()) { // for case, label:FLOW1
      init_visit(g, vis_queue, alg_data);
    } else {
      if (g.degree(vertex) == 1 && std::get<1>(alg_data).lower[vertex] != std::get<1>(alg_data).upper[vertex]) {
        std::get<1>(alg_data).lower[vertex] = std::get<1>(alg_data).upper[vertex] = ecc + 1;
        ++(std::get<2>(alg_data));
      }
    }

    return false;
  }

  friend inline bool operator>(const pluning_visitor &v1, const pluning_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const pluning_visitor &v1, const pluning_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  level_t ecc;
} __attribute__ ((packed));


template <typename graph_t, int k_num_sources>
void plun_single_degree_vertices(const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                 graph_t *const  graph,
                                 typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  typedef pluning_visitor<graph_t, k_num_sources> visitor_type;
  size_t num_pluned = 0;
  auto alg_data = std::forward_as_tuple(kbfs_vertex_data, ecc_vertex_data, num_pluned);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(graph, alg_data);
  // ------------------------------ Traversal ------------------------------ //
  {
    const double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); // init_visit -> queue
    MPI_Barrier(MPI_COMM_WORLD);
    const double time_end = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Pluning: " << time_end - time_start << "\t";
  }

  num_pluned = mpi_all_reduce(num_pluned, std::plus<level_t>(), MPI_COMM_WORLD);
  if (mpi_rank == 0) std::cout << "# pluned: " << num_pluned << std::endl;
}

// -------------------------------------------------------------------------------------------------------------- //
// compute_distance_score
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
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

template <typename graph_t, int k_num_sources>
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
template <typename graph_t, int k_num_sources, typename iterator_t>
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

template <typename graph_t, int k_num_sources>
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

template <typename graph_t, int k_num_sources, typename iterator_t>
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
      count_up(ecc_vertex_data.upper[*vitr], lower_count);
      count_up(ecc_vertex_data.lower[*vitr], upper_count);
    }
  }
}

template <typename graph_t, int k_num_sources>
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

} // namespace havoqgt
#endif //HAVOQGT_EXACT_ECCENTRICITY_HPP
