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
// find_highest_degree_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void find_highest_degree_vertices_helper(const graph_t *const graph, iterator_t vitr, iterator_t end,
                                         const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                         const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                         std::vector<size_t> &highest_degree_list,
                                         std::vector<uint64_t> &highest_degree_id_list)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level)
      continue; // skip unvisited vertices
    if (ecc_vertex_data.lower[*vitr] == ecc_vertex_data.upper[*vitr]) continue;

    const size_t degree = graph->degree(*vitr);
    if (highest_degree_list.size() < k_num_sources) {
      highest_degree_list.emplace_back(degree);
      highest_degree_id_list.emplace_back(graph->locator_to_label(*vitr));
    } else {
      const auto min_itr = std::min_element(highest_degree_list.begin(), highest_degree_list.end());
      const off_t min_pos = std::distance(highest_degree_list.begin(), min_itr);
      if (highest_degree_list[min_pos] < degree) {
        highest_degree_list[min_pos] = degree;
        highest_degree_id_list[min_pos] = graph->locator_to_label(*vitr);
      }
    }
  }
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> find_highest_degree_vertices(const graph_t *const graph,
                                                   const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                   const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  std::vector<size_t> global_highest_degree_list;
  std::vector<uint64_t> global_highest_degree_id_list;

  {
    std::vector<size_t> local_highest_degree_list;
    std::vector<uint64_t> local_highest_degree_id_list;
    find_highest_degree_vertices_helper<graph_t, k_num_sources>(graph,
                                                                graph->vertices_begin(),
                                                                graph->vertices_end(),
                                                                kbfs_vertex_data,
                                                                ecc_vertex_data,
                                                                local_highest_degree_list,
                                                                local_highest_degree_id_list);
    find_highest_degree_vertices_helper<graph_t, k_num_sources>(graph,
                                                                graph->controller_begin(),
                                                                graph->controller_end(),
                                                                kbfs_vertex_data,
                                                                ecc_vertex_data,
                                                                local_highest_degree_list,
                                                                local_highest_degree_id_list);
    assert(local_highest_degree_list.size() == global_highest_degree_id_list.size());

    havoqgt::mpi_all_gather(local_highest_degree_list, global_highest_degree_list, MPI_COMM_WORLD);
    havoqgt::mpi_all_gather(local_highest_degree_id_list, global_highest_degree_id_list, MPI_COMM_WORLD);
  }

  const size_t num_total_candidates = global_highest_degree_list.size();
  std::vector<uint64_t> selected_vertex_id_list;
  if (num_total_candidates > k_num_sources) {
    std::vector<std::pair<size_t, uint64_t>> table;
    table.resize(num_total_candidates);
    for (size_t i = 0; i < global_highest_degree_list.size(); ++i) {
      table[i] = std::make_pair(global_highest_degree_list[i], global_highest_degree_id_list[i]);
    }
    std::partial_sort(table.begin(), table.begin() + k_num_sources, table.end(),
                      [&](const std::pair<size_t, uint64_t> &a, const std::pair<size_t, uint64_t> &b) -> bool {
                        return (a.first > b.first);
                      });

    selected_vertex_id_list.resize(k_num_sources);
    for (size_t i = 0; i < k_num_sources; ++i) {
      selected_vertex_id_list[i] = table[i].second;
    }
  } else {
    selected_vertex_id_list = std::move(global_highest_degree_id_list);
  }

  return selected_vertex_id_list;
}

// -------------------------------------------------------------------------------------------------------------- //
// find_max_lower_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void find_max_lower_vertices_helper(const graph_t *const graph,
                                    iterator_t vitr, iterator_t end,
                                    const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                    const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                    std::vector<level_t> &level_list, std::vector<uint64_t> &vertex_id_list,
                                    const size_t max_num_sources)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level)
      continue; // skip unvisited vertices
    if (ecc_vertex_data.lower[*vitr] == ecc_vertex_data.upper[*vitr]) continue;

    const auto min_itr = std::min_element(kbfs_vertex_data.level[*vitr].begin(), kbfs_vertex_data.level[*vitr].end());
    const level_t min_level = *min_itr;
    if (level_list.size() < max_num_sources) {
      level_list.emplace_back(min_level);
      vertex_id_list.emplace_back(graph->locator_to_label(*vitr));
    } else {
      auto min_level_itr = std::min_element(level_list.begin(), level_list.end());
      const off_t min_pos = std::distance(level_list.begin(), min_level_itr);
      if (level_list[min_pos] < min_level) {
        level_list[min_pos] = min_level;
        vertex_id_list[min_pos] = graph->locator_to_label(*vitr);
      }
    }
  }
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> find_max_lower_vertices(const graph_t *const graph,
                                              const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                              const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                              const eecc_source_select_mode_tag::max_lower tag)
{
  std::vector<level_t> global_level_list;
  std::vector<uint64_t> global_vertex_id_list;

  {
    std::vector<level_t> local_level_list;
    std::vector<uint64_t> local_vertex_id_list;
    find_max_lower_vertices_helper<graph_t, k_num_sources>(graph, graph->vertices_begin(), graph->vertices_end(),
                                                           kbfs_vertex_data, ecc_vertex_data,
                                                           local_level_list, local_vertex_id_list, tag.max_sources);
    find_max_lower_vertices_helper<graph_t, k_num_sources>(graph, graph->controller_begin(), graph->controller_end(),
                                                           kbfs_vertex_data, ecc_vertex_data,
                                                           local_level_list, local_vertex_id_list, tag.max_sources);
    assert(local_level_list.size() == local_vertex_id_list.size());

    havoqgt::mpi_all_gather(local_level_list, global_level_list, MPI_COMM_WORLD);
    havoqgt::mpi_all_gather(local_vertex_id_list, global_vertex_id_list, MPI_COMM_WORLD);
  }

  const size_t num_total_candidates = global_level_list.size();
  std::vector<uint64_t> selected_vertex_id_list;
  if (num_total_candidates > tag.max_sources) {
    std::vector<std::pair<level_t, uint64_t>> table;
    table.resize(num_total_candidates);
    for (size_t i = 0; i < global_level_list.size(); ++i) {
      table[i] = std::make_pair(global_level_list[i], global_vertex_id_list[i]);
    }
    std::partial_sort(table.begin(), table.begin() + tag.max_sources, table.end(),
                      [&](const std::pair<level_t, uint64_t> &a,
                          const std::pair<level_t, uint64_t> &b) -> bool {
                        return (a.first > b.first);
                      });

    selected_vertex_id_list.resize(tag.max_sources);
    for (size_t i = 0; i < tag.max_sources; ++i) {
      selected_vertex_id_list[i] = table[i].second;
    }
  } else {
    selected_vertex_id_list = std::move(global_vertex_id_list);
  }

  return selected_vertex_id_list;
}

// -------------------------------------------------------------------------------------------------------------- //
// randomly_select_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
std::vector<uint64_t> randomly_select_source(graph_t *graph,
                                             const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                             const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector<uint64_t> local_vertex_id_list;
  std::vector<uint64_t> global_vertex_id_list;

  size_t num_non_resolved_vertices(0);
  for (auto itr = graph->vertices_begin(), end = graph->vertices_end(); itr != end; ++itr) {
    if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    if (ecc_vertex_data.lower[*itr] == ecc_vertex_data.upper[*itr]) continue;
    ++num_non_resolved_vertices;
  }

  for (auto itr = graph->controller_begin(), end = graph->controller_end(); itr != end; ++itr) {
    if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    if (ecc_vertex_data.lower[*itr] == ecc_vertex_data.upper[*itr]) continue;
    ++num_non_resolved_vertices;
  }

  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<uint64_t> dis(0, num_non_resolved_vertices - 1);

  // ----- randomly pick up vertices up to whichever smaller number of k sources or non resolved vertices ----- //
  while (local_vertex_id_list.size() < std::min(static_cast<size_t>(k_num_sources), num_non_resolved_vertices)) {
    for (auto itr = graph->vertices_begin(), end = graph->vertices_end(); itr != end; ++itr) {
      if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
      if (ecc_vertex_data.lower[*itr] == ecc_vertex_data.upper[*itr]) continue;
      if (dis(gen) >= k_num_sources) continue;
      uint64_t vid = graph->locator_to_label(*itr);
      if (std::find(local_vertex_id_list.begin(), local_vertex_id_list.end(), vid) != local_vertex_id_list.end())
        continue;
      local_vertex_id_list.emplace_back(graph->locator_to_label(*itr));
    }
    for (auto itr = graph->controller_begin(), end = graph->controller_end(); itr != end; ++itr) {
      if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
      if (ecc_vertex_data.lower[*itr] == ecc_vertex_data.upper[*itr]) continue;
      if (dis(gen) >= k_num_sources) continue;
      uint64_t vid = graph->locator_to_label(*itr);
      if (std::find(local_vertex_id_list.begin(), local_vertex_id_list.end(), vid) != local_vertex_id_list.end())
        continue;
      local_vertex_id_list.emplace_back(graph->locator_to_label(*itr));
    }
  }

  havoqgt::mpi_all_gather(local_vertex_id_list, global_vertex_id_list, MPI_COMM_WORLD);

  if (global_vertex_id_list.size() > k_num_sources) {
    std::mt19937_64 gen(global_vertex_id_list.size());
    std::shuffle(global_vertex_id_list.begin(), global_vertex_id_list.end(), gen);
    global_vertex_id_list.resize(k_num_sources);
  }

  return global_vertex_id_list;
}


// -------------------------------------------------------------------------------------------------------------- //
// select_2level_away_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
bool is_2_level_away(const typename graph_t::vertex_locator vertex,
                     const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                     const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  if (kbfs_vertex_data.level[vertex][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) return false;
  if (ecc_vertex_data.lower[vertex] == ecc_vertex_data.upper[vertex]) return false;

  for (level_t lv : kbfs_vertex_data.level[vertex]) {
    if (lv == 2) {
      return true;
    }
  }
  return false;
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_2_level_away_vertices(graph_t *graph,
                                                   const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                   const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector<uint64_t> local_vertex_id_list;
  std::vector<uint64_t> global_vertex_id_list;

  size_t num_non_resolved_vertices(0);
  for (auto itr = graph->vertices_begin(), end = graph->vertices_end(); itr != end; ++itr) {
    num_non_resolved_vertices += is_2_level_away<graph_t, k_num_sources>(*itr, kbfs_vertex_data, ecc_vertex_data);
  }
  for (auto itr = graph->controller_begin(), end = graph->controller_end(); itr != end; ++itr) {
    num_non_resolved_vertices += is_2_level_away<graph_t, k_num_sources>(*itr, kbfs_vertex_data, ecc_vertex_data);
  }

  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<uint64_t> dis(0, num_non_resolved_vertices - 1);

  // ---- randomly pick up vertices which are 2 levels away----- //
  while (local_vertex_id_list.size() < std::min(static_cast<size_t>(k_num_sources), num_non_resolved_vertices)) {
    for (auto itr = graph->vertices_begin(), end = graph->vertices_end(); itr != end; ++itr) {
      if (!is_2_level_away<graph_t, k_num_sources>(*itr, kbfs_vertex_data, ecc_vertex_data)) continue;
      if (dis(gen) >= k_num_sources) continue;
      uint64_t vid = graph->locator_to_label(*itr);
      if (std::find(local_vertex_id_list.begin(), local_vertex_id_list.end(), vid) != local_vertex_id_list.end())
        continue;
      local_vertex_id_list.emplace_back(graph->locator_to_label(*itr));
    }
    for (auto itr = graph->controller_begin(), end = graph->controller_end(); itr != end; ++itr) {
      if (!is_2_level_away<graph_t, k_num_sources>(*itr, kbfs_vertex_data, ecc_vertex_data)) continue;
      if (dis(gen) >= k_num_sources) continue;
      uint64_t vid = graph->locator_to_label(*itr);
      if (std::find(local_vertex_id_list.begin(), local_vertex_id_list.end(), vid) != local_vertex_id_list.end())
        continue;
      local_vertex_id_list.emplace_back(graph->locator_to_label(*itr));
    }
  }

  havoqgt::mpi_all_gather(local_vertex_id_list, global_vertex_id_list, MPI_COMM_WORLD);

  if (global_vertex_id_list.size() > k_num_sources) {
    std::mt19937_64 gen(global_vertex_id_list.size());
    std::shuffle(global_vertex_id_list.begin(), global_vertex_id_list.end(), gen);
    global_vertex_id_list.resize(k_num_sources);
  }

  return global_vertex_id_list;
}

// -------------------------------------------------------------------------------------------------------------- //
// compute_ecc_k_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void compute_ecc_k_source(iterator_t vitr, iterator_t end,
                          const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
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
template <typename graph_t, typename iterator_t>
void select_source_in_local_helper(const graph_t *const graph, iterator_t vitr, iterator_t end,
                                   const size_t max_size,
                                   std::function<bool(const uint64_t)> &is_candidate,
                                   std::function<bool(const uint64_t, const uint64_t)> &is_prior,
                                   std::vector<uint64_t> &vid_list)
{
  for (; vitr != end; ++vitr) {
    const uint64_t vid = graph->locator_to_label(*vitr);
    if (!is_candidate(vid)) continue;
    if (vid_list.size() < max_size) {
      vid_list.emplace_back(vid);
    } else {
      auto min_itr = std::min_element(vid_list.begin(), vid_list.end(),
                                      [&is_prior](const uint64_t x1, const uint64_t x2) -> bool {
                                        return !is_prior(x1, x2);
                                      });
      uint64_t &min_vid = *min_itr;
      if (is_prior(vid, min_vid)) {
        min_vid = vid; // kick out lowest priority element
      }
    }
  }
}

template <typename graph_t>
std::vector<uint64_t> select_source_in_local(const graph_t *const graph,
                                             const size_t max_size,
                                             std::function<bool(const uint64_t)> is_candidate,
                                             std::function<bool(const uint64_t, const uint64_t)> is_prior)
{
  std::vector<uint64_t> local_vid_list;
  select_source_in_local_helper(graph, graph->vertices_begin(), graph->vertices_end(),
                                max_size, is_candidate, is_prior, local_vid_list);
  select_source_in_local_helper(graph, graph->controller_begin(), graph->controller_end(),
                                max_size, is_candidate, is_prior, local_vid_list);

  return local_vid_list;
}

// -------------------------------------------------------------------------------------------------------------- //
// select_far_and_hdeg_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_far_and_hdeg_vertices_helper(graph_t *graph,
                                                          const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                          const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                                          const size_t max_num_sources, const bool prior_distance)
{
  auto local_vid_list = detail::select_source_in_local(graph, max_num_sources,
                                                       [&graph, &kbfs_vertex_data, &ecc_vertex_data](const uint64_t vid) -> bool {
                                                         auto locator = graph->label_to_locator(vid);
                                                         return ((kbfs_vertex_data.level[locator][0] !=
                                                                  kbfs_type<graph_t, k_num_sources>::unvisited_level) &&
                                                                 (ecc_vertex_data.lower[locator] !=
                                                                  ecc_vertex_data.upper[locator]));
                                                       },
                                                       [&graph, &kbfs_vertex_data, &ecc_vertex_data, &prior_distance](
                                                         const uint64_t vid1,
                                                         const uint64_t vid2) -> bool {
                                                         auto locator1 = graph->label_to_locator(vid1);
                                                         auto locator2 = graph->label_to_locator(vid2);

                                                         const size_t distance_score1 =
                                                           ecc_vertex_data.upper[locator1] - ecc_vertex_data.lower[locator1];
                                                         const size_t distance_score2 =
                                                           ecc_vertex_data.upper[locator2] - ecc_vertex_data.lower[locator2];
                                                         const size_t degree1 = graph->degree(locator1);
                                                         const size_t degree2 = graph->degree(locator2);

                                                         if (prior_distance) {
                                                           if (distance_score1 != distance_score2)
                                                             return (distance_score1 > distance_score2);
                                                           if (degree1 != degree2)
                                                             return (degree1 > degree2);
                                                         } else {
                                                           if (degree1 != degree2)
                                                             return (degree1 > degree2);
                                                           if (distance_score1 != distance_score2)
                                                             return (distance_score1 > distance_score2);
                                                         }

                                                         const uint64_t h1 = shifted_n_hash16(vid1, 64);
                                                         const uint64_t h2 = shifted_n_hash16(vid2, 64);
                                                         return (h1 > h2);
                                                       });

  std::vector<uint64_t> global_distance_score_list;
  {
    std::vector<uint64_t> local_distance_score_list;
    for (auto vid : local_vid_list) {
      auto locator = graph->label_to_locator(vid);
      local_distance_score_list.emplace_back(ecc_vertex_data.upper[locator] - ecc_vertex_data.lower[locator]);
    }
    havoqgt::mpi_all_gather(local_distance_score_list, global_distance_score_list, MPI_COMM_WORLD);
  }

  std::vector<size_t> global_degree_list;
  {
    std::vector<size_t> local_degree_list;
    for (auto vid : local_vid_list) {
      auto locator = graph->label_to_locator(vid);
      local_degree_list.emplace_back(graph->degree(locator));
    }
    havoqgt::mpi_all_gather(local_degree_list, global_degree_list, MPI_COMM_WORLD);
  }

  std::vector<uint64_t> global_vid_list;
  havoqgt::mpi_all_gather(local_vid_list, global_vid_list, MPI_COMM_WORLD);

  std::unordered_map<uint64_t, uint64_t> distance_score_table;
  std::unordered_map<uint64_t, size_t> degree_table;
  {
    for (size_t i = 0; i < global_vid_list.size(); ++i) {
      distance_score_table[global_vid_list[i]] = global_distance_score_list[i];
      degree_table[global_vid_list[i]] = global_degree_list[i];
    }
  }

  const size_t final_size = std::min(static_cast<size_t>(max_num_sources), global_vid_list.size());
  std::partial_sort(global_vid_list.begin(), global_vid_list.begin() + final_size, global_vid_list.end(),
                    [&distance_score_table, &degree_table, &prior_distance](const uint64_t vid1,
                                                                            const uint64_t vid2) -> bool {
                      const uint64_t distance_score1 = distance_score_table[vid1];
                      const uint64_t distance_score2 = distance_score_table[vid2];
                      const size_t degree1 = degree_table[vid1];
                      const size_t degree2 = degree_table[vid2];

                      if (prior_distance) {
                        if (distance_score1 != distance_score2)
                          return (distance_score1 > distance_score2);
                        if (degree1 != degree2)
                          return (degree1 > degree2);
                      } else {
                        if (degree1 != degree2)
                          return (degree1 > degree2);
                        if (distance_score1 != distance_score2)
                          return (distance_score1 > distance_score2);
                      }

                      const uint64_t h1 = shifted_n_hash16(vid1, 64);
                      const uint64_t h2 = shifted_n_hash16(vid2, 64);
                      return (h1 > h2);
                    });
  global_vid_list.resize(final_size);

  return global_vid_list;
};

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_far_and_hdeg_vertices(graph_t *graph,
                                                   const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                                   const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                                   eecc_source_select_mode_tag::far_and_hdeg tag)
{
  auto vid_list1 = select_far_and_hdeg_vertices_helper<graph_t, k_num_sources>(graph,
                                                                               kbfs_vertex_data, ecc_vertex_data,
                                                                               tag.max_num_far_sources, true);
  auto vid_list2 = select_far_and_hdeg_vertices_helper<graph_t, k_num_sources>(graph,
                                                                               kbfs_vertex_data, ecc_vertex_data,
                                                                               tag.max_num_hdeg_sources, false);

  // --- Merge selected sources removing duplicates --- //
  vid_list1.insert(vid_list1.end(), vid_list2.begin(), vid_list2.end());
  std::sort(vid_list1.begin(), vid_list1.end());
  auto end = std::unique(vid_list1.begin(), vid_list1.end());
  vid_list1.resize(std::distance(vid_list1.begin(), end));

  return vid_list1;
};


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


//template <typename level_t>
//std::unordered_map<std::pair<level_t, level_t>, size_t, pair_hash> bound_progress;

template <typename graph_t, int k_num_sources, typename iterator_t>
std::pair<size_t, size_t> bound_ecc(const graph_t *const graph,
                                    iterator_t vitr,
                                    iterator_t end,
                                    const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                    const std::vector<level_t> &k_source_ecc,
                                    typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
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
//      else {
//        const auto key = std::make_pair(current_lower, current_upper);
//        if (bound_progress.count(key) == 0) {
//          bound_progress[key] = 0;
//        }
//        ++bound_progress[key];
//      }
    }
    num_bounded_vertices += bounded;
    num_non_bounded_vertices += !bounded;
  }

  return std::make_pair(num_bounded_vertices, num_non_bounded_vertices);
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::random)
{
  return detail::randomly_select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::max_lower tag)
{
  return detail::find_max_lower_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data, tag);
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::high_degree)
{
  return detail::find_highest_degree_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::level2)
{
  return detail::select_2_level_away_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);
}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::far_and_hdeg tag)
{
  return detail::select_far_and_hdeg_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data, tag);
}

} // namespace detail


// -------------------------------------------------------------------------------------------------------------- //
// compute_eecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources>
size_t compute_eecc(graph_t *graph,
                    typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                    typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                    const std::vector<typename graph_t::vertex_locator> &source_locator_list)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  // ------------------------------ Compute ECC for k sources ------------------------------ //
  std::vector<level_t> k_source_ecc(source_locator_list.size(), 0);
  {
    const double time_start = MPI_Wtime();
    detail::compute_ecc_k_source<graph_t, k_num_sources>(graph->vertices_begin(), graph->vertices_end(),
                                                         kbfs_vertex_data, k_source_ecc);
    detail::compute_ecc_k_source<graph_t, k_num_sources>(graph->controller_begin(), graph->controller_end(),
                                                         kbfs_vertex_data, k_source_ecc);
    mpi_all_reduce_inplace(k_source_ecc, std::greater<level_t>(), MPI_COMM_WORLD);
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      if (source_locator_list[k].owner() == mpi_rank || source_locator_list[k].is_delegate()) {
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
      auto ret = detail::bound_ecc<graph_t, k_num_sources>(graph, graph->vertices_begin(), graph->vertices_end(),
                                                           kbfs_vertex_data, k_source_ecc, ecc_vertex_data);
      num_bounded_vertices += ret.first;
      num_non_bounded_vertices += ret.second;
    }

    {
      auto ret = detail::bound_ecc<graph_t, k_num_sources>(graph, graph->controller_begin(), graph->controller_end(),
                                                           kbfs_vertex_data, k_source_ecc, ecc_vertex_data);
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
template <typename graph_t, int k_num_sources, typename eecc_source_select_mode_tag_t>
std::vector<typename graph_t::vertex_locator>
select_source(graph_t *graph,
              typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
              typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
              eecc_source_select_mode_tag_t tag)
{
  auto source_id_list = detail::select_source_helper<graph_t, k_num_sources>(graph,
                                                                             kbfs_vertex_data,
                                                                             ecc_vertex_data,
                                                                             tag);
  std::vector<typename graph_t::vertex_locator> source_locator_list;
  for (auto id : source_id_list) {
    source_locator_list.emplace_back(graph->label_to_locator(id));
  }
  return source_locator_list;
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
  bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (std::get<0>(alg_data).level[vertex][0] == kbfs_t::unvisited_level) return false; // skip unvisited vertices
    if (std::get<1>(alg_data).lower[vertex] != std::get<1>(alg_data).upper[vertex]) return false;

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      if (!eitr.target().is_delegate()) // Don't send visitor to delegates
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
      assert(std::get<1>(alg_data).lower[vertex] == std::get<1>(alg_data).upper[vertex]);
      for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
        vis_queue->queue_visitor(pluning_visitor(eitr.target(), std::get<1>(alg_data).lower[vertex]));
      }
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
void plun_single_degree_vertices(graph_t *graph,
                                 typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
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
  if (mpi_rank == 0) std::cout << "Plunned: " << num_pluned << std::endl;
}


// -------------------------------------------------------------------------------------------------------------- //
// find_max_ecc_bound_from_neighbor
// -------------------------------------------------------------------------------------------------------------- //
template <typename Graph, int num_k_sources>
class ecc_bound_info_prop_visitor
{
 private:
  using kbfs_t = kbfs_type<Graph, num_k_sources>;

 public:
  typedef typename Graph::vertex_locator vertex_locator;

  ecc_bound_info_prop_visitor()
    : vertex(),
      lower(),
      upper() { }

  explicit ecc_bound_info_prop_visitor(vertex_locator _vertex)
    : vertex(_vertex),
      lower(),
      upper() { }

#pragma GCC diagnostic pop

  ecc_bound_info_prop_visitor(vertex_locator _vertex, level_t _lower, level_t _upper)
    : vertex(_vertex),
      lower(_lower),
      upper(_upper) { }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
  {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //

    if (std::get<0>(alg_data).level[vertex][0] == kbfs_t::unvisited_level) return false; // skip unvisited vertices

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(ecc_bound_info_prop_visitor(eitr.target(),
                                                           std::get<1>(alg_data).lower[vertex],
                                                           std::get<1>(alg_data).upper[vertex]));
    }

    return true; // trigger bcast to queue visitors from delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const
  {
    // -------------------------------------------------- //
    // This function applies sent data to the vertex (apply step)
    // -------------------------------------------------- //
    std::get<2>(alg_data).lower[vertex] = std::max((int)lower, std::get<2>(alg_data).lower[vertex] - 1);
    std::get<2>(alg_data).upper[vertex] = std::min((int)upper, std::get<2>(alg_data).upper[vertex] + 1);
    std::get<3>(alg_data)[vertex] |= (lower == upper);

    // return true to call pre_visit() at the master of delegated vertices required by visitor_queue.hpp 274
    // (label:FLOW2)
    return vertex.is_delegate();
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const
  {
    // return true, to bcast bitmap sent by label:FLOW2
    if (!vertex.get_bcast()) {
      return true;
    }

    pre_visit(alg_data);

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(ecc_bound_info_prop_visitor(eitr.target(),
                                                           std::get<1>(alg_data).lower[vertex],
                                                           std::get<1>(alg_data).upper[vertex]));
    }
    return false;
  }

  friend inline bool operator>(const ecc_bound_info_prop_visitor &v1, const ecc_bound_info_prop_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const ecc_bound_info_prop_visitor &v1, const ecc_bound_info_prop_visitor &v2)
  {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator vertex;
  level_t lower;
  level_t upper;
} __attribute__ ((packed));


template <typename graph_t, int k_num_sources>
void find_max_ecc_bound_from_neighbor(graph_t *graph,
                                      typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                      typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  typedef ecc_bound_info_prop_visitor<graph_t, k_num_sources> visitor_type;

  typename eecc_type<graph_t, k_num_sources>::vertex_data tmp_ecc(*graph);
  typename graph_t::template vertex_data<bool, std::allocator<bool>> flag(*graph);

  auto alg_data = std::forward_as_tuple(kbfs_vertex_data, ecc_vertex_data, tmp_ecc, flag);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(graph, alg_data);
  // ------------------------------ Traversal ------------------------------ //
  {
    const double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); // init_visit -> queue
    MPI_Barrier(MPI_COMM_WORLD);
    const double time_end = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Collect bound info: " << time_end - time_start << "\t";
  }

  size_t num_bounded(0);
  for (auto itr = graph->vertices_begin(), end = graph->vertices_end(); itr != end; ++itr) {
    if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    num_bounded += ((tmp_ecc.lower[*itr] == tmp_ecc.upper[*itr]) &&
                    (ecc_vertex_data.lower[*itr] != ecc_vertex_data.upper[*itr]));
  }
  for (auto itr = graph->controller_begin(), end = graph->controller_end(); itr != end; ++itr) {
    if (kbfs_vertex_data.level[*itr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    num_bounded += ((tmp_ecc.lower[*itr] == tmp_ecc.upper[*itr]) &&
                    (ecc_vertex_data.lower[*itr] != ecc_vertex_data.upper[*itr]));
  }

  num_bounded = mpi_all_reduce(num_bounded, std::plus<level_t>(), MPI_COMM_WORLD);
  if (mpi_rank == 0) std::cout << "Bounded: " << num_bounded << std::endl;
}

// -------------------------------------------------------------------------------------------------------------- //
// compute_distance_score
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void compute_distance_histgram_helper(iterator_t vitr, iterator_t end,
                                      typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                      typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                      const size_t max_distance,
                                      std::vector<size_t> &histgram)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level) continue;
    const uint64_t x = ecc_vertex_data.upper[*vitr] - ecc_vertex_data.lower[*vitr];
    if (x > max_distance) ++histgram[max_distance + 1];
    else ++histgram[x];
  }
}

template <typename graph_t, int k_num_sources>
std::vector<size_t> compute_distance_histgram(graph_t *graph,
                                              typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                              typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                              const size_t max_distance)
{
  std::vector<size_t> histgram(max_distance + 2, 0);
  compute_distance_histgram_helper<graph_t, k_num_sources>(graph->vertices_begin(), graph->vertices_end(),
                                                           kbfs_vertex_data, ecc_vertex_data, max_distance, histgram);
  compute_distance_histgram_helper<graph_t, k_num_sources>(graph->controller_begin(), graph->controller_end(),
                                                           kbfs_vertex_data, ecc_vertex_data, max_distance, histgram);

  mpi_all_reduce_inplace(histgram, std::plus<size_t>(), MPI_COMM_WORLD);

  return histgram;
}

} // namespace havoqgt
#endif //HAVOQGT_EXACT_ECCENTRICITY_HPP
