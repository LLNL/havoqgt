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

#include <havoqgt/environment.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>


namespace havoqgt
{

struct eecc_source_select_mode_tag
{
  struct rnd {};
  struct far {};
  struct hdeg {};
  struct lvl2 {};
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
// find_farthest_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename graph_t, int k_num_sources, typename iterator_t>
void find_farthest_vertices_helper(const graph_t *const graph,
                                   iterator_t vitr, iterator_t end,
                                   const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                   const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                   std::vector<level_t> &level_list, std::vector<uint64_t> &vertex_id_list)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == kbfs_type<graph_t, k_num_sources>::unvisited_level)
      continue; // skip unvisited vertices
    if (ecc_vertex_data.lower[*vitr] == ecc_vertex_data.upper[*vitr]) continue;

    const auto min_itr = std::min_element(kbfs_vertex_data.level[*vitr].begin(), kbfs_vertex_data.level[*vitr].end());
    const level_t min_level = *min_itr;
    if (level_list.size() < k_num_sources) {
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
std::vector<uint64_t> find_farthest_vertices(const graph_t *const graph,
                                             const typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                             const typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data)
{
  std::vector<level_t> global_level_list;
  std::vector<uint64_t> global_vertex_id_list;

  {
    std::vector<level_t> local_level_list;
    std::vector<uint64_t> local_vertex_id_list;
    find_farthest_vertices_helper<graph_t, k_num_sources>(graph, graph->vertices_begin(), graph->vertices_end(),
                                                          kbfs_vertex_data, ecc_vertex_data,
                                                          local_level_list, local_vertex_id_list);
    find_farthest_vertices_helper<graph_t, k_num_sources>(graph, graph->controller_begin(), graph->controller_end(),
                                                          kbfs_vertex_data, ecc_vertex_data,
                                                          local_level_list, local_vertex_id_list);
    assert(local_level_list.size() == local_vertex_id_list.size());

    havoqgt::mpi_all_gather(local_level_list, global_level_list, MPI_COMM_WORLD);
    havoqgt::mpi_all_gather(local_vertex_id_list, global_vertex_id_list, MPI_COMM_WORLD);
  }

  const size_t num_total_candidates = global_level_list.size();
  std::vector<uint64_t> selected_vertex_id_list;
  if (num_total_candidates > k_num_sources) {
    std::vector<std::pair<level_t, uint64_t>> table;
    table.resize(num_total_candidates);
    for (size_t i = 0; i < global_level_list.size(); ++i) {
      table[i] = std::make_pair(global_level_list[i], global_vertex_id_list[i]);
    }
    std::partial_sort(table.begin(), table.begin() + k_num_sources, table.end(),
                      [&](const std::pair<level_t, uint64_t> &a,
                          const std::pair<level_t, uint64_t> &b) -> bool {
                        return (a.first > b.first);
                      });

    selected_vertex_id_list.resize(k_num_sources);
    for (size_t i = 0; i < k_num_sources; ++i) {
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
                                           eecc_source_select_mode_tag::rnd)
{
  std::cout << "randomly select sources" << std::endl;
  return detail::randomly_select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);

}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::far)
{
  std::cout << "select farthest sources" << std::endl;
  return detail::find_farthest_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);

}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::hdeg)
{
  std::cout << "select highest degree sources" << std::endl;
  return detail::find_highest_degree_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);

}

template <typename graph_t, int k_num_sources>
std::vector<uint64_t> select_source_helper(graph_t *graph,
                                           typename kbfs_type<graph_t, k_num_sources>::vertex_data &kbfs_vertex_data,
                                           typename eecc_type<graph_t, k_num_sources>::vertex_data &ecc_vertex_data,
                                           eecc_source_select_mode_tag::lvl2)
{
  std::cout << "select 2 level away sources" << std::endl;
  return detail::select_2_level_away_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);

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
  auto source_id_list = detail::select_source_helper<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data, tag);
  std::vector<typename graph_t::vertex_locator> source_locator_list;
  for (auto id : source_id_list) {
    source_locator_list.emplace_back(graph->label_to_locator(id));
  }
  return source_locator_list;
}

} // namespace havoqgt
#endif //HAVOQGT_EXACT_ECCENTRICITY_HPP
