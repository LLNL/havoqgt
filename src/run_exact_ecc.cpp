//
// Created by Iwabuchi, Keita on 12/14/17.
//

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

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <chrono>
#include <random>
#include <assert.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/fixed_size_unordered_map.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>

// -------------------------------------------------------------------------------------------------------------- //
// Types
// -------------------------------------------------------------------------------------------------------------- //
using namespace havoqgt;

using segment_manager_t = havoqgt::distributed_db::segment_manager_type;
using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;


// -------------------------------------------------------------------------------------------------------------- //
// Data structure
// -------------------------------------------------------------------------------------------------------------- //
struct ecc_vertex_data_t
{
  using vertex_data_ecc_t = typename graph_t::vertex_data<level_t, std::allocator<level_t>>;

  explicit ecc_vertex_data_t(const graph_t &graph)
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


// -------------------------------------------------------------------------------------------------------------- //
// Parse command line
// -------------------------------------------------------------------------------------------------------------- //
void usage()
{
  if (havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -i <string>    - input graph base filename (required)\n"
              << " -b <string>    - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
              << " -s <int>:<int> - Colon separated k source vertices (required)\n"
              << " -h             - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char **argv, std::string &input_filename,
                    std::string &backup_filename,
                    std::vector<uint64_t> &source_id_list)
{
  if (havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:b:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's': {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_id_list.push_back(std::stoull(buf.c_str()));
        break;
      }
      case 'i':
        found_input_filename = true;
        input_filename = optarg;
        break;
      case 'b':
        backup_filename = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_input_filename || source_id_list.size() < k_num_sources) {
    usage();
    exit(-1);
  }
}

// -------------------------------------------------------------------------------------------------------------- //
// find_highest_degree_vertices
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
void find_highest_degree_vertices(const graph_t *const graph, iterator_t vitr, iterator_t end,
                                  const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                                  const ecc_vertex_data_t &ecc_vertex_data,
                                  std::vector<size_t> &highest_degree_list,
                                  std::vector<uint64_t> &highest_degree_id_list)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == k_unvisited_level) continue; // skip unvisited vertices
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

std::vector<uint64_t> select_source_2(const graph_t *const graph,
                                      const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                                      const ecc_vertex_data_t &ecc_vertex_data)
{
  std::vector<size_t> global_highest_degree_list;
  std::vector<uint64_t> global_highest_degree_id_list;

  {
    std::vector<size_t> local_highest_degree_list;
    std::vector<uint64_t> local_highest_degree_id_list;
    find_highest_degree_vertices(graph, graph->vertices_begin(), graph->vertices_end(),
                                 kbfs_vertex_data, ecc_vertex_data,
                                 local_highest_degree_list, local_highest_degree_id_list);
    find_highest_degree_vertices(graph, graph->controller_begin(), graph->controller_end(),
                                 kbfs_vertex_data, ecc_vertex_data,
                                 local_highest_degree_list, local_highest_degree_id_list);
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
// select_non_zero_degree_source
// -------------------------------------------------------------------------------------------------------------- //
void select_non_zero_degree_source(const graph_t *const graph, std::vector<uint64_t> &source_id_list)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  for (auto id_itr = source_id_list.begin(), end = source_id_list.end(); id_itr != end; ++id_itr) {
    uint64_t source_id = *id_itr;
    graph_t::vertex_locator source = graph->label_to_locator(source_id);
    uint64_t global_degree(0);
    while (true) {
      uint64_t local_degree = 0;
      source = graph->label_to_locator(source_id);
      if (source.is_delegate()) {
        break;
      }
      if (uint32_t(mpi_rank) == source.owner()) {
        local_degree = graph->degree(source);
      }
      global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
      const bool already_exist = (std::find(source_id_list.begin(), id_itr, source_id) != id_itr);
      if (global_degree == 0 || already_exist) ++source_id;
      else break;
    }
    if (uint32_t(mpi_rank) == source.owner()) {
      if (source_id != *id_itr) {
        std::cout << "\nVertex " << *id_itr << " has a degree of 0.   New source vertex = " << source_id << std::endl;
      } else {
        std::cout << "\nStarting vertex = " << source_id << std::endl;
      }
      std::cout << "delegate? = " << source.is_delegate() << std::endl;
      std::cout << "local_id = " << source.local_id() << std::endl;
      std::cout << "degree = " << graph->degree(source) << std::endl;
    }
    *id_itr = source_id;
  }
}

// -------------------------------------------------------------------------------------------------------------- //
// select_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
void find_farthest_vertices(const graph_t *const graph,
                            iterator_t vitr, iterator_t end,
                            const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                            const ecc_vertex_data_t &ecc_vertex_data,
                            std::vector<level_t> &level_list, std::vector<uint64_t> &vertex_id_list)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == k_unvisited_level) continue; // skip unvisited vertices
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

std::vector<uint64_t> select_source(const graph_t *const graph,
                                    const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                                    const ecc_vertex_data_t &ecc_vertex_data)
{
  std::vector<level_t> global_level_list;
  std::vector<uint64_t> global_vertex_id_list;

  {
    std::vector<level_t> local_level_list;
    std::vector<uint64_t> local_vertex_id_list;
    find_farthest_vertices(graph, graph->vertices_begin(), graph->vertices_end(),
                           kbfs_vertex_data, ecc_vertex_data,
                           local_level_list, local_vertex_id_list);
    find_farthest_vertices(graph, graph->controller_begin(), graph->controller_end(),
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
                      [&](const std::pair<level_t, uint64_t> &a, const std::pair<level_t, uint64_t> &b) -> bool {
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
std::vector<uint64_t> randomly_select_source(graph_t *graph,
                                             const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                                             const ecc_vertex_data_t &ecc_vertex_data,
                                             const size_t num_solved_vertices)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector<uint64_t> local_vertex_id_list;
  std::vector<uint64_t> global_vertex_id_list;


  std::random_device rd; /// not implemented in icc ?
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<uint64_t> dis(0, graph->max_local_vertex_id() - 1);

  for (size_t i = 0; i < k_num_sources; ++i) {
    const uint64_t start_id = dis(gen);
    uint64_t vertex_id = start_id;
    while (true) {
      auto vertex = graph->label_to_locator(vertex_id);
      if (kbfs_vertex_data.level[vertex][0] != k_unvisited_level &&
          ecc_vertex_data.lower[vertex] != ecc_vertex_data.upper[vertex] &&
          std::find(local_vertex_id_list.begin(), local_vertex_id_list.end(), vertex_id) != local_vertex_id_list.end()) {
        local_vertex_id_list.push_back(vertex_id);
        break;
      }
      vertex_id = (vertex_id + 1) % graph->max_local_vertex_id();
      if (start_id == vertex_id) goto NOT_FOUND;
    }
  }

NOT_FOUND:

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
template <typename iterator_t>
void compute_ecc_k_source(iterator_t vitr, iterator_t end,
                          const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                          std::vector<level_t> &k_source_ecc)
{
  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == k_unvisited_level) continue;
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      assert(kbfs_vertex_data.level[*vitr][k] != k_unvisited_level); // vertices must be reachable from all k sources
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


std::unordered_map<std::pair<level_t, level_t>, size_t, pair_hash> bound_progress;

template <typename iterator_t>
std::pair<size_t, size_t> bound_ecc(const graph_t *const graph,
                                    iterator_t vitr,
                                    iterator_t end,
                                    const kbfs_vertex_data_t<graph_t> &kbfs_vertex_data,
                                    const std::vector<level_t> &k_source_ecc,
                                    ecc_vertex_data_t &ecc_vertex_data)
{
  size_t num_bounded_vertices(0);
  size_t num_non_bounded_vertices(0);

  for (; vitr != end; ++vitr) {
    if (kbfs_vertex_data.level[*vitr][0] == k_unvisited_level) continue; // skip unvisited vertices

    level_t &current_lower = ecc_vertex_data.lower[*vitr];
    level_t &current_upper = ecc_vertex_data.upper[*vitr];
    if (current_lower == current_upper) continue; // Exact ecc has been already found

    bool bounded = false;
    for (size_t k = 0; k < k_source_ecc.size(); ++k) {
      assert(kbfs_vertex_data.level[*vitr][k] != k_unvisited_level);

      const level_t level = kbfs_vertex_data.level[*vitr][k];

      current_lower = std::max(current_lower, std::max(level, static_cast<uint16_t>(k_source_ecc[k] - level)));
      current_upper = std::min(current_upper, static_cast<uint16_t>(k_source_ecc[k] + level));

      if (current_lower == current_upper) {
        bounded = true;
        break;
      } else {
        const auto key = std::make_pair(current_lower, current_upper);
        if (bound_progress.count(key) == 0) {
          bound_progress[key] = 0;
        }
        ++bound_progress[key];
      }
    }
    num_bounded_vertices += bounded;
    num_non_bounded_vertices += !bounded;
  }

  return std::make_pair(num_bounded_vertices, num_non_bounded_vertices);
}


// -------------------------------------------------------------------------------------------------------------- //
// find_max_ecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
level_t find_max_ecc(iterator_t vitr, iterator_t end, ecc_vertex_data_t &ecc_vertex_data)
{
  level_t max_ecc = std::numeric_limits<level_t>::min();
  for (; vitr != end; ++vitr) {
    max_ecc = std::max(max_ecc, ecc_vertex_data.lower[*vitr]);
  }
  return max_ecc;
}

// -------------------------------------------------------------------------------------------------------------- //
// find_max_lower_upper
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
std::pair<level_t, level_t> find_max_lower_upper(iterator_t vitr, iterator_t end, ecc_vertex_data_t &ecc_vertex_data)
{
  level_t max_lower = std::numeric_limits<level_t>::min();
  level_t max_upper = std::numeric_limits<level_t>::min();
  for (; vitr != end; ++vitr) {
    if (ecc_vertex_data.upper[*vitr] == std::numeric_limits<level_t>::max()) continue; // skip unvisited vertices
    max_lower = std::max(max_lower, ecc_vertex_data.lower[*vitr]);
    max_upper = std::max(max_upper, ecc_vertex_data.upper[*vitr]);
  }

  return std::make_pair(max_lower, max_upper);
}


// -------------------------------------------------------------------------------------------------------------- //
// Main
// -------------------------------------------------------------------------------------------------------------- //
int main(int argc, char **argv)
{

  int mpi_rank(0), mpi_size(0);


  havoqgt::havoqgt_init(&argc, &argv);
  {
    // -------------------------------------------------------------------------------------------------------------- //
    //                                            Parse options & Prepare graph
    // -------------------------------------------------------------------------------------------------------------- //
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    havoqgt::get_environment();

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      havoqgt::get_environment().print();
      std::cout << "k_num_sources " << k_num_sources << std::endl;
      //print_system_info(false);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_filename;
    std::vector<uint64_t> source_id_list;

    parse_cmd_line(argc, argv, graph_input, backup_filename, source_id_list);
    MPI_Barrier(MPI_COMM_WORLD);

    if (backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_t *graph = ddb.get_segment_manager()->find<graph_t>("graph_obj").first;
    assert(graph != nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }

    select_non_zero_degree_source(graph, source_id_list);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Selected initial vertices" << std::endl;
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Compute diameter
    //-------------------------------------------------------------------------------------------------------------- //
    const double total_start_time = MPI_Wtime();
    kbfs_vertex_data_t<graph_t> kbfs_vertex_data(*graph);
    ecc_vertex_data_t ecc_vertex_data(*graph);
    if (mpi_rank == 0) std::cout << "Allocated vertex data: " << std::endl;

    // ------------------------------ Init vertex data ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      ecc_vertex_data.init();
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) std::cout << "Init ecc vertex data: " << time_end - time_start << std::endl;
    }

    // ------------------------------ Compute exact ECC ------------------------------ //
    size_t current_iteration_number(0);
    while (true) {
      if (mpi_rank == 0)
        std::cout << "\n==================== " << current_iteration_number << " ====================" << std::endl;

      // ------------------------------ Init KBFS ------------------------------ //
      std::vector<typename graph_t::vertex_locator> source_locator_list;
      {
        const double time_start = MPI_Wtime();
        kbfs_vertex_data.init();
        for (uint64_t source_id : source_id_list) {
          source_locator_list.push_back(graph->label_to_locator(source_id));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Init for KBFS: " << time_end - time_start << std::endl;
          std::cout << "# KBFS sources: " << source_locator_list.size() << std::endl;
          for (uint64_t source_id : source_id_list) std::cout << source_id << " ";
          std::cout << std::endl;
        }
      }

      // ------------------------------ KBFS ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        havoqgt::k_breadth_first_search_level_per_source(graph, kbfs_vertex_data, source_locator_list);
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();

        if (mpi_rank == 0) std::cout << "BFS Time: " << time_end - time_start << std::endl;
      }

      // ------------------------------ Compute ECC for k sources ------------------------------ //
      std::vector<level_t> k_source_ecc(source_id_list.size(), 0);
      size_t total_solved_vertices = 0;
      {
        const double time_start = MPI_Wtime();

        compute_ecc_k_source(graph->vertices_begin(), graph->vertices_end(), kbfs_vertex_data, k_source_ecc);
        compute_ecc_k_source(graph->controller_begin(), graph->controller_end(), kbfs_vertex_data, k_source_ecc);
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
          total_solved_vertices += k_source_ecc.size();
        }
      }

      // ------------------------------ Bound ECC ------------------------------ //
      size_t num_bounded_vertices(0);
      size_t num_non_bounded_vertices(0);
      bound_progress.clear();
      {
        const double time_start = MPI_Wtime();
        {
          auto ret = bound_ecc(graph,
                               graph->vertices_begin(),
                               graph->vertices_end(),
                               kbfs_vertex_data,
                               k_source_ecc,
                               ecc_vertex_data);
          num_bounded_vertices += ret.first;
          num_non_bounded_vertices += ret.second;
        }

        {
          auto ret = bound_ecc(graph,
                               graph->controller_begin(),
                               graph->controller_end(),
                               kbfs_vertex_data,
                               k_source_ecc,
                               ecc_vertex_data);
          num_bounded_vertices += ret.first;
          num_non_bounded_vertices += ret.second;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        size_t global_num_bounded_vertices(0);
        CHK_MPI(MPI_Reduce(&num_bounded_vertices,
                           &global_num_bounded_vertices,
                           1,
                           mpi_typeof(num_bounded_vertices),
                           MPI_SUM,
                           0,
                           MPI_COMM_WORLD));

        size_t global_num_non_bounded_vertices(0);
        CHK_MPI(MPI_Reduce(&num_non_bounded_vertices,
                           &global_num_non_bounded_vertices,
                           1,
                           mpi_typeof(num_non_bounded_vertices),
                           MPI_SUM,
                           0,
                           MPI_COMM_WORLD));

        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Update ECC bounds: " << time_end - time_start << std::endl;
          std::cout << "# bounded vertices: " << global_num_bounded_vertices << std::endl;
          std::cout << "# non-bounded vertices: " << global_num_non_bounded_vertices << std::endl;
          std::cout << "L\tU\tCount" << std::endl;
          for (auto itr : bound_progress) {
            std::cout << itr.first.first << " " << itr.first.second << " " << itr.second << std::endl;
          }
          total_solved_vertices += global_num_bounded_vertices;
        }
      }

      // ------------------------------ Check termination ------------------------------ //
      {
        const size_t global = mpi_all_reduce(num_non_bounded_vertices, std::logical_or<size_t>(), MPI_COMM_WORLD);
        if (!global) break;
        MPI_Barrier(MPI_COMM_WORLD); // Just in case
      }

      // ------------------------------ Select sources for next step ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        source_id_list = select_source(graph, kbfs_vertex_data, ecc_vertex_data);
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Select sources: " << time_end - time_start << std::endl;
        }
      }

      // ------------------------------ Diameter calculation termination test ------------------------------ //
      {
        const auto ret1 = find_max_lower_upper(graph->vertices_begin(), graph->vertices_end(), ecc_vertex_data);
        const auto ret2 = find_max_lower_upper(graph->controller_begin(), graph->controller_end(), ecc_vertex_data);

        level_t max_lower = std::max(ret1.first, ret2.first);
        level_t max_upper = std::max(ret1.second, ret2.second);

        max_lower = mpi_all_reduce(max_lower, std::greater<level_t>(), MPI_COMM_WORLD);
        max_upper = mpi_all_reduce(max_upper, std::greater<level_t>(), MPI_COMM_WORLD);

        if (mpi_rank == 0) std::cout << "terminate test " << max_lower << " : " << max_upper << std::endl;
      }
      ++current_iteration_number;
    } // End Exact ECC

    const double total_end_time = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Total execution time: " << total_end_time - total_start_time << std::endl;

    // ------------------------------ Compute diameter ------------------------------ //
    {
      // const double time_start = MPI_Wtime();
      const level_t max_ecc = std::max(find_max_ecc(graph->vertices_begin(), graph->vertices_end(), ecc_vertex_data),
                                       find_max_ecc(graph->controller_begin(),
                                                    graph->controller_end(),
                                                    ecc_vertex_data));
      level_t diameter;
      CHK_MPI(MPI_Reduce(&max_ecc, &diameter, 1, mpi_typeof(max_ecc), MPI_MAX, 0, MPI_COMM_WORLD));
      // const double time_end = MPI_Wtime();

      if (mpi_rank == 0) {
        // std::cout << "\nCompute diameter: " << time_end - time_start << std::endl;
        std::cout << "Diameter: " << diameter << std::endl;
      }
    } // End Compute diameter
  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}
