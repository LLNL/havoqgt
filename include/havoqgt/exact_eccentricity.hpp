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

#include <boost/container/flat_map.hpp>

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>
#include <havoqgt/kth_core_new.hpp>

#define TREE_FROM_ONLY_SOURCES 0
bool use_tk = false;
//bool use_new_max_u;
//bool use_skip_strategy;
bool use_soft_contribution_score = false;
bool use_ds_fix = false;
bool use_ds_adp = false;
bool use_hanging_tree = false;
bool use_ds_fix_middle_shel = false;

const int k_strategy_id_degree = 99;
const int k_strategy_id_ds = 999;

namespace havoqgt {
template <typename graph_allocator_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity_vertex_data {
 private:
  using graph_t = havoqgt::delegate_partitioned_graph<graph_allocator_t>;
  using ecc_t = typename graph_t::template vertex_data<level_t, std::allocator<level_t>>;
  using flag_t = typename graph_t::template vertex_data<uint8_t, std::allocator<uint8_t>>;
  using bool_t = typename graph_t::template vertex_data<bool, std::allocator<bool>>;
  using count_t = typename graph_t::template vertex_data<uint32_t, std::allocator<uint32_t>>;

 public:
  /// solved flags ///
  static constexpr uint8_t k_source = 1 << 0;
  static constexpr uint8_t k_bound = 1 << 1;
  static constexpr uint8_t k_tree = 1 << 2;

  explicit exact_eccentricity_vertex_data(const graph_t &graph)
      : m_graph(graph),
        m_lower(graph),
        m_upper(graph),
        m_just_solved_status(graph),
        m_num_unsolved_neighbors(graph),
        m_farthest_vertex(graph),
        m_ever_source_vertex(graph) {
#ifdef DEBUG
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
#endif
  }

  void init() {
    m_lower.reset(std::numeric_limits<level_t>::min());
    m_upper.reset(std::numeric_limits<level_t>::max());
    m_farthest_vertex.reset(false);
    m_ever_source_vertex.reset(false);
    reset_for_each_iteration();
  }

  void reset_for_each_iteration() {
    m_just_solved_status.reset(0);
    m_num_unsolved_neighbors.reset(0);
  }

  level_t &lower(const typename graph_t::vertex_locator vertex) {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return m_lower[vertex];
  }

  level_t &upper(const typename graph_t::vertex_locator vertex) {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return m_upper[vertex];
  }

  uint8_t &just_solved_status(const typename graph_t::vertex_locator vertex) {
#ifdef DEBUG
    assert(vertex.owner() == static_cast<uint32_t>(m_mpi_rank) || vertex.is_delegate());
#endif
    return m_just_solved_status[vertex];
  }

  uint32_t &num_unsolved_neighbors(const typename graph_t::vertex_locator vertex) {
    return m_num_unsolved_neighbors[vertex];
  }

  bool &farthest_vertex(const typename graph_t::vertex_locator vertex) {
    return m_farthest_vertex[vertex];
  }

  bool &ever_source_vertex(const typename graph_t::vertex_locator vertex) {
    return m_ever_source_vertex[vertex];
  }

 private:
  const graph_t &m_graph;
  ecc_t m_lower;
  ecc_t m_upper;
  flag_t m_just_solved_status;
  count_t m_num_unsolved_neighbors;
  bool_t m_farthest_vertex;
  bool_t m_ever_source_vertex;

#ifdef DEBUG
  int m_mpi_rank;
#endif
};

template <typename graph_allocator_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity {
 private:
  using graph_t = havoqgt::delegate_partitioned_graph<graph_allocator_t>;
  using vertex_locator_t = typename graph_t::vertex_locator;
  using source_score_function_t = std::function<uint64_t(const vertex_locator_t)>;
  using source_score_function_list_t = std::vector<source_score_function_t>;

  class hanging_tree_visitor;
  class take_pruning_visitor;
  class unsolved_visitor;
//  class articulation_visitor;
//  class single_source_visitor;

 public:
  using kbfs_t = k_breadth_first_search<graph_allocator_t, level_t, k_num_sources>;
  using kbfs_vertex_data_t = typename kbfs_t::vertex_data_t;
  using ecc_vertex_data_t = exact_eccentricity_vertex_data<graph_allocator_t, level_t, k_num_sources>;
  using kth_core_vertex_data_t = typename graph_t::template vertex_data<kth_core_data, std::allocator<kth_core_data>>;

  explicit exact_eccentricity(graph_t &graph, const std::set<int> &use_algorithm)
      : m_graph(graph),
        m_kbfs(graph),
        m_ecc_vertex_data(graph),
        m_source_info(),
        m_source_score_function_list(),
        m_progress_info(),
        m_2core_vertex_data(graph) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    use_tk = !!std::getenv("USE_TAKE");
    use_soft_contribution_score = !!std::getenv("USE_SOFT_CONT_SCORE");
    use_ds_fix = !!std::getenv("USE_DS_FIX");
    use_ds_adp = !!std::getenv("USE_DS_ADP");
    use_hanging_tree = !!std::getenv("USE_TREE");
    use_ds_fix_middle_shel = !!std::getenv("USE_DS_FIX_MS");
    if (mpi_rank == 0) {
      std::cout << "use_tk: " << use_tk << std::endl;
      std::cout << "use_soft_contribution_score: " << use_soft_contribution_score << std::endl;
      std::cout << "use_ds_fix  " << use_ds_fix << std::endl;
      std::cout << "use_ds_fix_middle_shel  " << use_ds_fix_middle_shel << std::endl;
      std::cout << "use_ds_adp  " << use_ds_adp << std::endl;
      std::cout << "use_hanging_tree  " << use_hanging_tree << std::endl;
    }

    if (use_ds_adp) set_strategy(use_algorithm);
  }

  void load_2core_info(const std::string path) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    if (!use_hanging_tree) {
      if (mpi_rank == 0)
        std::cout << "!!!! skip " << __FUNCTION__ << std::endl;
      return;
    }

    std::ifstream ifs(path);
    if (!ifs.is_open()) {
      std::cerr << "Can not open " << path << std::endl;
      std::abort();
    }
    uint64_t vid;
    bool alive;
    size_t height;
    size_t num_cut;
    uint64_t parent_id;
    while (ifs >> vid >> alive >> height >> num_cut >> parent_id) {
      auto locator = m_graph.label_to_locator(vid);
      if (locator.owner() == mpi_rank || locator.is_delegate()) {
        m_2core_vertex_data[locator].set_alive(alive);
        m_2core_vertex_data[locator].set_height(height);
        m_2core_vertex_data[locator].set_num_cut(num_cut);
        m_2core_vertex_data[locator].set_parent_id(parent_id);
      }
    }
  }

  // -------------------------------------------------------------------------------------------------------------- //
  // run exact ecc algorithm
  // -------------------------------------------------------------------------------------------------------------- //
  void run() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    m_ecc_vertex_data.init();

    while (true) {
      if (mpi_rank == 0) std::cout << "========== " << m_progress_info.iteration_no << " ==========" << std::endl;

      // -------------------- Select source -------------------- //
      {
        source_info_t new_source_info;
        const double time_start = MPI_Wtime();
        if (use_tk) {
          if (m_progress_info.iteration_no == 0) {
            new_source_info = select_initial_source_by_take_algorithm();
          } else {
            new_source_info = select_source_by_take_algorithm();
          }
        } else if (use_ds_fix) {
          if (m_progress_info.iteration_no == 0) {
            new_source_info = select_initial_source_by_degree();
          } else if (m_progress_info.num_unsolved <= k_num_sources) {
            new_source_info = select_left_sources(m_progress_info.num_unsolved);
          } else {
            new_source_info = select_source_by_double_sweep_fixed();
          }
          if (m_progress_info.num_unsolved > 0 && new_source_info.num_source() == 0) {
            new_source_info = select_left_sources(std::min(m_progress_info.num_unsolved, (size_t)k_num_sources));
          }
        } else if (use_ds_adp) {
          if (m_progress_info.iteration_no == 0) {
            new_source_info = select_initial_source_by_degree();
          } else if (m_progress_info.num_unsolved <= k_num_sources) {
            new_source_info = select_left_sources(m_progress_info.num_unsolved);
          } else {
            new_source_info = select_source_by_double_sweep_adptive();
          }
          if (m_progress_info.num_unsolved > 0 && new_source_info.num_source() == 0) {
            new_source_info = select_left_sources(std::min(m_progress_info.num_unsolved, (size_t)k_num_sources));
          }
        } else {
          if (mpi_comm_rank() == 0)
            std::cerr << "Invalid source selection algorithm" << std::endl;
          std::abort();
        }
        m_source_info = std::move(new_source_info);

        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Select sources took: " << time_end - time_start << std::endl;
          std::cout << "#souces: " << m_source_info.num_source() << std::endl;

          std::cout << "ID: ";
          for (auto v : m_source_info.source_list) std::cout << m_graph.locator_to_label(v) << " ";
          std::cout << std::endl;

          std::cout << "Strategy: ";
          for (auto s : m_source_info.strategy_list) std::cout << s << " ";
          std::cout << std::endl;

          std::cout << "Owner: ";
          for (auto v : m_source_info.source_list) std::cout << v.owner() << " ";
          std::cout << std::endl;
        }

        for (auto src : m_source_info.source_list) {
          if (src.owner() == mpi_comm_rank())
            m_ecc_vertex_data.ever_source_vertex(src) = true;
        }
      }

      // -------------------- Run KBFS -------------------- //
      {
        const double time_start = MPI_Wtime();
        m_kbfs.run(m_source_info.source_list);
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_comm_rank() == 0) {
          std::cout << "kBFS_took: " << time_end - time_start << std::endl;
        }
      }

      m_ecc_vertex_data.reset_for_each_iteration();

      // -------------------- Bounding algorithm -------------------- //
      {
        const double time_start = MPI_Wtime();
        compute_ecc_k_source();
        std::tie(m_progress_info.num_solved, m_progress_info.num_unsolved) = bound_ecc();
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "------------------------------------------------------" << std::endl;
          std::cout << "Bounding_algorithm_took:\t" << time_end - time_start << std::endl;
          std::cout << "#solved:\t" << m_progress_info.num_solved << std::endl;
          std::cout << "------------------------------------------------------" << std::endl;
        }
      }
      if (m_progress_info.num_unsolved == 0) return;

      // -------------------- Pruning -------------------- //
      if (std::getenv("USE_TAKE")) {
        const double time_start = MPI_Wtime();
        const size_t num_solved = solve_single_degree_vertices();
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Solve single degree vertices: " << time_end - time_start << std::endl;
          std::cout << "#solved_single_vertices: " << num_solved << std::endl;
          std::cout << "------------------------------------------------------" << std::endl;
        }

        m_progress_info.num_solved += num_solved;
        m_progress_info.num_unsolved -= num_solved;
      } else if (use_hanging_tree) {
        {
          const double time_start = MPI_Wtime();
          const size_t num_solved = solve_hanging_tree();
          const double time_end = MPI_Wtime();
          if (mpi_rank == 0) {
            std::cout << "Solve_hanging_tree_took:\t" << time_end - time_start << std::endl;
            std::cout << "#solved_by_tree:\t" << num_solved << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;
          }
          m_progress_info.num_solved += num_solved;
          m_progress_info.num_unsolved -= num_solved;
        }
        {
          const double time_start = MPI_Wtime();
          size_t num_solved;
          size_t num_unsolved;
          std::tie(num_solved, num_unsolved) = bound_ecc_for_leaf();
          MPI_Barrier(MPI_COMM_WORLD);
          const double time_end = MPI_Wtime();
          if (mpi_rank == 0) {
            std::cout << "Bounding_for_leaf_took:\t" << time_end - time_start << std::endl;
            std::cout << "#solved_by_leaf:\t  " << num_solved << std::endl;
            std::cout << "------------------------------------------------------" << std::endl;
            // std::cout << "#unsolved_by_leaf:\t" << num_unsolved << std::endl;
          }
          m_progress_info.num_solved += num_solved;
          m_progress_info.num_unsolved -= num_solved;
        }
      }

      if (mpi_rank == 0) {
        std::cout << "#unsolved:\t" << m_progress_info.num_unsolved << std::endl;
        std::cout << "------------------------------------------------------" << std::endl;
      }

      if (m_progress_info.num_unsolved == 0) return;

      // -------------------- unsolved neighbors -------------------- //
      if (use_ds_adp) {
        const double time_start = MPI_Wtime();
        count_unsolved_neighbors();
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Count unsolved neighbors took: " << time_end - time_start << std::endl;
          std::cout << "------------------------------------------------------" << std::endl;
        }
      }

      // -------------------- Progress report -------------------- //
      {
        auto histgram = compute_diff_score_histgram(4);
        if (mpi_rank == 0) {
          std::cout << "Diff score: ";
          for (auto n : histgram) std::cout << n << " ";
          std::cout << std::endl;
        }
        collect_unsolved_vertices_statistics();
      }
      ++m_progress_info.iteration_no;
    }
  }

// -------------------------------------------------------------------------------------------------------------- //
// find max ecc
// -------------------------------------------------------------------------------------------------------------- //
  uint16_t
  find_max_ecc() {
    uint16_t max_ecc = 0;
    for (auto vitr = m_graph.vertices_begin(), end = m_graph.vertices_end(); vitr != end; ++vitr)
      max_ecc = std::max(max_ecc, m_ecc_vertex_data.lower(*vitr));
    for (auto vitr = m_graph.controller_begin(), end = m_graph.controller_end(); vitr != end; ++vitr)
      max_ecc = std::max(max_ecc, m_ecc_vertex_data.lower(*vitr));

    max_ecc = mpi_all_reduce(max_ecc, std::greater<uint16_t>(), MPI_COMM_WORLD);

    return max_ecc;
  }

// -------------------------------------------------------------------------------------------------------------- //
// dump ecc
// -------------------------------------------------------------------------------------------------------------- //
  void dump_ecc(const std::string &file_name) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    std::ofstream ofs(file_name + "-" + std::to_string(mpi_rank));

    if (!ofs.is_open()) {
      std::cerr << "Can not open " << file_name << std::endl;
      return;
    }

    for (auto vitr = m_graph.vertices_begin(), end = m_graph.vertices_end(); vitr != end; ++vitr) {
      if (m_kbfs.vertex_data().visited_by(*vitr, 0))
        ofs << m_graph.locator_to_label(*vitr) << " " << m_ecc_vertex_data.lower(*vitr) << "\n";
    }

    for (auto vitr = m_graph.controller_begin(), end = m_graph.controller_end(); vitr != end; ++vitr) {
      if (m_kbfs.vertex_data().visited_by(*vitr, 0))
        ofs << m_graph.locator_to_label(*vitr) << " " << m_ecc_vertex_data.lower(*vitr) << "\n";
    }

    ofs.close();
  }

 private:

// -------------------------------------------------------------------------------------------------------------- //
// Data structures used in internal
// -------------------------------------------------------------------------------------------------------------- //
  class source_info_t {

   public:
    bool uniquely_add_source(vertex_locator_t source, int strategy) {
      auto ret = std::find(source_list.begin(), source_list.end(), source);
      if (ret != source_list.end()) return false;

      source_list.emplace_back(source);
      num_solved_list.emplace_back(0);
      strategy_list.emplace_back(strategy);
      ecc_list.emplace_back(std::numeric_limits<level_t>::min());

      return true;
    }

    size_t num_source() const {
#ifdef DEBUG
      assert(source_list.size() == num_solved_list.size()
                 && source_list.size() == strategy_list.size()
                 && source_list.size() == ecc_list.size());
#endif
      return source_list.size();
    }

    std::vector<vertex_locator_t> source_list;
    std::vector<size_t> num_solved_list;
    std::vector<int> strategy_list;
    std::vector<level_t> ecc_list;
  };

  struct progress_info_t {
    size_t iteration_no{0}; // Iteration No
    size_t num_solved{0}; // num solved veritces at the level; not total value
    size_t num_unsolved{0}; // num unsolved vertices
    size_t inner_iteration_cnt{0};
  };

// -------------------------------------------------------------------------------------------------------------- //
// source selection score functions
// -------------------------------------------------------------------------------------------------------------- //
  void set_strategy(const std::set<int> &use_algorithm) {
    if (use_algorithm.count(0))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return min_lower(vertex);
      });
    if (use_algorithm.count(1))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return degree_score(vertex);
      });
    if (use_algorithm.count(2))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return diff_score(vertex);
      });
    if (use_algorithm.count(3))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return level2_score(vertex);
      });
    if (use_algorithm.count(4))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return unsolved_neighbors_score(vertex);
      });
    if (use_algorithm.count(5))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return random_score(vertex);
      });
    if (use_algorithm.count(6))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return middle_shell_score(vertex);
      });
    if (use_algorithm.count(7))
      m_source_score_function_list.emplace_back([this](const vertex_locator_t vertex) -> uint64_t {
        return articulation_score(vertex);
      });
  }

  uint32_t hash_vertex_id(const uint64_t vid) {
    return static_cast<uint32_t>(detail::hash_nbits(vid, 32));
  }

  bool compare_vertex_by_random_hash(const uint64_t &vid_lhd, const uint64_t &vid_rhd) {
    const uint64_t h_lhd = hash_vertex_id(vid_lhd);
    const uint64_t h_rhd = hash_vertex_id(vid_rhd);
    return (h_lhd > h_rhd);
  }

  bool compare_vertex_by_random_hash(const vertex_locator_t &lhd, const vertex_locator_t &rhd) {
    const uint64_t vid_lhd = m_graph.locator_to_label(lhd);
    const uint64_t vid_rhd = m_graph.locator_to_label(rhd);
    return compare_vertex_by_random_hash(vid_lhd, vid_rhd); // To minimize replacement, smaller IDs have priority
  }

  uint64_t degree_score(const vertex_locator_t vertex) {
    return m_graph.degree(vertex);
  }

  uint64_t sum_level(const vertex_locator_t vertex) {
    uint64_t total(0);
    for (size_t k = 0; k < m_source_info.num_source(); ++k) {
      total += m_kbfs.vertex_data().level(vertex)[k];
    }
    return total;
  }

  uint64_t diff_score(const vertex_locator_t vertex) {
    return m_ecc_vertex_data.upper(vertex) - m_ecc_vertex_data.lower(vertex);
  }

  uint64_t level2_score(const vertex_locator_t vertex) {
    uint64_t total(0);
    for (size_t k = 0; k < m_source_info.num_source(); ++k) {
      total += (m_kbfs.vertex_data().level(vertex)[k] == 2);
    }
    return total;
  }

  uint64_t min_lower(const vertex_locator_t vertex) {
    // Note: lower numbers need to get higher priorities
    return std::numeric_limits<level_t>::max() - m_ecc_vertex_data.lower(vertex);
  }

  uint64_t max_upper(const vertex_locator_t vertex) {
    return m_ecc_vertex_data.upper(vertex);
  }

  uint64_t random_score(const vertex_locator_t vertex) {
    return hash_vertex_id(m_graph.locator_to_label(vertex));
  }

  uint64_t unsolved_neighbors_score(const vertex_locator_t vertex) {
    return m_ecc_vertex_data.num_unsolved_neighbors(vertex);
  }

  uint64_t articulation_score(const vertex_locator_t vertex) {
    return m_2core_vertex_data[vertex].get_num_cut();
  }

  uint64_t middle_shell_score(const vertex_locator_t vertex) {
    uint64_t total(0);
    for (size_t k = 0; k < m_source_info.num_source(); ++k) {
      if (m_kbfs.vertex_data().level(vertex)[k] <= m_source_info.ecc_list[k] / 2) {
        total += m_graph.degree(vertex) * m_kbfs.vertex_data().level(vertex)[k];
      }
    }
    return total;
  }

// -------------------------------------------------------------------------------------------------------------- //
// select source
// -------------------------------------------------------------------------------------------------------------- //
// ---------------------------------------- Takes's source selection algorithm ---------------------------------------- //
  source_info_t select_initial_source_by_take_algorithm() {
    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return (m_graph.degree(vertex) >= 1);
    };

    return select_source_by_take_algorithm(is_candidate);
  }

  source_info_t select_source_by_take_algorithm() {
    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return m_kbfs.vertex_data().visited_by(vertex, 0)
          && (m_ecc_vertex_data.lower(vertex) != m_ecc_vertex_data.upper(vertex));
    };
    return select_source_by_take_algorithm(is_candidate);
  };

  source_info_t select_source_by_take_algorithm(const std::function<bool(const vertex_locator_t)> &is_candidate) {

    std::vector<std::function<uint64_t(const vertex_locator_t)>> source_selection_strategy(2);

    if (m_progress_info.iteration_no % 2 == 0)
      source_selection_strategy[0] = [this](const vertex_locator_t vertex) -> uint64_t {
        return max_upper(vertex);
      };
    else
      source_selection_strategy[0] = [this](const vertex_locator_t vertex) -> uint64_t {
        return min_lower(vertex);
      };

    source_selection_strategy[1] = [this](const vertex_locator_t vertex) -> uint64_t {
      return degree_score(vertex);
    };

    auto source_candidate_list = select_source(k_num_sources, is_candidate, source_selection_strategy);
    source_info_t new_source_info;
    for (auto candidate : source_candidate_list) {
      new_source_info.uniquely_add_source(candidate, m_progress_info.iteration_no % 2);
    }
    return new_source_info;
  }

// ---------------------------------------- select_left_sources ---------------------------------------- //
  source_info_t select_left_sources(const size_t max_num_to_select) {
    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return m_kbfs.vertex_data().visited_by(vertex, 0)
          && (m_ecc_vertex_data.lower(vertex) != m_ecc_vertex_data.upper(vertex));
    };

    auto source_candidate_list = select_source(max_num_to_select,
                                               is_candidate,
                                               {[this](const vertex_locator_t vertex) -> uint64_t {
                                                 return degree_score(vertex);
                                               }});

    source_info_t new_source_info;
    for (auto candidate : source_candidate_list) {
      new_source_info.uniquely_add_source(candidate, k_strategy_id_degree);
    }

    return new_source_info;
  }

// ---------------------------------------- Initial source selection algorithm ---------------------------------------- //
  source_info_t select_initial_source_by_degree() {

    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return (m_graph.degree(vertex) >= 1);
    };

    auto source_candidate_list = select_source(k_num_sources,
                                               is_candidate,
                                               {
                                                   [this](const vertex_locator_t vertex) -> uint64_t {
                                                     return degree_score(vertex);
                                                   }
                                               });
    source_info_t new_source_info;
    for (auto candidate : source_candidate_list) {
      new_source_info.uniquely_add_source(candidate, k_strategy_id_degree);
    }
    return new_source_info;
  }

// ---------------------------------------- double sweep fixed source selection algorithhm ---------------------------------------- //
  source_info_t select_source_by_double_sweep_fixed() {
    // ---------- Select farthest vertices ---------- //
    const auto farthest_vertices = select_and_mark_farthest_vertices();
    source_info_t new_source_info;
    for (auto candidate : farthest_vertices) {
      new_source_info.uniquely_add_source(candidate, k_strategy_id_ds);
    }
    if (m_progress_info.num_unsolved < new_source_info.num_source()) {
      if (mpi_comm_rank() == 0) std::cerr << __func__ << "Logic error" << std::endl;
    }
    if (m_progress_info.num_unsolved == new_source_info.num_source() || k_num_sources == new_source_info.num_source()) {
      return new_source_info;
    }


    // ---------- Find rest of sources by the fixed strategy ---------- //
    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return m_kbfs.vertex_data().visited_by(vertex, 0)
          && (m_ecc_vertex_data.lower(vertex) != m_ecc_vertex_data.upper(vertex))
          && !m_ecc_vertex_data.farthest_vertex(vertex);
    };

    using func_list_t = std::vector<std::function<uint64_t(const vertex_locator_t)>>;
    func_list_t score_function_list;
    if (use_ds_fix_middle_shel) {
      score_function_list = func_list_t({
                                            [this](const vertex_locator_t vertex) -> uint64_t {
                                              return middle_shell_score(vertex);
                                            }});
    } else {
      score_function_list = func_list_t({
                                            [this](const vertex_locator_t vertex) -> uint64_t {
                                              return diff_score(vertex);
                                            },
                                            [this](const vertex_locator_t vertex) -> uint64_t {
                                              return min_lower(vertex);
                                            },
                                            [this](const vertex_locator_t vertex) -> uint64_t {
                                              return degree_score(vertex);
                                            }});
    }

    auto source_candidate_list = select_source(k_num_sources - new_source_info.num_source(),
                                               is_candidate,
                                               score_function_list);

    for (auto candidate : source_candidate_list) {
      new_source_info.uniquely_add_source(candidate, 0);
    }

    return new_source_info;
  }

// ---------------------------------------- For double sweep adaptive source selection algorithm ---------------------------------------- //
  source_info_t select_source_by_double_sweep_adptive() {
    // ---------- Select farthest vertices ---------- //
    const auto farthest_vertices = select_and_mark_farthest_vertices();
    source_info_t new_source_info;
    for (auto candidate : farthest_vertices) {
      new_source_info.uniquely_add_source(candidate, k_strategy_id_ds);
    }
    if (m_progress_info.num_unsolved < new_source_info.num_source()) {
      if (mpi_comm_rank() == 0) std::cerr << __func__ << "Logic error" << std::endl;
    }
    if (new_source_info.num_source() == m_progress_info.num_unsolved || new_source_info.num_source() == k_num_sources) {
      return new_source_info;
    }

    // ---------- Find rest of sources by adaptive method ---------- //
    const auto is_candidate = [this](const vertex_locator_t &vertex) -> bool {
      return m_kbfs.vertex_data().visited_by(vertex, 0)
          && (m_ecc_vertex_data.lower(vertex) != m_ecc_vertex_data.upper(vertex))
          && !m_ecc_vertex_data.farthest_vertex(vertex);
    };

    // ---------- Set contribution score ---------- //
    std::vector<size_t> strategy_contribution_score(m_source_score_function_list.size(), 0);
    if (m_progress_info.iteration_no > 1) {
      for (uint32_t i = 0; i < m_source_info.num_source(); ++i) {
        if (m_source_info.strategy_list[i] != k_strategy_id_degree
            && m_source_info.strategy_list[i] != k_strategy_id_ds)
          strategy_contribution_score[m_source_info.strategy_list[i]] += m_source_info.num_solved_list[i];
      }
      if (mpi_comm_rank() == 0) {
        std::cout << "Contribution score: ";
        for (auto n : strategy_contribution_score) std::cout << n << " ";
        std::cout << std::endl;
      }
    }

    select_source_by_contribution_score(is_candidate, strategy_contribution_score, new_source_info);

    return new_source_info;
  }

  void select_source_by_contribution_score(const std::function<bool(const vertex_locator_t)> &is_candidate,
                                           const std::vector<size_t> &strategy_contribution_score,
                                           source_info_t &new_source_info) {

    // Use a single strategy to select the remaining unsolved sources
    if (new_source_info.num_source() + m_progress_info.num_unsolved <= k_num_sources) {

      std::vector<vertex_locator_t> source_candidate_list = select_source(m_progress_info.num_unsolved, is_candidate,
                                                                          {
                                                                              [this](const vertex_locator_t vertex) -> uint64_t {
                                                                                return degree_score(vertex);
                                                                              }

                                                                          });
      // ----- Merge sources ----- //
      for (auto candidate : source_candidate_list) {
        new_source_info.uniquely_add_source(candidate, k_strategy_id_degree);
      }
      return;
    }

    size_t num_to_generate = k_num_sources - new_source_info.num_source();

    // ---------- Compute how many sources to be selected by each strategy using discrete_distribution ---------- //
    std::vector<size_t> wk_strategy_contribution_score(strategy_contribution_score);
    for (auto &n : wk_strategy_contribution_score) n += 1; // To avoid the case where all contirubtion scores are 0
    std::discrete_distribution<uint32_t> distribution(wk_strategy_contribution_score.begin(),
                                                      wk_strategy_contribution_score.end());

    std::vector<uint32_t> num_to_generate_by_strategy;
    if (num_to_generate > m_source_score_function_list.size()) {
      // use all strategies at least one time, initial value is 1
      num_to_generate_by_strategy.resize(m_source_score_function_list.size(), 1);
      num_to_generate -= m_source_score_function_list.size();
    } else {
      num_to_generate_by_strategy.resize(m_source_score_function_list.size(), 0);
    }

    std::mt19937 rnd(m_progress_info.num_solved); // seed can be any number but must be same among the all processes
    for (uint32_t i = 0; i < num_to_generate; ++i) {
      const uint32_t strategy_id = distribution(rnd);
      ++num_to_generate_by_strategy[strategy_id];
    }

    // ---------- Select sources ---------- //
    select_source_with_multiple_strategy(num_to_generate_by_strategy, is_candidate, new_source_info);
  }

  void select_source_with_multiple_strategy(const std::vector<uint32_t> num_to_generate_by_strategy,
                                            const std::function<bool(const vertex_locator_t)> &is_candidate,
                                            source_info_t &new_source_info) {

    for (uint32_t strategy_id = 0; strategy_id < num_to_generate_by_strategy.size(); ++strategy_id) {
      if (num_to_generate_by_strategy[strategy_id] == 0) continue;

      const size_t new_total_num_sources = new_source_info.num_source() + num_to_generate_by_strategy[strategy_id];
      const size_t num_to_generate
          = new_total_num_sources; // To get enough non-duplicated vertices, select more sources than actually needed
      std::vector<vertex_locator_t> source_candidate_list = select_source(num_to_generate, is_candidate,
                                                                          {m_source_score_function_list[strategy_id]});
      // ----- Merge sources ----- //
      for (auto candidate : source_candidate_list) {
        new_source_info.uniquely_add_source(candidate, strategy_id);
        if (new_source_info.num_source() == new_total_num_sources) break;
      }
    }
  }

// ---------------------------------------- General methods for source selection ---------------------------------------- //
  std::vector<vertex_locator_t>
  select_source(const size_t max_num_sources,
                const std::function<bool(const vertex_locator_t)> &is_candidate,
                const std::vector<std::function<uint64_t(const vertex_locator_t)>> &score_function_list) {

    std::vector<vertex_locator_t> local_candidate_list;
    select_source_in_local(max_num_sources, is_candidate,
                           score_function_list, local_candidate_list,
                           m_graph.vertices_begin(), m_graph.vertices_end());
    select_source_in_local(max_num_sources, is_candidate,
                           score_function_list, local_candidate_list,
                           m_graph.controller_begin(), m_graph.controller_end());

    return select_source_in_global(max_num_sources, score_function_list, local_candidate_list);
  }

  template <typename iterator_t>
  void select_source_in_local(const size_t max_num_sources,
                              const std::function<bool(const vertex_locator_t)> &is_candidate,
                              const std::vector<std::function<uint64_t(const vertex_locator_t)>> &score_function_list,
                              std::vector<vertex_locator_t> &selected_source_list,
                              iterator_t vitr, iterator_t end) {
    // Return true if 'lhd' has higher priority than 'rhd'
    auto is_prior = [this, &score_function_list](const vertex_locator_t &lhd, const vertex_locator_t &rhd) -> bool {
      for (auto &score_calculator : score_function_list) {
        const uint64_t score_lhd = score_calculator(lhd);
        const uint64_t score_rhd = score_calculator(rhd);
        if (score_lhd != score_rhd) return (score_lhd > score_rhd);
      }
      // --- Final tie braker--- //
      return compare_vertex_by_random_hash(lhd, rhd);
    };

    // To sort by ascending order
    auto less = [this, &is_prior](const vertex_locator_t &lhd, const vertex_locator_t &rhd) -> bool {
      return !is_prior(lhd, rhd);
    };

    for (; vitr != end; ++vitr) {
      auto locator = *vitr;
      if (!is_candidate(locator)) continue;
      if (selected_source_list.size() < max_num_sources) {
        selected_source_list.emplace_back(locator);
        if (selected_source_list.size() == max_num_sources)
          std::sort(selected_source_list.begin(), selected_source_list.end(), less);
      } else {
        if (is_prior(locator, selected_source_list[0])) {
          // kick out the lowest priority element
          selected_source_list[0] = locator;

          // Just move the lowest priority element at the begining
          std::partial_sort(selected_source_list.begin(), selected_source_list.begin() + 1,
                            selected_source_list.end(), less);
        }
      }
    }
  }

  std::vector<vertex_locator_t>
  select_source_in_global(const size_t max_num_sources,
                          const std::vector<std::function<uint64_t(const vertex_locator_t)>> &score_function_list,
                          const std::vector<vertex_locator_t> &local_candidate_list) {
    // -------------------- Construct and exchange score lists in global -------------------- //
    std::vector<std::vector<uint64_t>> global_score_matrix(score_function_list.size());
    for (size_t i = 0; i < score_function_list.size(); ++i) {
      std::vector<uint64_t> local_score_list;
      for (auto locator : local_candidate_list) {
        local_score_list.emplace_back(score_function_list[i](locator));
      }
      mpi_all_gather(local_score_list, global_score_matrix[i], MPI_COMM_WORLD);
    }

    std::vector<uint64_t> local_vid_list;
    for (auto locator : local_candidate_list) {
      local_vid_list.emplace_back(m_graph.locator_to_label(locator));
    }
    std::vector<uint64_t> global_vid_list;
    havoqgt::mpi_all_gather(local_vid_list, global_vid_list, MPI_COMM_WORLD);

    // -------------------- Select high priority sources based on the exchanged scores -------------------- //
    std::vector<std::unordered_map<uint64_t, uint64_t>> global_score_table_list(global_score_matrix.size());
    for (size_t i = 0; i < global_score_matrix.size(); ++i) {
      auto &global_score_table = global_score_table_list[i];
      auto &global_score_list = global_score_matrix[i];
      for (size_t j = 0; j < global_vid_list.size(); ++j) {
        global_score_table[global_vid_list[j]] = global_score_list[j];
      }
    }

    const size_t final_size = std::min(static_cast<size_t>(max_num_sources), global_vid_list.size());
    std::partial_sort(global_vid_list.begin(), global_vid_list.begin() + final_size, global_vid_list.end(),
                      [this, &global_score_table_list](const uint64_t lhd, const uint64_t rhd) -> bool {
                        for (auto &score_table : global_score_table_list) {
                          const uint64_t score1 = score_table[lhd];
                          const uint64_t score2 = score_table[rhd];
                          if (score1 != score2) return (score1 > score2);
                        }
                        return compare_vertex_by_random_hash(lhd, rhd);
                      });
    global_vid_list.resize(final_size);

    std::vector<vertex_locator_t> global_selected_source_list;
    for (uint64_t vid : global_vid_list) {
      global_selected_source_list.emplace_back(m_graph.label_to_locator(vid));
    }

    return global_selected_source_list;
  }

// -------------------------------------------------------------------------------------------------------------- //
// compute_ecc_k_source
// -------------------------------------------------------------------------------------------------------------- //
  template <typename iterator_t>
  void compute_ecc_k_source_helper(std::vector<level_t> &k_source_ecc, iterator_t vitr, iterator_t end) {
    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue;
      for (size_t k = 0; k < k_source_ecc.size(); ++k) {
        k_source_ecc[k] = std::max(m_kbfs.vertex_data().level(*vitr)[k], k_source_ecc[k]);
      }
    }
  }

  void compute_ecc_k_source() {
    m_source_info.ecc_list.resize(m_source_info.num_source(), std::numeric_limits<level_t>::min());
    compute_ecc_k_source_helper(m_source_info.ecc_list, m_graph.vertices_begin(), m_graph.vertices_end());
    compute_ecc_k_source_helper(m_source_info.ecc_list, m_graph.controller_begin(), m_graph.controller_end());

    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    mpi_all_reduce_inplace(m_source_info.ecc_list, std::greater<level_t>(), MPI_COMM_WORLD);
    for (size_t k = 0; k < m_source_info.ecc_list.size(); ++k) {
      auto &source = m_source_info.source_list[k];
      if (source.owner() == static_cast<uint32_t>(mpi_rank) || source.is_delegate()) {
        m_ecc_vertex_data.lower(source) = m_ecc_vertex_data.upper(source) = m_source_info.ecc_list[k];
        m_ecc_vertex_data.just_solved_status(source) = ecc_vertex_data_t::k_source;
      }
    }

    if (mpi_rank == 0) {
      std::cout << "ECC for k sources: ";
      for (level_t ecc : m_source_info.ecc_list) std::cout << ecc << " ";
      std::cout << std::endl;
    }
  }

// -------------------------------------------------------------------------------------------------------------- //
// bounding algorithm
// -------------------------------------------------------------------------------------------------------------- //
  std::pair<size_t, size_t> bound_ecc() {
    const auto ret_vrtx = bound_ecc_helper(m_graph.vertices_begin(), m_graph.vertices_end());
    const auto ret_ctrl = bound_ecc_helper(m_graph.controller_begin(), m_graph.controller_end());

    size_t num_solved = ret_vrtx.first + ret_ctrl.first;
    size_t num_unsolved = ret_vrtx.second + ret_ctrl.second;

    num_solved = mpi_all_reduce(num_solved, std::plus<size_t>(), MPI_COMM_WORLD);
    num_unsolved = mpi_all_reduce(num_unsolved, std::plus<size_t>(), MPI_COMM_WORLD);
    mpi_all_reduce_inplace(m_source_info.num_solved_list, std::plus<size_t>(), MPI_COMM_WORLD);

    return std::make_pair(num_solved, num_unsolved);
  }

  template <typename iterator_t>
  std::pair<size_t, size_t> bound_ecc_helper(iterator_t vitr, iterator_t end) {
    size_t num_solved(0);
    size_t num_unsolved(0);

    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue; // skip unvisited vertices

      level_t &lower = m_ecc_vertex_data.lower(*vitr);
      level_t &upper = m_ecc_vertex_data.upper(*vitr);
      if (lower == upper) continue; // Exact ecc has been already found

      const level_t old_lower = lower;
      const level_t old_upper = upper;
      for (size_t k = 0; k < m_source_info.num_source(); ++k) {

        const level_t level = m_kbfs.vertex_data().level(*vitr)[k];
        const level_t ecc = m_source_info.ecc_list[k];

        lower = compute_lower(lower, level, ecc);
        upper = compute_upper(upper, level, ecc);

        if (lower == upper) {
          break;
        }
      }

      if (use_soft_contribution_score) {
        if (old_lower < lower || upper < old_upper) {
          for (size_t k = 0; k < m_source_info.num_source(); ++k) {
            const level_t level = m_kbfs.vertex_data().level(*vitr)[k];
            const level_t ecc = m_source_info.ecc_list[k];
            if (lower == compute_lower_candidate(level, ecc)) {
              ++m_source_info.num_solved_list[k];
            } else if (upper == compute_upper_candidate(level, ecc)) {
              ++m_source_info.num_solved_list[k];
            }
          }
        }
      } else {
        if (lower == upper) {
          m_ecc_vertex_data.just_solved_status(*vitr) = ecc_vertex_data_t::k_bound;

          // ----- Compute contribution score ----- //
          for (size_t k = 0; k < m_source_info.num_source(); ++k) {
            const level_t level = m_kbfs.vertex_data().level(*vitr)[k];
            const level_t ecc = m_source_info.ecc_list[k];

            if (old_lower != lower && lower == compute_lower_candidate(level, ecc)) {
              ++m_source_info.num_solved_list[k];
            } else if (old_upper != upper && upper == compute_upper_candidate(level, ecc)) {
              ++m_source_info.num_solved_list[k];
            }
          }
        }
      }

      if (lower == upper) {
        ++num_solved;
      } else {
        ++num_unsolved;
      }

    }

    return std::make_pair(num_solved, num_unsolved);
  }

  level_t compute_lower_candidate(const level_t distance, const level_t ecc) {
    return std::max(distance, static_cast<uint16_t>(ecc - distance));
  }

  level_t compute_lower(const level_t current_lower, const level_t distance, const level_t ecc) {
    return std::max(current_lower, compute_lower_candidate(distance, ecc));
  }

  level_t compute_upper_candidate(const level_t distance, const level_t ecc) {
    return static_cast<uint16_t>(distance + ecc);
  }

  level_t compute_upper(const level_t current_upper, const level_t distance, const level_t ecc) {
    return std::min(current_upper, compute_upper_candidate(distance, ecc));
  }

// -------------------------------------------------------------------------------------------------------------- //
// bound_ecc_for_leaf
// -------------------------------------------------------------------------------------------------------------- //
  std::pair<size_t, size_t> bound_ecc_for_leaf() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    std::vector<uint16_t> source_height(m_source_info.num_source(), 0);
    for (size_t k = 0; k < m_source_info.num_source(); ++k) {
      const auto src = m_source_info.source_list[k];
      if (src.owner() == mpi_rank) {
        if (m_ecc_vertex_data.lower(src) > m_2core_vertex_data[src].get_height())
          source_height[k] = m_2core_vertex_data[src].get_height();
      }
    }

    mpi_all_reduce_inplace(source_height, std::greater<uint16_t>(), MPI_COMM_WORLD);

    const auto ret_vrtx = bound_ecc_for_leaf_helper(source_height, m_graph.vertices_begin(), m_graph.vertices_end());
    const auto
        ret_ctrl = bound_ecc_for_leaf_helper(source_height, m_graph.controller_begin(), m_graph.controller_end());

    size_t num_solved = ret_vrtx.first + ret_ctrl.first;
    size_t num_unsolved = ret_vrtx.second + ret_ctrl.second;

    num_solved = mpi_all_reduce(num_solved, std::plus<size_t>(), MPI_COMM_WORLD);
    num_unsolved = mpi_all_reduce(num_unsolved, std::plus<size_t>(), MPI_COMM_WORLD);

    return std::make_pair(num_solved, num_unsolved);
  }

  template <typename iterator_t>
  std::pair<size_t, size_t> bound_ecc_for_leaf_helper(const std::vector<uint16_t> &height,
                                                      iterator_t vitr, iterator_t end) {
    size_t num_solved(0);
    size_t num_unsolved(0);

    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue; // skip unvisited vertices

      level_t &lower = m_ecc_vertex_data.lower(*vitr);
      level_t &upper = m_ecc_vertex_data.upper(*vitr);
      if (lower == upper) continue; // Exact ecc has been already found

      for (size_t k = 0; k < m_source_info.num_source(); ++k) {
        if (height[k] == 0) continue;

        const level_t level = m_kbfs.vertex_data().level(*vitr)[k] + height[k];
        const level_t ecc = m_source_info.ecc_list[k] + height[k];

        lower = compute_lower(lower, level, ecc);
        upper = compute_upper(upper, level, ecc);

        if (lower == upper) {
          break;
        }
      }

      if (lower == upper) {
        ++num_solved;
        m_ecc_vertex_data.just_solved_status(*vitr) = ecc_vertex_data_t::k_bound;
      } else {
        ++num_unsolved;
      }

    }

    return std::make_pair(num_solved, num_unsolved);
  }

// -------------------------------------------------------------------------------------------------------------- //
// solve_single_degree_vertices
// -------------------------------------------------------------------------------------------------------------- //
  size_t solve_single_degree_vertices() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    size_t local_num_pruned = 0;
    auto alg_data = std::forward_as_tuple(m_ecc_vertex_data, m_2core_vertex_data, local_num_pruned);
    auto vq = create_visitor_queue<take_pruning_visitor, havoqgt::detail::visitor_priority_queue>(&m_graph,
                                                                                                  alg_data);
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD);

    const size_t global_num_pruned = mpi_all_reduce(local_num_pruned, std::plus<level_t>(), MPI_COMM_WORLD);

    return global_num_pruned;
  }

  size_t solve_hanging_tree() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    size_t local_num_pruned = 0;
    auto alg_data = std::forward_as_tuple(m_kbfs.vertex_data(),
                                          m_ecc_vertex_data,
                                          m_2core_vertex_data,
                                          local_num_pruned);
    auto vq = create_visitor_queue<hanging_tree_visitor, havoqgt::detail::visitor_priority_queue>(&m_graph,
                                                                                                  alg_data);
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD);

    const size_t global_num_pruned = mpi_all_reduce(local_num_pruned, std::plus<level_t>(), MPI_COMM_WORLD);

    return global_num_pruned;
  }

// -------------------------------------------------------------------------------------------------------------- //
// count_unsolved_neighbors
// -------------------------------------------------------------------------------------------------------------- //
  void count_unsolved_neighbors() {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    auto alg_data = std::forward_as_tuple(m_kbfs.vertex_data(), m_ecc_vertex_data);
    auto vq = create_visitor_queue<unsolved_visitor, havoqgt::detail::visitor_priority_queue>(&m_graph, alg_data);
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD);
  }

// -------------------------------------------------------------------------------------------------------------- //
// mark_farthest_vertices
// -------------------------------------------------------------------------------------------------------------- //
//  void mark_farthest_vertices() {
//    auto count_farthest = [&](const vertex_locator_t &vertex) {
//      for (size_t i = 0; i < m_source_info.num_source(); ++i) {
//        if (m_kbfs.vertex_data().level(vertex)[i] == m_source_info.ecc_list[i])
//          m_ecc_vertex_data.farthest_vertex(vertex) = true;
//      }
//    };
//
//    for (auto itr = m_graph.vertices_begin(), end = m_graph.vertices_end();
//         itr != end; ++itr) {
//      count_farthest(*itr);
//    }
//    for (auto itr = m_graph.controller_begin(), end = m_graph.controller_end();
//         itr != end; ++itr) {
//      count_farthest(*itr);
//    }
//  }

  std::vector<typename graph_t::vertex_locator> select_and_mark_farthest_vertices() {
    const uint64_t nan = std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> farthest_vertices(m_source_info.num_source(), nan);

    // Find single farthest vertex for each source
    {
      auto find_farthest_vertices = [&](const vertex_locator_t &vertex) {
        for (size_t k = 0; k < m_source_info.num_source(); ++k) {
          if (m_kbfs.vertex_data().level(vertex)[k] != m_source_info.ecc_list[k]) continue; // not a farthest vertex

          m_ecc_vertex_data.farthest_vertex(vertex) = true;

          if (farthest_vertices[k] == nan) {
            farthest_vertices[k] = m_graph.locator_to_label(vertex);
          } else {
            if (hash_vertex_id(farthest_vertices[k]) > hash_vertex_id(m_graph.locator_to_label(vertex))) {
              farthest_vertices[k] = m_graph.locator_to_label(vertex);
            }
          }
        }
      };

      for (auto itr = m_graph.vertices_begin(), end = m_graph.vertices_end();
           itr != end; ++itr) {
        find_farthest_vertices(*itr);
      }
      for (auto itr = m_graph.controller_begin(), end = m_graph.controller_end();
           itr != end; ++itr) {
        find_farthest_vertices(*itr);
      }
    }

    std::vector<uint64_t> never_source_farthest_vertices;
    {
      for (size_t k = 0; k < m_source_info.num_source(); ++k) {
        std::vector<uint64_t> global_farthest_vertices;
        mpi_all_gather(farthest_vertices[k], global_farthest_vertices, MPI_COMM_WORLD);

        std::sort(global_farthest_vertices.begin(),
                  global_farthest_vertices.end(),
                  [&](const uint64_t &lhs, const uint64_t &rhs) -> bool {
                    if (lhs == nan || rhs == nan)
                      return std::less<uint64_t>()(lhs, rhs);
                    else
                      return std::less<uint64_t>()(hash_vertex_id(lhs), hash_vertex_id(rhs));
                  });
        const uint64_t farthest_vertex = global_farthest_vertices[0];
        if (farthest_vertex == nan) continue;

        unsigned char selected_before = 0; // mpi can't handle bool
        auto locator = m_graph.label_to_locator(farthest_vertex);
        if (mpi_comm_rank() == locator.owner()) {
          selected_before = m_ecc_vertex_data.ever_source_vertex(locator);
        }
        selected_before = mpi_all_reduce(selected_before, std::greater<unsigned char>(), MPI_COMM_WORLD);
        if (!selected_before) never_source_farthest_vertices.emplace_back(farthest_vertex);
      }
    }

    // Remove duplicated vertices
    std::sort(never_source_farthest_vertices.begin(), never_source_farthest_vertices.end());
    auto end = std::unique(never_source_farthest_vertices.begin(), never_source_farthest_vertices.end());
    never_source_farthest_vertices.resize(std::distance(never_source_farthest_vertices.begin(), end));

    std::vector<typename graph_t::vertex_locator> return_list;
    for (auto id : never_source_farthest_vertices) {
      return_list.emplace_back(m_graph.label_to_locator(id));
    }

    if (mpi_comm_rank() == 0) {
      std::cout << "#farthest_vertices\t" << return_list.size() << std::endl;
    }

    return return_list;
  }

// -------------------------------------------------------------------------------------------------------------- //
// compute_diff_score
// -------------------------------------------------------------------------------------------------------------- //
  template <typename iterator_t>
  void compute_diff_score_histgram_helper(const size_t max_diff, std::vector<size_t> &histgram,
                                          iterator_t vitr, iterator_t end) {
    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue;
      const uint64_t diff = m_ecc_vertex_data.upper(*vitr) - m_ecc_vertex_data.lower(*vitr);
      assert(diff >= 0);
      if (diff > max_diff) ++histgram[max_diff + 1];
      else ++histgram[diff];
    }
  }

  std::vector<size_t> compute_diff_score_histgram(const size_t max_diff) {
    // values larger than max_diff are gathered to the last element, i.e., histgram[max_diff + 1]
    std::vector<size_t> histgram(max_diff + 2, 0);
    compute_diff_score_histgram_helper(max_diff, histgram, m_graph.vertices_begin(), m_graph.vertices_end());
    compute_diff_score_histgram_helper(max_diff, histgram, m_graph.controller_begin(), m_graph.controller_end());

    mpi_all_reduce_inplace(histgram, std::plus<size_t>(), MPI_COMM_WORLD);

    return histgram;
  }

// -------------------------------------------------------------------------------------------------------------- //
// collect information about unsolved vertices
// -------------------------------------------------------------------------------------------------------------- //
  template <typename iterator_t>
  void collect_unsolved_vertices_statistics_helper(std::map<size_t, size_t> &degree_count,
                                                   std::map<size_t, size_t> &level_count,
                                                   std::map<size_t, size_t> &lower_count,
                                                   std::map<size_t, size_t> &upper_count,
                                                   iterator_t vitr, iterator_t end) {
    auto count_up = [](const size_t key, std::map<size_t, size_t> &table) {
      if (table.count(key) == 0) table[key] = 0;
      ++table[key];
    };

    for (; vitr != end; ++vitr) {
      if (!m_kbfs.vertex_data().visited_by(*vitr, 0)) continue;
      if (m_ecc_vertex_data.upper(*vitr) == m_ecc_vertex_data.lower(*vitr)) continue;
      {
        size_t degree = m_graph.degree(*vitr);
        if (degree > 10) {
          degree = std::pow(10, static_cast<uint64_t>(std::log10(degree)));
        }
        count_up(degree, degree_count);
      }
      {
        level_t average_level = 0;
        for (size_t k = 0; k < m_source_info.num_source(); ++k) {
          average_level += m_kbfs.vertex_data().level(*vitr)[k];
        }
        average_level /= m_source_info.num_source();
        count_up(average_level, level_count);
      }
      {
        count_up(m_ecc_vertex_data.lower(*vitr), lower_count);
        count_up(m_ecc_vertex_data.upper(*vitr), upper_count);
      }
    }
  }

  void collect_unsolved_vertices_statistics() {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    std::map<size_t, size_t> degree_count;
    std::map<size_t, size_t> level_count;
    std::map<size_t, size_t> lower_count;
    std::map<size_t, size_t> upper_count;

    collect_unsolved_vertices_statistics_helper(degree_count, level_count, lower_count, upper_count,
                                                m_graph.vertices_begin(), m_graph.vertices_end());
    collect_unsolved_vertices_statistics_helper(degree_count, level_count, lower_count, upper_count,
                                                m_graph.controller_begin(), m_graph.controller_end());

    auto all_gather = [](std::map<size_t, size_t> &table) {
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

    auto print = [](std::map<size_t, size_t> &table) {
      for (auto elem : table) {
        std::cout << elem.first << " " << elem.second << ", ";
      }
      std::cout << std::endl;
    };

    if (mpi_rank == 0) {
      std::cout << "degree: ";
      print(degree_count);
      std::cout << "level: ";
      print(level_count);
      std::cout << "lower: ";
      print(lower_count);
      std::cout << "upper: ";
      print(upper_count);
    }
  }

  graph_t &m_graph;
  kbfs_t m_kbfs;
  ecc_vertex_data_t m_ecc_vertex_data;
  source_info_t m_source_info;
  source_score_function_list_t m_source_score_function_list;
  progress_info_t m_progress_info;
  kth_core_vertex_data_t m_2core_vertex_data;
  std::vector<uint32_t> m_strategy_num_to_use;
// std::vector<size_t> m_contribution_score;
};

template <typename graph_allocator_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity<graph_allocator_t, level_t, k_num_sources>::take_pruning_visitor {
 private:
  enum index {
    ecc_data = 0,
    core_info = 1,
    count_num_pruned = 2
  };

 public:
  take_pruning_visitor()
      : vertex(),
        ecc() {}

  explicit take_pruning_visitor(vertex_locator_t _vertex)
      : vertex(_vertex),
        ecc() {}

#pragma GCC diagnostic pop

  take_pruning_visitor(vertex_locator_t _vertex, level_t _cee)
      : vertex(_vertex),
        ecc(_cee) {}

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (std::get<index::ecc_data>(alg_data).just_solved_status(vertex) != ecc_vertex_data_t::k_source) return false;
    if (std::get<index::core_info>(alg_data)[vertex].get_num_cut() == 0) return false;

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      if (!eitr.target().is_delegate()) // Don't send visitor to delegates because they are obviously not degree 1 vertices
        vis_queue->queue_visitor(take_pruning_visitor(eitr.target(),
                                                      std::get<index::ecc_data>(alg_data).lower(vertex)));
    }

    return true; // trigger bcast from masters of delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const {
    return true;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    if (vertex.get_bcast()) { // for case, label:FLOW1
      for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
        if (!eitr.target().is_delegate()) // Don't send visitor to delegates because they are obviously not degree 1 vertices
          vis_queue->queue_visitor(take_pruning_visitor(eitr.target(),
                                                        std::get<index::ecc_data>(alg_data).lower(vertex)));
      }
    } else {
      // ----- Update ecc of 1 degree vertices ----- //
      if (g.degree(vertex) == 1 &&
          std::get<index::ecc_data>(alg_data).lower(vertex) != std::get<index::ecc_data>(alg_data).upper(vertex)) {
        std::get<index::ecc_data>(alg_data).lower(vertex) = std::get<index::ecc_data>(alg_data).upper(vertex) = ecc + 1;
        ++(std::get<index::count_num_pruned>(alg_data));
      }
    }

    return false;
  }

  friend inline bool operator>(const take_pruning_visitor &v1, const take_pruning_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const take_pruning_visitor &v1, const take_pruning_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator_t vertex;
  level_t ecc;
} __attribute__ ((packed));

template <typename graph_allocator_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity<graph_allocator_t, level_t, k_num_sources>::unsolved_visitor {
 private:
  enum index {
    kbfs_data = 0,
    ecc_data = 1
  };

 public:
  unsolved_visitor()
      : vertex() {}

  explicit unsolved_visitor(vertex_locator_t _vertex)
      : vertex(_vertex) {}

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
    if (std::get<index::kbfs_data>(alg_data).visited_by(vertex, 0))
      return false; // skip unvisited vertices

    if (std::get<index::ecc_data>(alg_data).lower(vertex) == std::get<index::ecc_data>(alg_data).upper(vertex))
      return false; // skip solved vertices

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(unsolved_visitor(eitr.target()));
    }

    return true; // trigger bcast from masters of delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const {
    return true;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    if (vertex.get_bcast()) { // for visitors triggered by label:FLOW1
      for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
        vis_queue->queue_visitor(unsolved_visitor(eitr.target()));
      }
    } else {
      ++std::get<index::ecc_data>(alg_data).num_unsolved_neighbors(vertex);
    }

    return false;
  }

  friend inline bool operator>(const unsolved_visitor &v1, const unsolved_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const unsolved_visitor &v1, const unsolved_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator_t vertex;
} __attribute__ ((packed));

template <typename graph_allocator_t, typename level_t, uint32_t k_num_sources>
class exact_eccentricity<graph_allocator_t, level_t, k_num_sources>::hanging_tree_visitor {
 private:
  enum index {
    kbfs_data = 0,
    ecc_data = 1,
    k_core_data = 2,
    num_solved = 3,
  };

 public:
  hanging_tree_visitor()
      : vertex() {}

  explicit hanging_tree_visitor(const vertex_locator_t _vertex)
      : vertex(_vertex) {}

  explicit hanging_tree_visitor(const vertex_locator_t _vertex, const uint16_t _ecc, const uint16_t _height)
      : vertex(_vertex),
        ecc(_ecc),
        height(_height) {}

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    // -------------------------------------------------- //
    // This function issues visitors for neighbors (scatter step)
    // -------------------------------------------------- //
#if TREE_FROM_ONLY_SOURCES
    if (std::get<index::ecc_data>(alg_data).just_solved_status(vertex) != ecc_vertex_data_t::k_source)
      return false; // skip non source vertices
#else
    if (!std::get<index::ecc_data>(alg_data).just_solved_status(vertex))
      return false; // skip non just solved vertices
#endif

    if (std::get<index::k_core_data>(alg_data)[vertex].get_num_cut() == 0)
      return false; // skip as it does not have trees

    if (std::get<index::k_core_data>(alg_data)[vertex].get_height()
        >= std::get<index::ecc_data>(alg_data).lower(vertex))
      return false; // skip vertices who have squal or higher trees than their ecc values

    if (std::get<index::k_core_data>(alg_data)[vertex].get_height() == 0) {
      std::cerr << "std::get<index::k_core_data>(alg_data)[vertex].get_height() == 0" << std::endl;
      std::abort();
    }

    for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
      vis_queue->queue_visitor(hanging_tree_visitor(eitr.target(),
                                                    std::get<index::ecc_data>(alg_data).lower(vertex) + 1,
                                                    std::get<index::k_core_data>(alg_data)[vertex].get_height() - 1));
    }

    return true; // trigger bcast from masters of delegates (label:FLOW1)
  }

  template <typename AlgData>
  bool pre_visit(AlgData &alg_data) const {
    if (std::get<index::k_core_data>(alg_data)[vertex].get_alive())
      return false; // skip alive vertices

    if (height < std::get<index::k_core_data>(alg_data)[vertex].get_height())
      return false; // only goes to the leaf side

    if (std::get<index::ecc_data>(alg_data).upper(vertex) < ecc)
      return false;

    return true;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(graph_t &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const {
    // -- Only dead vertices and delegate vertices should be here -- //
    auto my_height = std::get<index::k_core_data>(alg_data)[vertex].get_height();

    if (vertex.get_bcast()) {
      for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
        vis_queue->queue_visitor(hanging_tree_visitor(eitr.target(),
                                                      std::get<index::ecc_data>(alg_data).lower(vertex) + 1,
                                                      my_height - 1));
      }
      return false;
    }

    // All alive vertices in here MUST BE vertex.get_bcast() == true
    if (std::get<index::k_core_data>(alg_data)[vertex].get_alive()) {
      std::cerr << __FUNCTION__ << " I'm alive but not get_bcast()" << std::endl;
      std::abort();
    }

    if (height < my_height) { // Logic error
      std::cerr << height << " < " << my_height << std::endl;
      std::abort();
    }

    if (std::get<index::ecc_data>(alg_data).upper(vertex) < ecc) {
      std::cerr << std::get<index::ecc_data>(alg_data).upper(vertex) << " < " << ecc << std::endl;
      std::abort(); // Logic error
    }

    if (std::get<index::ecc_data>(alg_data).lower(vertex) != std::get<index::ecc_data>(alg_data).upper(vertex)) {

      std::get<index::ecc_data>(alg_data).lower(vertex) = std::get<index::ecc_data>(alg_data).upper(vertex) = ecc;
      std::get<index::ecc_data>(alg_data).just_solved_status(vertex) = ecc_vertex_data_t::k_tree;
      ++std::get<index::num_solved>(alg_data);

      if (my_height == 0) return false; // we are at a leaf

      for (auto eitr = g.edges_begin(vertex), end = g.edges_end(vertex); eitr != end; ++eitr) {
        vis_queue->queue_visitor(hanging_tree_visitor(eitr.target(),
                                                      ecc + static_cast<uint16_t>(1),
                                                      my_height - static_cast<uint16_t>(1)));
      }
    }

    return false;
  }

  friend inline bool operator>(const hanging_tree_visitor &v1, const hanging_tree_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  friend inline bool operator<(const hanging_tree_visitor &v1, const hanging_tree_visitor &v2) {
    return v1.vertex < v2.vertex; // or source?
  }

  vertex_locator_t vertex;
  uint16_t ecc;
  uint16_t height; // we only go to leaf side
} __attribute__ ((packed));

} // namespace havoqgt

#endif //HAVOQGT_EXACT_ECCENTRICITY_HPP
