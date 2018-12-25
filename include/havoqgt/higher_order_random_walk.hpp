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

#ifndef HAVOQGT_HIGHER_ORDER_RANDOM_WALK_HPP
#define HAVOQGT_HIGHER_ORDER_RANDOM_WALK_HPP

#define WRITE_LOG 0
#if WRITE_LOG
std::ofstream *ofs_log;
#define LOG_OUT (*ofs_log)
//#define LOG_OUT std::cerr
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <tuple>
#include <cassert>
#include <chrono>
#include <sstream>
#include <unordered_map>

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace ho_rw {

template <typename graph_type, int k_num_histories>
class algo_data_type {
 private:
  using vertex_locator = typename graph_type::vertex_locator;
  using count_table_type = typename graph_type::template vertex_data<uint64_t, std::allocator < uint64_t>>;
  using history_select_prob_table_type = std::array<int, k_num_histories>;

 public:
  using history_type = std::array<uint64_t, k_num_histories>;

  explicit algo_data_type(graph_type *const graph,
                          const bool personalized,
                          const vertex_locator &seed,
                          const std::size_t num_local_walkers,
                          const uint64_t die_rate,
                          const uint64_t max_walk_length,
                          const bool warp_to_seed,
                          const history_select_prob_table_type &history_select_prob_table)
      : m_graph(graph),
        m_personalized(personalized),
        m_seed_locator(seed),
        m_num_local_walkers(num_local_walkers),
        m_die_rate(die_rate),
        m_max_walk_length(max_walk_length),
        m_num_launched_walkers(0),
        m_rnd_gen((uint32_t)std::chrono::system_clock::now().time_since_epoch().count()),
        m_warp_to_seed(warp_to_seed),
        m_num_dead_table(*graph),
        m_num_visits_table(*graph),
        m_history_select_prob_table(history_select_prob_table) {
    m_num_dead_table.reset(0);
    m_num_visits_table.reset(0);
  }

  bool rw_launcher(const vertex_locator vertex) {
    if (m_graph->degree(vertex) == 0) return false;

    // Personalized model
    if (m_personalized) {
      if (m_seed_locator == vertex && m_num_launched_walkers < m_num_local_walkers) {
        ++m_num_launched_walkers;
        return true;
      }
    } else {
      std::uniform_int_distribution<uint64_t> dis(0, m_graph->num_local_vertices() - 1);
      if (dis(m_rnd_gen) < m_num_local_walkers) {
        ++m_num_launched_walkers;
        return true;
      }
    }

    return false;
  }

  history_type new_history() {
    history_type history;
    for (int i = 0; i < k_num_histories; ++i) history[i] = m_graph->locator_to_label(warp_machine());
    return history;
  }

  bool russian_roulette() {
    std::uniform_int_distribution<int> dis(0, 99);
    return (dis(m_rnd_gen) < m_die_rate);
  }

  vertex_locator warp_machine() {
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->max_global_vertex_id() - 1);
    return m_graph->label_to_locator(dis(m_rnd_gen));
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex) {
    assert(m_graph->degree(vertex) > 0);
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->degree(vertex) - 1);
    const uint64_t offset = dis(m_rnd_gen);
    auto edge = m_graph->edges_begin(vertex) + offset;
    return edge.target();
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex, const history_type &history) {
    for (int i = 0; i < k_num_histories; ++i) {
      for (auto eitr = m_graph->edges_begin(vertex), end = m_graph->edges_end(vertex); eitr != end; ++eitr) {
        if (m_graph->locator_to_label(eitr.target()) != history[i]) continue;

        std::uniform_int_distribution<int> dis(0, 99);
        if (dis(m_rnd_gen) < m_history_select_prob_table[i]) return eitr.target();
      }
    }
    return neighbor_roulette(vertex);
  }

  uint64_t max_walk_length() const {
    return m_max_walk_length;
  }

  void increment_num_dead(const vertex_locator vertex) {
    ++m_num_dead_table[vertex];
  }

  uint64_t num_dead(const vertex_locator vertex) const {
    return m_num_dead_table[vertex];
  }

  void increment_num_visits(const vertex_locator vertex) {
    ++m_num_visits_table[vertex];
  }

  uint64_t num_visits(const vertex_locator vertex) const {
    return m_num_visits_table[vertex];
  }

  uint64_t num_launched_walkers() const {
    return m_num_launched_walkers;
  }

  bool personalized() const {
    return m_personalized;
  }

  vertex_locator seed() const {
    return m_seed_locator;
  }

  graph_type *graph() {
    return m_graph;
  }

 private:
  graph_type *m_graph;
  bool m_personalized;
  vertex_locator m_seed_locator;
  uint64_t m_num_local_walkers;
  uint64_t m_die_rate;
  uint64_t m_max_walk_length;
  uint64_t m_num_launched_walkers;
  std::mt19937 m_rnd_gen;
  bool m_warp_to_seed;
  count_table_type m_num_dead_table;
  count_table_type m_num_visits_table;
  history_select_prob_table_type m_history_select_prob_table;
};

template <typename graph_type, int k_num_histories>
class rw_visitor {
 private:
  using history_type = std::array<uint64_t, k_num_histories>;

 public:
  using vertex_locator = typename graph_type::vertex_locator;

  rw_visitor()
      : vertex(),
        walk_length(0),
        history() {}

  explicit rw_visitor(vertex_locator _vertex)
      : vertex(_vertex), walk_length(0), id(-1), history() {}

  rw_visitor(vertex_locator _vertex, int _walk_length, int _id, const history_type &_history)
      : vertex(_vertex), walk_length(_walk_length), id(_id), history(_history) {}

  template <typename visitor_queue_handle, typename alg_data_type>
  bool init_visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    bool launched = false;
    if (alg_data.personalized()) {
      while (alg_data.rw_launcher(vertex)) {
#if WRITE_LOG
        LOG_OUT << "Personalized start from " << g.locator_to_label(vertex) << std::endl;
#endif
        vis_queue->queue_visitor(rw_visitor(vertex,
                                            1,
                                            g.locator_to_label(vertex),
                                            alg_data.new_history()));
        launched = true;
      }
    } else {
      if (alg_data.rw_launcher(vertex)) {
#if WRITE_LOG
        LOG_OUT << "Start from " << g.locator_to_label(vertex) << std::endl;
#endif
        vis_queue->queue_visitor(rw_visitor(vertex,
                                            1,
                                            g.locator_to_label(vertex),
                                            alg_data.new_history()));
        launched = true;
      }
    }
    return launched;
  }

  template <typename alg_data_type>
  bool pre_visit(alg_data_type &alg_data) const {
    if (walk_length == alg_data.max_walk_length()) {
#if WRITE_LOG
      LOG_OUT << "Has reached visit limit " << serialize(*(alg_data.graph())) << "\n";
#endif
      return false;
    }

    return true;
  }

  template <typename visitor_queue_handle, typename alg_data_type>
  bool visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    alg_data.increment_num_visits(vertex);

    if (alg_data.russian_roulette()) { // Die at here and increment the die score
      alg_data.increment_num_dead(vertex);
#if WRITE_LOG
      LOG_OUT << "Increment die score at " << serialize(g) << "\n";
#endif
      return false;
    }

    if (g.degree(vertex) == 0) {
      auto target = alg_data.warp_machine();
      rw_visitor new_visitor = rw_visitor(target, walk_length + 1, id, append_history(g.locator_to_label(vertex)));
#if WRITE_LOG
      LOG_OUT << "Dead end (out degree 0) at " << serialize(g) << ". Restart from " << new_visitor.serialize(g)
              << "\n";
#endif
      vis_queue->queue_visitor(new_visitor);
      return true;
    }

    auto target = alg_data.neighbor_roulette(vertex, history);
    rw_visitor new_visitor = rw_visitor(target, walk_length + 1, id, append_history(g.locator_to_label(vertex)));
#if WRITE_LOG
    LOG_OUT << "Visit " << serialize(g) << " -> " << new_visitor.serialize(g) << "\n";
#endif

    vis_queue->queue_visitor(new_visitor);
    return true;
  }

  friend inline
  bool operator>(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex != v2.vertex) && !(v1.vertex < v2.vertex);
  }

  friend inline bool operator<(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex < v2.vertex);
  }

  history_type append_history(const uint64_t new_vertex) const {
    history_type new_history;

    for (int i = 0; i < k_num_histories - 1; ++i) new_history[i] = history[i + 1];
    new_history[k_num_histories - 1] = new_vertex;

    return new_history;
  }

  std::string serialize(const graph_type &g) const {
    std::stringstream ss;
    ss << "id " << id << ", lenght " << walk_length << ", vertex " << g.locator_to_label(vertex) << ", history";
    for (int i = 0; i < k_num_histories; ++i) {
      ss << " [" << i << "] " << history[i];
    }
    return ss.str();
  }

  vertex_locator vertex;
  int walk_length;
  int id;
  history_type history;
};

template <typename graph_type, int k_num_histories>
void run_rw(graph_type *const graph, algo_data_type<graph_type, k_num_histories> &algo_data) {

#if WRITE_LOG
  ofs_log = new std::ofstream(std::string("./log_") + std::to_string(mpi_rank));
#endif
  auto vq = havoqgt::create_visitor_queue<rw_visitor<graph_type, k_num_histories>,
                                          havoqgt::detail::visitor_priority_queue>(graph, algo_data);

  MPI_Barrier(MPI_COMM_WORLD);
  const double start_time = MPI_Wtime();
  vq.init_visitor_traversal();
  MPI_Barrier(MPI_COMM_WORLD);
  const double end_time = MPI_Wtime();
#if WRITE_LOG
  ofs_log->close();
    delete ofs_log;
#endif
  if (havoqgt::mpi_comm_rank() == 0) {
    std::cout << "Execution time: " << end_time - start_time << std::endl;
  }
}

template <typename graph_type>
std::vector<std::pair<double, uint64_t>>
compute_global_top_scores(const graph_type *const graph,
                          const uint64_t num_top,
                          const std::function<double(typename graph_type::vertex_locator)> &score_function) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  // Find the top vertices
  std::map<double, uint64_t> top_vertices;
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    const auto score = score_function(*vitr);
    // std::cout << graph->locator_to_label(*vitr) << " >> " << score << std::endl;
    if (top_vertices.size() < num_top && 0 < score) {
      top_vertices.emplace(score, graph->locator_to_label(*vitr));
    } else if (top_vertices.begin()->first < score) {
      top_vertices.emplace(score, graph->locator_to_label(*vitr));
      top_vertices.erase(top_vertices.begin());
    }
  }

  // Separate the array to use MPI
  std::vector<double> local_top_score;
  std::vector<uint64_t> local_top_id;
  for (const auto elem : top_vertices) {
    local_top_score.emplace_back(elem.first);
    local_top_id.emplace_back(elem.second);
  }

  std::vector<double> global_top_score;
  std::vector<uint64_t> global_top_id;
  havoqgt::mpi_all_gather(local_top_score, global_top_score, MPI_COMM_WORLD);
  havoqgt::mpi_all_gather(local_top_id, global_top_id, MPI_COMM_WORLD);

  std::vector<std::pair<double, uint64_t>> global_top;
  for (uint64_t i = 0; i < std::min((std::size_t)num_top, (std::size_t)global_top_score.size()); ++i) {
    global_top.emplace_back(std::make_pair(global_top_score.at(i), global_top_id.at(i)));
  }
  std::sort(global_top.begin(), global_top.end(),
            [](const std::pair<double, uint64_t> &lh, std::pair<double, uint64_t> &rh) -> bool {
              return (lh.first > rh.first);
            });
  return global_top;
}

template <typename graph_type, typename algo_data_type>
std::pair<std::vector<std::pair < double, uint64_t>>,
std::vector<std::pair<double, uint64_t>>>
compute_top_scores(const graph_type *const graph,
                   const algo_data_type &algo_data,
                   const std::size_t num_top) {
  using vertex_locator = typename graph_type::vertex_locator;

  uint64_t local_total_visits = 0;
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    local_total_visits += algo_data.num_visits(*vitr);
  }
  const uint64_t
      global_total_visits = havoqgt::mpi_all_reduce(local_total_visits, std::plus<uint64_t>(), MPI_COMM_WORLD);

  uint64_t local_total_dead = 0;
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    local_total_dead += algo_data.num_dead(*vitr);
  }
  const uint64_t global_total_dead = havoqgt::mpi_all_reduce(local_total_dead, std::plus<uint64_t>(), MPI_COMM_WORLD);

  const uint64_t global_num_walkers = havoqgt::mpi_all_reduce(algo_data.num_launched_walkers(),
                                                              std::plus<uint64_t>(),
                                                              MPI_COMM_WORLD);

  if (havoqgt::mpi_comm_rank() == 0) {
    std::cout << "#launched wk: " << global_num_walkers << std::endl;
    std::cout << "sum of visits: " << global_total_visits << std::endl;
    std::cout << "sum of dead: " << global_total_dead << std::endl;
  }

  const auto visit_score_func = [&algo_data, &global_total_visits](vertex_locator vertex) -> double {
    return static_cast<double>(algo_data.num_visits(vertex)) / global_total_visits;
  };
  const auto top_visit_scores = compute_global_top_scores(graph, num_top, visit_score_func);

  const auto dead_score_func = [&algo_data, &global_total_dead](vertex_locator vertex) -> double {
    return static_cast<double>(algo_data.num_dead(vertex)) / global_total_dead;
  };
  const auto top_dead_scores = compute_global_top_scores(graph, num_top, dead_score_func);

  return std::make_pair(top_visit_scores, top_dead_scores);
}

double compute_recall(const std::vector<uint64_t> &top_vertices,
                      const uint64_t target_community_id,
                      const std::size_t community_size,
                      const std::unordered_map<uint64_t, uint64_t> &true_community_table) {
  assert(top_vertices.size() <= community_size);

  if (community_size == 0) return 0;

  std::size_t count_true_positives = 0;
  for (const auto vid : top_vertices) {
    if (true_community_table.at(vid) == target_community_id) {
      ++count_true_positives;
    }
  }

  return (double)count_true_positives / top_vertices.size();
}

}

#endif //HAVOQGT_HIGHER_ORDER_RANDOM_WALK_HPP
