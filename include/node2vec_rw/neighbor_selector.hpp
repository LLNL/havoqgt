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


#ifndef HAVOQGT_NEIGHBOR_SELECTOR_HPP
#define HAVOQGT_NEIGHBOR_SELECTOR_HPP

#include <node2vec_rw/alias_method.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <node2vec_rw/k_bfs_level_diff_table.hpp>

namespace node2vec_rw {

namespace {
template <typename graph_allocator_type>
using dp_graph = havoqgt::delegate_partitioned_graph<graph_allocator_type>;

template <typename graph_allocator_type>
using dp_vertex_locator = typename dp_graph<graph_allocator_type>::vertex_locator;
}

template <typename graph_allocator_type>
bool check_wedge_closer(const dp_vertex_locator<graph_allocator_type> &previous,
                        const dp_vertex_locator<graph_allocator_type> &current,
                        const dp_vertex_locator<graph_allocator_type> &next,
                        const k_bfs_level_diff_table<graph_allocator_type> &k_bfs_level_diff_table) {
  for (std::size_t i = 0; i < k_bfs_level_diff_table.k_size(current, previous); ++i) {
    if (std::abs(k_bfs_level_diff_table.get(current, previous, i) - k_bfs_level_diff_table.get(current, next, i)) == 2) {
      return true;
    }
  }
  return false;
}

template <typename graph_allocator_type, typename edge_weight_type>
double compute_transitioning_probability(const dp_vertex_locator<graph_allocator_type> &previous_vertex,
                                         const dp_vertex_locator<graph_allocator_type> &current_vertex,
                                         const dp_vertex_locator<graph_allocator_type> &target_vertex,
                                         const edge_weight_type &edge_weight,
                                         const k_bfs_level_diff_table<graph_allocator_type> &k_bfs_level_diff_table,
                                         const double p,
                                         const double q,
                                         const bool initial_step) {
  double probability;

  if (initial_step) {
    probability = edge_weight;
  } else {
    if (target_vertex == previous_vertex) {
      probability = edge_weight * p; // TODO: check direction
    } else if (check_wedge_closer(previous_vertex, current_vertex, target_vertex, k_bfs_level_diff_table)) {
      probability = edge_weight;
    } else {
      probability = edge_weight * q;
    }
  }

  return probability;
}

template <typename graph_allocator_type, typename edge_weight_data_type>
class neighbor_selector_with_rejection_sampling {

 private:
  using graph_type = dp_graph<graph_allocator_type>;
  using vertex_type = dp_vertex_locator<graph_allocator_type>;

 public:
  neighbor_selector_with_rejection_sampling(const graph_type &graph,
                                            const edge_weight_data_type &edge_weight_data,
                                            const k_bfs_level_diff_table<graph_allocator_type> &k_bfs_level_diff_table,
                                            const double p,
                                            const double q,
                                            const bool small_edge_weight_variance = false)
      : m_p(p),
        m_q(q),
        m_graph(graph),
        m_edge_weight_data(edge_weight_data),
        m_k_bfs_level_diff_table(k_bfs_level_diff_table),
        m_proposal_distribution_table(),
        m_proposal_probability_density_table() {
    priv_set_proposal_distribution();
  }

  /// \brief Select a next vertex using the Rejection Sampling
  /// If initial_step is true, the value of previous_vertex is not used
  template <typename random_generator_t>
  vertex_type select(const vertex_type &previous_vertex,
                     const vertex_type &current_vertex,
                     const bool initial_step,
                     random_generator_t *random_generator) const {
    assert(m_proposal_distribution_table.count(current_vertex));
    const auto &proposal_distribution = m_proposal_distribution_table.at(current_vertex);

    assert(m_proposal_probability_density_table.count(current_vertex));
    const auto &proposal_probability_densities = m_proposal_probability_density_table.at(current_vertex);

    std::size_t count_rejections = 0;
    while (true) {
      // Sample an edge randomly with the proposal distribution
      const auto sample_offset = proposal_distribution(*random_generator);
      assert(sample_offset < m_graph.degree(current_vertex));
      assert(sample_offset < proposal_probability_densities.size());
      const double p = proposal_probability_densities.at(sample_offset);

      // Calculate the transitioning probability to the picked up edge
      assert(sample_offset < m_graph.degree(current_vertex));
      const auto sample_edge_itr = m_graph.edges_begin(current_vertex) + sample_offset;
      const auto sample_vertex = sample_edge_itr.target();
      const auto sample_probability = compute_transitioning_probability(previous_vertex,
                                                                        current_vertex,
                                                                        sample_vertex,
                                                                        m_edge_weight_data[sample_edge_itr],
                                                                        m_k_bfs_level_diff_table,
                                                                        m_p,
                                                                        m_q,
                                                                        initial_step);

      std::uniform_real_distribution<> uniform_real_distribution(std::numeric_limits<double>::epsilon(), 1.0);
      const auto u = uniform_real_distribution(*random_generator);
      const bool accept = ((double)u < (double)sample_probability / p);
      if (accept) {
        return sample_vertex;
      }
      ++count_rejections;
    }

    assert(false);
  }

 private:
  struct vertex_hash {
    std::size_t operator()(const vertex_type &v) const noexcept {
      return hash<decltype(v.local_id())>{}(v.local_id());
    }
  };

  void priv_set_proposal_distribution() {

    // Allocate space first
    for (auto v_itr = m_graph.vertices_begin(), v_end = m_graph.vertices_end(); v_itr != v_end; ++v_itr) {
      const auto vertex = *v_itr;
      if (m_graph.degree(vertex) == 0) continue;

      m_proposal_distribution_table.emplace(vertex, typename decltype(m_proposal_distribution_table)::mapped_type());
      m_proposal_probability_density_table.emplace(vertex,
                                                   typename decltype(m_proposal_probability_density_table)::mapped_type());
    }

    const double max_factor = std::max((double)1.0, (double)std::max(1.0f / m_p, 1.0f / m_q));

    for (auto v_itr = m_graph.vertices_begin(), v_end = m_graph.vertices_end(); v_itr != v_end; ++v_itr) {
      const auto vertex = *v_itr;

      if (m_graph.degree(vertex) == 0) continue;

      std::vector<double> proposal_probability_density;
      proposal_probability_density.reserve(m_graph.degree(vertex));
      for (auto e_itr = m_graph.edges_begin(vertex), e_end = m_graph.edges_end(vertex); e_itr != e_end; ++e_itr) {
        proposal_probability_density.emplace_back(max_factor * m_edge_weight_data[e_itr]);
      }

      assert(m_proposal_distribution_table.count(vertex));
      m_proposal_distribution_table.at(vertex).reset_table(proposal_probability_density.begin(), proposal_probability_density.end());

      assert(m_proposal_probability_density_table.count(vertex));
      m_proposal_probability_density_table.at(vertex) = std::move(proposal_probability_density);
    }
  }

  const graph_type &m_graph;
  const edge_weight_data_type& m_edge_weight_data;
  const k_bfs_level_diff_table<graph_allocator_type> &m_k_bfs_level_diff_table;
  const double m_p{1.0}; // return parameter
  const double m_q{1.0}; // DFS parameter
  std::unordered_map<vertex_type, alias_method<uint64_t, double>, vertex_hash> m_proposal_distribution_table;
  std::unordered_map<vertex_type, std::vector<double>, vertex_hash> m_proposal_probability_density_table;
};

#if 0
template <typename graph_type>
class neighbor_selector_with_rejection_sampling_fixed_propose_dist {

 private:
  using vertex_t = typename graph_type::vertex_type;

 public:
  template <typename level_type>
  neighbor_selector_with_rejection_sampling_fixed_propose_dist(const double p,
                                                               const double q,
                                                               const graph_type &graph,
                                                               const kbfs::level_table_type<level_type> &kbfs_level_table)
      : m_p(p),
        m_q(q),
        m_graph(graph),
        m_neighbor_selector_helper(p, q, graph, kbfs_level_table),
        m_max_probability_table() {
    priv_find_max_probability();
  }

  /// \brief Select a next vertex using the Rejection Sampling
  template <typename random_generator_t>
  std::pair<vertex_t, std::size_t> select(const vertex_t &previous_vertex,
                                          const vertex_t &current_vertex,
                                          const std::size_t step,
                                          random_generator_t *random_generator) const {
    const double max_probability = m_max_probability_table.at(current_vertex);
    std::size_t count_rejections = 0;
    while (true) {
      // Sample an edge randomly with the uniform distribution
      std::uniform_int_distribution<uint64_t> uniform_int_distribution(0, m_graph.degree(current_vertex) - 1);
      const auto sample_offset = uniform_int_distribution(*random_generator);

      // Calculate the transitioning probability to the picked up edge
      assert(sample_offset < m_graph.degree(current_vertex));
      const auto sample_edge_itr = m_graph.edges_begin(current_vertex) + sample_offset;
      const auto sample_vertex = *sample_edge_itr;
      const double sample_probability
          = m_neighbor_selector_helper.compute_transitioning_probability(previous_vertex,
                                                                         current_vertex,
                                                                         sample_vertex,
                                                                         m_graph.edge_weight(sample_edge_itr),
                                                                         (step == 0));

      std::uniform_real_distribution<> uniform_real_distribution(std::numeric_limits<double>::epsilon(), 1.0);
      const double u = uniform_real_distribution(*random_generator);
      const bool accept = (u < sample_probability / max_probability);
      if (accept) {
        return std::make_pair(sample_vertex, count_rejections);
      }
      ++count_rejections;
    }
    assert(false);
  }

 private:
  void priv_find_max_probability() {

    // Allocate space first
    for (auto v_itr = m_graph.vertices_begin(), v_end = m_graph.vertices_end(); v_itr != v_end; ++v_itr) {
      const auto vid = *v_itr;
      m_max_probability_table.emplace(vid, (double)0.0);
    }

    const double max_factor = std::max((double)1.0, (double)std::max(1.0f / m_p, 1.0f / m_q));

    OMP_DIRECTIVE(parallel)
    {
      const auto vertex_range = utility::cal_range(0,
                                                   m_graph.num_vertices() - 1,
                                                   omp::get_thread_num(),
                                                   omp::get_num_threads());
      for (auto offset = vertex_range.first; offset <= vertex_range.second; ++offset) {
        const auto vid = *(m_graph.vertices_begin() + offset);

        double max_probability = 0.0;
        for (auto e_itr = m_graph.edges_begin(vid), e_end = m_graph.edges_end(vid); e_itr != e_end; ++e_itr) {
          max_probability = std::max(max_probability, max_factor * m_graph.edge_weight(e_itr));
        }
        m_max_probability_table.at(vid) = max_probability;
      }
    }
  }

  const double m_p{1.0}; // return parameter
  const double m_q{1.0}; // DFS parameter
  const graph_type &m_graph;
  const neighbor_selector_helper<graph_type> m_neighbor_selector_helper;
  std::unordered_map<uint64_t, double> m_max_probability_table;
};
#endif
}

#endif //HAVOQGT_NEIGHBOR_SELECTOR_HPP
