// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED
#define HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED

#include <limits>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>

namespace havoqgt {

template <typename Graph>
class cc_visitor {
 private:
  static constexpr uint64_t k_dummy_cc = std::numeric_limits<uint64_t>::max();

 public:
  typedef typename Graph::vertex_locator vertex_locator;
  cc_visitor() {}
  cc_visitor(vertex_locator _vertex, uint64_t _cc)
      : vertex(_vertex), m_cc(_cc) {}

  cc_visitor(vertex_locator _vertex) : vertex(_vertex), m_cc(k_dummy_cc) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    auto& v_predicate = std::get<2>(alg_data);
    if (!v_predicate(vertex)) {
      return false;
    }
    auto& graph   = std::get<0>(alg_data);
    auto& cc_data = std::get<1>(alg_data);
    auto cc = (m_cc == k_dummy_cc) ? graph.locator_to_label(vertex) : m_cc;

    if (cc < cc_data[vertex]) {
      cc_data[vertex] = cc;
      return true;
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    auto& graph       = std::get<0>(alg_data);
    auto& cc_data     = std::get<1>(alg_data);
    auto& v_predicate = std::get<2>(alg_data);
    auto& e_predicate = std::get<3>(alg_data);
    if (!v_predicate(vertex)) {
      return false;
    }

    auto cc = (m_cc == k_dummy_cc) ? graph.locator_to_label(vertex) : m_cc;

    if (cc <= cc_data[vertex]) {
      cc_data[vertex] = cc;
      for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        auto neighbor = eitr.target();
        if (e_predicate(eitr)) {
          if (cc < graph.locator_to_label(neighbor)) {
            cc_visitor new_visitor(neighbor, cc);
            vis_queue->queue_visitor(new_visitor);
          }
        }
      }
      return true;
    }
    return false;
  }

  friend inline bool operator>(const cc_visitor& v1, const cc_visitor& v2) {
    if (v2.m_cc < v1.m_cc) {
      return true;
    } else if (v1.m_cc < v2.m_cc) {
      return false;
    }
    if (v1.vertex == v2.vertex) return false;
    return !(v1.vertex < v2.vertex);
  }

  vertex_locator vertex;
  uint64_t m_cc;
};

template <typename TGraph, typename CCData>
void connected_components(TGraph* g, CCData& cc_data) {
  auto v_predicate = [](const auto& vi) { return true; };
  auto e_predicate = [](const auto& e) { return true; };
  connected_components(g, cc_data, v_predicate, e_predicate);
}

template <typename TGraph, typename CCData, typename VFunction,
          typename EFunction>
void connected_components(TGraph* g, CCData& cc_data, VFunction& v_predicate,
                          EFunction& e_predicate) {
  typedef cc_visitor<TGraph> visitor_type;

  for (auto vitr = g->vertices_begin(); vitr != g->vertices_end(); ++vitr) {
    cc_data[*vitr] = g->locator_to_label(*vitr);
  }
  for (auto citr = g->controller_begin(); citr != g->controller_end(); ++citr) {
    cc_data[*citr] = g->locator_to_label(*citr);
  }

  auto alg_data = std::forward_as_tuple(*g, cc_data, v_predicate, e_predicate);
  auto vq = create_visitor_queue<visitor_type, detail::visitor_priority_queue>(
      g, alg_data);
  vq.init_visitor_traversal();
}

}  // end namespace havoqgt

#endif  // HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED
