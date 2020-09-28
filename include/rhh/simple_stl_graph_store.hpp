// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_INCLUDE_RHH_SIMPLE_STL_GRAPH_STORE_HPP_
#define HAVOQGT_INCLUDE_RHH_SIMPLE_STL_GRAPH_STORE_HPP_

#include <boost/container/scoped_allocator.hpp>
#include <boost/container/vector.hpp>
#include <boost/unordered_map.hpp>
#include <rhh/hash.hpp>

namespace rhh {

template <typename _vertex_type, typename _vertex_value_type,
          typename _edge_value_type,
          typename _allocator_type = std::allocator<std::byte>>
class simple_stl_graph_store {
 public:
  using vertex_type       = _vertex_type;
  using vertex_value_type = _vertex_value_type;
  using edge_value_type   = _edge_value_type;
  using allocator_type    = _allocator_type;

 private:
  using edge_list_allocator_t =
      typename std::allocator_traits<allocator_type>::template rebind_alloc<
          std::pair<const vertex_type, edge_value_type>>;

  using edge_list_t = boost::unordered_multimap<
      vertex_type, edge_value_type, rhh::hash<vertex_type>,
      std::equal_to<vertex_type>,
      boost::container::scoped_allocator_adaptor<edge_list_allocator_t>>;

  using adjacencly_list_allocator_t =
      typename std::allocator_traits<allocator_type>::template rebind_alloc<
          std::pair<const vertex_type, edge_list_t>>;
  using adjacencly_list_t = boost::unordered_map<
      vertex_type, edge_list_t, rhh::hash<vertex_type>,
      std::equal_to<vertex_type>,
      boost::container::scoped_allocator_adaptor<adjacencly_list_allocator_t>>;

  using vertex_table_allocator_t =
      typename std::allocator_traits<allocator_type>::template rebind_alloc<
          std::pair<const vertex_type, vertex_value_type>>;
  using vertex_table_t =
      boost::unordered_map<vertex_type, vertex_value_type,
                           rhh::hash<vertex_type>, std::equal_to<vertex_type>,
                           vertex_table_allocator_t>;

 public:
  using vertex_iterator = typename adjacencly_list_t::const_iterator;
  using edge_iterator   = typename edge_list_t::const_iterator;

  explicit simple_stl_graph_store(
      const allocator_type &alloc = allocator_type())
      : m_adj_list(alloc), m_vertex_table(alloc) {}

  simple_stl_graph_store(const simple_stl_graph_store &)     = default;
  simple_stl_graph_store(simple_stl_graph_store &&) noexcept = default;

  ~simple_stl_graph_store() = default;

  simple_stl_graph_store &operator=(const simple_stl_graph_store &) = default;
  simple_stl_graph_store &operator=(simple_stl_graph_store &&) noexcept =
      default;

  bool has_vertex(const vertex_type &vertex) const {
    return m_vertex_table.count(vertex) > 0;
  }

  bool insert_vertex(
      const vertex_type &      vertex,
      const vertex_value_type &vertex_value = vertex_value_type()) {
    if (has_vertex(vertex)) return false;
    m_vertex_table.emplace(vertex, vertex_value);
    return true;
  }

  bool insert_vertex(vertex_type &&    vertex,
                     vertex_value_type vertex_value = vertex_value_type()) {
    if (has_vertex(vertex)) return false;
    m_vertex_table.emplace(std::move(vertex), std::move(vertex_value));
    return true;
  }

  void insert_edge(const vertex_type &    source_vertex,
                   const vertex_type &    target_vertex,
                   const edge_value_type &edge_value = edge_value_type()) {
    m_adj_list[source_vertex].emplace(target_vertex, edge_value);
  }

  void insert_edge(vertex_type &&source_vertex, vertex_type &&target_vertex,
                   edge_value_type edge_value = edge_value_type()) {
    m_adj_list[std::move(source_vertex)].emplace(std::move(target_vertex), std::move(edge_value));
  }

  std::size_t num_vertices() const { return m_adj_list.size(); }

  std::size_t degree(const vertex_type &vertex) const {
    if (m_adj_list.count(vertex) == 0) return 0;
    return m_adj_list.at(vertex).size();
  }

  vertex_value_type &vertex_value(const vertex_type &vertex) {
    return m_vertex_table.at(vertex);
  }

  const vertex_value_type &vertex_value(const vertex_type &vertex) const {
    return m_vertex_table.at(vertex);
  }

  edge_value_type &edge_value(const vertex_type &source,
                              const vertex_type &target) {
    return m_vertex_table.at(source).at(target);
  }

  const edge_value_type &edge_value(const vertex_type &source,
                                    const vertex_type &target) const {
    return m_vertex_table.at(source).at(target);
  }

  vertex_iterator vertices_begin() const { return m_adj_list.begin(); }

  vertex_iterator vertices_end() const { return m_adj_list.end(); }

  edge_iterator edges_begin(
      const vertex_iterator &source_vertex_iterator) const {
    return source_vertex_iterator->second.begin();
  }

  edge_iterator edges_end(const vertex_iterator &source_vertex_iterator) const {
    return source_vertex_iterator->second.end();
  }

 private:
  adjacencly_list_t m_adj_list;
  vertex_table_t    m_vertex_table;
};

}  // namespace rhh

#endif  // HAVOQGT_INCLUDE_RHH_SIMPLE_STL_GRAPH_STORE_HPP_
