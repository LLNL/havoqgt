#ifndef GRAPHSTORE_BASELINE_HPP
#define GRAPHSTORE_BASELINE_HPP

#include <boost/unordered_map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>


namespace graphstore {

#define ERASE_WITH_SWAP 1

template <typename _vertex_type,
          typename _vertex_property_data_type,
          typename _edge_property_data_type,
          typename _segment_manager_type>
class graphstore_baseline
{
public:
  using vertex_type                 = _vertex_type;
  using vertex_property_data_type   = _vertex_property_data_type;
  using edge_property_data_type     = _edge_property_data_type;
  using segment_manager_type        = _segment_manager_type;

  /// iterator
  class vertex_iterator;
  class adjacent_edge_iterator;

private:
  using size_type             = size_t;
  template <typename t1, typename t2>
  using pair_type             = graphstore::utility::packed_pair<t1, t2>;
  using allocator_type        = boost::interprocess::allocator<void, segment_manager_type>;


  /// adjacency list
  using edge_vec_element_type = pair_type<vertex_type, edge_property_data_type>;
  using vec_allocator_type    = boost::interprocess::allocator<edge_vec_element_type, segment_manager_type>;
  using edge_vec_type         = boost::interprocess::vector<edge_vec_element_type, vec_allocator_type>;

  /// vertex table
  using map_value_type        = pair_type<vertex_property_data_type, edge_vec_type>;
  using map_element_type      = std::pair<const vertex_type, map_value_type>;
  using map_allocator_type    = boost::interprocess::allocator<map_element_type, segment_manager_type>;
  using map_table_type        = boost::unordered_map<vertex_type,
                                                     map_value_type,
                                                     boost::hash<vertex_type>,
                                                     std::equal_to<vertex_type>,
                                                     map_allocator_type>;

public:

  explicit graphstore_baseline(segment_manager_type* segment_manager) :
    m_allocator(segment_manager),
    m_map_table(m_allocator),
    m_num_edges(0)
  {}


  /// -------- Lookup -------- ///
  std::pair<vertex_iterator, vertex_iterator> vertices()
  {
    return std::make_pair(vertices_begin(), vertices_end());
  }

  vertex_iterator vertices_begin()
  {
    return vertex_iterator(m_map_table.begin());
  }

  vertex_iterator vertices_end()
  {
    return vertex_iterator(m_map_table.end());
  }

  adjacent_edge_iterator adjacent_edge_begin(const vertex_type& src_vrt)
  {
    const auto itr = m_map_table.find(src_vrt);
    edge_vec_type& edge_vec = itr->second.second;
    return adjacent_edge_iterator(edge_vec.begin());
  }

  adjacent_edge_iterator adjacent_edge_end(const vertex_type& src_vrt)
  {
    return adjacent_edge_iterator(m_map_table.find(src_vrt)->second.second.end());
  }

  // Returns the degree of a vertex.
  size_type degree(const vertex_type& vertex)
  {
    const auto value = m_map_table.find(vertex);
    if (value == m_map_table.end()) {
      return 0;
    } else {
      edge_vec_type& edge_vec = value->second.second;
      return edge_vec.size();
    }
  }

  /// -------- Modifiers ------- ////
  bool insert_edge(const vertex_type& src, const vertex_type& trg, const edge_property_data_type& edge_property)
  {

    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) { // new vertex
      const vertex_property_data_type dummy_vp = false;
      map_value_type map_value(dummy_vp,
                               edge_vec_type(1, edge_vec_element_type(trg, edge_property), m_allocator));
      m_map_table.insert(map_element_type(src, map_value));
    } else {
      edge_vec_type& edge_vec = value->second.second;

      /// --- find duplicated edge --- ///
      for (const auto edge : edge_vec) {
        if (edge.first ==  trg) {
          return false; /// found duplicated edge
        }
      }

      /// --- since there is no empty spase, add it at end --- ///
      edge_vec.push_back(edge_vec_element_type(trg, edge_property));
    }

    ++m_num_edges;
    return true;
  }

  /// allows duplicated insertion
  bool insert_edge_dup(const vertex_type& src, const vertex_type& trg, const edge_property_data_type& edge_property)
  {

    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) { // new vertex
      const vertex_property_data_type dummy_vp = false;
      map_value_type map_value(dummy_vp,
                               edge_vec_type(1, edge_vec_element_type(trg, edge_property), m_allocator));
      m_map_table.insert(map_element_type(src, map_value));
    } else {
      edge_vec_type& edge_vec = value->second.second;
      /// --- just simply add it at end --- ///
      edge_vec.push_back(edge_vec_element_type(trg, edge_property));
    }

    ++m_num_edges;
    return true;
  }

  bool erase_edge(const vertex_type& src, const vertex_type& trg)
  {
    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) {
      return false;
    }

    edge_vec_type& edge_vec = value->second.second;
    for (auto itr = edge_vec.begin(), end = edge_vec.end(); itr != end; ++itr) {
      if (itr->first ==  trg) {
#if ERASE_WITH_SWAP
        using std::swap;
        swap(*itr, *(end-1));
        edge_vec.erase(end-1, end);
#else
        edge_vec.erase(itr, itr+1);
#endif
        if (edge_vec.size() == 0) {
          m_map_table.erase(value);
        }
        --m_num_edges;
        return true;
      }
    }

    return false;
  }

  /// this fuction deletes duplicated edges
  size_type erase_edge_dup(const vertex_type& src, const vertex_type& trg)
  {
    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) {
      return 0;
    }

    size_t count = 0;
    edge_vec_type& edge_vec = value->second.second;
    while (true) {
      bool erased_edge = false;
      for (auto itr = edge_vec.begin(), end = edge_vec.end(); itr != end; ++itr) {
        if (itr->first ==  trg) {
#if ERASE_WITH_SWAP
        using std::swap;
        swap(*itr, *(end-1));
        edge_vec.erase(end-1, end);
#else
        edge_vec.erase(itr, itr+1);
#endif
          erased_edge = true;
          ++count;
          break;
        }
      }
      if (!erased_edge) break;
    }

    if (edge_vec.size() == 0) {
      m_map_table.erase(value);
    }

    m_num_edges -= count;
    return count;
  }

  inline vertex_property_data_type& vertex_property_data(const vertex_type& vertex)
  {
    return m_map_table.find(vertex)->second.first;
  }

  // Returns the associated edge data given a source and the sources edge.
  edge_property_data_type& edge_property_data(const vertex_type& src, const vertex_type& trg) {
    auto value = m_map_table.find(src);
    edge_vec_type& edge_vec = value->second.second;

    if (value == m_map_table.end()) {
      std::cout << src << ":" << trg;
      std::cout << " Bad edge: vertex lookup failed. ";  exit(-1);
    }

    for (auto itr = edge_vec.begin(); itr != edge_vec.end(); itr++) {
      if (itr->first == trg) {
        return itr->second;
      }
    }
    std::cout << src << ":" << trg;
    std::cout << " Bad edge: could not find destination. ";  exit(-1);
    return edge_vec.end()->second;
  }


  void opt()
  { }

  inline size_type num_edges() const
  {
    return m_num_edges;
  }

  void clear()
  { }

  void print_status(const int level) const
  { }

  void fprint_all_elements(std::ofstream& of)
  { }

  allocator_type m_allocator;
  map_table_type m_map_table;
  size_type m_num_edges;
};

}

#include <havoqgt/graphstore/baseline/baseline_vertex_iterator.hpp>
#include <havoqgt/graphstore/baseline/baseline_adjacent_edge_iterator.hpp>

#endif // GRAPHSTORE_BASELINE_HPP

