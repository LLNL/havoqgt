#ifndef GRAPHSTORE_BASELINE_HPP
#define GRAPHSTORE_BASELINE_HPP

#include <boost/unordered_map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>


namespace graphstore {


template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class graphstore_baseline
{
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
  /// iterator
  class vertex_iterator;
  class adjacent_edge_iterator;


public:

  explicit graphstore_baseline(segment_manager_type* segment_manager) :
    m_allocator(segment_manager),
    m_map_table(m_allocator)
  {}


  /// -------- Lookup -------- ///
  vertex_iterator vertices_begin()
  {
    return vertex_iterator(this);
  }

  static vertex_iterator vertices_end()
  {
    return vertex_iterator();
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

    return true;
  }

  size_type erase_edge(vertex_type& src, vertex_type& trg)
  {
    auto value = m_map_table.find(src);
    if (value == m_map_table.end()) {
      return 0;
    }

    size_t count = 0;
    edge_vec_type& edge_vec = value->second.second;
    for (auto itr = edge_vec.begin(), end = edge_vec.end(); itr != end; ++itr) {
      if (itr->first ==  trg) {
        edge_vec.erase(itr, itr+1);
        ++count;
      }
    }

    if (edge_vec.size() == 0) {
      m_map_table.erase(value);
    }

    return count;
  }

  inline vertex_property_data_type& vertex_property_data(const vertex_type& vertex)
  {
    return m_map_table.find(vertex)->second.first;
  }

//  typename edge_vec_type::iterator adjacencylist(const vertex_type& src)
//  {
//    const auto itr = m_map_table.find(src);
//    edge_vec_type& edge_vec = itr->second.second;
//    return edge_vec.begin();
//  }

//  typename edge_vec_type::iterator adjacencylist_end(const vertex_type& src)
//  {
//    const auto itr = m_map_table.find(src);
//    edge_vec_type& edge_vec = itr->second.second;
//    return edge_vec.end();
//  }


  void opt()
  {

  }

  void clear()
  {

  }

  void print_status(const int level) const
  {

  }

  void fprint_all_elements(std::ofstream& of)
  {
  }

  allocator_type m_allocator;
  map_table_type m_map_table;
};


template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class graphstore_baseline<vertex_type,
                          vertex_property_data_type,
                          edge_property_data_type,
                          segment_manager_type>::vertex_iterator
{
 private:
  using graphstore_type       = graphstore_baseline<vertex_type,
                                                    vertex_property_data_type,
                                                    edge_property_data_type,
                                                    segment_manager_type>;
  using self_type             = graphstore_type::vertex_iterator;
  using table_vertex_iterator = typename graphstore_type::map_table_type::iterator;

 public:

  vertex_iterator () :
    m_iterator(graphstore_type::map_table_type::end())
  { }


  vertex_iterator(graphstore_type* gstore) :
    m_iterator(gstore->m_map_table.begin())
  { }

  void swap(self_type &other) noexcept
  {
    using std::swap;
    swap(m_iterator, other.m_iterator);
  }

  self_type &operator++ () // Pre-increment
  {
    find_next_value();
    return *this;
  }

  self_type operator++ (int) // Post-increment
  {
    self_type tmp(*this);
    find_next_value();
    return tmp;
  }

  // two-way comparison: v.begin() == v.cbegin() and vice versa
  bool operator == (const self_type &rhs) const
  {
    return is_equal(rhs);
  }

  bool operator != (const self_type &rhs) const
  {
    return !is_equal(rhs);
  }

  const vertex_type& source_vertex()
  {
    return m_iterator.first;
  }

  const vertex_property_data_type& property_data()
  {
    return m_iterator.second.first;
  }


private:

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_iterator == rhs.m_iterator);
  }

  inline void find_next_value()
  {
    ++m_iterator;
  }

  table_vertex_iterator m_iterator;
};


template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class graphstore_baseline<vertex_type,
                          vertex_property_data_type,
                          edge_property_data_type,
                          segment_manager_type>::adjacent_edge_iterator
{
 private:
  using graphstore_type    = graphstore_baseline<vertex_type,
                                                     vertex_property_data_type,
                                                     edge_property_data_type,
                                                     segment_manager_type>;
  using self_type          = graphstore_type::adjacent_edge_iterator;
  using edge_iterator_type = typename graphstore_type::edge_vec_type::iterator;

 public:

  explicit adjacent_edge_iterator (edge_iterator_type iterator) :
    m_iterator(iterator)
  { }

  void swap(self_type &other) noexcept
  {
    using std::swap;
    swap(m_iterator, other.m_iterator);
  }

  self_type &operator++ () // Pre-increment
  {
    find_next_value();
    return *this;
  }

  self_type operator++ (int) // Post-increment
  {
    self_type tmp(*this);
    find_next_value();
    return tmp;
  }

  // two-way comparison: v.begin() == v.cbegin() and vice versa
  bool operator == (const self_type &rhs) const
  {
    return is_equal(rhs);
  }

  bool operator != (const self_type &rhs) const
  {
    return !is_equal(rhs);
  }

  const vertex_type& target_vertex()
  {
    return m_iterator->first;
  }

  const edge_property_data_type& property_data()
  {
    return m_iterator->second;
  }


private:

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_iterator == rhs.m_iterator);
  }

  inline void find_next_value()
  {
    ++m_iterator;
  }

  edge_iterator_type m_iterator;
};

}
#endif // GRAPHSTORE_BASELINE_HPP

