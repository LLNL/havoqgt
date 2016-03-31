#ifndef BASELINE_MAP_ADJACENT_EDGE_ITERATOR_HPP
#define BASELINE_MAP_ADJACENT_EDGE_ITERATOR_HPP

#include <havoqgt/graphstore/baseline/baseline_map.hpp>

namespace graphstore {

template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class graphstore_baseline_map<vertex_type,
                              vertex_property_data_type,
                              edge_property_data_type,
                              segment_manager_type>::adjacent_edge_iterator
{
 private:
  using graphstore_type    = graphstore_baseline_map<vertex_type,
                                                     vertex_property_data_type,
                                                     edge_property_data_type,
                                                     segment_manager_type>;
  using self_type          = graphstore_type::adjacent_edge_iterator;
  using edge_iterator_type = typename graphstore_type::edge_map_table_type::iterator;

 public:

  explicit adjacent_edge_iterator (edge_iterator_type&& iterator) :
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

  edge_property_data_type& property_data()
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

#endif // BASELINE_MAP_ADJACENT_EDGE_ITERATOR_HPP

