#ifndef DEGAWARERHH_ADJACENT_EDGE_ITERATOR_HPP
#define DEGAWARERHH_ADJACENT_EDGE_ITERATOR_HPP

#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>

namespace graphstore {

template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type,
          size_t middle_high_degree_threshold>
class degawarerhh<vertex_type,
                  vertex_property_data_type,
                  edge_property_data_type,
                  segment_manager_type,
                  middle_high_degree_threshold>::adjacent_edge_iterator
{

 private:
  using graphstore_type            = degawarerhh<vertex_type,
                                                            vertex_property_data_type,
                                                            edge_property_data_type,
                                                            segment_manager_type,
                                                            middle_high_degree_threshold>;
  using self_type                  = graphstore_type::adjacent_edge_iterator;
  using low_deg_edge_iterator_type = typename graphstore_type::low_deg_table_type::value_iterator;
  using mh_deg_edge_iterator_type  = typename graphstore_type::mh_edge_chunk_type::whole_iterator;


 public:

  adjacent_edge_iterator () :
    m_low_itr(),
    m_mh_itr()
  { }


  adjacent_edge_iterator (graphstore_type* gstore, const vertex_type& src_vrt) :
    m_low_itr(gstore->m_low_degree_table->find(src_vrt)),
    m_mh_itr(mh_adjacent_edge_begin(gstore, src_vrt))
  { }


  void swap(self_type &other) noexcept
  {
    using std::swap;
    swap(m_low_itr, other.m_low_itr);
    swap(m_mh_itr, other.m_mh_itr);
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

  /// TODO: handle an error when m_mh_itr.is_end() == true
  const vertex_type& target_vertex()
  {
    if (!m_low_itr.is_end()) {
      return m_low_itr->second;
    }
    return m_mh_itr->key;
  }

  /// TODO: handle an error when m_mh_itr.is_end() == true
  const edge_property_data_type& property_data()
  {
    if (!m_low_itr.is_end()) {
      return m_low_itr->third;
    }
    return m_mh_itr->value;
  }


 private:

  inline static mh_deg_edge_iterator_type mh_adjacent_edge_begin (graphstore_type* graphstore, const vertex_type& src_vrt)
  {
    const auto itr_matrix = graphstore->m_mh_degree_table->find(src_vrt);
    if (!itr_matrix.is_end()) {
      typename graphstore_type::mh_edge_chunk_type* adj_list = itr_matrix->second;
      return adj_list->begin();
    } else {
      return graphstore_type::mh_edge_chunk_type::end();
    }
  }

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_low_itr == rhs.m_low_itr) && (m_mh_itr == rhs.m_mh_itr);
  }

  inline void find_next_value()
  {
    if (!m_low_itr.is_end()) {
      ++m_low_itr;
    } else if (!m_mh_itr.is_end()) {
      ++m_mh_itr;
    }
  }

  low_deg_edge_iterator_type m_low_itr;
  mh_deg_edge_iterator_type m_mh_itr;
};

}
#endif // DEGAWARERHH_ADJACENT_EDGE_ITERATOR_HPP

