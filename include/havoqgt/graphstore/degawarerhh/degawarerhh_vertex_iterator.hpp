#ifndef DEGAWARERHH_VERTEX_ITERATOR_HPP
#define DEGAWARERHH_VERTEX_ITERATOR_HPP

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
                  middle_high_degree_threshold>::vertex_iterator
{

 private:
  using graphstore_type             = degawarerhh<vertex_type,
                                                  vertex_property_data_type,
                                                  edge_property_data_type,
                                                  segment_manager_type,
                                                  middle_high_degree_threshold>;
  using self_type                   = graphstore_type::vertex_iterator;
  using ldeg_table_iterator_type = typename graphstore_type::ldeg_table_type::whole_iterator;
  using mhdeg_table_iterator_type  = typename graphstore_type::mhdeg_table_type::whole_iterator;

 public:

  /// ---- Constructors ----
  /// initialize to an 'end iterator'
  vertex_iterator() :
  m_low_itr(ldeg_table_iterator_type::end()),
  m_mh_itr(mhdeg_table_iterator_type::end())
  { }

  vertex_iterator (const ldeg_table_iterator_type& low_itr,
                   const mhdeg_table_iterator_type& mh_itr) :
    m_low_itr(low_itr),
    m_mh_itr(mh_itr)
  { }

  vertex_iterator (ldeg_table_iterator_type&& low_itr,
                   mhdeg_table_iterator_type&& mh_itr) :
    m_low_itr(std::move(low_itr)),
    m_mh_itr(std::move(mh_itr))
  { }

  /// Copy constructor
  vertex_iterator (const vertex_iterator& other) :
    m_low_itr(other.m_low_itr),
    m_mh_itr(other.m_mh_itr)
  { }

  /// Move constructor
  vertex_iterator (vertex_iterator&& other) :
    m_low_itr(std::move(other.m_low_itr)),
    m_mh_itr(std::move(other.m_mh_itr))
  { }

  /// Copy assignment operators
  vertex_iterator& operator=(const vertex_iterator& other)
  {
    m_low_itr = other.m_low_itr;
    m_mh_itr  = other.m_mh_itr;
    return *this;
  }

  /// Move assignment operators
  vertex_iterator& operator=(vertex_iterator&& other)
  {
    m_low_itr = std::move(other.m_low_itr);
    m_mh_itr  = std::move(other.m_mh_itr);
    return *this;
  }

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
  const vertex_type& source_vertex()
  {
    if (!m_low_itr.is_end()) {
      return m_low_itr->key;
    }
    return m_mh_itr->key;
  }

  /// TODO: handle an error when m_mh_itr.is_end() == true
  vertex_property_data_type& property_data()
  {
    if (!m_low_itr.is_end()) {
      return m_low_itr->value.first;
    }
    return m_mh_itr->value.first;
  }

//  static adjacent_edge_iterator end()
//  {
//    return adjacent_edge_iterator(ldeg_table_iterator_type::end(),
//                                  mhdeg_table_iterator_type::end());
//  }

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_low_itr == rhs.m_low_itr) && (m_mh_itr == rhs.m_mh_itr);
  }


 private:

  inline void find_next_value()
  {
    if (!m_low_itr.is_end()) {
      ++m_low_itr;
      while (!m_low_itr.is_end() && !m_low_itr->value.fourth)
        ++m_low_itr;
      return;
    }
    if (!m_mh_itr.is_end()) {
      ++m_mh_itr;
      return;
    }
  }

  ldeg_table_iterator_type m_low_itr;
  mhdeg_table_iterator_type m_mh_itr;
};



template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class degawarerhh<vertex_type,
                  vertex_property_data_type,
                  edge_property_data_type,
                  segment_manager_type,
                  1>::vertex_iterator
{

 private:
  using graphstore_type             = degawarerhh<vertex_type,
                                                  vertex_property_data_type,
                                                  edge_property_data_type,
                                                  segment_manager_type,
                                                  1>;
  using self_type                 = graphstore_type::vertex_iterator;
  using ldeg_table_iterator_type  = typename graphstore_type::ldeg_table_type::whole_iterator;
  using mhdeg_table_iterator_type = typename graphstore_type::mhdeg_table_type::whole_iterator;

 public:

  /// ---- Constructors ----
  /// initialize to an 'end iterator'
  vertex_iterator() :
    m_mh_itr(mhdeg_table_iterator_type::end())
  { }

  vertex_iterator (const mhdeg_table_iterator_type& mh_itr) :
    m_mh_itr(mh_itr)
  { }

  vertex_iterator (mhdeg_table_iterator_type&& mh_itr) :
    m_mh_itr(std::move(mh_itr))
  { }

  /// Copy constructor
  vertex_iterator (const vertex_iterator& other) :
    m_mh_itr(other.m_mh_itr)
  { }

  /// Move constructor
  vertex_iterator (vertex_iterator&& other) :
    m_mh_itr(std::move(other.m_mh_itr))
  { }

  /// Copy assignment operators
  vertex_iterator& operator=(const vertex_iterator& other)
  {
    m_mh_itr  = other.m_mh_itr;
    return *this;
  }

  /// Move assignment operators
  vertex_iterator& operator=(vertex_iterator&& other)
  {
    m_mh_itr  = std::move(other.m_mh_itr);
    return *this;
  }

  void swap(self_type &other) noexcept
  {
    using std::swap;
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
  const vertex_type& source_vertex()
  {
    return m_mh_itr->key;
  }

  /// TODO: handle an error when m_mh_itr.is_end() == true
  vertex_property_data_type& property_data()
  {
    return m_mh_itr->value.first;
  }

  static adjacent_edge_iterator end()
  {
    return adjacent_edge_iterator(mhdeg_table_iterator_type::end());
  }

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_mh_itr == rhs.m_mh_itr);
  }


 private:

  inline void find_next_value()
  {
    if (!m_mh_itr.is_end()) {
      ++m_mh_itr;
      return;
    }
  }

  mhdeg_table_iterator_type m_mh_itr;
};
}
#endif // DEGAWARERHH_VERTEX_ITERATOR_HPP

