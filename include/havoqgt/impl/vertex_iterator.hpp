// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::vertex_iterator
    : public std::iterator<std::input_iterator_tag, vertex_locator, ptrdiff_t,
                           const vertex_locator* const, const vertex_locator&> {
 public:
  vertex_iterator()
    : m_ptr_graph(NULL) {};
  vertex_iterator& operator++();
  vertex_iterator operator++(int);

  vertex_iterator& operator--();
  vertex_iterator operator--(int);

  bool is_equal(const vertex_iterator& x) const;

  friend bool operator==(const vertex_iterator& x,
                         const vertex_iterator& y) { return x.is_equal(y); }

  friend bool operator!=(const vertex_iterator& x,
                         const vertex_iterator& y) { return !(x.is_equal(y)); }


  const vertex_locator& operator*()        const { return m_locator; }
  const vertex_locator* const operator->() const { return &m_locator; }

 private:
  friend class delegate_partitioned_graph;

  vertex_iterator(uint64_t index, const delegate_partitioned_graph*  pgraph);
  vertex_iterator(uint64_t index, const delegate_partitioned_graph*  pgraph, bool is_delegate); 
  void update_locator();

  const delegate_partitioned_graph*  m_ptr_graph;
  uint64_t                                m_owned_vert_index; // also, for delegate vertices
  vertex_locator                          m_locator;
  bool					  m_is_delegate; 
};



////////////////////////////////////////////////////////////////////////////////
//                             Vertex Iterator                                //
////////////////////////////////////////////////////////////////////////////////

template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
vertex_iterator(uint64_t index, const delegate_partitioned_graph<SegmentManager>*  pgraph)
  : m_ptr_graph(pgraph)
  , m_owned_vert_index(index) 
  , m_is_delegate(false) {
  update_locator();
}

template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
vertex_iterator(uint64_t index, const delegate_partitioned_graph<SegmentManager>*  pgraph, bool is_delegate) 
  : m_ptr_graph(pgraph)
  , m_owned_vert_index(index) 
  , m_is_delegate(is_delegate) {
  update_locator();
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator&
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator++() {
  ++m_owned_vert_index;
  update_locator();
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator++(int) {
  vertex_iterator to_return = *this;
  ++m_owned_vert_index;
  update_locator();
  return to_return;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator&
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator--() {
  assert(m_owned_vert_index > 0);
  --m_owned_vert_index;
  update_locator();
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator--(int) {
  assert(m_owned_vert_index > 0);
  vertex_iterator to_return = *this;
  --m_owned_vert_index;
  update_locator();
  return to_return;
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x) const {
  assert(m_ptr_graph        == x.m_ptr_graph);
  return m_owned_vert_index == x.m_owned_vert_index;
}

template <typename SegmentManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& y) {
  return x.is_equal(y);

}

template <typename SegmentManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegmentManager>
inline void
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
update_locator() {
  if (m_is_delegate) {
    if(m_owned_vert_index < m_ptr_graph->m_delegate_info.size()) {
      uint32_t owner = m_owned_vert_index % m_ptr_graph->m_mpi_size;
      m_locator = vertex_locator(true, m_owned_vert_index, owner);
    }  
  } else {
    for(; m_owned_vert_index < m_ptr_graph->m_owned_info.size()
          && m_ptr_graph->m_owned_info[m_owned_vert_index].is_delegate == true;
          ++ m_owned_vert_index);
    if(m_owned_vert_index < m_ptr_graph->m_owned_info.size()) {
      assert(m_ptr_graph->m_owned_info[m_owned_vert_index].is_delegate == false);
      uint32_t owner = m_ptr_graph->m_mpi_rank;
      m_locator = vertex_locator(false, m_owned_vert_index, owner);
    }
  }
}

}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_
