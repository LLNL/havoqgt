// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_
#define HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::edge_iterator {
 public:
  edge_iterator()
    : m_ptr_graph(NULL) {};

  edge_iterator& operator++();
  edge_iterator operator++(int);

  edge_iterator& operator+=(uint64_t offset);
  edge_iterator operator+(uint64_t offset);

  bool is_equal(const edge_iterator& x) const;

  friend bool operator==(const edge_iterator& x,
                         const edge_iterator& y) {return x.is_equal(y); }

  friend bool operator!=(const edge_iterator& x,
                         const edge_iterator& y) {return !(x.is_equal(y)); }

  vertex_locator source() const { return m_source; }
  vertex_locator target() const;
  //edge_data_type edge_data() const;  

 protected:
  friend class delegate_partitioned_graph;
  template <typename T1, typename T2> friend class edge_data;
  edge_iterator(vertex_locator source, uint64_t edge_offset,
                const delegate_partitioned_graph* const pgraph);

  vertex_locator                          m_source;
  uint64_t                                m_edge_offset;
  const delegate_partitioned_graph* const m_ptr_graph;
};


////////////////////////////////////////////////////////////////////////////////
//                               Edge Iterator                                //
////////////////////////////////////////////////////////////////////////////////
/**
 * \class delegate_partitioned_graph::edge_iterator
 * \details Put details here for class
 */
/**
 * @
 */
template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::edge_iterator::
edge_iterator(vertex_locator source,
              uint64_t edge_offset,
              const delegate_partitioned_graph* const pgraph)
  : m_source(source)
  , m_edge_offset(edge_offset)
  , m_ptr_graph(pgraph) { }

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator&
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator++() {
  ++m_edge_offset;
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator++(int) {
  edge_iterator to_return = *this;
  ++m_edge_offset;
  return to_return;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator&
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator+=(const uint64_t offset) {
  m_edge_offset += offset;
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator+(const uint64_t offset) {
  edge_iterator to_return = *this;
  to_return += offset;
  return to_return;
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::edge_iterator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x) const {
    assert(m_source      == x.m_source);
    assert(m_ptr_graph   == x.m_ptr_graph);
    return m_edge_offset == x.m_edge_offset;
}

template <typename SegmentManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& y) {
  return x.is_equal(y);

}

template <typename SegmentManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_locator
delegate_partitioned_graph<SegmentManager>::edge_iterator::target() const {
  if(m_source.is_delegate()) {
    assert(m_edge_offset < m_ptr_graph->m_delegate_targets_size);
    assert(m_ptr_graph->m_delegate_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
    return m_ptr_graph->m_delegate_targets[m_edge_offset];
  }
  assert(m_edge_offset < m_ptr_graph->m_owned_targets_size);
  assert(m_ptr_graph->m_owned_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
  return m_ptr_graph->m_owned_targets[m_edge_offset];
}

//template <typename SegmentManager>
//inline
//typename delegate_partitioned_graph<SegmentManager>::edge_data_type
//delegate_partitioned_graph<SegmentManager>::edge_iterator::edge_data() const {
//  if(m_source.is_delegate()) {
//    assert(m_edge_offset < m_ptr_graph->m_delegate_targets_size);
//    assert(m_ptr_graph->m_delegate_targets[m_edge_offset].m_owner_dest <
//          m_ptr_graph->m_mpi_size);
//    return m_ptr_graph->m_edge_data.m_delegate_edge_data[m_edge_offset]; 
//  }
//  assert(m_edge_offset < m_ptr_graph->m_owned_targets_size);
//  assert(m_ptr_graph->m_owned_targets[m_edge_offset].m_owner_dest <
//          m_ptr_graph->m_mpi_size);
//  return m_ptr_graph->m_edge_data.m_owned_edge_data[m_edge_offset];
//`}

}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_
