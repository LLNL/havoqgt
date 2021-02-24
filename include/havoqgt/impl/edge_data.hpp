// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_
#define HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename Allocator>
template <typename T, typename AllocatorOther>
class delegate_partitioned_graph<Allocator>::edge_data {
 public:
   typedef typename bip::vector< T, other_allocator<AllocatorOther, T>>::iterator iterator;
  typedef T value_type;

  T& operator[](const edge_iterator& itr) {
    if (itr.m_source.is_delegate()) {
      assert(itr.m_edge_offset < m_delegate_edge_data.size());
      return m_delegate_edge_data[itr.m_edge_offset];
    }
    assert(itr.m_edge_offset < m_owned_edge_data.size());
    return m_owned_edge_data[itr.m_edge_offset];
  }
  const T& operator[](const edge_iterator& itr) const {
    if (itr.m_source.is_delegate()) {
      assert(itr.m_edge_offset < m_delegate_edge_data.size());
      return m_delegate_edge_data[itr.m_edge_offset];
    }
    assert(itr.m_edge_offset < m_owned_edge_data.size());
    return m_owned_edge_data[itr.m_edge_offset];
  }

  void reset(const T& r) {
    for (size_t i = 0; i < m_owned_edge_data.size(); ++i) {
      m_owned_edge_data[i] = r;
    }
    for (size_t i = 0; i < m_delegate_edge_data.size(); ++i) {
      m_delegate_edge_data[i] = r;
    }
  }

  iterator delegate_begin() { return m_delegate_edge_data.begin(); }
  iterator delegate_end() { return m_delegate_edge_data.end(); }
  iterator owned_begin() { return m_owned_edge_data.begin(); }
  iterator owned_end() { return m_owned_edge_data.end(); }

  // private:
  friend class delegate_partitioned_graph;
  edge_data(const delegate_partitioned_graph& dpg, AllocatorOther allocate = AllocatorOther() )
    : m_owned_edge_data(allocate)
    , m_delegate_edge_data(allocate) {
    //m_owned_edge_data.resize(dpg.m_owned_targets_size);
    //m_delegate_edge_data.resize(dpg.m_delegate_targets_size);
    resize(dpg);
  }

  edge_data(AllocatorOther allocate = AllocatorOther() )
    : m_owned_edge_data(allocate)
    , m_delegate_edge_data(allocate) {}

  void resize(const delegate_partitioned_graph& dpg) {
    m_owned_edge_data.resize(dpg.m_owned_targets_size);
    m_delegate_edge_data.resize(dpg.m_delegate_targets_size);
  }

  // edge_data(uint64_t owned_size, uint64_t delegate_size,
  //     SegManagerOther* sm);
  // edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init,
  //     SegManagerOther* sm);

// private:
 protected:
  bip::vector< T, other_allocator<AllocatorOther, T>>      m_owned_edge_data;
  bip::vector< T, other_allocator<AllocatorOther, T>>      m_delegate_edge_data;
};

////////////////////////////////////////////////////////////////////////////////
//                                edge_data                                 //
////////////////////////////////////////////////////////////////////////////////
// template <typename SegmentManager>
// template<typename T, typename Allocator>
// delegate_partitioned_graph<SegmentManager>::edge_data<T,Allocator>::
// edge_data(const delegate_partitioned_graph& dpg, Allocator allocate =
// Allocator() )
//   : m_owned_edge_data(allocate)
//   , m_delegate_edge_data(allocate) {
//   m_owned_vert_data.resize(dpg.m_owned_targets.size());
//   m_delegate_data.resize(dpg.m_delegate_targets.size());
//   }
//
//
// template <typename SegmentManager>
// template<typename T, typename Allocator>
// delegate_partitioned_graph<SegmentManager>::edge_data<T,Allocator>::
// edge_data(uint64_t owned_size, uint64_t delegate_size, SegManagerOther* sm)
//   : m_owned_edge_data(sm->template get_allocator<T>())
//   , m_delegate_edge_data(sm->template get_allocator<T>()) {
//   m_owned_edge_data.resize(owned_size);
//   m_delegate_edge_data.resize(delegate_size);
//   }
//
// template <typename SegmentManager>
// template<typename T, typename Allocator>
// delegate_partitioned_graph<SegmentManager>::edge_data<T, Allocator>::
// edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init,
// SegManagerOther* sm)
//   : m_owned_edge_data(owned_size, init, sm->template get_allocator<T>())
//   , m_delegate_edge_data(delegate_size, init, sm->template
//   get_allocator<T>()) { }

// template <typename SegmentManager>
// template<typename T, typename Allocator>
// T&
// delegate_partitioned_graph<SegmentManager>::edge_data<T, Allocator>::
// operator[](const edge_iterator& itr) {
//   if(itr.m_source.is_delegate()) {
//     assert(itr.m_edge_offset < m_delegate_edge_data.size());
//     return m_delegate_edge_data[itr.m_edge_offset];
//   }
//   assert(itr.m_edge_offset < m_owned_edge_data.size());
//   return m_owned_edge_data[itr.m_edge_offset];
// }
//
// template <typename SegmentManager>
// template<typename T, typename Allocator>
// const T&
// delegate_partitioned_graph<SegmentManager>::edge_data<T,
// Allocator>::operator[](const edge_iterator& itr) const {
//   if(itr.m_source.is_delegate()) {
//     assert(itr.m_edge_offset < m_delegate_edge_data.size());
//     return m_delegate_edge_data[itr.m_edge_offset];
//   }
//   assert(itr.m_edge_offset < m_owned_edge_data.size());
//   return m_owned_edge_data[itr.m_edge_offset];
// }

}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_
