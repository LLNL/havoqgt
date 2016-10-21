/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Re-written by Steven Feldman <feldman12@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_
#define HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {
namespace mpi {

template <typename SegementManager>
template <typename T, typename Allocator>
class delegate_partitioned_graph<SegementManager>::edge_data {
 public:
   typedef typename bip::vector< T, Allocator >::iterator iterator;
  typedef T value_type;

  edge_data() {}

  T&       operator[] (const edge_iterator& itr) {
    if(itr.m_source.is_delegate()) {
      assert(itr.m_edge_offset < m_delegate_edge_data.size());
      return m_delegate_edge_data[itr.m_edge_offset];
    }
    assert(itr.m_edge_offset < m_owned_edge_data.size());
    return m_owned_edge_data[itr.m_edge_offset];
  }
  const T& operator[] (const edge_iterator& itr) const {
    if(itr.m_source.is_delegate()) {
      assert(itr.m_edge_offset < m_delegate_edge_data.size());
      return m_delegate_edge_data[itr.m_edge_offset];
    }
    assert(itr.m_edge_offset < m_owned_edge_data.size());
    return m_owned_edge_data[itr.m_edge_offset];
  }

  void reset(const T& r) {
    for(size_t i=0; i<m_owned_edge_data.size(); ++i) {
      m_owned_edge_data[i] = r;
    }
    for(size_t i=0; i<m_delegate_edge_data.size(); ++i) {
      m_delegate_edge_data[i] = r;
    }
  }

  iterator delegate_begin() { return m_delegate_edge_data.begin(); }
  iterator delegate_end()   { return m_delegate_edge_data.end(); }
  iterator owned_begin()    { return m_owned_edge_data.begin(); }
  iterator owned_end()      { return m_owned_edge_data.end(); }

//private:
  friend class delegate_partitioned_graph;
  edge_data(const delegate_partitioned_graph& dpg, Allocator allocate = Allocator() )
    : m_owned_edge_data(allocate)
    , m_delegate_edge_data(allocate) {
    //m_owned_edge_data.resize(dpg.m_owned_targets_size);
    //m_delegate_edge_data.resize(dpg.m_delegate_targets_size);
    resize(dpg);
  }

  edge_data(Allocator allocate = Allocator() )
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
  bip::vector< T, Allocator >      m_owned_edge_data;
  bip::vector< T, Allocator >      m_delegate_edge_data;
};

////////////////////////////////////////////////////////////////////////////////
//                                edge_data                                 //
////////////////////////////////////////////////////////////////////////////////
// template <typename SegmentManager>
// template<typename T, typename Allocator>
// delegate_partitioned_graph<SegmentManager>::edge_data<T,Allocator>::
// edge_data(const delegate_partitioned_graph& dpg, Allocator allocate = Allocator() )
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
// edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
//   : m_owned_edge_data(owned_size, init, sm->template get_allocator<T>())
//   , m_delegate_edge_data(delegate_size, init, sm->template get_allocator<T>()) { }

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
// delegate_partitioned_graph<SegmentManager>::edge_data<T, Allocator>::operator[](const edge_iterator& itr) const {
//   if(itr.m_source.is_delegate()) {
//     assert(itr.m_edge_offset < m_delegate_edge_data.size());
//     return m_delegate_edge_data[itr.m_edge_offset];
//   }
//   assert(itr.m_edge_offset < m_owned_edge_data.size());
//   return m_owned_edge_data[itr.m_edge_offset];
// }

}  // mpi
}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_EDGE_DATA_HPP_
