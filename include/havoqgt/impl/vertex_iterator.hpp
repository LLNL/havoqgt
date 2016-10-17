
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

#ifndef HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {
namespace mpi {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::vertex_iterator
    : public std::iterator<std::input_iterator_tag, vertex_locator, ptrdiff_t,
                           const vertex_locator* const, const vertex_locator&> {
 public:
  vertex_iterator()
    : m_ptr_graph(NULL) {};
  vertex_iterator& operator++();
  vertex_iterator operator++(int);

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

}  // mpi
}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERTEX_ITERATOR_HPP_
