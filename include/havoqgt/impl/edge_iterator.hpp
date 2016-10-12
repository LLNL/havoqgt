
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

#ifndef HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_
#define HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {
namespace mpi {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::edge_iterator {
 public:
  edge_iterator()
    : m_ptr_graph(NULL) {};

  edge_iterator& operator++();
  edge_iterator operator++(int);

  bool is_equal(const edge_iterator& x) const;

  friend bool operator==(const edge_iterator& x,
                         const edge_iterator& y) {return x.is_equal(y); }

  friend bool operator!=(const edge_iterator& x,
                         const edge_iterator& y) {return !(x.is_equal(y)); }

  vertex_locator source() const { return m_source; }
  vertex_locator target() const;
  edge_data_type edge_data() const;  

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

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_data_type
delegate_partitioned_graph<SegmentManager>::edge_iterator::edge_data() const {
  if(m_source.is_delegate()) {
    assert(m_edge_offset < m_ptr_graph->m_delegate_targets_size);
    assert(m_ptr_graph->m_delegate_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
    //return m_ptr_graph->m_delegate_targets[m_edge_offset];
    return m_ptr_graph->m_edge_data.m_delegate_edge_data[m_edge_offset]; 
  }
  assert(m_edge_offset < m_ptr_graph->m_owned_targets_size);
  assert(m_ptr_graph->m_owned_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
  //return m_ptr_graph->m_owned_targets[m_edge_offset];
  return m_ptr_graph->m_edge_data.m_owned_edge_data[m_edge_offset];
}

}  // mpi
}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_EDGE_ITERATOR_HPP_
