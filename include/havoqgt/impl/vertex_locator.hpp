
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

#ifndef HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::vertex_locator {
 public:
  vertex_locator() {
    m_is_delegate  = 0;
    m_is_bcast     = 0;
    m_is_intercept = 0;
    #pragma GCC diagnostic ignored "-Woverflow"   /// NOTE:  is there a better way to clean these overflows?
    m_owner_dest   = std::numeric_limits<uint32_t>::max();
    m_local_id     = std::numeric_limits<uint64_t>::max();
    #pragma GCC diagnostic pop
  }

  bool is_valid() const {
    delegate_partitioned_graph<SegementManager>::vertex_locator conv;
    #pragma GCC diagnostic ignored "-Woverflow"   /// NOTE:  is there a better way to clean these overflows?
    conv.m_local_id = std::numeric_limits<uint64_t>::max();
    conv.m_owner_dest = std::numeric_limits<uint64_t>::max();
    #pragma GCC diagnostic pop


    return (m_local_id != conv.m_local_id || m_owner_dest != conv.m_owner_dest);
  }
  bool is_delegate_master() const {
     return (is_delegate() && ((m_local_id % havoqgt_env()->world_comm().size())
                               == havoqgt_env()->world_comm().rank()));
  }

  bool is_delegate() const { return m_is_delegate == 1;}
  uint32_t owner() const { return m_owner_dest; }
  void set_dest(uint32_t dest) { m_owner_dest = dest; assert(m_owner_dest == dest);}
  uint64_t local_id() const { return m_local_id;}
  bool is_equal(const vertex_locator x) const;
  uint32_t get_bcast() const { return m_is_bcast ; }
  void set_bcast(uint32_t bcast) { m_is_bcast = bcast; }
  bool is_intercept() const { return m_is_intercept == 1;}
  void set_intercept(bool intercept) { m_is_intercept = intercept; }

  friend bool operator==(const vertex_locator& x,
                         const vertex_locator& y) {return x.is_equal(y); }
  friend bool operator<(const vertex_locator& x,
                         const vertex_locator& y) {
    if (x.m_is_delegate == y.m_is_delegate) {
      if (x.m_owner_dest == y.m_owner_dest) {
        return x.m_local_id < y.m_local_id;
      }
      else {
        return x.m_owner_dest < y.m_owner_dest;
      }

    } else {
      return x.m_is_delegate < y.m_is_delegate;
    }
  }


  friend bool operator!=(const vertex_locator& x,
                         const vertex_locator& y) {return !(x.is_equal(y)); }

 private:
  friend class delegate_partitioned_graph;
  unsigned int m_is_delegate  : 1;

  unsigned int m_is_bcast     : 1;
  unsigned int m_is_intercept : 1;
  unsigned int m_owner_dest   : 20;
  uint64_t     m_local_id     : 39;

  vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest);
} __attribute__ ((packed)) ;


///////////////////////////////////////////////////////////////////////////////
//                           Vertex Locator                                  //
///////////////////////////////////////////////////////////////////////////////
/**
 * @class  delegate_partitioned_graph::vertex_locator
 * @details Here are some very important details.
 */
/**
 *
 */
template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vertex_locator::
vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest) {
  m_is_bcast     = 0;
  m_is_intercept = 0;

  if (is_delegate) {
    m_is_delegate = true;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    if(!(m_is_delegate == true
        && m_local_id    == local_id
        && m_owner_dest  == owner_dest)) { std::cerr << "ERROR:  vertex_locator()" << std::endl; exit(-1);}
  } else {
    m_is_delegate = false;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    if(!(m_is_delegate == false
        && m_owner_dest  == owner_dest
        && m_local_id    == local_id)) { std::cerr << "ERROR:  vertex_locator()" << std::endl; exit(-1);}

  }
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::vertex_locator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::vertex_locator x) const {
  return m_is_delegate  == x.m_is_delegate
      && m_is_bcast     == x.m_is_bcast
      && m_is_intercept == x.m_is_intercept
      && m_owner_dest   == x.m_owner_dest
      && m_local_id     == x.m_local_id;
}


}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_
