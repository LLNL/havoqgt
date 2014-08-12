
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

#ifndef HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {
namespace mpi {

template <typename SegementManager>
template <typename T, typename SegManagerOther>
class delegate_partitioned_graph<SegementManager>::vertex_data {
 public:
  typedef T value_type;
  vertex_data() {}

  T&       operator[] (const vertex_locator& locator);
  const T& operator[] (const vertex_locator& locator) const;

  void reset(const T& r) {
    for(size_t i=0; i<m_owned_vert_data.size(); ++i) {
      m_owned_vert_data[i] = r;
    }
    for(size_t i=0; i<m_delegate_data.size(); ++i) {
      m_delegate_data[i] = r;
    }
  }

  void all_reduce() {
    std::vector<T> tmp_in(m_delegate_data.begin(), m_delegate_data.end());
    std::vector<T> tmp_out(tmp_in.size(), 0);
    mpi_all_reduce(tmp_in, tmp_out, std::plus<T>(), MPI_COMM_WORLD);
    std::copy(tmp_out.begin(), tmp_out.end(), m_delegate_data.begin());
  }

  vertex_data(uint64_t owned_data_size, uint64_t delegate_size,
    SegManagerOther* segment_manager);
  vertex_data(uint64_t owned_data_size, uint64_t delegate_size,
    const T& init, SegManagerOther* segment_manager);

 private:
  bip::vector<T, bip::allocator<T, SegManagerOther>  > m_owned_vert_data;
  bip::vector<T, bip::allocator<T, SegManagerOther>  > m_delegate_data;
};

////////////////////////////////////////////////////////////////////////////////
//                                vertex_data                                 //
////////////////////////////////////////////////////////////////////////////////
template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::vertex_data<T,SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, SegManagerOther* sm)
  : m_owned_vert_data(sm->template get_allocator<T>())
  , m_delegate_data(sm->template get_allocator<T>()) {
  m_owned_vert_data.resize(owned_data_size);
  m_delegate_data.resize(delegate_size);
  }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
  : m_owned_vert_data(owned_data_size, init, sm->template get_allocator<T>())
  , m_delegate_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
T&
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::
operator[](const vertex_locator& locator) {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
const T&
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::operator[](const vertex_locator& locator) const {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

}  // mpi
}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_
