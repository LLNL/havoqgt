// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_VERT_INFO_HPP_
#define HAVOQGT_MPI_IMPL_VERT_INFO_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::vert_info {
 public:
  vert_info(bool in_is_delegate, uint64_t in_delegate_id,
            uint64_t in_low_csr_idx);

  uint32_t is_delegate :  1;
  uint32_t delegate_id : 24;
  uint64_t low_csr_idx : 39;

  friend bool operator==(const vert_info& x, const vert_info& y){
    return  (x.is_delegate == y.is_delegate) &&
            (x.delegate_id == y.delegate_id) &&
            (x.low_csr_idx == y.low_csr_idx);
  }

  friend bool operator!=(const vert_info& x, const vert_info& y){
    return !(x == y);
  }
};

////////////////////////////////////////////////////////////////////////////////
//                                vert_info                                   //
////////////////////////////////////////////////////////////////////////////////

template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vert_info::
vert_info(bool in_is_delegate, uint64_t in_delegate_id, uint64_t in_low_csr_idx)
  : is_delegate(in_is_delegate)
  , delegate_id(in_delegate_id)
  , low_csr_idx(in_low_csr_idx) {
  assert(is_delegate == in_is_delegate);
  assert(delegate_id == in_delegate_id);
  assert(low_csr_idx == in_low_csr_idx);
  assert(sizeof(vert_info) == 8);
}

}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERT_INFO_HPP_
