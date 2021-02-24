// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_DETAIL_ITERATOR_HPP_INCLUDED
#define HAVOQGT_MPI_DETAIL_ITERATOR_HPP_INCLUDED

#include <havoqgt/mpi.hpp>

namespace havoqgt { namespace detail {

/**
 * Tests if all processes' iterator range is empty
 * @param itr begin iterator
 * @param itr_end end iterator
 * @param comm MPI communicator
 * @return true, if all processes have empty range, else false
 */
template <typename Iterator>
bool global_iterator_range_empty(Iterator itr, Iterator itr_end, MPI_Comm comm) {
  uint32_t my_unfinished = itr != itr_end;
  uint32_t ranks_unfinished = mpi_all_reduce(my_unfinished, 
                                             std::plus<uint32_t>(), 
                                             comm);
  return ranks_unfinished == 0;
}

}}

#endif //HAVOQGT_MPI_DETAIL_ITERATOR_HPP_INCLUDED