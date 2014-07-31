
/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Steven Feldman <feldman12@llnl.gov>.
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

#ifndef _HAVOQGT_CACHE_UTILI_HPP
#define _HAVOQGT_CACHE_UTILI_HPP

// #include <boost/interprocess/managed_mapped_file.hpp>
// #include <boost/interprocess/allocators/allocator.hpp>



template<typename Vector>
void advise_vector_rand(Vector vec) {
  void * addr = vec.data();
  size_t length = vec.size() * sizeof(vec[0]);
  assert(madvise(addr, length, MADV_RANDOM) == 0);
}


template<typename Vector>
void flush_advise_vector_dont_need(Vector vec) {
  void * addr = vec.data();
  size_t length = vec.size() * sizeof(vec[0]);
  assert(msync(addr, length, MS_SYNC) == 0);
  assert(madvise(addr, length, MADV_DONTNEED) == 0);
}


template<typename Vector>
void flush_vector(Vector vec) {
  void * addr = vec.data();
  size_t length = vec.size() * sizeof(vec[0]);
  assert(msync(addr, length, MS_SYNC) == 0);
}

template<typename Vector>
void flush_advise_vector(Vector vec) {
  void * addr = vec.data();
  size_t length = vec.size() * sizeof(vec[0]);
  assert(msync(addr, length, MS_SYNC) == 0);
  assert(madvise(addr, length, MADV_DONTNEED) == 0);
  assert(madvise(addr, length, MADV_RANDOM) == 0);
}


// template<typename mapped_t>
// void custom_flush(mapped_t * mapped) {

//   mapped->flush();

//   boost::interprocess::mapped_region::advice_types advise;
//   advise = boost::interprocess::mapped_region::advice_types::advice_dontneed;
//   assert(mapped->advise(advise));

//   advise = boost::interprocess::mapped_region::advice_types::advice_random;
//   assert(mapped->advise(advise));

// }



// template <typename SegmentManager>
// void
// delegate_partitioned_graph<SegmentManager>::
// try_flush(MPI_Comm comm) {
// #if 1
//   static uint64_t check_id = 0;

//   if ((m_mpi_rank  % 24) == (check_id++) % 24) {

//     uint32_t dirty_kb;
//     {
//       FILE *pipe;
//       pipe = popen("grep Dirty /proc/meminfo | awk '{print $2}'", "r" );
//       fscanf(pipe, "%u", &dirty_kb);
//       pclose(pipe);
//     }
//     const uint32_t dirty_threshold_kb = DIRTY_THRESHOLD_GB * 1000000;

//     if (dirty_kb > dirty_threshold_kb) {
//       m_flush_func();
//     }
//   }


// #else
//   bool do_flush;
//   if (m_mpi_rank == 0) {
//     uint32_t dirty_kb;
//     {
//       FILE *pipe;
//       pipe = popen("grep Dirty /proc/meminfo | awk '{print $2}'", "r" );
//       fscanf(pipe, "%u", &dirty_kb);
//       pclose(pipe);
//     }
//     const uint32_t dirty_threshold_kb = DIRTY_THRESHOLD_GB * 1000000;
//     do_flush = (dirty_kb > dirty_threshold_kb);
//   }

//   MPI_Bcast(&do_flush, 1, mpi_typeof(do_flush), 0, comm);


//   if (do_flush) {
//     m_flush_func();
//   }
// #endif
// }
#endif  // _HAVOQGT_CACHE_UTILI_HPP

