/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. 
 * Written by Roger Pearce <rpearce@llnl.gov>. 
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


#ifndef HAVOQGT_OMP_OMP_HPP_INCLUDED
#define HAVOQGT_OMP_OMP_HPP_INCLUDED


#include <omp.h>
#include <detail/omp.hpp>
#include <assert.h>


namespace havoqgt { namespace omp {

/** 
 * Initializes the parallel environment for havoqgt.
 * Typically called early in main(). 
 */
inline void init_environment() {
  // Start parallel region and get thread information.
  #pragma omp parallel
  {
    detail::__tls_thread_num() = omp_get_thread_num();
    detail::__num_threads() = omp_get_num_threads();
    assert(detail::__tls_thread_num() == omp_get_thread_num());
    assert(detail::__num_threads() == omp_get_num_threads());
  }
}

/**
 * Returns the threads identifier.
 */
inline int thread_num() {
  return detail::__tls_thread_num();
}

/**
 * Returns the number of threads that havoqgt uses.
 */
inline int num_threads() {
  return detail::__num_threads();
}

/**
 * Asserts that current execution is in a sequential region.
 */
inline void assert_sequential() {
  assert(omp_in_parallel() == 0);
}

/**
 * Asserts that current execution is in a parallel region.
 */
inline void assert_parallel() {
  assert(omp_in_parallel() == 1);
  assert(thread_num() == omp_get_thread_num());
  assert(num_threads() == omp_get_num_threads());
}


} } //end havoqgt::omp


#endif //HAVOQGT_OMP_OMP_HPP_INCLUDED