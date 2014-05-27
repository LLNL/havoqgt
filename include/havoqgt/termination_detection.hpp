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


#ifndef HAVOQGT_OMP_TERMINATION_DETECTION_HPP_INCLUDED
#define HAVOQGT_OMP_TERMINATION_DETECTION_HPP_INCLUDED

#include <omp.hpp>
#include <stdint.h>

namespace havoqgt { namespace omp {

/**
 * @todo Replace thread local storage with pointers to private memory.  
 *       Currently only supports having a single termination detection at a time.
 */

class termination_detection {
public:
  termination_detection() 
  {
    reset();
  }

  void reset() 
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      tp_sent()      = 0;
      tp_completed() = 0;
    }
    m_global_sent      = 0;
    m_global_completed = 0;
    m_threads_pending = 0;
  }

  void inc_sent(uint64_t inc=1) { check_new_pending(); tp_sent() += inc; }
  void inc_completed(uint64_t inc=1) { check_new_pending(); tp_completed() += inc;}

  bool test_for_termination() 
  {
    bool clear_pending = false;
    if(tp_sent() > 0) {
      clear_pending = true;
      #pragma omp atomic
        m_global_sent += tp_sent();
      tp_sent() = 0;
    }

    if(tp_completed() > 0) {
      clear_pending = true;
      #pragma omp atomic
        m_global_completed += tp_completed();
      tp_completed() = 0;
    }

    if(clear_pending) {
      #pragma omp atomic
      --m_threads_pending;
    }


    /*// first test
    __sync_synchronize();
    // this is a hack.  Really need to just make sure every thread has tested at lest once.
    if(m_global_completed == 0) return false;
    uint64_t first_grab = m_global_completed;
    if(m_global_completed != m_global_sent) return false;

    // start second test
    __sync_synchronize();
    if(m_global_completed == m_global_sent) {
      if(first_grab == m_global_completed) return true;
    }
    return false;
    */
    __sync_synchronize();
    return (m_global_completed == m_global_sent && m_threads_pending == 0 && m_global_completed > 0);
  }

  uint64_t global_completed() const { return m_global_completed; }

private:
  uint64_t& tp_sent() {
    static __thread uint64_t _tp_sent;
    return _tp_sent;
  }
  uint64_t& tp_completed() {
    static __thread uint64_t _tp_completed;
    return _tp_completed;
  }

  void check_new_pending() {
    if(tp_sent() == 0 && tp_completed() == 0) {
      #pragma omp atomic
      ++m_threads_pending;
    }
  }
  uint64_t m_global_sent      __attribute__ ((aligned (64)));
  uint64_t m_global_completed __attribute__ ((aligned (64)));
  uint64_t m_threads_pending __attribute__ ((aligned (64)));
};



} } //end havoqgt::omp

#endif //HAVOQGT_OMP_TERMINATION_DETECTION_HPP_INCLUDED
