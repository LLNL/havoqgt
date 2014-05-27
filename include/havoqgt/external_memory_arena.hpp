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


#ifndef HAVOQGT_OMP_EXTERNAL_MEMORY_ARENA_HPP_INCLUDED
#define HAVOQGT_OMP_EXTERNAL_MEMORY_ARENA_HPP_INCLUDED

#include <omp.hpp>
#include <vector>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

namespace havoqgt { namespace omp { 

class external_memory_arena
{
public:
  typedef boost::interprocess::managed_mapped_file::segment_manager segment_manager;

  external_memory_arena(const char* fname, uint64_t fsize_bytes)
    : m_vec_ptr( num_threads ())
  {
    namespace bip = boost::interprocess;
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      std::stringstream ssfname;
      ssfname << fname << "_" << thread_num();
      m_vec_ptr[thread_num()] = new bip::managed_mapped_file(bip::create_only, ssfname.str().c_str(), fsize_bytes);
    }
  }
  external_memory_arena(const char* fname)
    : m_vec_ptr( num_threads ())
  {
    namespace bip = boost::interprocess;
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      std::stringstream ssfname;
      ssfname << fname << "_" << thread_num();
      m_vec_ptr[thread_num()] = new bip::managed_mapped_file(bip::open_only, ssfname.str().c_str());
    }
  }
  ~external_memory_arena() 
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      delete m_vec_ptr[thread_num()];
    }
  }
  template <typename T>
  struct allocator {
    typedef typename boost::interprocess::allocator<T,boost::interprocess::managed_mapped_file::segment_manager> type;
  };

  template<typename T>
  typename allocator<T>::type make_allocator() {
    return typename allocator<T>::type(get_sm());
  }

  boost::interprocess::managed_mapped_file::segment_manager* get_sm()
  {
    return m_vec_ptr[thread_num()]->get_segment_manager();
  }

  static void remove(const char* fname)
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      std::stringstream ssfname;
      ssfname << fname << "_" << thread_num();
      ::remove(ssfname.str().c_str());
    }
  }

private:
  /**
   * Not Copyable
   */
  external_memory_arena(const external_memory_arena&);

  std::vector<boost::interprocess::managed_mapped_file*> m_vec_ptr;
};

}}

#endif //HAVOQGT_OMP_EXTERNAL_MEMORY_ARENA_HPP_INCLUDED