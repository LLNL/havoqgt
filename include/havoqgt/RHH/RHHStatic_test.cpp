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


/// How to use
/// $ g++ std=c++11 RHHStatic_test.cpp
/// $ ./a.out

#include <iostream>
#include <cstdint>
#include <cassert>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/adaptive_pool.hpp>

#include "RHHUtility.hpp"
#include "RHHMgrStatic.hpp"

#define DEBUG(msg) do { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << msg << std::endl; } while (0)
#define DEBUG2(x) do  { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << #x << " =\t" << x << std::endl; } while (0)
#define DISP_VAR(x) do  { std::cout << #x << " =\t" << x << std::endl; } while (0)


int main (void)
{
  
  boost::interprocess::managed_shared_memory segment(create_only,
                                "MySharedMemory",  //segment name
                                65536);
  
  RHH::AllocatorsHolder *holder = new RHH::AllocatorsHolder(segment.get_segment_manager());
  RHH::RHHMgrStatic<uint64_t, RHH::NoValueType> *rhh = new RHH::RHHMgrStatic<uint64_t, RHH::NoValueType>();

  std::cout << "-------insert--------\n";
  for (uint64_t i = 0; i < 64; i++) {
    std::cout << i << ": ";
    rhh->insert_uniquely_static(holder, i, NULL, i);
  }

  
//  std::cout << "-------allocate chained RHHStatic--------\n";
//  RHH::RHHStatic<uint64_t, unsigned char, 256> *rhh2 = new RHH::RHHStatic<uint64_t, unsigned char, 256>();
//  rhh2->m_next_ = rhh;
//  std::cout << "-------insert--------\n";
//  for (uint64_t i = 0; i < 64; i++) {
//    std::cout << i << ": ";
//    insertion_result(rhh2->insert_uniquely(i, 0));
//  }
//  std::cout << "--------erase-------\n";
//  for (uint64_t i = 0; i < 64; i++) {
//    std::cout << i << ": ";
//    erasion_result(rhh2->erase(i));
//  }
//  
//  std::cout << "-------allocate chained RHHStatic--------\n";
//  RHH::RHHStatic<uint64_t, unsigned char, 256> *rhh3 = new RHH::RHHStatic<uint64_t, unsigned char, 256>();
//  rhh3->m_next_ = rhh2;
//  std::cout << "------insert---------\n";
//  for (uint64_t i = 0; i < 64; i++) {
//    std::cout << i << ": ";
//    insertion_result(rhh3->insert_uniquely(i, 0));
//  }
//  std::cout << "--------erase-------\n";
//  for (uint64_t i = 0; i < 64; i++) {
//    std::cout << i << ": ";
//    erasion_result(rhh3->erase(i));
//  }
  
  delete rhh;
//  delete rhh2;
//  delete rhh3;
  
  return 0;
}