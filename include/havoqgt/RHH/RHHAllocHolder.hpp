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
#ifndef RHHH_RHHALLOCHOLDER_HPP_INCLUDED
#define RHHH_RHHALLOCHOLDER_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include "RHHCommon.hpp"
#include "RHHStatic.hpp"

namespace RHH {

  namespace bip = boost::interprocess;
  typedef bip::managed_mapped_file  mapped_t;
  typedef mapped_t::segment_manager segment_manager_t;

  /// Memory size (Byte) = capacity * (  1  /// property block
  ///                                  + 8  /// key block
  ///                                  + 1) /// value block
  ///                                  + 8   /// next pointer
  enum {
    capacityNormalArray1 = 2ULL,
    capacityNormalArray2 = 4ULL,
    capacityNormalArray3 = 8ULL,
    capacityRHHStatic_1 = 8ULL,
    capacityRHHStatic_2 = 16ULL,
    capacityRHHStatic_3 = 32ULL,
    capacityRHHStatic_4 = 64ULL,
    capacityRHHStatic_5 = 128ULL,
    capacityRHHStatic_6 = 256ULL,
    capacityRHHStatic_7 = 512ULL,
    capacityRHHStatic_8 = 1024ULL,
    capacityRHHStatic_9 = 2048ULL,
    capacityRHHStatic_10 = 4096ULL,
    capacityRHHStatic_11 = 8192ULL,
    capacityRHHStatic_12 = 16384ULL,
    capacityRHHStatic_13 = 32768ULL,
    capacityRHHStatic_14 = 65536ULL,
    capacityRHHStatic_15 = 131072ULL,
    capacityRHHStatic_16 = 262144ULL,
    capacityRHHStatic_17 = 524288ULL
  };

  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_1> RHHStaticNoVal_1;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_2> RHHStaticNoVal_2;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_3> RHHStaticNoVal_3;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_4> RHHStaticNoVal_4;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_5> RHHStaticNoVal_5;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_6> RHHStaticNoVal_6;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_7> RHHStaticNoVal_7;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_8> RHHStaticNoVal_8;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_9> RHHStaticNoVal_9;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_10> RHHStaticNoVal_10;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_11> RHHStaticNoVal_11;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_12> RHHStaticNoVal_12;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_13> RHHStaticNoVal_13;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_14> RHHStaticNoVal_14;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_15> RHHStaticNoVal_15;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_16> RHHStaticNoVal_16;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_17> RHHStaticNoVal_17;

#if 0
  /// chuk size = 4MB
  typedef bip::node_allocator<RHHStaticNoVal_1, segment_manager_t,  47662>  Allocator_RHH_NoVal_1;
  typedef bip::node_allocator<RHHStaticNoVal_2, segment_manager_t,  24966>  Allocator_RHH_NoVal_2;
  typedef bip::node_allocator<RHHStaticNoVal_3, segment_manager_t,  12787>  Allocator_RHH_NoVal_3;
  typedef bip::node_allocator<RHHStaticNoVal_4, segment_manager_t,   6472>  Allocator_RHH_NoVal_4;
  typedef bip::node_allocator<RHHStaticNoVal_5, segment_manager_t,   3256>  Allocator_RHH_NoVal_5;
  typedef bip::node_allocator<RHHStaticNoVal_6, segment_manager_t,   1633>  Allocator_RHH_NoVal_6;
  typedef bip::node_allocator<RHHStaticNoVal_7, segment_manager_t,   817>  Allocator_RHH_NoVal_7;
  typedef bip::node_allocator<RHHStaticNoVal_8, segment_manager_t,   409>  Allocator_RHH_NoVal_8;
  typedef bip::node_allocator<RHHStaticNoVal_9, segment_manager_t,   204>  Allocator_RHH_NoVal_9;
  typedef bip::node_allocator<RHHStaticNoVal_10, segment_manager_t,  102>  Allocator_RHH_NoVal_10;
  typedef bip::node_allocator<RHHStaticNoVal_11, segment_manager_t,  51>  Allocator_RHH_NoVal_11;
  typedef bip::node_allocator<RHHStaticNoVal_12, segment_manager_t,  25>  Allocator_RHH_NoVal_12;
  typedef bip::node_allocator<RHHStaticNoVal_13, segment_manager_t,  12>  Allocator_RHH_NoVal_13;
  typedef bip::node_allocator<RHHStaticNoVal_14, segment_manager_t,  6>  Allocator_RHH_NoVal_14;
  typedef bip::node_allocator<RHHStaticNoVal_15, segment_manager_t,  3>  Allocator_RHH_NoVal_15;
  typedef bip::node_allocator<RHHStaticNoVal_16, segment_manager_t,  2>  Allocator_RHH_NoVal_16;
  typedef bip::node_allocator<RHHStaticNoVal_17, segment_manager_t,  1>  Allocator_RHH_NoVal_17;
  typedef bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_normalarray_t;
  typedef bip::allocator<unsigned char, segment_manager_t> allocator_raw_t;
#else
  /// chuk size = 4KB
  typedef bip::node_allocator<RHHStaticNoVal_1, segment_manager_t,  47>  Allocator_RHH_NoVal_1;
  typedef bip::node_allocator<RHHStaticNoVal_2, segment_manager_t,  24>  Allocator_RHH_NoVal_2;
  typedef bip::node_allocator<RHHStaticNoVal_3, segment_manager_t,  12>  Allocator_RHH_NoVal_3;
  typedef bip::node_allocator<RHHStaticNoVal_4, segment_manager_t,   6>  Allocator_RHH_NoVal_4;
  typedef bip::node_allocator<RHHStaticNoVal_5, segment_manager_t,   3>  Allocator_RHH_NoVal_5;
  typedef bip::node_allocator<RHHStaticNoVal_6, segment_manager_t,   1>  Allocator_RHH_NoVal_6;
  typedef bip::node_allocator<RHHStaticNoVal_7, segment_manager_t,   1>  Allocator_RHH_NoVal_7;
  typedef bip::node_allocator<RHHStaticNoVal_8, segment_manager_t,   1>  Allocator_RHH_NoVal_8;
  typedef bip::node_allocator<RHHStaticNoVal_9, segment_manager_t,   1>  Allocator_RHH_NoVal_9;
  typedef bip::node_allocator<RHHStaticNoVal_10, segment_manager_t,  1>  Allocator_RHH_NoVal_10;
  typedef bip::node_allocator<RHHStaticNoVal_11, segment_manager_t,  1>  Allocator_RHH_NoVal_11;
  typedef bip::node_allocator<RHHStaticNoVal_12, segment_manager_t,  1>  Allocator_RHH_NoVal_12;
  typedef bip::node_allocator<RHHStaticNoVal_13, segment_manager_t,  1>  Allocator_RHH_NoVal_13;
  typedef bip::node_allocator<RHHStaticNoVal_14, segment_manager_t,  1>  Allocator_RHH_NoVal_14;
  typedef bip::node_allocator<RHHStaticNoVal_15, segment_manager_t,  1>  Allocator_RHH_NoVal_15;
  typedef bip::node_allocator<RHHStaticNoVal_16, segment_manager_t,  1>  Allocator_RHH_NoVal_16;
  typedef bip::node_allocator<RHHStaticNoVal_17, segment_manager_t,  1>  Allocator_RHH_NoVal_17;
  typedef bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_normalarray_t;
  typedef bip::allocator<unsigned char, segment_manager_t> allocator_raw_t;
#endif
  ///  =========================================================================== ///
  ///                             Allocator Holder
  ///  =========================================================================== ///
  class AllocatorsHolder
  {
  public:
    explicit AllocatorsHolder(segment_manager_t* segment_manager)
    : allocator_rhh_noval_1(segment_manager)
    , allocator_rhh_noval_2(segment_manager)
    , allocator_rhh_noval_3(segment_manager)
    , allocator_rhh_noval_4(segment_manager)
    , allocator_rhh_noval_5(segment_manager)
    , allocator_rhh_noval_6(segment_manager)
    , allocator_rhh_noval_7(segment_manager)
    , allocator_rhh_noval_8(segment_manager)
    , allocator_rhh_noval_9(segment_manager)
    , allocator_rhh_noval_10(segment_manager)
    , allocator_rhh_noval_11(segment_manager)
    , allocator_rhh_noval_12(segment_manager)
    , allocator_rhh_noval_13(segment_manager)
    , allocator_rhh_noval_14(segment_manager)
    , allocator_rhh_noval_15(segment_manager)
    , allocator_rhh_noval_16(segment_manager)
    , allocator_rhh_noval_17(segment_manager)
    , allocator_normalarray(segment_manager)
    , allocator_raw(segment_manager)
    { }


  /// size = capacity * (1 + 8 + 1) + 8
  /// probedistance = 1 byte, key = 8 byte, value block = 1 byte
  Allocator_RHH_NoVal_1 allocator_rhh_noval_1;
  Allocator_RHH_NoVal_2 allocator_rhh_noval_2;
  Allocator_RHH_NoVal_3 allocator_rhh_noval_3;
  Allocator_RHH_NoVal_4 allocator_rhh_noval_4;
  Allocator_RHH_NoVal_5 allocator_rhh_noval_5;
  Allocator_RHH_NoVal_6 allocator_rhh_noval_6;
  Allocator_RHH_NoVal_7 allocator_rhh_noval_7;
  Allocator_RHH_NoVal_8 allocator_rhh_noval_8;
  Allocator_RHH_NoVal_9 allocator_rhh_noval_9;
  Allocator_RHH_NoVal_10 allocator_rhh_noval_10;
  Allocator_RHH_NoVal_11 allocator_rhh_noval_11;
  Allocator_RHH_NoVal_12 allocator_rhh_noval_12;
  Allocator_RHH_NoVal_13 allocator_rhh_noval_13;
  Allocator_RHH_NoVal_14 allocator_rhh_noval_14;
  Allocator_RHH_NoVal_15 allocator_rhh_noval_15;
  Allocator_RHH_NoVal_16 allocator_rhh_noval_16;
  Allocator_RHH_NoVal_17 allocator_rhh_noval_17;
  allocator_normalarray_t allocator_normalarray;
  allocator_raw_t allocator_raw;

  private:
    // static inline uint64_t cal_next_highest_power_of_2(uint64_t x)
    // {
    //   --x;
    //   x != x >> 1ULL;
    //   x != x >> 2ULL;
    //   x != x >> 4ULL;
    //   x != x >> 8ULL;
    //   x != x >> 16ULL;
    //   x != x >> 32ULL;
    //   return ++x;
    // }

  };
}

#endif