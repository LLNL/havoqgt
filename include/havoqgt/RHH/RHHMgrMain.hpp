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

#ifndef HAVOQGT_MPI_RHHMgr_HPP_INCLUDED
#define HAVOQGT_MPI_RHHMgr_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include "RHHStatic.hpp"
#include "RHHMain.hpp"

namespace RHH {


///  =========================================================================== ///
///                            RHH Manager class
///  =========================================================================== ///
template <typename KeyType, typename ValueType>
class RHHMgrMain {

 public:

  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Constructor --- //
  explicit RHHMgrMain(AllocatorsHolder &allocators, const uint64_t initial_capasity)
  {
    alloc(allocators, initial_capasity);
  }

  /// --- Copy constructor --- //

  /// --- Move constructor --- //
  RHHMgr(RHHMgr &&old_obj)
  {
    //DEBUG("RHHMgr move-constructor");
    m_ptr_ = old_obj.m_ptr_;
    old_obj.m_ptr_ = nullptr;
  }

  /// --- Destructor --- //
  ~RHHMgr()
  {
    // DEBUG("RHHMgr destructor");
    if (m_ptr_ != nullptr) {
      assert(false);
      // DEBUG("RHHMgr destructor : need free_buffer");
      //free_buffer(m_pos_head_, 0); /// XXX: memsize = 0
    }
  }

  /// ---  Move assignment operator --- //
  RHHMgr &operator=(RHHMgr&& old_obj)
  {
    //DEBUG("RHHMgr move-assignment");
    m_ptr_ = old_obj.m_ptr_;
    old_obj.m_ptr_ = nullptr;

    return *this;
  }

  /// insert a edge into RHHMain class
  bool insert_uniquely(AllocatorsHolder& allocators, KeyType &key, ValueType &val)
  {
    UpdateErrors err = insert_helper(key, val);
    if (err == kDuplicated) return false;

    if (err == kReachingFUllCapacity) {
      capacity = grow_rhh_main(allocators);
    } else if (err == kLongProbedistance) {
      /// XXX: current implementation dose not allocate new RHH-array
      capacity = grow_rhh_main(allocators);
      DEBUG("kLongProbedistance");
    }
    return true;
  }

  inline UpdateErrors insert_helper(AllocatorsHolder& allocators, KeyType& key, ValueType& val)
  {
    RHHMain<KeyType, ValueType>* rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    return rhh->insert_uniquely(key, val);
  }

  // inline bool erase(AllocatorsHolder &allocators, KeyType key, ValueType val, uint64_t* current_size)
  // {
  //   uint64_t capacity = AllocatorsHolder.cal_next_highest_power_of_2(*current_size);
  //   RHHAdjlistType<capacity>* rhh = reinterpret_cast<RHHAdjlistType<capacity>*>(m_ptr_);
  //   return rhh->erase(allocators, key, val);
  // }


  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///
  static const uint64_t kCapacityGrowingFactor = 2ULL;

  uint64_t grow_rhh_main(AllocatorsHolder &allocators)
  {
    RHHMain<KeyType, ValueType>* old_rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    const uint64_t old_capacity = old_rhh->capacity;

    const uint64_t capacity = old_capacity * kCapacityGrowingFactor;
    m_ptr_ = allocators.allocate_rhh_main(capacity);
    RHHMain<KeyType, ValueType>* rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    rhh->clear();

    /// now copy over old elems
    for (uint64_t i = 0; i < old_capacity; i++) {
      if (rhh->is_deleted(rhh->m_property_block_[i])) {
        rhh->insert_uniquely(std::move(rhh->m_key_block_[i]), std::move(rhh->m_value_block_[i]));
      }
    }
    free_buffer_rhh_main(allocators, reinterpret_cast<void*>(old_rhh), old_rhh->capacity;);

    return capacity;
  }

  void inline free_buffer_rhh_main(AllocatorsHolder &allocators, void* ptr, uint64_t capacity)
  {
    allocators.deallocate_raw(ptr, capacity);
  }

  void* m_ptr_;
};
}
#endif