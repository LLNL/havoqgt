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

#ifndef HAVOQGT_MPI_RHHMGRSTATIC_HPP_INCLUDED
#define HAVOQGT_MPI_RHHMGRSTATIC_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include "RHHCommon.hpp"
#include "RHHAllocHolder.hpp"
#include "RHHStatic.hpp"

 namespace RHH {

  ///  =========================================================================== ///
  ///                            RHH Manager class
  ///  =========================================================================== ///
  template <typename KeyType, typename ValueType>
  class RHHMgrStatic {

  public:

    ///  ------------------------------------------------------ ///
    ///              Constructor / Destructor
    ///  ------------------------------------------------------ ///
    /// --- Constructor --- //
    explicit RHHMgrStatic(AllocatorsHolder& allocators)
    {
      m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_1.allocate(1).get());
      RHHStaticNoVal1* new_rhh = reinterpret_cast<RHHStaticNoVal1*>(m_ptr_);
      new_rhh->clear(false);
    }

    /// --- Copy constructor --- //

    /// --- Move constructor --- //
    RHHMgrStatic(RHHMgrStatic &&old_obj)
    {
      //DEBUG("RHHMgr move-constructor");
      m_ptr_ = old_obj.m_ptr_;
      old_obj.m_ptr_ = nullptr;
    }

    /// --- Destructor --- //
    ~RHHMgrStatic()
    {
      // DEBUG("RHHMgr destructor");
      if (m_ptr_ != nullptr) {
        DEBUG("Have to deallcate");
      }
    }

    /// ---  Move assignment operator --- //
    RHHMgrStatic &operator=(RHHMgrStatic&& old_obj)
    {
      //DEBUG("RHHMgr move-assignment");
      m_ptr_ = old_obj.m_ptr_;
      old_obj.m_ptr_ = nullptr;

      return *this;
    }


    ///  ------------------------------------------------------ ///
    ///              Public function
    ///  ------------------------------------------------------ ///
    /// insert a key into RHHStatic class
    bool insert_uniquely_static(AllocatorsHolder& allocators, KeyType& key, ValueType& val, uint64_t current_size)
    {
      uint64_t current_capacity;
      UpdateErrors err = insert_helper_static(key, val, current_size+1, &current_capacity);
      if (err == kDuplicated) return false;
      ++current_size;

      if (err == kLongProbedistance) {
        /// XXX: current implementation dose not allocate new RHH-array
        ///DEBUG("kLongProbedistance");
        grow_rhh_static(allocators, current_capacity, true);

      } else  if (current_size % current_capacity == 0) {
        /// full-capacity
        ///DEBUG("full-capacity");
        grow_rhh_static(allocators, current_capacity, false);
      }

      return true;
    }

    UpdateErrors insert_helper_static(KeyType& key, ValueType& val, const uint64_t require_capacity, uint64_t *current_capacity)
    {
      UpdateErrors err;

      if (require_capacity <= capacityRHHStatic1) {
        RHHStaticNoVal1* rhh = reinterpret_cast<RHHStaticNoVal1*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic1;

      } else if (require_capacity <= capacityRHHStatic2) {
        RHHStaticNoVal2* rhh = reinterpret_cast<RHHStaticNoVal2*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic2;

      }  else if (require_capacity <= capacityRHHStatic3) {
        RHHStaticNoVal3* rhh = reinterpret_cast<RHHStaticNoVal3*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic3;

      }  else if (require_capacity <= capacityRHHStatic4) {
        RHHStaticNoVal4* rhh = reinterpret_cast<RHHStaticNoVal4*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic4;

      }  else if (require_capacity <= capacityRHHStatic5) {
        RHHStaticNoVal5* rhh = reinterpret_cast<RHHStaticNoVal5*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic5;

      } else {
        RHHStaticNoVal6* rhh = reinterpret_cast<RHHStaticNoVal6*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = capacityRHHStatic6;
      }

      return err;
    }

    // inline bool erase(AllocatorsHolder &allocators, KeyType key, uint64_t* current_size)
    // {
    //   uint64_t capacity = AllocatorsHolder.cal_next_highest_power_of_2(*current_size);
    //   RHHAdjlistType<capacity>* rhh = reinterpret_cast<RHHAdjlistType<capacity>*>(m_ptr_);
    //   return rhh->erase(allocators, key);
    // }


    ///  ------------------------------------------------------ ///
    ///              Private Member Functions
    ///  ------------------------------------------------------ ///

    /// XXX: should use template function
#define ALLOCATE_AND_MOVE(OLD_C, NEW_C) do{\
    RHHStaticNoVal##OLD_C* old_rhh = reinterpret_cast<RHHStaticNoVal##OLD_C*>(m_ptr_);\
    m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##NEW_C.allocate(1).get());\
    RHHStaticNoVal##NEW_C* new_rhh = reinterpret_cast<RHHStaticNoVal##NEW_C*>(m_ptr_);\
    new_rhh->clear(false);\
    const uint64_t old_capacity = old_rhh->capacity();\
    for (uint64_t i = 0; i < old_capacity; ++i) {\
      if (old_rhh->is_valid(i)) {\
        new_rhh->insert_uniquely(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);\
      }\
    }\
    allocators.deallocate_rhh_static(reinterpret_cast<void*>(old_rhh), capacityRHHStatic##OLD_C);\
  }while(0)

  uint64_t grow_rhh_static(AllocatorsHolder &allocators, uint64_t current_capacity, bool has_long_probedistance)
  {
    if (current_capacity <= capacityRHHStatic1) {
      ALLOCATE_AND_MOVE(1, 2);
    } else if (current_capacity <= capacityRHHStatic2){
      ALLOCATE_AND_MOVE(2, 3);
    } else if (current_capacity <= capacityRHHStatic3){
      ALLOCATE_AND_MOVE(3, 4);
    } else if (current_capacity <= capacityRHHStatic4){
      ALLOCATE_AND_MOVE(4, 5);
    } else if (current_capacity <= capacityRHHStatic5){
      ALLOCATE_AND_MOVE(5, 6);
    } else {
      /// allocate chained rhh
      RHHStaticNoVal6* old_rhh = reinterpret_cast<RHHStaticNoVal6*>(m_ptr_);
      m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_6.allocate(1).get());
      RHHStaticNoVal6* new_rhh = reinterpret_cast<RHHStaticNoVal6*>(m_ptr_);
      new_rhh->clear(false);
      new_rhh->m_next_ = old_rhh;

      if (has_long_probedistance) {
        /// move longprobedisrance element to new array
        const uint64_t old_capacity = old_rhh->capacity();
        for (uint64_t i = 0; i < old_capacity; ++i) {
          if (old_rhh->is_valid(i) && old_rhh->is_longprobedistance(i)) {
            new_rhh->insert_uniquely(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);
            old_rhh->delete_key(i);
          }
        }
      }

    }

  }

  void inline free_buffer_rhh_static(AllocatorsHolder &allocators, void* ptr, uint64_t capacity)
  {
    allocators.deallocate_rhh_static(ptr, capacity);
  }

  void* m_ptr_;
};
}

#endif