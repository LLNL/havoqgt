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
#include "RHHUtility.hpp"
#include "RHHStatic.hpp"

namespace RHH {
  
  ///  =========================================================================== ///
  ///                            RHH Manager class
  ///  =========================================================================== ///
  template <typename KeyType, typename ValueType>
  class RHHMgrStatic {
    
  public:
    
    template<C>
    typedef RHHStatic<KeyType, ValueType, C, PropertyBlockType> RHHAdjlistType;
    
    
    ///  ------------------------------------------------------ ///
    ///              Constructor / Destructor
    ///  ------------------------------------------------------ ///
    /// --- Constructor --- //
    explicit RHHMgrStatic(AllocatorsHolder &allocators)
    {
      m_ptr_ = reinterpret_cast<void*>(allocator_rhh_2.allocate(1).get());
      AllocatorsHolder::RHHStaticNoVal2* new_rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal2*>(m_ptr_);
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
        assert(false);
        // DEBUG("RHHMgr destructor : need free_buffer");
        //free_buffer(m_pos_head_, 0); /// XXX: memsize = 0
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
      UpdateErrors err = insert_helper_static(key, val, current_size, &current_capacity);
      if (err == kDuplicated) return false;
      
      if (++current_size > current_capacity / 100ULL * 90ULL) {
        /// near-capacity
        grow_rhh_static(allocators, current_size);
      } else if (err == kLongProbedistance) {
        /// XXX: current implementation dose not allocate new RHH-array
        grow_rhh_static(allocators, current_size);
        DEBUG("kLongProbedistance");
      }
      return true;
    }
    
    UpdateErrors insert_helper_static(KeyType& key, ValueType& val, const uint64_t current_size, uint64_t *current_capacity)
    {
      UpdateErrors err;
      
      if (current_size <= AllocatorsHolder::capacityRHHStatic1) {
        AllocatorsHolder::RHHStaticNoVal* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal1*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic1;
      } else if (current_size <= AllocatorsHolder::capacity2) {
        AllocatorsHolder::RHHStaticNoVa2* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal2*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic2;
      }  else if (current_size <= AllocatorsHolder::capacity3) {
        AllocatorsHolder::RHHStaticNoVal3* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal3*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic3;
      }  else if (current_size <= AllocatorsHolder::capacity4) {
        AllocatorsHolder::RHHStaticNoVal4* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal4*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic4;
      }  else if (current_size <= AllocatorsHolder::capacity5) {
        AllocatorsHolder::RHHStaticNoVal5* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal5*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic5;
      } else {
        AllocatorsHolder::RHHStaticNoVal6* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal6*>(m_ptr_);
        err = rhh->insert_uniquely(key, val);
        *current_capacity = AllocatorsHolder::capacityRHHStatic6;
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
#define ALLOCATE_AND_MOVE(OLD_C, NEW_C) do {
    AllocatorsHolder::RHHStaticNoVal ## OLD_C* old_rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal ## OLD_C*>(m_ptr_);
    m_ptr_ = reinterpret_cast<void*>(allocator_rhh_ ## NEW_C.allocate(1).get());
    AllocatorsHolder::RHHStaticNoVal ## NEW_C* new_rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal ## NEW_C*>(m_ptr_);
    
    new_rhh->clear(false);
    const uint64_t old_capacity = old_rhh->capacity;
    for (uint64_t i = 0; i < old_capacity; ++i) {
      if (!old_rhh->is_deleted(old_rhh->m_property_block_[i])) {
        new_rhh->insert_uniquely(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);
      }
    }
  } while(0)
    
    uint64_t grow_rhh_static(AllocatorsHolder &allocators, uint64_t current_size)
  {
    if (current_size <= AllocatorsHolder::capacityRHHStatic1) {
      ALLOCATE_AND_MOVE(1, 2);
    } else if (current_size <= AllocatorsHolder::capacityRHHStatic2){
      ALLOCATE_AND_MOVE(2, 3);
    } else if (current_size <= AllocatorsHolder::capacityRHHStatic3){
      ALLOCATE_AND_MOVE(3, 4);
    } else if (current_size <= AllocatorsHolder::capacityRHHStatic4){
      ALLOCATE_AND_MOVE(4, 5);
    } else if (current_size <= AllocatorsHolder::capacityRHHStatic5){
      ALLOCATE_AND_MOVE(5, 6);
    } else {
      /// allocate chained rhh
      AllocatorsHolder::RHHStaticNoVal6* old_rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal6*>(m_ptr_);
      m_ptr_ = reinterpret_cast<void*>(allocator_rhh_6.allocate(1).get());
      AllocatorsHolder::RHHStaticNoVal6* new_rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal6*>(m_ptr_);
      new_rhh->m_next_ = old_rhh;
      new_rhh->clear(false);
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