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
#include <fstream>
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
      // DEBUG("RHHMgrStatic");
      m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_1.allocate(1).get());
      RHHStaticNoVal1* new_rhh = reinterpret_cast<RHHStaticNoVal1*>(m_ptr_);
      new_rhh->m_next_ = nullptr;
      new_rhh->clear(false);
    }

    /// --- Copy constructor --- //

    /// --- Move constructor --- //
    RHHMgrStatic(RHHMgrStatic &&old_obj)
    {
      DEBUG("RHHMgr move-constructor");
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
      DEBUG("RHHMgr move-assignment");
      m_ptr_ = old_obj.m_ptr_;
      old_obj.m_ptr_ = nullptr;

      return *this;
    }


    ///  ------------------------------------------------------ ///
    ///              Public function
    ///  ------------------------------------------------------ ///
    /// insert a key into RHHStatic class
    bool insert_uniquely(AllocatorsHolder& allocators, KeyType& key, ValueType& val, uint64_t current_size)
    {

      UpdateErrors err = insert_uniquely_helper(allocators, key, val, current_size);
      if (err == kDuplicated) return false;
      ++current_size;

      if (err == kLongProbedistance) {
        /// XXX: current implementation dose not allocate new RHH-array
        grow_rhh_static(allocators, current_size, true);
      }

      return true;
    }


#define INSERT_AND_CHECK_CAPACITY(OLD_C, NEW_C) \
    do { \
        RHHStaticNoVal##OLD_C* rhh = reinterpret_cast<RHHStaticNoVal##OLD_C*>(m_ptr_);\
        if (!rhh->try_unique_key_insertion(key, val))\
          return kDuplicated;\
        if (current_size == capacityRHHStatic##OLD_C) {\
          grow_rhh_static(allocators, current_size, false);\
          RHHStaticNoVal##NEW_C* rhh_new = reinterpret_cast<RHHStaticNoVal##NEW_C*>(m_ptr_);\
          return rhh_new->insert_uniquely(key, val);\
        }\
        return rhh->insert_uniquely(key, val);\
    } while (0)

    UpdateErrors insert_uniquely_helper(AllocatorsHolder &allocators, KeyType& key, ValueType& val, const uint64_t current_size)
    {
      if (current_size <= capacityRHHStatic1) {
        INSERT_AND_CHECK_CAPACITY(1, 2);

      } else if (current_size <= capacityRHHStatic2) {
        INSERT_AND_CHECK_CAPACITY(2, 3);

      }  else if (current_size <= capacityRHHStatic3) {
        INSERT_AND_CHECK_CAPACITY(3, 4);

      }  else if (current_size <= capacityRHHStatic4) {
        INSERT_AND_CHECK_CAPACITY(4, 5);

      }  else if (current_size <= capacityRHHStatic5) {
        INSERT_AND_CHECK_CAPACITY(5, 6);

      } else {
        RHHStaticNoVal6* rhh = reinterpret_cast<RHHStaticNoVal6*>(m_ptr_);
        if (!rhh->try_unique_key_insertion(key, val))
          return kDuplicated;
        if (current_size % capacityRHHStatic6 == 0) {
          grow_rhh_static(allocators, current_size, false);
          RHHStaticNoVal6* rhh_new = reinterpret_cast<RHHStaticNoVal6*>(m_ptr_);
          return rhh_new->insert_uniquely(key, val);
        }
        return rhh->insert_uniquely(key, val);

      }
      return kSucceed;
    }

    // inline bool erase(AllocatorsHolder &allocators, KeyType key, uint64_t current_size)
    // {
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
    new_rhh->m_next_ = nullptr;\
    new_rhh->clear(false);\
    const uint64_t old_capacity = old_rhh->capacity();\
    for (uint64_t i = 0; i < old_capacity; ++i) {\
      if (old_rhh->is_valid(i)) {\
        new_rhh->insert_uniquely(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);\
      }\
    }\
    free_buffer_rhh_static(allocators, reinterpret_cast<void*>(old_rhh), capacityRHHStatic##OLD_C);\
  }while(0)

  void grow_rhh_static(AllocatorsHolder &allocators, uint64_t current_capacity, bool has_long_probedistance)
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
      /// FIXME: deallocation error!
    if (capacity == capacityRHHStatic1) {
      allocators.allocator_rhh_noval_1.deallocate(bip::offset_ptr<RHHStaticNoVal1>(reinterpret_cast<RHHStaticNoVal1*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic2) {
      allocators.allocator_rhh_noval_2.deallocate(bip::offset_ptr<RHHStaticNoVal2>(reinterpret_cast<RHHStaticNoVal2*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic3) {
      allocators.allocator_rhh_noval_3.deallocate(bip::offset_ptr<RHHStaticNoVal3>(reinterpret_cast<RHHStaticNoVal3*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic4) {
      allocators.allocator_rhh_noval_4.deallocate(bip::offset_ptr<RHHStaticNoVal4>(reinterpret_cast<RHHStaticNoVal4*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic5) {
      allocators.allocator_rhh_noval_5.deallocate(bip::offset_ptr<RHHStaticNoVal5>(reinterpret_cast<RHHStaticNoVal5*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic6) {
      allocators.allocator_rhh_noval_6.deallocate(bip::offset_ptr<RHHStaticNoVal6>(reinterpret_cast<RHHStaticNoVal6*>(ptr)), 1);
    } else {
      assert(false);
    }
  }


  /// XXX: should use template function
#define DISP_KEYS(C, PRFX, OF) \
  do{ \
    RHHStaticNoVal##C* rhh = reinterpret_cast<RHHStaticNoVal##C*>(m_ptr_); \
    rhh->disp_keys(PRFX, OF); \
  } while (0)

  void disp_keys(uint64_t current_capacity, std::string prefix, std::ofstream& output_file)
  {
    if (current_capacity <= capacityRHHStatic1) {
      DISP_KEYS(1, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic2){
      DISP_KEYS(2, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic3){
      DISP_KEYS(3, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic4){
      DISP_KEYS(4, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic5){
      DISP_KEYS(5, prefix, output_file);
    } else {
      DISP_KEYS(6, prefix, output_file);
    }
  }


  void* m_ptr_;
};
}

#endif