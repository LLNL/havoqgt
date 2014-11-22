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
      RHHStaticNoVal_1* new_rhh = reinterpret_cast<RHHStaticNoVal_1*>(m_ptr_);
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
        RHHStaticNoVal_##OLD_C* rhh = reinterpret_cast<RHHStaticNoVal_##OLD_C*>(m_ptr_);\
        if (!rhh->try_unique_key_insertion(key, val))\
          return kDuplicated;\
        if (current_size == capacityRHHStatic_##OLD_C) {\
          grow_rhh_static(allocators, current_size, false);\
          RHHStaticNoVal_##NEW_C* rhh_new = reinterpret_cast<RHHStaticNoVal_##NEW_C*>(m_ptr_);\
          return rhh_new->insert_uniquely(key, val);\
        }\
        return rhh->insert_uniquely(key, val);\
    } while (0)

    UpdateErrors insert_uniquely_helper(AllocatorsHolder &allocators, KeyType& key, ValueType& val, const uint64_t current_size)
    {
      if (current_size <= capacityRHHStatic_1) {
        INSERT_AND_CHECK_CAPACITY(1, 2);

      } else if (current_size <= capacityRHHStatic_2) {
        INSERT_AND_CHECK_CAPACITY(2, 3);

      } else if (current_size <= capacityRHHStatic_3) {
        INSERT_AND_CHECK_CAPACITY(3, 4);

      } else if (current_size <= capacityRHHStatic_4) {
        INSERT_AND_CHECK_CAPACITY(4, 5);

      } else if (current_size <= capacityRHHStatic_5) {
        INSERT_AND_CHECK_CAPACITY(5, 6);

      } else if (current_size <= capacityRHHStatic_6) {
        INSERT_AND_CHECK_CAPACITY(6, 7);

      } else if (current_size <= capacityRHHStatic_7) {
        INSERT_AND_CHECK_CAPACITY(7, 8);

      } else if (current_size <= capacityRHHStatic_8) {
        INSERT_AND_CHECK_CAPACITY(8, 9);

      } else if (current_size <= capacityRHHStatic_9) {
        INSERT_AND_CHECK_CAPACITY(9, 10);

      } else if (current_size <= capacityRHHStatic_10) {
        INSERT_AND_CHECK_CAPACITY(10, 11);

      } else if (current_size <= capacityRHHStatic_11) {
        INSERT_AND_CHECK_CAPACITY(11, 12);

      } else if (current_size <= capacityRHHStatic_12) {
        INSERT_AND_CHECK_CAPACITY(12, 13);

      } else if (current_size <= capacityRHHStatic_13) {
        INSERT_AND_CHECK_CAPACITY(13, 14);

      } else if (current_size <= capacityRHHStatic_14) {
        INSERT_AND_CHECK_CAPACITY(14, 15);

      } else if (current_size <= capacityRHHStatic_15) {
        INSERT_AND_CHECK_CAPACITY(15, 16);

      } else if (current_size <= capacityRHHStatic_16) {
        INSERT_AND_CHECK_CAPACITY(16, 17);

      } else {
        RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
        if (!rhh->try_unique_key_insertion(key, val))
          return kDuplicated;
        if (current_size % (uint64_t)(capacityRHHStatic_17 * 0.9f) == 0) {
          grow_rhh_static(allocators, current_size, false);
          RHHStaticNoVal_17* rhh_new = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
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
    RHHStaticNoVal_##OLD_C* old_rhh = reinterpret_cast<RHHStaticNoVal_##OLD_C*>(m_ptr_);\
    m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##NEW_C.allocate(1).get());\
    RHHStaticNoVal_##NEW_C* new_rhh = reinterpret_cast<RHHStaticNoVal_##NEW_C*>(m_ptr_);\
    new_rhh->m_next_ = nullptr;\
    new_rhh->clear(false);\
    const uint64_t old_capacity = old_rhh->capacity();\
    for (uint64_t i = 0; i < old_capacity; ++i) {\
      if (old_rhh->is_valid(i)) {\
        new_rhh->insert_uniquely(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);\
      }\
    }\
    free_buffer_rhh_static(allocators, reinterpret_cast<void*>(old_rhh), capacityRHHStatic_##OLD_C);\
  }while(0)

  void grow_rhh_static(AllocatorsHolder &allocators, uint64_t current_capacity, bool has_long_probedistance)
  {
    if (current_capacity <= capacityRHHStatic_1) {
      ALLOCATE_AND_MOVE(1, 2);
    } else if (current_capacity <= capacityRHHStatic_2){
      ALLOCATE_AND_MOVE(2, 3);
    } else if (current_capacity <= capacityRHHStatic_3){
      ALLOCATE_AND_MOVE(3, 4);
    } else if (current_capacity <= capacityRHHStatic_4){
      ALLOCATE_AND_MOVE(4, 5);
    } else if (current_capacity <= capacityRHHStatic_5){
      ALLOCATE_AND_MOVE(5, 6);
    } else if (current_capacity <= capacityRHHStatic_6) {
        ALLOCATE_AND_MOVE(6, 7);
    } else if (current_capacity <= capacityRHHStatic_7) {
      ALLOCATE_AND_MOVE(7, 8);
    } else if (current_capacity <= capacityRHHStatic_8) {
      ALLOCATE_AND_MOVE(8, 9);
    } else if (current_capacity <= capacityRHHStatic_9) {
      ALLOCATE_AND_MOVE(9, 10);
    } else if (current_capacity <= capacityRHHStatic_10) {
      ALLOCATE_AND_MOVE(10, 11);
    } else if (current_capacity <= capacityRHHStatic_11) {
      ALLOCATE_AND_MOVE(11, 12);
    } else if (current_capacity <= capacityRHHStatic_12) {
      ALLOCATE_AND_MOVE(12, 13);
    } else if (current_capacity <= capacityRHHStatic_13) {
      ALLOCATE_AND_MOVE(13, 14);
    } else if (current_capacity <= capacityRHHStatic_14) {
      ALLOCATE_AND_MOVE(14, 15);
    } else if (current_capacity <= capacityRHHStatic_15) {
      ALLOCATE_AND_MOVE(15, 16);
    } else if (current_capacity <= capacityRHHStatic_16) {
      ALLOCATE_AND_MOVE(16, 17);
    } else {

      /// allocate chained rhh
      RHHStaticNoVal_17* old_rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_17.allocate(1).get());
      RHHStaticNoVal_17* new_rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
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
    if (capacity == capacityRHHStatic_1) {
      allocators.allocator_rhh_noval_1.deallocate(bip::offset_ptr<RHHStaticNoVal_1>(reinterpret_cast<RHHStaticNoVal_1*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_2) {
      allocators.allocator_rhh_noval_2.deallocate(bip::offset_ptr<RHHStaticNoVal_2>(reinterpret_cast<RHHStaticNoVal_2*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_3) {
      allocators.allocator_rhh_noval_3.deallocate(bip::offset_ptr<RHHStaticNoVal_3>(reinterpret_cast<RHHStaticNoVal_3*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_4) {
      allocators.allocator_rhh_noval_4.deallocate(bip::offset_ptr<RHHStaticNoVal_4>(reinterpret_cast<RHHStaticNoVal_4*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_5) {
      allocators.allocator_rhh_noval_5.deallocate(bip::offset_ptr<RHHStaticNoVal_5>(reinterpret_cast<RHHStaticNoVal_5*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_6) {
      allocators.allocator_rhh_noval_6.deallocate(bip::offset_ptr<RHHStaticNoVal_6>(reinterpret_cast<RHHStaticNoVal_6*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_7) {
      allocators.allocator_rhh_noval_7.deallocate(bip::offset_ptr<RHHStaticNoVal_7>(reinterpret_cast<RHHStaticNoVal_7*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_8) {
      allocators.allocator_rhh_noval_8.deallocate(bip::offset_ptr<RHHStaticNoVal_8>(reinterpret_cast<RHHStaticNoVal_8*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_9) {
      allocators.allocator_rhh_noval_9.deallocate(bip::offset_ptr<RHHStaticNoVal_9>(reinterpret_cast<RHHStaticNoVal_9*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_10) {
      allocators.allocator_rhh_noval_10.deallocate(bip::offset_ptr<RHHStaticNoVal_10>(reinterpret_cast<RHHStaticNoVal_10*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_11) {
      allocators.allocator_rhh_noval_11.deallocate(bip::offset_ptr<RHHStaticNoVal_11>(reinterpret_cast<RHHStaticNoVal_11*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_12) {
      allocators.allocator_rhh_noval_12.deallocate(bip::offset_ptr<RHHStaticNoVal_12>(reinterpret_cast<RHHStaticNoVal_12*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_13) {
      allocators.allocator_rhh_noval_13.deallocate(bip::offset_ptr<RHHStaticNoVal_13>(reinterpret_cast<RHHStaticNoVal_13*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_14) {
      allocators.allocator_rhh_noval_14.deallocate(bip::offset_ptr<RHHStaticNoVal_14>(reinterpret_cast<RHHStaticNoVal_14*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_15) {
      allocators.allocator_rhh_noval_15.deallocate(bip::offset_ptr<RHHStaticNoVal_15>(reinterpret_cast<RHHStaticNoVal_15*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_16) {
      allocators.allocator_rhh_noval_16.deallocate(bip::offset_ptr<RHHStaticNoVal_16>(reinterpret_cast<RHHStaticNoVal_16*>(ptr)), 1);
    } else if (capacity == capacityRHHStatic_17) {
      allocators.allocator_rhh_noval_17.deallocate(bip::offset_ptr<RHHStaticNoVal_17>(reinterpret_cast<RHHStaticNoVal_17*>(ptr)), 1);
    } else {
      assert(false);
    }
  }


/// XXX: should use template function
#define DISP_KEYS(C, PRFX, OF) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->disp_keys(PRFX, OF); \
  } while (0)

  void disp_keys(uint64_t current_capacity, std::string prefix, std::ofstream& output_file)
  {
    if (current_capacity <= capacityRHHStatic_1) {
      DISP_KEYS(1, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_2){
      DISP_KEYS(2, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_3){
      DISP_KEYS(3, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_4){
      DISP_KEYS(4, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_5){
      DISP_KEYS(5, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_6){
      DISP_KEYS(6, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_7){
      DISP_KEYS(7, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_8){
      DISP_KEYS(8, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_9){
      DISP_KEYS(9, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_10){
      DISP_KEYS(10, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_11){
      DISP_KEYS(11, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_12){
      DISP_KEYS(12, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_13){
      DISP_KEYS(13, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_14){
      DISP_KEYS(14, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_15){
      DISP_KEYS(15, prefix, output_file);
    } else if (current_capacity <= capacityRHHStatic_16){
      DISP_KEYS(16, prefix, output_file);
    } else {
      DISP_KEYS(17, prefix, output_file);
    }
  }

  uint64_t cal_depth(uint64_t current_capacity)
  {
    if (current_capacity <= capacityRHHStatic_16){
      return 0;
    } else {
      RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      return rhh->cal_depth();
    }
  }

/// XXX: should use template function
#define DISP_PROBEDISTANCE(C, OF) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->disp_probedistance(OF); \
  } while (0)

  void disp_probedistance(uint64_t current_capacity, std::ofstream& output_file)
  {
    if (current_capacity <= capacityRHHStatic_1) {
      DISP_PROBEDISTANCE(1, output_file);
    } else if (current_capacity <= capacityRHHStatic_2){
      DISP_PROBEDISTANCE(2, output_file);
    } else if (current_capacity <= capacityRHHStatic_3){
      DISP_PROBEDISTANCE(3, output_file);
    } else if (current_capacity <= capacityRHHStatic_4){
      DISP_PROBEDISTANCE(4, output_file);
    } else if (current_capacity <= capacityRHHStatic_5){
      DISP_PROBEDISTANCE(5, output_file);
    } else if (current_capacity <= capacityRHHStatic_6){
      DISP_PROBEDISTANCE(6, output_file);
    } else if (current_capacity <= capacityRHHStatic_7){
      DISP_PROBEDISTANCE(7, output_file);
    } else if (current_capacity <= capacityRHHStatic_8){
      DISP_PROBEDISTANCE(8, output_file);
    } else if (current_capacity <= capacityRHHStatic_9){
      DISP_PROBEDISTANCE(9, output_file);
    } else if (current_capacity <= capacityRHHStatic_10){
      DISP_PROBEDISTANCE(10, output_file);
    } else if (current_capacity <= capacityRHHStatic_11){
      DISP_PROBEDISTANCE(11, output_file);
    } else if (current_capacity <= capacityRHHStatic_12){
      DISP_PROBEDISTANCE(12, output_file);
    } else if (current_capacity <= capacityRHHStatic_13){
      DISP_PROBEDISTANCE(13, output_file);
    } else if (current_capacity <= capacityRHHStatic_14){
      DISP_PROBEDISTANCE(14, output_file);
    } else if (current_capacity <= capacityRHHStatic_15){
      DISP_PROBEDISTANCE(15, output_file);
    } else if (current_capacity <= capacityRHHStatic_16){
      DISP_PROBEDISTANCE(16, output_file);
    } else {
      DISP_PROBEDISTANCE(17, output_file);
    }
  }

  void* m_ptr_;
};
}

#endif
