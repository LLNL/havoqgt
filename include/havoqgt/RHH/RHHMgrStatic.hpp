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
      m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_2.allocate(1).get());
      RHHStaticNoVal_2* new_rhh = reinterpret_cast<RHHStaticNoVal_2*>(m_ptr_);
      new_rhh->m_next_ = nullptr;
      new_rhh->clear(false);
    }

    /// --- Copy constructor --- //
    RHHMgrStatic(const RHHMgrStatic &old_obj)
    {
      assert(false);
    }

    /// --- Move constructor --- //
    RHHMgrStatic(RHHMgrStatic &&old_obj)
    {
      assert(false);
      m_ptr_ = old_obj.m_ptr_;
      old_obj.m_ptr_ = nullptr;
    }

    /// --- Destructor --- //
    ~RHHMgrStatic()
    {
      assert(m_ptr_ != nullptr);
    }

    /// ---  Move assignment operator --- //
    RHHMgrStatic &operator=(RHHMgrStatic&& old_obj)
    {
      assert(false);
      m_ptr_ = old_obj.m_ptr_;
      old_obj.m_ptr_ = nullptr;

      return *this;
    }


    ///  ------------------------------------------------------ ///
    ///              Public function
    ///  ------------------------------------------------------ ///
    /// insert a key into RHHStatic class
    inline bool insert_uniquely(AllocatorsHolder& allocators, KeyType& key, ValueType& val, uint64_t current_size)
    {
      UpdateErrors err = insert_uniquely_helper(allocators, key, val, current_size);
      return (err != kDuplicated);
    }

    inline void free(AllocatorsHolder &allocators, uint64_t size)
    {
      free_buffer_rhh_static(allocators, size);
    }

    inline bool delete_key(AllocatorsHolder& allocators, KeyType& key, uint64_t current_size)
    {
      return delete_key_helper(allocators, key, current_size);
    }


    template<typename CurRHH, typename CurRHHAllocator, uint64_t CurCapacity,
             typename SrkRHH, typename SrkRHHAllocator, uint64_t SrkCapacity>
    bool delete_key_and_shrink(AllocatorsHolder &allocators, CurRHHAllocator cur_rhh_allocator, SrkRHHAllocator srk_rhh_allocator,
                              const KeyType& key, uint64_t current_size)
    {
      /// --- Delete target element --- ///
      CurRHH* cur_rhh = reinterpret_cast<CurRHH*>(m_ptr_);
      if (!cur_rhh->delete_key(key)) {
        return false;
      }
      --current_size;

      /// --- If current size is less than previous size capacity, shrink RHH --- ///
      if (current_size <= SrkCapacity * kFullCalacityFactor) {
        SrkRHH* srk_rhh = allocate_rhh_static<SrkRHHAllocator, SrkRHH>(srk_rhh_allocator);
        insert_rhh_elements_and_free<CurRHHAllocator, CurRHH, CurCapacity>(allocators, cur_rhh_allocator, cur_rhh, current_size);
      }
      return true;
    }

    template<typename Allocator, typename RHH>
    inline RHH *allocate_rhh_static(Allocator& allocator)
    {
      m_ptr_ = reinterpret_cast<void*>(allocator.allocate(1).get());
      RHH* rhh = reinterpret_cast<RHH*>(m_ptr_);
      rhh->m_next_ = nullptr;
      rhh->clear();
      return rhh;
    }

    template<typename Allocator, typename RHH>
    inline void deallocate_rhh_static(Allocator& allocator, RHH *rhh)
    {
      allocator.deallocate(bip::offset_ptr<RHH>(rhh), 1);
    }

    template<typename RHHAllocator, typename RHH, uint64_t Capacity>
    inline void insert_rhh_elements_and_free(AllocatorsHolder& allocators, RHHAllocator& rhh_allocator, RHH& rhh, const uint64_t current_size)
    {
      while (rhh) {
        RHH* next_rhh = rhh->m_next_;
        for (int i = 0; i < Capacity; ++i) {
          if (rhh->is_valid(i))
            insert_helper(allocators, rhh->m_key_block_[i], rhh->m_value_block_[i], current_size);
        }
        deallocate_rhh_static<RHHAllocator, RHH>(rhh_allocator, rhh);
        rhh = next_rhh;
      }
    }

#define DELETE_KEY_AND_SHRINK(SRK_C, CUR_C) \
    do { \
        RHHStaticNoVal_##CUR_C* cur_rhh = reinterpret_cast<RHHStaticNoVal_##CUR_C*>(m_ptr_);\
        if (!cur_rhh->delete_key(key)) {\
          return false;\
        }\
        --current_size;\
        if (current_size <= capacityRHHStatic_##SRK_C * kFullCalacityFactor) {\
          m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##SRK_C.allocate(1).get());\
          RHHStaticNoVal_##SRK_C* prv_rhh = reinterpret_cast<RHHStaticNoVal_##SRK_C*>(m_ptr_);\
          prv_rhh->m_next_ = nullptr;\
          prv_rhh->clear(false);\
          while (1) {\
            RHHStaticNoVal_##CUR_C* cur_next_rhh = cur_rhh->m_next_;\
            for (int i = 0; i < capacityRHHStatic_##CUR_C; ++i) {\
              if (cur_rhh->is_valid(i)) {\
                insert_helper(allocators, cur_rhh->m_key_block_[i], cur_rhh->m_value_block_[i], current_size);\
              }\
            }\
            allocators.allocator_rhh_noval_##CUR_C.deallocate(bip::offset_ptr<RHHStaticNoVal_##CUR_C>(cur_rhh), 1);\
            if (cur_next_rhh == nullptr) {\
              break;\
            }\
            cur_rhh = cur_next_rhh;\
          }\
        }\
        return true;\
    } while (0)

    /// TODO: use binary tree ?
    bool delete_key_helper(AllocatorsHolder &allocators, const KeyType& key, uint64_t current_size)
    {
      if (current_size <= capacityRHHStatic_1 * kFullCalacityFactor) {
        assert(false);
        RHHStaticNoVal_1* rhh = reinterpret_cast<RHHStaticNoVal_1*>(m_ptr_);
        return rhh->delete_key(key);

      } else if (current_size <= capacityRHHStatic_2 * kFullCalacityFactor) {
        /// NOTE: If current_size is equal or less than capacityRHHStatic_1,
        /// we use normal array.
        RHHStaticNoVal_2* rhh = reinterpret_cast<RHHStaticNoVal_2*>(m_ptr_);
        return rhh->delete_key(key);

      } else if (current_size <= capacityRHHStatic_3 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(2, 3);

      } else if (current_size <= capacityRHHStatic_4 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(3, 4);

      } else if (current_size <= capacityRHHStatic_5 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(4, 5);

      } else if (current_size <= capacityRHHStatic_6 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(5, 6);

      } else if (current_size <= capacityRHHStatic_7 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(6, 7);

      } else if (current_size <= capacityRHHStatic_8 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(7, 8);

      } else if (current_size <= capacityRHHStatic_9 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(8, 9);

      } else if (current_size <= capacityRHHStatic_10 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(9, 10);

      } else if (current_size <= capacityRHHStatic_11 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(10, 11);

      } else if (current_size <= capacityRHHStatic_12 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(11, 12);

      } else if (current_size <= capacityRHHStatic_13 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(12, 13);

      } else if (current_size <= capacityRHHStatic_14 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(13, 14);

      } else if (current_size <= capacityRHHStatic_15 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(14, 15);

      } else if (current_size <= capacityRHHStatic_16 * kFullCalacityFactor) {
        DELETE_KEY_AND_SHRINK(15, 16);

      } else {
        DELETE_KEY_AND_SHRINK(16, 17);

      }

    }

    void get_elemnts_array(const uint64_t current_size, KeyType *key_array, ValueType *val_array)
    {
      if (current_size <= capacityRHHStatic_1 * kFullCalacityFactor) {
        RHHStaticNoVal_1* rhh = reinterpret_cast<RHHStaticNoVal_1*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_2 * kFullCalacityFactor) {
        RHHStaticNoVal_2* rhh = reinterpret_cast<RHHStaticNoVal_2*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_3 * kFullCalacityFactor) {
        RHHStaticNoVal_3* rhh = reinterpret_cast<RHHStaticNoVal_3*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_4 * kFullCalacityFactor) {
        RHHStaticNoVal_4* rhh = reinterpret_cast<RHHStaticNoVal_4*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_5 * kFullCalacityFactor) {
        RHHStaticNoVal_5* rhh = reinterpret_cast<RHHStaticNoVal_5*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_6 * kFullCalacityFactor) {
        RHHStaticNoVal_6* rhh = reinterpret_cast<RHHStaticNoVal_6*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_7 * kFullCalacityFactor) {
        RHHStaticNoVal_7* rhh = reinterpret_cast<RHHStaticNoVal_7*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_8 * kFullCalacityFactor) {
        RHHStaticNoVal_8* rhh = reinterpret_cast<RHHStaticNoVal_8*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_9 * kFullCalacityFactor) {
        RHHStaticNoVal_9* rhh = reinterpret_cast<RHHStaticNoVal_9*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_10 * kFullCalacityFactor) {
        RHHStaticNoVal_10* rhh = reinterpret_cast<RHHStaticNoVal_10*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_11 * kFullCalacityFactor) {
        RHHStaticNoVal_11* rhh = reinterpret_cast<RHHStaticNoVal_11*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_12 * kFullCalacityFactor) {
        RHHStaticNoVal_12* rhh = reinterpret_cast<RHHStaticNoVal_12*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_13 * kFullCalacityFactor) {
        RHHStaticNoVal_13* rhh = reinterpret_cast<RHHStaticNoVal_13*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_14 * kFullCalacityFactor) {
        RHHStaticNoVal_14* rhh = reinterpret_cast<RHHStaticNoVal_14*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_15 * kFullCalacityFactor) {
        RHHStaticNoVal_15* rhh = reinterpret_cast<RHHStaticNoVal_15*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else if (current_size <= capacityRHHStatic_16 * kFullCalacityFactor) {
        RHHStaticNoVal_16* rhh = reinterpret_cast<RHHStaticNoVal_16*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);

      } else {
        RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
        rhh->get_elemnts_array(key_array, val_array);
      }

    }


/// XXX: should use template function
#define FPRINT_KEYS(C, PRFX, OF) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->fprint_keys(PRFX, OF); \
  } while (0)

  void fprint_keys(uint64_t size, std::string prefix, std::ofstream& output_file)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      FPRINT_KEYS(1, prefix, output_file);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor){
      FPRINT_KEYS(2, prefix, output_file);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor){
      FPRINT_KEYS(3, prefix, output_file);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor){
      FPRINT_KEYS(4, prefix, output_file);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor){
      FPRINT_KEYS(5, prefix, output_file);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor){
      FPRINT_KEYS(6, prefix, output_file);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor){
      FPRINT_KEYS(7, prefix, output_file);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor){
      FPRINT_KEYS(8, prefix, output_file);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor){
      FPRINT_KEYS(9, prefix, output_file);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor){
      FPRINT_KEYS(10, prefix, output_file);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor){
      FPRINT_KEYS(11, prefix, output_file);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor){
      FPRINT_KEYS(12, prefix, output_file);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor){
      FPRINT_KEYS(13, prefix, output_file);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor){
      FPRINT_KEYS(14, prefix, output_file);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor){
      FPRINT_KEYS(15, prefix, output_file);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor){
      FPRINT_KEYS(16, prefix, output_file);
    } else {
      FPRINT_KEYS(17, prefix, output_file);
    }
  }

/// XXX: should use template function
#define DISP_ELEMS(C) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->disp_elems(0); \
  } while (0)

  void disp_elems(uint64_t size)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      DISP_ELEMS(1);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor){
      DISP_ELEMS(2);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor){
      DISP_ELEMS(3);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor){
      DISP_ELEMS(4);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor){
      DISP_ELEMS(5);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor){
      DISP_ELEMS(6);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor){
      DISP_ELEMS(7);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor){
      DISP_ELEMS(8);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor){
      DISP_ELEMS(9);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor){
      DISP_ELEMS(10);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor){
      DISP_ELEMS(11);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor){
      DISP_ELEMS(12);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor){
      DISP_ELEMS(13);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor){
      DISP_ELEMS(14);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor){
      DISP_ELEMS(15);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor){
      DISP_ELEMS(16);
    } else {
      DISP_ELEMS(17);
    }
  }

  uint64_t cal_depth(uint64_t size)
  {
    if (size <= capacityRHHStatic_16){
      return 0;
    } else {
      RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      return rhh->cal_depth();
    }
  }

/// XXX: should use template function
#define FPRINT_PROBEDISTANCE(C, OF) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->fprint_probedistance(OF); \
  } while (0)

  void fprint_probedistance(uint64_t size, std::ofstream& output_file)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      FPRINT_PROBEDISTANCE(1, output_file);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(2, output_file);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(3, output_file);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(4, output_file);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(5, output_file);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(6, output_file);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(7, output_file);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(8, output_file);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(9, output_file);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(10, output_file);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(11, output_file);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(12, output_file);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(13, output_file);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(14, output_file);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(15, output_file);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor){
      FPRINT_PROBEDISTANCE(16, output_file);
    } else {
      FPRINT_PROBEDISTANCE(17, output_file);
    }
  }


private:
  ///  ------------------------------------------------------ ///
  ///              Private
  ///  ------------------------------------------------------ ///

  // ---------  static variables ------------ ///
  static constexpr double kFullCalacityFactor  = 0.9f; /// Must be bigger than 0.5f


    ///  ------------------------------------------------------ ///
    ///              Private Member Functions
    ///  ------------------------------------------------------ ///

#define INSERT_AND_CHECK_CAPACITY(CUR_C, NEW_C) \
    do { \
        RHHStaticNoVal_##CUR_C* rhh = reinterpret_cast<RHHStaticNoVal_##CUR_C*>(m_ptr_);\
        if (!rhh->try_unique_key_insertion(key, val)) {\
          return kDuplicated;\
        }\
        if (current_size+1 > capacityRHHStatic_##CUR_C * kFullCalacityFactor) {\
          grow_rhh_static(allocators, current_size, false);\
          RHHStaticNoVal_##NEW_C* rhh_new = reinterpret_cast<RHHStaticNoVal_##NEW_C*>(m_ptr_);\
          return rhh_new->insert(key, val);\
        }\
        const UpdateErrors err = rhh->insert(key, val);\
        if (err == kLongProbedistance) {\
          m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##CUR_C.allocate(1).get());\
          RHHStaticNoVal_##CUR_C* new_rhh = reinterpret_cast<RHHStaticNoVal_##CUR_C*>(m_ptr_);\
          new_rhh->clear(false);\
          new_rhh->m_next_ = rhh;\
          for(int i = 0; i < capacityRHHStatic_##CUR_C; ++i) {\
            if (rhh->is_valid(i) && rhh->is_longprobedistance(i)) {\
              const UpdateErrors err2 = new_rhh->insert(rhh->m_key_block_[i], rhh->m_value_block_[i]);\
              assert(err2 != kLongProbedistance);\
              rhh->delete_key(i);\
            }\
          }\
        }\
    } while (0)

    /// TODO: use binary tree ?
    UpdateErrors insert_uniquely_helper(AllocatorsHolder &allocators, KeyType& key, ValueType& val, const uint64_t current_size)
    {
      if (current_size <= capacityRHHStatic_1 * kFullCalacityFactor) {
        assert(false);
        INSERT_AND_CHECK_CAPACITY(1, 2);

      } else if (current_size <= capacityRHHStatic_2 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(2, 3);

      } else if (current_size <= capacityRHHStatic_3 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(3, 4);

      } else if (current_size <= capacityRHHStatic_4 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(4, 5);

      } else if (current_size <= capacityRHHStatic_5 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(5, 6);

      } else if (current_size <= capacityRHHStatic_6 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(6, 7);

      } else if (current_size <= capacityRHHStatic_7 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(7, 8);

      } else if (current_size <= capacityRHHStatic_8 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(8, 9);

      } else if (current_size <= capacityRHHStatic_9 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(9, 10);

      } else if (current_size <= capacityRHHStatic_10 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(10, 11);

      } else if (current_size <= capacityRHHStatic_11 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(11, 12);

      } else if (current_size <= capacityRHHStatic_12 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(12, 13);

      } else if (current_size <= capacityRHHStatic_13 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(13, 14);

      } else if (current_size <= capacityRHHStatic_14 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(14, 15);

      } else if (current_size <= capacityRHHStatic_15 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(15, 16);

      } else if (current_size <= capacityRHHStatic_16 * kFullCalacityFactor) {
        INSERT_AND_CHECK_CAPACITY(16, 17);

      } else {
        RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
        if (!rhh->try_unique_key_insertion(key, val))
          return kDuplicated;
        if (current_size % (uint64_t)(capacityRHHStatic_17 * kFullCalacityFactor) == 0) {
          grow_rhh_static(allocators, current_size, false);
          RHHStaticNoVal_17* rhh_new = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
          return rhh_new->insert(key, val);
        }
        const UpdateErrors err = rhh->insert(key, val);
        if (err == kLongProbedistance) {
          m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_17.allocate(1).get());
          RHHStaticNoVal_17* new_rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
          new_rhh->clear(false);
          new_rhh->m_next_ = rhh;
          for(int i = 0; i < capacityRHHStatic_17; ++i) {
            if (rhh->is_valid(i) && rhh->is_longprobedistance(i)) {
              const UpdateErrors err2 = new_rhh->insert(rhh->m_key_block_[i], rhh->m_value_block_[i]);
              assert(err2 != kLongProbedistance);
              rhh->delete_key(i);
            }
          }
        }
      }
      return kSucceed;
    }


#define INSERT_WITHOUT_CAPACITY_CHECKING(CUR_C, NEW_C) \
    do { \
        RHHStaticNoVal_##CUR_C* rhh = reinterpret_cast<RHHStaticNoVal_##CUR_C*>(m_ptr_);\
        UpdateErrors err = rhh->insert(key, val);\
        if (err == kLongProbedistance) {\
          m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##CUR_C.allocate(1).get());\
          RHHStaticNoVal_##CUR_C* new_rhh = reinterpret_cast<RHHStaticNoVal_##CUR_C*>(m_ptr_);\
          new_rhh->clear(false);\
          new_rhh->m_next_ = rhh;\
          for(int i = 0; i < capacityRHHStatic_##CUR_C; ++i) {\
            if (rhh->is_valid(i) && rhh->is_longprobedistance(i)) {\
              const UpdateErrors err = new_rhh->insert(rhh->m_key_block_[i], rhh->m_value_block_[i]);\
              assert(err != kLongProbedistance);\
              rhh->delete_key(i);\
            }\
          }\
        }\
    } while (0)

    /// TODO: use binary tree ?
    UpdateErrors insert_helper(AllocatorsHolder &allocators, KeyType& key, ValueType& val, const uint64_t current_size)
    {
      if (current_size <= capacityRHHStatic_1 * kFullCalacityFactor) {
        assert(false);
        INSERT_WITHOUT_CAPACITY_CHECKING(1, 2);

      } else if (current_size <= capacityRHHStatic_2 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(2, 3);

      } else if (current_size <= capacityRHHStatic_3 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(3, 4);

      } else if (current_size <= capacityRHHStatic_4 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(4, 5);

      } else if (current_size <= capacityRHHStatic_5 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(5, 6);

      } else if (current_size <= capacityRHHStatic_6 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(6, 7);

      } else if (current_size <= capacityRHHStatic_7 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(7, 8);

      } else if (current_size <= capacityRHHStatic_8 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(8, 9);

      } else if (current_size <= capacityRHHStatic_9 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(9, 10);

      } else if (current_size <= capacityRHHStatic_10 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(10, 11);

      } else if (current_size <= capacityRHHStatic_11 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(11, 12);

      } else if (current_size <= capacityRHHStatic_12 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(12, 13);

      } else if (current_size <= capacityRHHStatic_13 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(13, 14);

      } else if (current_size <= capacityRHHStatic_14 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(14, 15);

      } else if (current_size <= capacityRHHStatic_15 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(15, 16);

      } else if (current_size <= capacityRHHStatic_16 * kFullCalacityFactor) {
        INSERT_WITHOUT_CAPACITY_CHECKING(16, 17);

      } else {
        RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
        const UpdateErrors err = rhh->insert(key, val);
        if (err == kLongProbedistance) {
          m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_17.allocate(1).get());
          RHHStaticNoVal_17* new_rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
          new_rhh->clear(false);
          new_rhh->m_next_ = rhh;
          for(int i = 0; i < capacityRHHStatic_17; ++i) {
            if (rhh->is_valid(i) && rhh->is_longprobedistance(i)) {
              const UpdateErrors err2 = new_rhh->insert(rhh->m_key_block_[i], rhh->m_value_block_[i]);
              assert(err2 != kLongProbedistance);
              rhh->delete_key(i);
            }
          }
        }
      }
      return kSucceed;
    }
    /// XXX: should use template function
  #define ALLOCATE_AND_MOVE(OLD_C, NEW_C) do{\
    RHHStaticNoVal_##OLD_C* old_rhh = reinterpret_cast<RHHStaticNoVal_##OLD_C*>(m_ptr_);\
    m_ptr_ = reinterpret_cast<void*>(allocators.allocator_rhh_noval_##NEW_C.allocate(1).get());\
    RHHStaticNoVal_##NEW_C* new_rhh = reinterpret_cast<RHHStaticNoVal_##NEW_C*>(m_ptr_);\
    new_rhh->m_next_ = nullptr;\
    new_rhh->clear(false);\
    while (1) {\
      RHHStaticNoVal_##OLD_C* old_next_rhh = old_rhh->m_next_;\
      for (uint64_t i = 0; i < capacityRHHStatic_##OLD_C; ++i) {\
        if (old_rhh->is_valid(i)) {\
          new_rhh->insert(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);\
        }\
      }\
      allocators.allocator_rhh_noval_##OLD_C.deallocate(bip::offset_ptr<RHHStaticNoVal_##OLD_C>(old_rhh), 1);\
      if (old_next_rhh == nullptr) {\
        break;\
      }\
      old_rhh = old_next_rhh;\
    }\
  }while(0)

  /// TODO: use binary tree ?
  void grow_rhh_static(AllocatorsHolder &allocators, uint64_t size, bool has_long_probedistance)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      assert(false);
      ALLOCATE_AND_MOVE(1, 2);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor){
      ALLOCATE_AND_MOVE(2, 3);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor){
      ALLOCATE_AND_MOVE(3, 4);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor){
      ALLOCATE_AND_MOVE(4, 5);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor){
      ALLOCATE_AND_MOVE(5, 6);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(6, 7);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(7, 8);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(8, 9);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(9, 10);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(10, 11);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(11, 12);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(12, 13);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(13, 14);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(14, 15);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor) {
      ALLOCATE_AND_MOVE(15, 16);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor) {
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
            new_rhh->insert(old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);
            old_rhh->delete_key(i);
          }
        }
      }

    }

  }

  /// TODO: use binary tree ?
  void free_buffer_rhh_static(AllocatorsHolder &allocators, uint64_t size)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_1.deallocate(bip::offset_ptr<RHHStaticNoVal_1>(reinterpret_cast<RHHStaticNoVal_1*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_2.deallocate(bip::offset_ptr<RHHStaticNoVal_2>(reinterpret_cast<RHHStaticNoVal_2*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_3.deallocate(bip::offset_ptr<RHHStaticNoVal_3>(reinterpret_cast<RHHStaticNoVal_3*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_4.deallocate(bip::offset_ptr<RHHStaticNoVal_4>(reinterpret_cast<RHHStaticNoVal_4*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_5.deallocate(bip::offset_ptr<RHHStaticNoVal_5>(reinterpret_cast<RHHStaticNoVal_5*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_6.deallocate(bip::offset_ptr<RHHStaticNoVal_6>(reinterpret_cast<RHHStaticNoVal_6*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_7.deallocate(bip::offset_ptr<RHHStaticNoVal_7>(reinterpret_cast<RHHStaticNoVal_7*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_8.deallocate(bip::offset_ptr<RHHStaticNoVal_8>(reinterpret_cast<RHHStaticNoVal_8*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_9.deallocate(bip::offset_ptr<RHHStaticNoVal_9>(reinterpret_cast<RHHStaticNoVal_9*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_10.deallocate(bip::offset_ptr<RHHStaticNoVal_10>(reinterpret_cast<RHHStaticNoVal_10*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_11.deallocate(bip::offset_ptr<RHHStaticNoVal_11>(reinterpret_cast<RHHStaticNoVal_11*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_12.deallocate(bip::offset_ptr<RHHStaticNoVal_12>(reinterpret_cast<RHHStaticNoVal_12*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_13.deallocate(bip::offset_ptr<RHHStaticNoVal_13>(reinterpret_cast<RHHStaticNoVal_13*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_14.deallocate(bip::offset_ptr<RHHStaticNoVal_14>(reinterpret_cast<RHHStaticNoVal_14*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_15.deallocate(bip::offset_ptr<RHHStaticNoVal_15>(reinterpret_cast<RHHStaticNoVal_15*>(m_ptr_)), 1);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor) {
      allocators.allocator_rhh_noval_16.deallocate(bip::offset_ptr<RHHStaticNoVal_16>(reinterpret_cast<RHHStaticNoVal_16*>(m_ptr_)), 1);
    } else {

      RHHStaticNoVal_17 *rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      while (1) {
        void* next = rhh->m_next_;
        allocators.allocator_rhh_noval_17.deallocate(bip::offset_ptr<RHHStaticNoVal_17>(rhh), 1);
        if (next == nullptr) {
          break;
        }
        rhh = reinterpret_cast<RHHStaticNoVal_17*>(next);
      }

    }
  }


  void* m_ptr_;
};
}

#endif
