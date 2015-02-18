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
    allocate_rhh_static<Allocator_RHH_NoVal_2, RHHStaticNoVal_2>(allocators.allocator_rhh_noval_2);
  }

  /// --- Explicitly Deleted Copy and Move Functions -- ///
  ~RHHMgrStatic() =delete;
  RHHMgrStatic(const RHHMgrStatic&) =delete;
  RHHMgrStatic& operator=(const RHHMgrStatic&) =delete;
  RHHMgrStatic(const RHHMgrStatic&&) =delete;
  RHHMgrStatic& operator=(const RHHMgrStatic&&) =delete;

  ///  ------------------------------------------------------ ///
  ///              Public function
  ///  ------------------------------------------------------ ///
  /// insert a key into RHHStatic class
  inline bool insert_uniquely(AllocatorsHolder& allocators, KeyType& key, ValueType& val, const uint64_t current_size)
  {
    UpdateErrors err = insert_uniquely_capacity_selector(allocators, key, val, current_size);
    return (err != kDuplicated);
  }

  inline void free(AllocatorsHolder &allocators, uint64_t size)
  {
    dealocate_rhh_static_capacity_selector(allocators, size);
  }

  inline bool delete_key(AllocatorsHolder& allocators, KeyType& key, uint64_t current_size)
  {
    return delete_key_capacity_selector(allocators, key, current_size);
  }

  inline void get_elemnts_array(const uint64_t current_size, KeyType *key_array, ValueType *val_array)
  {
    get_elemnts_array_capacity_selector(current_size, key_array, val_array);
  }

  bool get_next_key(int64_t *current_key_pos, const uint64_t current_size, KeyType *key)
  {
    return get_next_key_select_capacity(current_key_pos, current_size, key);
  }

  bool get_next_key_select_capacity(int64_t *current_key_pos, const uint64_t size, KeyType *key)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      RHHStaticNoVal_1 *rhh = reinterpret_cast<RHHStaticNoVal_1*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor) {
      RHHStaticNoVal_2 *rhh = reinterpret_cast<RHHStaticNoVal_2*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor) {
      RHHStaticNoVal_3 *rhh = reinterpret_cast<RHHStaticNoVal_3*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor) {
      RHHStaticNoVal_4 *rhh = reinterpret_cast<RHHStaticNoVal_4*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor) {
      RHHStaticNoVal_5 *rhh = reinterpret_cast<RHHStaticNoVal_5*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor) {
      RHHStaticNoVal_6 *rhh = reinterpret_cast<RHHStaticNoVal_6*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor) {
      RHHStaticNoVal_7 *rhh = reinterpret_cast<RHHStaticNoVal_7*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor) {
      RHHStaticNoVal_8 *rhh = reinterpret_cast<RHHStaticNoVal_8*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor) {
      RHHStaticNoVal_9 *rhh = reinterpret_cast<RHHStaticNoVal_9*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor) {
      RHHStaticNoVal_10 *rhh = reinterpret_cast<RHHStaticNoVal_10*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor) {
      RHHStaticNoVal_11 *rhh = reinterpret_cast<RHHStaticNoVal_11*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor) {
      RHHStaticNoVal_12 *rhh = reinterpret_cast<RHHStaticNoVal_12*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor) {
      RHHStaticNoVal_13 *rhh = reinterpret_cast<RHHStaticNoVal_13*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor) {
      RHHStaticNoVal_14 *rhh = reinterpret_cast<RHHStaticNoVal_14*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor) {
      RHHStaticNoVal_15 *rhh = reinterpret_cast<RHHStaticNoVal_15*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor) {
      RHHStaticNoVal_16 *rhh = reinterpret_cast<RHHStaticNoVal_16*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

    } else {
      RHHStaticNoVal_17 *rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      return rhh->get_next_key(current_key_pos, key);

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

// XXX: should use template function
#define CAL_AVERAGE_PROBEDISTANCE(C) \
  do{ \
    RHHStaticNoVal_##C* rhh = reinterpret_cast<RHHStaticNoVal_##C*>(m_ptr_); \
    rhh->cal_average_probedistance(total, n); \
  } while (0)

  void cal_average_probedistance(uint64_t size, uint64_t* total, uint64_t* n)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      CAL_AVERAGE_PROBEDISTANCE(1);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(2);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(3);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(4);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(5);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(6);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(7);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(8);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(9);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(10);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(11);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(12);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(13);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(14);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(15);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor){
      CAL_AVERAGE_PROBEDISTANCE(16);
    } else {
      CAL_AVERAGE_PROBEDISTANCE(17);
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
  template<typename CurRHHAllocator, typename CurRHH, uint64_t CurCapacity,
  typename GrwRHHAllocator, typename GrwRHH, uint64_t GrwCapacity>
  inline UpdateErrors insert_with_grow(CurRHHAllocator& cur_rhh_allocator, GrwRHHAllocator& grw_rhh_allocator,
                                        const KeyType& key, const ValueType& val, const uint64_t current_size)
  {
    CurRHH* cur_rhh = reinterpret_cast<CurRHH*>(m_ptr_);
    if (!cur_rhh->try_unique_key_insertion(key, val))  return kDuplicated;

    /// --- grow rhh  --- ///
    if (current_size+1 > CurCapacity * kFullCalacityFactor) {
      change_capacity<CurRHHAllocator, CurRHH, CurCapacity, GrwRHHAllocator, GrwRHH, GrwCapacity>(cur_rhh_allocator, grw_rhh_allocator);
      return insert_without_grow<GrwRHHAllocator, GrwRHH, GrwCapacity>(grw_rhh_allocator, key, val);
    } else {
      return insert_without_grow<CurRHHAllocator, CurRHH, CurCapacity>(cur_rhh_allocator, key, val);
    }
  }

  template<typename RHHAllocator, typename RHH, uint64_t Capacity>
  inline UpdateErrors insert_without_grow(RHHAllocator& rhh_allocator, const KeyType& key, const ValueType& val)
  {
    RHH* rhh = reinterpret_cast<RHH*>(m_ptr_);
    KeyType key_long_prbdst;
    ValueType val_long_prbdst;
    const UpdateErrors err = rhh->insert(key, val, &key_long_prbdst, &val_long_prbdst);
    if (err == kLongProbedistance) {
      move_longprobedistance_element<RHHAllocator, RHH, Capacity>(rhh_allocator, rhh, key_long_prbdst, val_long_prbdst);
    }
    return err;
  }

  template<typename RHHAllocator, typename RHH, uint64_t Capacity>
  inline void move_longprobedistance_element(RHHAllocator& rhh_allocator, RHH* rhh, KeyType& key_long_prbdst, ValueType& val_long_prbdst)
  {
    RHH* new_rhh = allocate_rhh_static<RHHAllocator, RHH>(rhh_allocator);
    new_rhh->m_next_ = rhh;
    new_rhh->insert(key_long_prbdst, val_long_prbdst);
  }

  template<typename CurRHHAllocator, typename CurRHH, uint64_t CurCapacity,
  typename SrkRHHAllocator, typename SrkRHH, uint64_t SrkCapacity>
  inline bool delete_key_and_shrink(AllocatorsHolder &allocators, CurRHHAllocator& cur_rhh_allocator, SrkRHHAllocator& srk_rhh_allocator,
    const KeyType& key, const uint64_t current_size)
  {
    /// --- Delete target element --- ///
    CurRHH* cur_rhh = reinterpret_cast<CurRHH*>(m_ptr_);
    if (!cur_rhh->delete_key(key)) return false;

    /// --- If current size is less than previous size capacity, shrink RHH --- ///
    if (current_size - 1 <= SrkCapacity * kFullCalacityFactor) {
      change_capacity<CurRHHAllocator, CurRHH, CurCapacity, SrkRHHAllocator, SrkRHH, SrkCapacity>
                     (cur_rhh_allocator, srk_rhh_allocator);
    }
    return true;
  }

  template<typename CurRHHAllocator, typename CurRHH, uint64_t CurCapacity,
           typename NewRHHAllocator, typename NewRHH, uint64_t NewCapacity>
  inline NewRHH *change_capacity(CurRHHAllocator& cur_rhh_allocator, NewRHHAllocator& new_rhh_allocator)
  {
    CurRHH* cur_rhh = reinterpret_cast<CurRHH*>(m_ptr_);
    NewRHH* new_rhh = allocate_rhh_static<NewRHHAllocator, NewRHH>(new_rhh_allocator);
    extract_old_rhh_elements<CurRHHAllocator, CurRHH, CurCapacity, NewRHHAllocator, NewRHH, NewCapacity>
                            (cur_rhh_allocator, new_rhh_allocator, cur_rhh);
    return new_rhh;
  }

  template<typename OldRHHAllocator, typename OldRHH, uint64_t OldCapacity,
           typename CurRHHAllocator, typename CurRHH, uint64_t CurCapacity>
  inline void extract_old_rhh_elements(OldRHHAllocator& old_rhh_allocator, CurRHHAllocator& cur_rhh_allocator, OldRHH* old_rhh)
  {
    while (old_rhh) {
      OldRHH* next_old_rhh = old_rhh->m_next_;
      for (int i = 0; i < OldCapacity; ++i) {
        if (old_rhh->is_valid(i))
          insert_without_grow<CurRHHAllocator, CurRHH, CurCapacity>(cur_rhh_allocator, old_rhh->m_key_block_[i], old_rhh->m_value_block_[i]);
      }
      deallocate_rhh_static<OldRHHAllocator, OldRHH>(old_rhh_allocator, old_rhh);
      old_rhh = next_old_rhh;
    }
  }

  /// ------------- Allocators, deallocators  ---------- ///

  template<typename RHHAllocator, typename RHH>
  inline RHH *allocate_rhh_static(RHHAllocator& allocator)
  {
    m_ptr_ = reinterpret_cast<void*>(allocator.allocate(1).get());
    RHH* rhh = reinterpret_cast<RHH*>(m_ptr_);
    rhh->m_next_ = nullptr;
    rhh->clear();
    return rhh;
  }

  template<typename RHHAllocator, typename RHH>
  inline void deallocate_chained_rhh_static(RHHAllocator& allocator, void *ptr)
  {
    RHH *rhh = reinterpret_cast<RHH*>(ptr);
    while (rhh) {
      RHH* next = rhh->m_next_;
      allocator.deallocate(bip::offset_ptr<RHH>(rhh), 1);
      rhh = next;
      /// TODO: call deallocate_free_chunks() periodically
    }
  }

  template<typename RHHAllocator, typename RHH>
  inline void deallocate_rhh_static(RHHAllocator& allocator, void *ptr)
  {
    RHH *rhh = reinterpret_cast<RHH*>(ptr);
    allocator.deallocate(bip::offset_ptr<RHH>(rhh), 1);
    /// TODO: call deallocate_free_chunks() periodically
  }


  /// ------------- Capacity selectors  ---------- ///

#define INSERT_WITH_GROW(CUR_C, GRW_C) \
  do {\
    return insert_with_grow<Allocator_RHH_NoVal_ ## CUR_C,\
    RHHStaticNoVal_ ## CUR_C,\
    capacityRHHStatic_ ## CUR_C,\
    Allocator_RHH_NoVal_ ## GRW_C,\
    RHHStaticNoVal_ ## GRW_C,\
    capacityRHHStatic_ ## GRW_C>\
    (allocators.allocator_rhh_noval_ ## CUR_C,\
     allocators.allocator_rhh_noval_ ## GRW_C,\
     key,\
     val,\
     current_size);\
  } while (0)

  /// TODO: use binary tree ?
  UpdateErrors insert_uniquely_capacity_selector(AllocatorsHolder &allocators, KeyType& key, ValueType& val, const uint64_t current_size)
  {
    if (current_size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      assert(false);
      INSERT_WITH_GROW(1, 2);

    } else if (current_size <= capacityRHHStatic_2 * kFullCalacityFactor) {
      INSERT_WITH_GROW(2, 3);

    } else if (current_size <= capacityRHHStatic_3 * kFullCalacityFactor) {
      INSERT_WITH_GROW(3, 4);

    } else if (current_size <= capacityRHHStatic_4 * kFullCalacityFactor) {
      INSERT_WITH_GROW(4, 5);

    } else if (current_size <= capacityRHHStatic_5 * kFullCalacityFactor) {
      INSERT_WITH_GROW(5, 6);

    } else if (current_size <= capacityRHHStatic_6 * kFullCalacityFactor) {
      INSERT_WITH_GROW(6, 7);

    } else if (current_size <= capacityRHHStatic_7 * kFullCalacityFactor) {
      INSERT_WITH_GROW(7, 8);

    } else if (current_size <= capacityRHHStatic_8 * kFullCalacityFactor) {
      INSERT_WITH_GROW(8, 9);

    } else if (current_size <= capacityRHHStatic_9 * kFullCalacityFactor) {
      INSERT_WITH_GROW(9, 10);

    } else if (current_size <= capacityRHHStatic_10 * kFullCalacityFactor) {
      INSERT_WITH_GROW(10, 11);

    } else if (current_size <= capacityRHHStatic_11 * kFullCalacityFactor) {
      INSERT_WITH_GROW(11, 12);

    } else if (current_size <= capacityRHHStatic_12 * kFullCalacityFactor) {
      INSERT_WITH_GROW(12, 13);

    } else if (current_size <= capacityRHHStatic_13 * kFullCalacityFactor) {
      INSERT_WITH_GROW(13, 14);

    } else if (current_size <= capacityRHHStatic_14 * kFullCalacityFactor) {
      INSERT_WITH_GROW(14, 15);

    } else if (current_size <= capacityRHHStatic_15 * kFullCalacityFactor) {
      INSERT_WITH_GROW(15, 16);

    } else if (current_size <= capacityRHHStatic_16 * kFullCalacityFactor) {
      INSERT_WITH_GROW(16, 17);

    } else {
      RHHStaticNoVal_17* rhh = reinterpret_cast<RHHStaticNoVal_17*>(m_ptr_);
      if (!rhh->try_unique_key_insertion(key, val)) return kDuplicated;
      if (current_size % (uint64_t)(capacityRHHStatic_17 * kFullCalacityFactor) == 0) {
        RHHStaticNoVal_17* rhh_new = allocate_rhh_static<Allocator_RHH_NoVal_17, RHHStaticNoVal_17>(allocators.allocator_rhh_noval_17);
        rhh_new->m_next_ = rhh;
        return rhh_new->insert(key, val);
      }
      KeyType key_long_prbdst;
      ValueType val_long_prbdst;
      const UpdateErrors err = rhh->insert(key, val, &key_long_prbdst, &val_long_prbdst);
      if (err == kLongProbedistance) {
        move_longprobedistance_element<Allocator_RHH_NoVal_17, RHHStaticNoVal_17, capacityRHHStatic_17>
                                      (allocators.allocator_rhh_noval_17, rhh, key_long_prbdst, val_long_prbdst);
      }
      return err;
    }

  }


#define DELETE_KEY_AND_SHRINK(SRK_C, CUR_C)\
    do {\
      return delete_key_and_shrink<Allocator_RHH_NoVal_ ## CUR_C,\
      RHHStaticNoVal_ ## CUR_C,\
      capacityRHHStatic_ ## CUR_C,\
      Allocator_RHH_NoVal_ ## SRK_C,\
      RHHStaticNoVal_ ## SRK_C,\
      capacityRHHStatic_ ## SRK_C>\
      (allocators,\
       allocators.allocator_rhh_noval_ ## CUR_C,\
       allocators.allocator_rhh_noval_ ## SRK_C,\
       key,\
       current_size);\
    } while (0)

  /// TODO: use binary tree ?
    bool delete_key_capacity_selector(AllocatorsHolder &allocators, const KeyType& key, uint64_t current_size)
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


  /// TODO: use binary tree ?
  void dealocate_rhh_static_capacity_selector(AllocatorsHolder &allocators, uint64_t size)
  {
    if (size <= capacityRHHStatic_1 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_1, RHHStaticNoVal_1>(allocators.allocator_rhh_noval_1, m_ptr_);
    } else if (size <= capacityRHHStatic_2 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_2, RHHStaticNoVal_2>(allocators.allocator_rhh_noval_2, m_ptr_);
    } else if (size <= capacityRHHStatic_3 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_3, RHHStaticNoVal_3>(allocators.allocator_rhh_noval_3, m_ptr_);
    } else if (size <= capacityRHHStatic_4 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_4, RHHStaticNoVal_4>(allocators.allocator_rhh_noval_4, m_ptr_);
    } else if (size <= capacityRHHStatic_5 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_5, RHHStaticNoVal_5>(allocators.allocator_rhh_noval_5, m_ptr_);
    } else if (size <= capacityRHHStatic_6 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_6, RHHStaticNoVal_6>(allocators.allocator_rhh_noval_6, m_ptr_);
    } else if (size <= capacityRHHStatic_7 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_7, RHHStaticNoVal_7>(allocators.allocator_rhh_noval_7, m_ptr_);
    } else if (size <= capacityRHHStatic_8 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_8, RHHStaticNoVal_8>(allocators.allocator_rhh_noval_8, m_ptr_);
    } else if (size <= capacityRHHStatic_9 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_9, RHHStaticNoVal_9>(allocators.allocator_rhh_noval_9, m_ptr_);
    } else if (size <= capacityRHHStatic_10 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_10, RHHStaticNoVal_10>(allocators.allocator_rhh_noval_10, m_ptr_);
    } else if (size <= capacityRHHStatic_11 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_11, RHHStaticNoVal_11>(allocators.allocator_rhh_noval_11, m_ptr_);
    } else if (size <= capacityRHHStatic_12 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_12, RHHStaticNoVal_12>(allocators.allocator_rhh_noval_12, m_ptr_);
    } else if (size <= capacityRHHStatic_13 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_13, RHHStaticNoVal_13>(allocators.allocator_rhh_noval_13, m_ptr_);
    } else if (size <= capacityRHHStatic_14 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_14, RHHStaticNoVal_14>(allocators.allocator_rhh_noval_14, m_ptr_);
    } else if (size <= capacityRHHStatic_15 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_15, RHHStaticNoVal_15>(allocators.allocator_rhh_noval_15, m_ptr_);
    } else if (size <= capacityRHHStatic_16 * kFullCalacityFactor) {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_16, RHHStaticNoVal_16>(allocators.allocator_rhh_noval_16, m_ptr_);
    } else {
      deallocate_chained_rhh_static<Allocator_RHH_NoVal_17, RHHStaticNoVal_17>(allocators.allocator_rhh_noval_17, m_ptr_);
    }
  }

  void get_elemnts_array_capacity_selector(const uint64_t current_size, KeyType *key_array, ValueType *val_array)
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


    void* m_ptr_;
  };
}

#endif
