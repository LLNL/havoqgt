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

#define DEBUG(msg) do { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << msg << std::endl; } while (0)
#define DEBUG2(x) do  { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << #x << " =\t" << x << std::endl; } while (0)
#define DISP_VAR(x) do  { std::cout << #x << " =\t" << x << std::endl; } while (0)

namespace havoqgt {
namespace mpi {

namespace bip = boost::interprocess;

///  =========================================================================== ///
///                             Allocator Holder
///  =========================================================================== ///
class AllocatorsHolder
{
  typedef bip::managed_mapped_file  mapped_t;
  typedef mapped_t::segment_manager segment_manager_t;

  enum {
    capacity1 = 8,
    capacity2 = 16,
    capacity3 = 32,
    capacity4 = 64,
    capacity5 = 128,
    capacity6 = 256,
    max_capacity = 256,
  }

  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity1> RHHStaticNoVal1;
  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity2> RHHStaticNoVal2;
  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity3> RHHStaticNoVal3;
  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity4> RHHStaticNoVal4;
  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity5> RHHStaticNoVal5;
  typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity6> RHHStaticNoVal6;

  /// capacity  memsize # nodes per chunk
  /// 1 11  372.4
  /// 2 20  204.8
  /// 4 38  107.8
  /// 8 74  55.4
  /// 16  146 28.1
  /// 32  290 14.1
  /// 64  578 7.1
  /// 128 1154  3.5
  /// 256 2306  1.8
  /// 512 4610  0.9
  bip::node_allocator<RHHStaticNoVal1, segment_manager_t,  55>  allocator_rhh_noval_1;
  bip::node_allocator<RHHStaticNoVal2, segment_manager_t,  28>  allocator_rhh_noval_2;
  bip::node_allocator<RHHStaticNoVal3, segment_manager_t,  14>  allocator_rhh_noval_3;
  bip::node_allocator<RHHStaticNoVal4, segment_manager_t,   7>  allocator_rhh_noval_4;
  bip::node_allocator<RHHStaticNoVal5, segment_manager_t,   3>  allocator_rhh_noval_5;
  bip::node_allocator<RHHStaticNoVal6, segment_manager_t,   1>  allocator_rhh_noval_6;
  bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_key_array;
  bip::allocator<unsigned char, segment_manager_t> allocator_raw;

  explicit AllocatorsHolder(segment_manager_t segment_manager)
  : allocator_rhh_noval_1(segment_manager)
  , allocator_rhh_noval_2(segment_manager)
  , allocator_rhh_noval_3(segment_manager)
  , allocator_rhh_noval_4(segment_manager)
  , allocator_rhh_noval_5(segment_manager)
  , allocator_rhh_noval_6(segment_manager)
  , allocator_key_array(segment_manager)
  , allocator_raw(segment_manager)
  { }

  void* allocate_rhh_static(const uint64_t required_capacity)
  {
    if (required_capacity < capacity1) {
      return reinterpret_cast<void*>(allocator_key_array.allocate(required_capacity).get());
    } else if (required_capacity <= capacity1){
      return reinterpret_cast<void*>(allocator_rhh_1.allocate(1).get());
    } else if (required_capacity <= capacity2){
      return reinterpret_cast<void*>(allocator_rhh_2.allocate(1).get());
    } else if (required_capacity <= capacity3){
      return reinterpret_cast<void*>(allocator_rhh_3.allocate(1).get());
    } else if (required_capacity <= capacity4){
      return reinterpret_cast<void*>(allocator_rhh_4.allocate(1).get());
    } else if (required_capacity <= capacity5){
      return reinterpret_cast<void*>(allocator_rhh_5.allocate(1).get());
    } else if (required_capacity <= capacity6){
      return reinterpret_cast<void*>(allocator_rhh_6.allocate(1).get());
    } else {
      assert(false);
    }
  }

  /// this function allocate rhh which has 64bit of property array, key block and value block
  void* allocate_rhh_main(const uint64_t required_capacity)
  {
      const uint64_t memsize = RHHMain<KeyType, ValueType>.cal_memsize(required_capacity);
      return reinterpret_cast<void*>(allocator_raw.allocate(memsize).get());
  }

  void deallocate_rhh_static(void* ptr, uint64_t capacity)
  {
    if (required_capacity < capacity1) {
      allocator_value.deallocate(ptr, required_capacity);
    } else if (required_capacity <= capacity1){
      allocator_rhh_1.deallocate(ptr, 1);
    } else if (required_capacity <= capacity2){
      allocator_rhh_2.deallocate(ptr, 1);
    } else if (required_capacity <= capacity3){
      allocator_rhh_3.deallocate(ptr, 1);
    } else if (required_capacity <= capacity4){
      allocator_rhh_4.deallocate(ptr, 1);
    } else if (required_capacity <= capacity5){
      allocator_rhh_5.deallocate(ptr, 1);
    } else if (required_capacity <= capacity6){
      allocator_rhh_6.deallocate(ptr, 1);
    } else {
      assert(false);
    }
  }

  void allocate_rhh_main(void* ptr, uint64_t capacity)
  {
    allocator_raw.deallocate(ptr, capacity);
  }

  static inline uint64_t cal_capacity(uint64_t size)
  {
    return cal_next_highest_power_of_2(size);
  }

 private:
  static inline uint64_t cal_next_highest_power_of_2(uint64_t x)
  {
    --x;
    x != x >> 1ULL;
    x != x >> 2ULL;
    x != x >> 4ULL;
    x != x >> 8ULL;
    x != x >> 16ULL;
    x != x >> 32ULL;
    return ++x;
  }

};


///  =========================================================================== ///
///                            RHH Manager class
///  =========================================================================== ///
template <typename KeyType, typename ValueType>
class RHHMgr {

 public:

  template<C>
  typedef RHHStatic<KeyType, ValueType, C, PropertyBlockType> RHHAdjlistType;
  typedef void* NoValueType;

  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Constructor --- //
  explicit RHHMgr(AllocatorsHolder &allocators, const uint64_t initial_capasity)
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

  inline UpdateErrors insert_helper(AllocatorsHolder& allocators, KeyType& key, ValueType& val, const uint64_t current_size)
  {
    RHHMain<KeyType, ValueType>* rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    return rhh->insert_uniquely(key, val);
  }

  /// insert a key into RHHStatic class
  bool insert_uniquely(AllocatorsHolder& allocators, KeyType& key, uint64_t current_size)
  {
    UpdateErrors err = insert_helper(key, current_size);
    if (err == kDuplicated) return false;

    uint64_t new_capacity;
    if (err == kReachingFUllCapacity) {
      new_capacity = grow(allocators, current_size);
    } else if (err == kLongProbedistance) {
      /// XXX: current implementation dose not allocate new RHH-array
      new_capacity = grow(allocators, current_size);
      DEBUG("kLongProbedistance");
    }
    return true;
  }

  UpdateErrors insert_helper(KeyType& key, const uint64_t current_size)
  {
    if (current_size < AllocatorsHolder::capacity1) {
      return insert_simple_array(key);
    } else if (current_size <= AllocatorsHolder::capacity1) {
      AllocatorsHolder::RHHStaticNoVal1* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal1*>(m_ptr_);
      return rhh->insert_uniquely(key);
    }  else if (current_size <= AllocatorsHolder::capacity2) {
      AllocatorsHolder::RHHStaticNoVal2* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal2*>(m_ptr_);
      return rhh->insert_uniquely(key);
    }  else if (current_size <= AllocatorsHolder::capacity3) {
      AllocatorsHolder::RHHStaticNoVal3* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal3*>(m_ptr_);
      return rhh->insert_uniquely(key);
    }  else if (current_size <= AllocatorsHolder::capacity4) {
      AllocatorsHolder::RHHStaticNoVal4* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal4*>(m_ptr_);
      return rhh->insert_uniquely(key);
    }  else if (current_size <= AllocatorsHolder::capacity5) {
      AllocatorsHolder::RHHStaticNoVal5* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal5*>(m_ptr_);
      return rhh->insert_uniquely(key);
    }  else if (current_size <= AllocatorsHolder::capacity6) {
      AllocatorsHolder::RHHStaticNoVal6* rhh = reinterpret_cast<AllocatorsHolder::RHHStaticNoVal6*>(m_ptr_);
      return rhh->insert_uniquely(key);
    } else {
      assert(1);
    }
  }


  // inline bool erase(AllocatorsHolder &allocators, KeyType key, ValueType val, uint64_t* current_size)
  // {
  //   uint64_t capacity = AllocatorsHolder.cal_next_highest_power_of_2(*current_size);
  //   RHHAdjlistType<capacity>* rhh = reinterpret_cast<RHHAdjlistType<capacity>*>(m_ptr_);
  //   return rhh->erase(allocators, key, val);
  // }

  // inline bool erase(AllocatorsHolder &allocators, KeyType key, uint64_t* current_size)
  // {
  //   uint64_t capacity = AllocatorsHolder.cal_next_highest_power_of_2(*current_size);
  //   RHHAdjlistType<capacity>* rhh = reinterpret_cast<RHHAdjlistType<capacity>*>(m_ptr_);
  //   return rhh->erase(allocators, key);
  // }


  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///
  static const uint64_t kCapacityGrowingFactor = 2ULL;

  uint64_t grow_rhh_main(AllocatorsHolder &allocators)
  {
    RHHMain<KeyType, ValueType>* old_rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);

    uint64_t capacity = old_rhh->capacity; * kCapacityGrowingFactor;
    alloc(allocators, capacity);
    /// now copy over old elems
    RHHMain<KeyType, ValueType>* rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    rhh->move_elems_from(old_rhh);
    free_buffer_rhh_main(allocators, reinterpret_cast<void*>(old_rhh), old_rhh->capacity;);

    return capacity;
  }

  void allocate_rhh_main(AllocatorsHolder &allocators, const uint64_t capacity)
  {
    m_ptr_ = allocators.allocate_rhh_main(capacity);
    RHHMain<KeyType, ValueType>* rhh = reinterpret_cast<RHHMain<KeyType, ValueType>*>(m_ptr_);
    rhh->reset_property_block();
  }

  uint64_t grow_rhh_static(AllocatorsHolder &allocators)
  {
    /// Depends on capacity, expand current rhh array OR allocate new rhh array and make chain.
    ///
  }

  void allocate_rhh_static(AllocatorsHolder &allocators, const uint64_t capacity)
  {
    ///
    ///
  }

  void inline free_buffer_rhh_main(AllocatorsHolder &allocators, void* ptr, uint64_t capacity)
  {
    allocators.deallocate_raw(ptr, capacity);
  }

  void inline free_buffer_rhh_static(AllocatorsHolder &allocators, void* ptr, uint64_t capacity)
  {
    allocators.deallocate_rhh_static(ptr, capacity);
  }


  void* m_ptr_;
};

} /// namespace mpi
} /// namespace havoqgt

#endif