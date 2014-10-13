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

#ifndef HAVOQGT_MPI_RHHMAIN_HPP_INCLUDED
#define HAVOQGT_MPI_RHHMAIN_HPP_INCLUDED
#include "RHHCommon.hpp"
#include "RHHAllocHolder.hpp"
#include "RHHMgrStatic.hpp"

 namespace RHH {

template <typename KeyType, typename ValueType>
  class RHHMain {

  public:

  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///
  /// --- Constructor --- //
    RHHMain(AllocatorsHolder &allocators, const uint64_t initial_capasity)
    : m_num_elems_(0)
    , m_capacity_(initial_capasity)
    , m_ptr_(nullptr)
    {
      DEBUG("RHHMain constructor");
      allocate_rhh_main(allocators);
    }

  /// --- Copy constructor --- //
  /// --- Move constructor --- //

  /// --- Destructor --- //
    ~RHHMain()
    {
     DEBUG("RHHMain destructor");
   }

  /// ---  Move assignment operator --- //


  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///
   bool insert_uniquely(AllocatorsHolder& allocators, KeyType key, ValueType val)
   {

    ValueWrapperType value;
    value.value = val;
    const uint64_t mask = cal_mask();
    ProbeDistanceType dist = 0;
    const int64_t pos_key = find_key(key, &dist, mask);
    if (pos_key == kInvaridIndex) {
        /// 1. new vertex
      insert_directly_with_growing(allocators, std::move(key), std::move(value));
      return true;
    }

    const uint64_t size = extract_size(property(pos_key));
    if (size == 1) {
      if (is_equal_element(key, val, pos_key)) return false;
        /// 2. convert data structure to array model
      ValueType* array = reinterpret_cast<ValueType*>(allocators.allocator_normalarray.allocate(1).get());
      array[0] = m_value_block_[pos_key].value;
      array[1] = val;
      m_value_block_[pos_key].value_array = array;
      NormalArrayBitmapType btmp = 0x03;
      set_bitmap(pos_key, btmp);

    } else if (size <= capacityNormalArray3) {
        /// 3. we already have a normal-array
      NormalArrayBitmapType btmp = extract_bitmap(property(pos_key));
      ValueType* array = m_value_block_[pos_key].value_array;
      uint64_t capacity = cal_next_pow2(size);

        /// Check duplicated elements
      for (int i=0; i < capacity; ++i) {
        if (((btmp >> i) & 0xfe) && (array[i] == val)) return false;
      }

      if (size == capacityNormalArray3) {
          /// 3-1. convert data strucure normal-array to RHHStatic
        ValueType* array = m_value_block_[pos_key].value_array;
        ValueWrapperType value_wrapper;
        new(&value_wrapper.adj_list) RHHAdjalistType(allocators);
        uint64_t array_size = 0;
        unsigned char dmy = 0;
        array_size += value_wrapper.adj_list.insert_uniquely(allocators, val, dmy, array_size);
        for (uint64_t i = 0; i < capacityNormalArray3; ++i) {
          array_size += value_wrapper.adj_list.insert_uniquely(allocators, array[i], dmy, array_size);
        }
        allocators.allocator_normalarray.deallocate(bip::offset_ptr<ValueType>(array), capacity);
        set_size(pos_key, array_size);

        return true;
      }

      if (size == capacity) {
          /// grow array
        const int new_capacity = capacity * kCapacityGrowingFactor;
        ValueType* new_array = reinterpret_cast<ValueType*>(allocators.allocator_normalarray.allocate(new_capacity).get());
        memcpy(new_array, array, sizeof(ValueType)*capacity);
        allocators.allocator_normalarray.deallocate(bip::offset_ptr<ValueType>(array), capacity);
        m_value_block_[pos_key].value_array = new_array;
        array = new_array;

      }

        /// 3-2. insert value into normal array
      const int pos = ffs(btmp) - 1;
      array[pos] = val;
      btmp |= (1<<pos);
      set_bitmap(pos_key, btmp);
      set_size(pos_key, size+1);

    } else {
      RHHAdjalistType& rhh_adj_list = m_value_block_[pos_key].adj_list;
      unsigned char dmy = 0;
      bool err = rhh_adj_list.insert_uniquely(allocators, value.value, dmy, size);
      if (err) set_size(pos_key, size+1);

      return err;
    }

    return true;
  }

    // bool erase(KeyType &key)
    // {
    //   ProbeDistanceType dist = 0;
    //   const int64_t pos = find_key(key, &dist, cal_mask());
    //   if (pos == kInvaridIndex) {
    //     return false;
    //   }

    //   delete_key(pos);
    //   --m_num_elems_;
    //   return true;
    // }

  uint64_t size()
  {
    return m_num_elems_;
  }

  // void disp_elems()
  // {

  //   for (uint64_t i = 0; i < m_capacity_; ++i) {
  //     PropertyBlockType prop = property(i);
  //     if (prop == kClearedValue || is_deleted(prop))
  //       constitute;

  //     const uint64_t size = extract_size(prop);
  //     if (size == 1) {
  //       std::cout << m_key_block_[i] << "\t" << m_value_block_[i] << std::endl;
  //     } else if (size <= capacityNormalArray3) {
  //       NormalArrayBitmapType btmp = extract_bitmap(property(pos_key));
  //       ValueType* array = m_value_block_[pos_key].value_array;
  //       uint64_t capacity = cal_next_pow2(size);
  //       for (int j=0; j < capacity; ++j) {
  //         if (((btmp >> j) & 0xfe) && (array[j] == val))
  //           std::cout << m_key_block_[j] << "\t" << array[i] << std::endl;;
  //       }
  //     } else {

  //     }

  //   }

  // }

private:
  ///  ------------------------------------------------------ ///
  ///              Private
  ///  ------------------------------------------------------ ///

  /// ---------  Typedefs and Enums ------------ ///
  typedef RHHMgrStatic<ValueType, NoValueType> RHHAdjalistType;
  typedef uint64_t PropertyBlockType;
  typedef int64_t ProbeDistanceType;
  typedef uint64_t HashType;
  typedef unsigned char NormalArrayBitmapType;

  union ValueWrapperType {
    ValueType value;
    ValueType* value_array;
    RHHAdjalistType adj_list;

    ValueWrapperType()
    {
    //DEBUG("union constructor");
    }

    ValueWrapperType(ValueWrapperType &&old_obj)
    {
      value = old_obj.value;
      old_obj.value = 0;
    }

    ~ValueWrapperType()
    {
    //DEBUG("union destructor");
    }

    ValueWrapperType &operator=(ValueWrapperType &&old_obj)
    {
      value = old_obj.value;
      old_obj.value = 0;
      return *this;
    }

  };

  // ---------  static variables ------------ ///
  /// tombstone (1), size (48), bitmap (8), probe distance (7)
  static const PropertyBlockType kTombstoneMask               = 0x8000000000000000LL; /// mask value to mark as deleted
  static const PropertyBlockType kPropertyDirectInserted      = 0x0000000000008000LL; ///
  static const PropertyBlockType kProbedistanceMask           = 0x000000000000007FLL; /// mask value to extract probe distance
  static const PropertyBlockType kClearProbedistanceMask      = 0xFFFFFFFFFFFFFF80LL; /// mask value to clear probe distance
  static const ProbeDistanceType kLongProbedistanceThreshold  = 128LL;
  //static const PropertyBlockType kSizeMask                    = 0x7FFFFFFFFFFF8000LL; /// mask value to extract size
  static const PropertyBlockType kClearSize                   = 0x8000000000007FFFLL; /// mask value to clear size
  static const PropertyBlockType kClearedValue                = 0x7FFFFFFFFFFFFFFFLL; /// value repsents cleared space
  static const int64_t kInvaridIndex                          = -1LL;
  static constexpr double kFullCalacityFactor                 = 0.9;

  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///

  /// ------ Private member functions: algorithm core ----- ///
  /// XXX: we assume that key can use static_cast to HashType
  inline HashType hash_key(KeyType& key)
  {
    return static_cast<HashType>(key);
  }

  inline uint64_t cal_mask() const
  {
    return m_capacity_ - 1ULL;
  }

  inline int64_t cal_desired_pos(HashType hash, const uint64_t mask)
  {
    return hash & mask;
  }

  inline ProbeDistanceType cal_probedistance(HashType hash, const int64_t slot_index, const uint64_t mask)
  {
    return ((slot_index + (m_capacity_) - cal_desired_pos(hash, mask)) & mask);
  }

  inline PropertyBlockType& property(const int64_t ix)
  {
    return m_property_block_[ix];
  }

  inline PropertyBlockType property(const int64_t ix) const
  {
    return const_cast<RHHMain*>(this)->property(ix);
  }

  inline PropertyBlockType update_probedistance(const PropertyBlockType prop, const ProbeDistanceType probedistance)
  {
    return ((prop & kClearProbedistanceMask) | (probedistance & kProbedistanceMask));
  }

  // inline void delete_elem(const int64_t positon)
  // {
  //   property(positon) |= kTombstoneMask;
  // }

  // inline void delete_key(const int64_t positon)
  // {
  //   property(positon) |= kTombstoneMask;
  // }

  inline static bool is_deleted(const PropertyBlockType prop)
  {
    return (prop & kTombstoneMask) == kTombstoneMask;
  }

  inline static ProbeDistanceType extract_probedistance(PropertyBlockType prop)
  {
    return prop & kProbedistanceMask;
  }

  inline static ProbeDistanceType extract_size(PropertyBlockType prop)
  {
    return (prop >> 15ULL) & 0x0000ffffffffffffLL;
  }

  inline static unsigned char extract_bitmap(PropertyBlockType prop)
  {
    return (prop >> 7ULL) & 0xffULL;
  }

  inline static void set_bit(const int64_t pos, int64_t ix)
  {
    property(pos) = property(pos) | (0x40ULL << ix);
  }

  inline void set_bitmap(const int64_t pos, NormalArrayBitmapType btmp)
  {
    property(pos) |= (0x40ULL << btmp);
  }

  inline void set_size(const int64_t pos, const uint64_t sz)
  {
    property(pos) = (property(pos) & kClearSize) | (sz << 15ULL);
  }

  uint64_t grow_rhh_main(AllocatorsHolder &allocators)
  {
    const uint64_t old_capacity = m_capacity_;
    unsigned char* old_ptr = m_ptr_;
    PropertyBlockType* old_property_block = reinterpret_cast<PropertyBlockType*>(m_ptr_);
    KeyType* old_key_block                = reinterpret_cast<KeyType*>(&m_property_block_[old_capacity]);
    ValueWrapperType* old_value_block     = reinterpret_cast<ValueWrapperType*>(&m_key_block_[old_capacity]);

    m_capacity_ = old_capacity * kCapacityGrowingFactor;
    allocate_rhh_main(allocators);

    /// now copy over old elems
    for (uint64_t i = 0; i < old_capacity; ++i) {
      const PropertyBlockType old_prop = old_property_block[i];
      if (old_prop != kClearedValue && !is_deleted(old_prop)) {
        insert_helper(std::move(old_key_block[i]), std::move(old_value_block[i]), old_prop);
      }
    }

    free_buffer_rhh_main(allocators, old_ptr, old_capacity);

  }

  void allocate_rhh_main(AllocatorsHolder &allocators)
  {
    uint64_t mem_size = 0;
    mem_size += sizeof(PropertyBlockType) * m_capacity_;
    mem_size += sizeof(KeyType) * m_capacity_;
    mem_size += sizeof(ValueWrapperType) * m_capacity_;

    m_ptr_ = reinterpret_cast<unsigned char *>(allocators.allocator_raw.allocate(mem_size).get());
    m_property_block_ = reinterpret_cast<PropertyBlockType*>(m_ptr_);
    m_key_block_      = reinterpret_cast<KeyType*>(&m_property_block_[m_capacity_]);
    m_value_block_    = reinterpret_cast<ValueWrapperType*>(&m_key_block_[m_capacity_]);

    for (uint64_t i = 0; i < m_capacity_; ++i) {
      m_property_block_[i] = kClearedValue;
    }
  }

  void inline free_buffer_rhh_main(AllocatorsHolder &allocators, unsigned char* ptr, uint64_t capacity)
  {
    // uint64_t mem_size = 0;
    // mem_size += sizeof(PropertyBlockType) * capacity;
    // mem_size += sizeof(KeyType) * capacity;
    // mem_size += sizeof(ValueType) * capacity;

    allocators.allocator_raw.deallocate(bip::offset_ptr<unsigned char>(ptr), 0);
  }

  /// ------ Private member functions: insertion ----- ///
  inline void insert_directly_with_growing(AllocatorsHolder &allocators, KeyType&& key, ValueWrapperType&& value)
  {
    ++m_num_elems_;
    UpdateErrors err = insert_helper(std::move(key), std::move(value), kPropertyDirectInserted);

    if (m_num_elems_ >= m_capacity_*kFullCalacityFactor || err == kLongProbedistance) {
      grow_rhh_main(allocators);
    }
  }

  UpdateErrors insert_helper(KeyType &&key, ValueWrapperType &&value, PropertyBlockType prop)
  {
    const uint64_t mask = cal_mask();
    int64_t pos = cal_desired_pos(hash_key(key), mask);
    ProbeDistanceType dist = 0;

    while(true) {

      PropertyBlockType existing_elem_property = property(pos);
      if(existing_elem_property == kClearedValue)
      {
        return construct(pos, std::move(key), std::move(value), prop, dist);
      }

      /// If the existing elem has probed less than or "equal to" us, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if(is_deleted(existing_elem_property))
        {
          return construct(pos, std::move(key), std::move(value), prop, dist);
        }
        m_property_block_[pos] = update_probedistance(prop, dist);
        std::swap(key, m_key_block_[pos]);
        std::swap(value, m_value_block_[pos]);
        dist = extract_probedistance(existing_elem_property);
        prop = existing_elem_property;
      }

      pos = (pos+1) & mask;
      ++dist;
    }
  }

  inline UpdateErrors construct(const int64_t ix, KeyType&& key, ValueWrapperType&& val, const PropertyBlockType prop, const ProbeDistanceType dist)
  {
    m_property_block_[ix] = update_probedistance(prop, dist);
    m_key_block_[ix] = std::move(key);
    m_value_block_[ix] = std::move(val);
    if (dist >= kLongProbedistanceThreshold) return kLongProbedistance;
    else return kSucceed;
  }


  /// ------ Private member functions: search ----- ///

  int64_t find_key(KeyType& key, ProbeDistanceType* dist, const uint64_t mask)
  {
    const HashType hash = hash_key(key);
    int64_t pos = cal_desired_pos(hash, mask);

    while(true) {
      ProbeDistanceType existing_elem_property = property(pos);
      if (existing_elem_property == kClearedValue) { /// free space is found
        return kInvaridIndex;
      } else if (*dist > extract_probedistance(existing_elem_property)) {
        return kInvaridIndex;
      } else if (!is_deleted(existing_elem_property) && m_key_block_[pos] == key) {
        /// found !
        return pos;
      }
      pos = (pos+1) & mask;
      *dist = *dist + 1;
    }
  }

  /// ------ Private member functions: utility ----- ///
  inline bool is_equal_element(const KeyType& key, const ValueType& value, const int64_t pos) const
  {
    if (is_deleted(property(pos)) || m_key_block_[pos] != key || m_value_block_[pos].value != value)
      return false;
    else
      return true;
  }


  inline unsigned char cal_next_pow2(unsigned char n)
  {
    --n;
    n |= (n >> 1);
    n |= (n >> 2);
    n |= (n >> 4);
    ++n;
    return n;
  }

  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
  uint64_t m_num_elems_;
  uint64_t m_capacity_;
  unsigned char* m_ptr_;
  PropertyBlockType* m_property_block_;
  KeyType* m_key_block_;
  ValueWrapperType* m_value_block_;

};

} /// namespace RHH

#endif
