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

namespace RHH {

template <typename KeyType, typename ValueType>
class RHHMain {

  public:

  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Constructor --- //
  RHHMain()
  {
    DEBUG("RHHMain constructor");
  }

  /// --- Copy constructor --- //

  /// --- Move constructor --- //
  RHHMain(RHHMain &&old_obj)
  {
    DEBUG("RHHMain move-constructor");
    // m_property_block_ = old_obj.m_property_block_;
    // old_obj.m_property_block_ = NULL;
  }

  /// --- Destructor --- //
  ~RHHMain()
  {
    DEBUG("RHHMain destructor");
  }

  /// ---  Move assignment operator --- //
  RHHMain &operator=(RHHMain&& old_obj)
  {
    DEBUG("RHHMain move-assignment");
    // m_pos_head_ = old_obj.m_pos_head_;
    // old_obj.m_pos_head_ = nullptr;
    return *this;
  }


  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///
  UpdateErrors insert_uniquely(AllocatorsHolder& allocators, KeyType key, ValueType val)
  {
    // std::cout << key << "," << val << std::endl; //D
    // disp_elements(); //D

    ValueWrapperType value;
    value.value = val;
    const uint64_t mask = cal_mask();
    ProbeDistanceType dist = 0;
    const int64_t pos_key = find_key(key, &dist, mask);
    if (pos_key == kInvaridIndex) {
      /// 1. new vertex
      return insert_directly(std::move(key), std::move(value));
    }

    const uint64_t size = extract_size(property(pos_key));
    if (size == 1) {
      /// FIXME: check duplication
      /// convert data structure
      ValueType* array = allocators->allocate_normal_array(2);
      array[0] = m_value_type_[pos_key].value;
      array[1] = val;
      m_value_type_[pos_key].value_array = array;
    } else if (size < AllocatorsHolder::capacityNormalArray3) {
      /// FIXME: check duplication
      ValueType* array = m_value_type_[pos_key].value_array;
      const int pos = ffs(extract_bitmap(property(pos_key))) - 1;
      array[pos] = val;
      set_bit(pos_key, pos);
    } else if (size == AllocatorsHolder::capacityNormalArray3) {
          /// 4. convert data strucure of value from direct-insertion model to adj-list-insertion model
      //DEBUG("Convert start, current initial alloc size = 2");
      ValueWrapperType value_wrapper;
      new(&value_wrapper.adj_list) RHHMgr(kInitialCapacityAdjlist, false);
      value_wrapper.adj_list.insert_uniquely(allocators, std::move(value.value), extract_capacity(property(pos_key)));
      for (uint64_t k = 0; k < kDirectInsertionThreshold; ++k) {
        value_wrapper.adj_list.insert_uniquely(allocators, std::move(m_value_type_[moved_pos_list[k]].value), extract_capacity(property(pos_key)));
        delete_elem(moved_pos_list[k]);
      }
    /// insert adjacency-list
    construct(pos_key, std::move(key), std::move(value_wrapper), mask);

    } else {
      RHHMgr<ValueType, NoValueType>& rhh_adj_list = m_value_type_[pos_key].adj_list;
      return rhh_adj_list.insert_uniquely(allocators, std::move(value.value), NULL, extract_capacity(property(pos_key)));
    }

    if (is_adj_list(property(pos_key))) {
      /// 2. non-low-degree
      RHHMgr<ValueType, NoValueType>& rhh_adj_list = m_value_type_[pos_key].adj_list;
      return rhh_adj_list.insert_uniquely(allocators, std::move(value.value), extract_capacity(property(pos_key)));
    }

    /// Check whether we already have a same edge
    if (get_pos_equal_element_from_here(key, value.value, pos_key, dist, mask) != kInvaridIndex) {
      return kDuplicated;
    }

    int64_t moved_pos_list[kDirectInsertionThreshold];
    for (int64_t pos = pos_key, count_key = 0;; pos = (pos+1) & mask, ++dist) {
      PropertyBlockType existing_elem_property = property(pos);
      if (existing_elem_property == kClearedValue) {
        /// 3. low-degree, we found a space to insert
        construct(pos, std::move(key), std::move(value), false, mask);
        return kSucceed;
      } else if (dist > extract_probedistance(existing_elem_property)) {
        /// 3. low-degree
        return insert_directly(std::move(key), std::move(value));
      }else if (!is_deleted(existing_elem_property) && m_key_block_[pos] == key) {
        moved_pos_list[count_key] = pos;
        if (++count_key >= kDirectInsertionThreshold) break;
      }
    }



    return kSucceed;
  }

  bool erase(KeyType &key)
  {
    ProbeDistanceType dist = 0;
    const int64_t pos = find_key(key, &dist, cal_mask());
    if (pos == kInvaridIndex) {
      return false;
    }

    delete_key(pos);
    --m_num_elems_;
    return true;
  }

  uint64_t size()
  {
    return m_m_num_elems__;
  }

  void inline clear()
  {
    for (uint64_t k = 0; k < m_capacity_; ++k) {
      m_property_block_[k] = kClearedValue;
    }
  }

  void move_elems_from(void *ptr, uint64_t capacity)
  {

  }


 private:
  ///  ------------------------------------------------------ ///
  ///              Private
  ///  ------------------------------------------------------ ///

  /// ---------  Typedefs and Enums ------------ ///
  typedef RHHMgr<KeyType, RHHMgr::NoValueType> RHHAdjalistType;
  typedef uint64_t PropertyBlockType;
  typedef int64_t ProbeDistanceType;
  typedef uint64_t HashType;

  template<Cap>
  union ValueWrapperType {
    ValueType value;
    ValueType* value_array;
    RHHMgr<KeyType, ValueType> adj_list;

    ValueWrapperType() {
      //DEBUG("union constructor");
    }
    ValueWrapperType(ValueWrapperType &&old_obj)
    {
      //DEBUG("union move-constructor");
      value = old_obj.value;
      old_obj.value = 0;
    }
    ~ValueWrapperType()
    {
      //DEBUG("union destructor");
    }
    ValueWrapperType &operator=(ValueWrapperType &&old_obj)
    {
     // DEBUG("union move-assigment");
      value = old_obj.value;
      old_obj.value = 0;
      return *this;
    }
  };

  // ---------  static variables ------------ ///
  static const uint64_t kDirectInsertionThreshold = 1ULL; /// NOTE: We assume that this value is small (may be less than 10)

  /// tombstone (1), size (48), bitmap (8), probe distance (7)
  static const PropertyBlockType kTombstoneMask     = 0x8000000000000000LL; /// mask value to mark as deleted
  static const PropertyBlockType kPropertyDirectInserted = 0x0000000000001000LL;
  //static const PropertyBlockType kAdjacencylistMask = 0x4000000000000000LL; /// mask value represents whether value is adjacency-list(RHH) or just a value
  static const PropertyBlockType kProbedistanceMask = 0x0000000000000FFFLL; /// mask value to extract probe distance
  static const PropertyBlockType kClearProbedistanceMask = 0xFFFFFFFFFFFFF000LL; /// mask value to clear probe distance
  static const ProbeDistanceType kLongProbedistanceThreshold = 257LL;
  static const PropertyBlockType kCapacityMask      = 0x7FFFFFFFFFFFF000LL; /// mask value to extract probe distance
  static const PropertyBlockType kClearCapacity     = 0x8000000000000FFFLL; /// mask value to extract probe distance

  static const PropertyBlockType kClearedValue      = 0x7FFFFFFFFFFFFFFFLL; /// value repsents cleared space
  static const int64_t kInvaridIndex = -1LL;


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

  /// FIXME
  // inline PropertyBlockType cal_property(const PropertyBlockType prop, HashType hash, const int64_t slot_index, const uint64_t mask)
  // {
  //   return  | cal_probedistance(hash, slot_index, mask);
  // }

  inline PropertyBlockType cal_property(const PropertyBlockType prop, const ProbeDistanceType probedistance)
  {
    return ((prop & kClearProbedistanceMask) | (probedistance & kProbedistanceMask));
  }

  inline void delete_elem(const int64_t positon)
  {
    if (is_adj_list(property(positon))) {
      m_value_type_[positon].adj_list.free(allocator);
      //m_value_type_[positon].adj_list.~RHHMain();
    }
    property(positon) |= kTombstoneMask;
  }

  inline void delete_key(const int64_t positon)
  {
    property(positon) |= kTombstoneMask;
  }

  inline static bool is_deleted(const PropertyBlockType prop)
  {
    return (prop & kTombstoneMask) == kTombstoneMask;
  }

  inline static bool is_adj_list(const PropertyBlockType prop)
  {
    /// if current capacity is "MORE THAN" kDirectInsertionThreshold,
    /// cuurent value is adj-list.
    return (prop & kCapacityMask) > kDirectInsertionThreshold;
  }

  inline static ProbeDistanceType extract_probedistance(PropertyBlockType prop)
  {
    return prop & kProbedistanceMask;
  }

  inline static ProbeDistanceType extract_size(PropertyBlockType prop)
  {
    return prop & kCapacityMask;
  }

  inline static void set_capacity(const int64_t positon, const uint64_t cap)
  {
    property(positon) = (property(positon) & kClearCapacity) | (cap << 12UUL);
  }

  /// ------ Private member functions: insertion ----- ///
  inline UpdateErrors insert_directly(KeyType&& key, ValueWrapperType&& value)
  {
    ++m_num_elems_;
    UpdateErrors err = insert_helper(std::move(key), std::move(value), kPropertyDirectInserted, cal_mask());

    if (m_num_elems_ >= m_capacity_*kResizeFactor) {
      return kReachingFUllCapacity;
    } else {
      return err;
    }

  }

  UpdateErrors insert_helper(KeyType &&key, ValueWrapperType &&value, PropertyBlockType prop, const uint64_t mask)
  {
    int64_t pos = cal_desired_pos(hash_key(key), mask);
    ProbeDistanceType dist = 0;

    while(true) {
      PropertyBlockType existing_elem_property = property(pos);

      if(existing_elem_property == kClearedValue)
      {
        return construct(pos, std::move(key), std::move(value), prop, mask);
      }

      /// If the existing elem has probed less than or "equal to" us, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if(is_deleted(existing_elem_property))
        {
          return construct(pos, std::move(key), std::move(value), prop, mask);
        }
        m_property_block_[pos] = cal_property(prop, dist);
        std::swap(key, m_key_block_[pos]);
        std::swap(value, m_value_type_[pos]);
        dist = extract_probedistance(existing_elem_property);
        prop = propexisting_elem_property;
      }

      pos = (pos+1) & mask;
      ++dist;
    }
  }

  inline UpdateErrors construct(const int64_t ix, const HashType hash, KeyType&& key, ValueWrapperType&& val, const PropertyBlockType prop, const uint64_t mask)
  {
    ProbeDistanceType dist = cal_probedistance(hash, ix, mask);
    m_property_block_[ix] = cal_property(prop, dist);
    m_key_block_[ix] = std::move(key);
    m_value_type_[ix] = std::move(val);
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
  inline void get_all_elem_key_list(KeyType pos_list[])
  {
    uint64_t capacity = m_capacity_;
    uint64_t m_num_elems_ = m_num_elems_;
    for (int64_t pos=0, cout=0; pos < capacity && cout < m_num_elems_; ++pos) {
      PropertyBlockType prop = property(pos);
      if (is_deleted(prop) || prop == kClearedValue) continue;
      pos_list[cout] = m_key_block_[pos];
      ++cout;
    }
  }

  inline int64_t get_pos_equal_element_from_here(const KeyType& key, const ValueType& value, int64_t pos, ProbeDistanceType dist, const uint64_t mask) const
  {
    /// Check whether we have a same key-value pair.
    PropertyBlockType existing_elem_property = property(pos);
    while(true) {
      if (existing_elem_property==kClearedValue || dist > extract_probedistance(existing_elem_property))
        return kInvaridIndex;
      if (is_equal_element(key, value, pos, existing_elem_property))
        return pos;
      pos=(pos+1)&mask;
      ++dist;
      existing_elem_property = property(pos);
    }
  }

  inline bool is_equal_element(const KeyType& key, const ValueType& value, const int64_t pos, const PropertyBlockType prop) const
  {
    if (is_deleted(prop) || m_key_block_[pos] != key || m_value_type_[pos].value != value)
      return false;
    else
      return true;
  }


  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
  uint64_t m_num_elems_;
  uint64_t m_capacity_;
  PropertyBlockType* m_property_block_;
  KeyType* m_key_block_;
  ValueWrapperType* m_value_type_;

};

} /// namespace mpi
} /// namespace havoqgt

#endif
