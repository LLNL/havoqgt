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

#ifndef HAVOQGT_MPI_RHHSTATIC_HPP_INCLUDED
#define HAVOQGT_MPI_RHHSTATIC_HPP_INCLUDED

#include <fstream>
#include <havoqgt/detail/hash.hpp>
#include "RHHCommon.hpp"

namespace RHH {

  /// --------------------------------------------------------------------------------------- ///
  ///                                RHH static class
  /// --------------------------------------------------------------------------------------- ///
template <typename KeyType, typename ValueType, uint64_t Capacity>
class RHHStatic {

 public:

  /// ---------  Typedefs and Enums ------------ ///
  typedef RHHStatic<KeyType, ValueType, Capacity> RHHStaticType;
  typedef unsigned char PropertyBlockType;
  typedef unsigned char ProbeDistanceType;


  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Explicitly Deleted Functions -- ///
  RHHStatic() =delete;
  ~RHHStatic() =delete;
  RHHStatic(const RHHStatic&) =delete;
  RHHStatic& operator=(const RHHStatic&) =delete;
  RHHStatic(const RHHStatic&&) =delete;
  RHHStatic& operator=(RHHStatic&& old_obj) =delete;

  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///
  inline bool try_unique_key_insertion (const KeyType& key, const ValueType& val) const
  {
    const int64_t pos_key = find_key(key, true);
    return (pos_key == kInvaridIndex);
  }

  inline UpdateErrors insert(KeyType key, ValueType val, KeyType* key_long_prbdst = nullptr, ValueType* val_long_prbdst = nullptr)
  {
    return insert_helper(std::move(key), std::move(val), key_long_prbdst, val_long_prbdst);
  }

  bool delete_key(const KeyType& key)
  {
    ProbeDistanceType dist = 0;
    const int64_t pos = find_key(key, false);
    if (pos != kInvaridIndex) {
      delete_elem(pos);
      return true;
    }

    if (m_next_ != nullptr) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      return next_rhh->delete_key(key);
    }

    return false;
  }

  inline void delete_elem(const int64_t positon)
  {
    property(positon) |= kTombstoneMask;
  }

  inline void clear_elem(const int64_t positon)
  {
    property(positon) = kEmptyValue;
  }

  inline void clear(const bool is_clear_recursively = false)
  {
    for (uint64_t i = 0; i < Capacity; ++i) {
      clear_elem(i);
    }
    if (is_clear_recursively && m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      next_rhh->clear(true);
    }
  }

  inline static uint64_t capacity()
  {
    return Capacity;
  }

  inline bool is_valid(uint64_t pos) const
  {
    return (!is_empty(property(pos)) && !is_deleted(property(pos)));
  }

  inline bool is_longprobedistance(const int64_t positon)
  {
    return extract_probedistance(property(positon)) >= kLongProbedistanceThreshold;
  }

  inline void get_elemnts_array(KeyType *key_array, ValueType *val_array)
  {
    uint64_t pos = 0;
    for (uint64_t i = 0; i < Capacity; ++i) {
      if (is_valid(i)) {
        key_array[pos] = m_key_block_[i];
        val_array[pos] = m_value_block_[i];
        ++pos;
      }
    }
  }

  /// test code
  inline bool get_next_key(int64_t* current_key_pos, KeyType* key)
  {
    for (; *current_key_pos < Capacity; *current_key_pos = *current_key_pos + 1) {
      if (is_valid(*current_key_pos)) {
        *key = m_key_block_[*current_key_pos];
        *current_key_pos = *current_key_pos + 1;
        return true;
      }
    }
    return false;
  }

  uint64_t cal_depth()
  {
    if (m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      return next_rhh->cal_depth() + 1ULL;
    }
    return 1;
  }

  void disp_elems(const uint64_t level)
  {
    // std::cout << ">> --------- disp ---------" << std::endl;
    for (uint64_t i = 0; i < Capacity; ++i) {
      std::cout << "[" << level << " : " << i << " ]" << static_cast<int>(m_property_block_[i]) << " : " << m_key_block_[i] << " : " << static_cast<int>(m_value_block_[i]) << std::endl;
    }
    if (m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      next_rhh->disp_elems(level+1ULL);
    }
    // std::cout << "<< ------------------" << std::endl;
  }

  void fprint_keys(std::string prefix, std::ofstream& output_file)
  {
    for (uint64_t i = 0; i < Capacity; ++i) {
      if (is_valid(i)) {
        output_file << prefix << "\t" << m_key_block_[i] << std::endl;
      }
    }
    if (m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      next_rhh->fprint_keys(prefix, output_file);
    }
  }

  void fprint_probedistance(std::ofstream& output_file)
  {
    for (uint64_t i = 0; i < Capacity; ++i) {
      if (!is_empty(property(i))) {
        uint64_t d = extract_probedistance(property(i));
        output_file << Capacity << " " << d << std::endl;
      }
    }
    if (m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      next_rhh->fprint_probedistance(output_file);
    }
  }

  void cal_average_probedistance(uint64_t* sum , uint64_t* n)
  {
    for (uint64_t i = 0; i < Capacity; ++i) {
      if (!is_empty(property(i))) {
        *sum += extract_probedistance(property(i));
        ++(*n);
      }
    }

    if (m_next_) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      next_rhh->cal_average_probedistance(sum, n);
    }
  }

private:
  /// ---------  Typedefs and Enums ------------ ///
  typedef uint64_t HashType;
  static const PropertyBlockType kTombstoneMask     = 0x80; /// mask value to mark as deleted
  static const PropertyBlockType kProbedistanceMask = 0x7F; /// mask value to extract probe distance
  static const PropertyBlockType kEmptyValue        = 0x7F; /// value repsents cleared space
  static const int64_t kInvaridIndex                = -1LL;
  static const uint64_t kMask = Capacity - 1ULL;
  static const ProbeDistanceType kLongProbedistanceThreshold  = 30;  /// must be less than max value of probedirance

  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///
  /// ------ Private member functions: algorithm core ----- ///
  inline HashType hash_key(const KeyType& key) const
  {
#if 0
    return static_cast<HashType>(key);
#else
    /// The below hash function can work very good for 'sparse' graphs
    // return static_cast<HashType>(havoqgt::detail::hash32(static_cast<uint32_t>(key)));
    using namespace havoqgt::detail;
    return static_cast<HashType>(static_cast<uint64_t>(hash32(key>>32ULL)) << 32ULL | hash32(key));
#endif
  }

  inline int64_t cal_desired_pos(HashType hash) const
  {
    return hash & kMask;
  }

  inline ProbeDistanceType cal_probedistance(HashType hash, const int64_t slot_index) const
  {
    return ((slot_index + (Capacity) - cal_desired_pos(hash)) & kMask);
  }

  inline PropertyBlockType& property(const int64_t ix)
  {
    return m_property_block_[ix];
  }

  inline PropertyBlockType property(const int64_t ix) const
  {
    return const_cast<RHHStaticType*>(this)->property(ix);
  }

  inline PropertyBlockType cal_property(HashType hash, const int64_t slot_index) const
  {
    return cal_probedistance(hash, slot_index);
  }

  inline PropertyBlockType cal_property(const ProbeDistanceType probedistance) const
  {
    return probedistance;
  }

  inline bool is_deleted(const PropertyBlockType prop) const
  {
    return (prop & kTombstoneMask);
  }

  inline bool is_empty(const PropertyBlockType prop) const
  {
    return (prop == kEmptyValue);
  }


  inline ProbeDistanceType extract_probedistance(const PropertyBlockType prop) const
  {
    return prop & kProbedistanceMask;
  }

  inline void insert_into_tempolayspace(KeyType &&key, ValueType &&val)
  {
    for (int64_t i = 0; i < Capacity; ++i) {
      const PropertyBlockType elem_property = property(i);
      if (!is_valid(i)) {
        construct(i, kLongProbedistanceThreshold, std::move(key), std::move(val));
        return;
      }
    }
    assert(false);
  }

  UpdateErrors insert_helper(KeyType &&key, ValueType &&val, KeyType* key_long_prbdst, ValueType* val_long_prbdst)
  {
    int64_t pos = cal_desired_pos(hash_key(key));
    ProbeDistanceType dist = 0;
    UpdateErrors err = kSucceed;

    while(true) {

      PropertyBlockType existing_elem_property = property(pos);

      if(is_empty(existing_elem_property))
      {
        if (dist >= kLongProbedistanceThreshold) {
          *key_long_prbdst = key;
          *val_long_prbdst = val;
          //insert_into_tempolayspace(std::move(key), std::move(val));
          return kLongProbedistance;
        }
        construct(pos, dist, std::move(key), std::move(val));
        break;
      }

      /// If the existing elem has probed equal or less than new, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if (dist >= kLongProbedistanceThreshold) {
          *key_long_prbdst = key;
          *val_long_prbdst = val;
          // insert_into_tempolayspace(std::move(key), std::move(val));
          return kLongProbedistance;
        }
        if(is_deleted(existing_elem_property))
        {
          construct(pos, dist, std::move(key), std::move(val));
          break;
        }
        m_property_block_[pos] = dist;
        std::swap(key, m_key_block_[pos]);
        std::swap(val, m_value_block_[pos]);
        dist = extract_probedistance(existing_elem_property);
      }

      pos = (pos+1) & kMask;
      ++dist;
    }

    return kSucceed;
  }

  inline void construct(const int64_t ix, const ProbeDistanceType dist, KeyType&& key, ValueType&& val)
  {
    m_property_block_[ix] = static_cast<PropertyBlockType>(dist);
    m_key_block_[ix] = std::move(key);
    m_value_block_[ix] = std::move(val);
  }

  /// ------ Private member functions: search ----- ///
  int64_t find_key(const KeyType& key, const bool is_check_recursively) const
  {

    ProbeDistanceType dist = 0;
    const HashType hash = hash_key(key);
    int64_t pos = cal_desired_pos(hash);

    while(true) {

      ProbeDistanceType existing_elem_property = property(pos);
      if (is_empty(existing_elem_property)) { /// free space is found
        break;
      } else if (dist > extract_probedistance(existing_elem_property)) {
        break;
      } else if (!is_deleted(existing_elem_property) && m_key_block_[pos] == key) {
        /// found !
        return pos;
      }
      pos = (pos+1) & kMask;
      ++dist;
    }

    /// Find a key from chained RHH
    if (is_check_recursively && m_next_ != nullptr) {
      RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
      return next_rhh->find_key(key, true);
    }

    return kInvaridIndex;
  }

  /// ------ Private member functions: utility ----- ///


  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
 public: /// TODO: use friend class
  RHHStaticType* m_next_;
  PropertyBlockType m_property_block_[Capacity];
  KeyType m_key_block_[Capacity];
  ValueType m_value_block_[Capacity];

};
}

#endif
