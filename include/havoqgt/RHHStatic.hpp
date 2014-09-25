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

namespace havoqgt {
namespace mpi {

  enum UpdateErrors {
    kSucceed,
    kDuplicated,
    kReachingFUllCapacity,
    kLongProbedistance
  };

namespace bip = boost::interprocess;

/// --------------------------------------------------------------------------------------- ///
///                                RHH static class
/// --------------------------------------------------------------------------------------- ///
template <typename KeyType, uint64_t Capacity>
class RHHStatic {

  public:
  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Constructor --- //
  RHHStatic()
  : m_next_(nullptr)
  {
    DEBUG("RHHStatic constructor");
  }

  /// --- Copy constructor --- //

  /// --- Move constructor --- //
  // XXX: ???
  RHHStatic(RHHStatic &&old_obj)
  {
    DEBUG("RHHStatic move-constructor");
    // m_property_block_ = old_obj.m_property_block_;
    // old_obj.m_property_block_ = NULL;
  }

  /// --- Destructor --- //
  ~RHHStatic()
  {
    DEBUG("RHHStatic destructor");
  }

  /// ---  Move assignment operator --- //
  // XXX: ???
  RHHStatic &operator=(RHHStatic&& old_obj)
  {
    DEBUG("RHHStatic move-assignment");
    // m_pos_head_ = old_obj.m_pos_head_;
    // old_obj.m_pos_head_ = nullptr;
    return *this;
  }


  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///
  UpdateErrors bool insert_uniquely(KeyType&& key)
  {
    ProbeDistanceType dist = 0;
    const uint64_t mask = cal_mask(); // sinse this function is called by parent RHH, we should calvurate mask at here
    const int64_t pos_same_key = find_key(key, &dist, mask);

    if (pos_same_key != kInvaridIndex) {
      return kDuplicated;
    }

    return insert_directly(std::move(key));
  }

  bool erase(KeyType &key)
  {
    ProbeDistanceType dist = 0;
    const int64_t pos = lookup_first_position(key, &dist, cal_mask());
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


 private:
  /// ---------  Typedefs and Enums ------------ ///
  typedef unsigned char PropertyBlockType;
  typedef unsigned char ProbeDistanceType;
  typedef uint64_t HashType;


  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///
  static const PropertyBlockType kTombstoneMask     = 0x80; /// mask value to mark as deleted
  static const PropertyBlockType kProbedistanceMask = 0x7F; /// mask value to extract probe distance
  static const PropertyBlockType kClearedValue      = 0x7F; /// value repsents cleared space
  static const ProbeDistanceType kLongProbedistanceThreshold = 32LL;
  static const int64_t kInvaridIndex = -1LL;


  /// ------ Private member functions: algorithm core ----- ///
  /// XXX: we assume that key can use static_cast to HashType
  inline HashType hash_key(KeyType& key)
  {
    return static_cast<HashType>(key);
  }

  inline uint64_t cal_mask() const
  {
    return Capacity - 1ULL;
  }

  inline int64_t cal_desired_pos(HashType hash, const uint64_t mask)
  {
    return hash & mask;
  }

  inline ProbeDistanceType cal_probedistance(HashType hash, const int64_t slot_index, const uint64_t mask)
  {
    return ((slot_index + (Capacity) - cal_desired_pos(hash, mask)) & mask);
  }

  inline PropertyBlockType& property(const int64_t ix)
  {
    return m_property_block_[ix];
  }

  inline PropertyBlockType property(const int64_t ix) const
  {
    return const_cast<RHHStatic*>(this)->property(ix);
  }

  inline PropertyBlockType cal_property(HashType hash, const int64_t slot_index, const uint64_t mask)
  {
    return cal_probedistance(hash, slot_index, mask);
  }

  inline PropertyBlockType cal_property(const ProbeDistanceType probedistance)
  {
    return probedistance;
  }

  inline void delete_elem(const int64_t positon)
  {
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

  inline static ProbeDistanceType extract_probedistance(PropertyBlockType prop)
  {
    return prop & kProbedistanceMask;
  }

  UpdateErrors insert_helper(KeyType &&key, const uint64_t mask)
  {
    int64_t pos = cal_desired_pos(hash_key(key), mask);
    ProbeDistanceType dist = 0;
    UpdateErrors err;
    while(true) {
      if (dist >= kLongProbedistanceThreshold) {
        err = kLongProbedistance
        break;
      }
      if (dist >= Capacity) {
        err = kReachingFUllCapacity;
        break;
      }

      PropertyBlockType existing_elem_property = property(pos);

      if(existing_elem_property == kClearedValue)
      {
        return construct(pos, std::move(key), mask);
      }

      /// If the existing elem has probed less than or "equal to" us, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if(is_deleted(existing_elem_property))
        {
          return construct(pos, std::move(key), mask);
        }
        m_property_block_[pos] = dist;
        std::swap(key, m_key_block_[pos]);
        dist = extract_probedistance(existing_elem_property);
      }

      pos = (pos+1) & mask;
      ++dist;
    }

    if (m_next_ != nullptr)  {
      RHHStatic<KeyType, ValueType, Capacity> *next_rhh = reinterpret_cast<RHHStatic<KeyType, ValueType, Capacity> *>(m_next_);
      return next_rhh->insert_helper(key, mask);
    }

    return err;
  }

  inline void construct(const int64_t ix, const HashType hash, KeyType&& key, const uint64_t mask)
  {
    m_property_block_[ix] = cal_probedistance(hash, ix, mask);
    m_key_block_[ix] = std::move(key);
  }


  /// ------ Private member functions: search ----- ///
  int64_t find_key(KeyType& key, ProbeDistanceType* dist, const uint64_t mask)
  {
    const HashType hash = hash_key(key);
    int64_t pos = cal_desired_pos(hash, mask);

    while(true) {
      ProbeDistanceType existing_elem_property = property(pos);
      if (existing_elem_property == kClearedValue) { /// free space is found
        break;
      } else if (*dist > extract_probedistance(existing_elem_property)) {
        break;
      } else if (!is_deleted(existing_elem_property) && m_key_block_[pos] == key) {
        /// found !
        return pos;
      }
      pos = (pos+1) & mask;
      *dist = *dist + 1;
    }

    /// Find a key from chained RHH
    if (m_next_ != nullptr)  {
      RHHStatic<KeyType, ValueType, Capacity> *next_rhh = reinterpret_cast<RHHStatic<KeyType, ValueType, Capacity> *>(m_next_);
      return next_rhh->find_key(key, dist, mask);
    }

    return kInvaridIndex;
  }

  /// ------ Private member functions: utility ----- ///


  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
  uint64_t m_num_elems_;
  RHHStatic* m_next_;
  PropertyBlockType m_property_block_[Capacity];
  KeyType m_key_block_[Capacity];
  /// ValueWrapperType m_value_type_[Capacity];

};
} /// namespace mpi
} /// namespace havoqgt

#endif
