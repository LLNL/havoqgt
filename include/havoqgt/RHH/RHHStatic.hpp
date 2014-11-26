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

    /// --- Constructor --- //
    RHHStatic()
    : m_next_(nullptr)
    {
      clear(false);
    }

    /// --- Copy constructor --- //

    /// --- Move constructor --- //
    RHHStatic(RHHStatic&& old_obj)
    {
      assert(false); ///
    }

    /// --- Destructor --- //
    ~RHHStatic()
    {

    }

    /// ---  Move assignment operator --- //
    RHHStatic &operator=(RHHStatic&& old_obj)
    {
      assert(false);
      return *this;
    }


    ///  ------------------------------------------------------ ///
    ///              Public Member Functions
    ///  ------------------------------------------------------ ///
    inline bool try_unique_key_insertion (KeyType& key, ValueType& val)
    {
      const int64_t pos_key = find_key(key, true);
      if (pos_key == kInvaridIndex)
        return true;
      else
        return false;
    }

    inline UpdateErrors insert_uniquely(KeyType key, ValueType val)
    {
      return insert_helper(std::move(key), std::move(val));
    }


    // UpdateErrors insert_uniquely(KeyType key, ValueType val)
    // {
    //   // if (m_next_ != nullptr)
    //   //   std::cout << "!!" << std::hex << m_next_ << std::endl; //D

    //   const int64_t pos_key = find_key(key, true);
    //   if (pos_key != kInvaridIndex) {
    //     return kDuplicated;
    //   }

    //   return insert_helper(std::move(key), std::move(val));
    // }

    bool erase(KeyType key)
    {
      ProbeDistanceType dist = 0;
      const int64_t pos = find_key(key, false);
      if (pos != kInvaridIndex) {
        delete_key(pos);
        return true;
      }

      if (m_next_ != nullptr) {
        RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
        return next_rhh->erase(key);
      }

      return false;
    }

    inline void clear(bool is_clear_recursively)
    {
      for (uint64_t i = 0; i < Capacity; ++i) {
        m_property_block_[i] = kEmptyValue;
      }
      if (is_clear_recursively && m_next_) {
        RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
        next_rhh->clear(true);
      }
    }

    inline void delete_key(const int64_t positon)
    {
      property(positon) |= kTombstoneMask;
    }

    inline static uint64_t capacity()
    {
      return Capacity;
    }

    inline bool is_valid(uint64_t pos) const
    {
      return ( (property(pos) != kEmptyValue) && !is_deleted(property(pos)) );
    }

    inline bool is_longprobedistance(const int64_t positon)
    {
      return extract_probedistance(property(positon)) >= kLongProbedistanceThreshold;
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
        std::cout << "[" << i+Capacity*level << "] " << static_cast<int>(m_property_block_[i]) << " : " << m_key_block_[i] << " : " << static_cast<int>(m_value_block_[i]) << std::endl;
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
        if (is_valid(i)) {
          uint64_t d = extract_probedistance(property(i));
          output_file << Capacity << " " << d << std::endl;
        }
      }
      if (m_next_) {
        RHHStaticType *next_rhh = reinterpret_cast<RHHStaticType *>(m_next_);
        next_rhh->fprint_probedistance(output_file);
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

    static const ProbeDistanceType kLongProbedistanceThreshold  = 120;  /// must be less than max value of probedirance

    ///  ------------------------------------------------------ ///
    ///              Private Member Functions
    ///  ------------------------------------------------------ ///
    /// ------ Private member functions: algorithm core ----- ///
    inline HashType hash_key(KeyType& key)
    {
#if 1
      return static_cast<HashType>(key);
#else
      std::hash<KeyType> hash_fn;
      return hash_fn(key);
#endif
    }

    inline int64_t cal_desired_pos(HashType hash)
    {
      return hash & kMask;
    }

    inline ProbeDistanceType cal_probedistance(HashType hash, const int64_t slot_index)
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

    inline PropertyBlockType cal_property(HashType hash, const int64_t slot_index)
    {
      return cal_probedistance(hash, slot_index);
    }

    inline PropertyBlockType cal_property(const ProbeDistanceType probedistance)
    {
      return probedistance;
    }

    inline static bool is_deleted(const PropertyBlockType prop)
    {
      return (prop & kTombstoneMask) == kTombstoneMask;
    }

    inline static ProbeDistanceType extract_probedistance(PropertyBlockType prop)
    {
      return prop & kProbedistanceMask;
    }

    UpdateErrors insert_helper(KeyType &&key, ValueType &&val)
    {
      int64_t pos = cal_desired_pos(hash_key(key));
      ProbeDistanceType dist = 0;
      UpdateErrors err = kSucceed;

      while(true) {

        PropertyBlockType existing_elem_property = property(pos);

        if(existing_elem_property == kEmptyValue)
        {
          if (dist >= kLongProbedistanceThreshold) {
            assert(false);
            err = kLongProbedistance;
            dist = kLongProbedistanceThreshold;
          }
          construct(pos, dist, std::move(key), std::move(val));
          break;
        }

        /// If the existing elem has probed equal or less than new, then swap places with existing
        /// elem, and keep going to find another slot for that elem.
        if (extract_probedistance(existing_elem_property) <= dist)
        {
          if (dist >= kLongProbedistanceThreshold) {
            assert(false);
            err = kLongProbedistance;
            dist = kLongProbedistanceThreshold;
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

      return err;
    }

    inline void construct(const int64_t ix, const ProbeDistanceType dist, KeyType&& key, ValueType&& val)
    {
      m_property_block_[ix] = static_cast<PropertyBlockType>(dist);
      m_key_block_[ix] = std::move(key);
      m_value_block_[ix] = std::move(val);
    }

    /// ------ Private member functions: search ----- ///
    int64_t find_key(KeyType& key, bool is_check_recursively)
    {

      ProbeDistanceType dist = 0;
      const HashType hash = hash_key(key);
      int64_t pos = cal_desired_pos(hash);

      while(true) {

        ProbeDistanceType existing_elem_property = property(pos);
        if (existing_elem_property == kEmptyValue) { /// free space is found
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
        //std::cout << "* "; //D
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
