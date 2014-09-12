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

#ifndef HAVOQGT_MPI_RHHAdjacencyMatrix_HPP_INCLUDED
#define HAVOQGT_MPI_RHHAdjacencyMatrix_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>

#define DEBUG(msg) do { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << msg << std::endl; } while (0)
#define DEBUG2(x) do  { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << #x << " =\t" << x << std::endl; } while (0)
#define DISP_VAR(x) do  { std::cout << #x << " =\t" << x << std::endl; } while (0)

#define USE_SEPARATE_HASH_ARRAY 0

namespace havoqgt {
namespace mpi {

namespace bip = boost::interprocess;

/// Note: Since we use 0 to indicate that the elem has never been used at all,
///       key 0 and 1 use same hash value (1).
///       This class is supporting duplicated-key,
///       however not supporting duplicated-element.
template <typename ElemType, typename SegmentAllocator>


class RHHAdjacencyMatrix {

  // template<typename T>grow
  // using SegmentAllocator = bip::allocator<T, SegementManager>;

 public:

  ///
  /// Iterator
  ///
  // class elem_iterator : public std::iterator<std::input_iterator_tag, elem, ptrdiff_t, const elem* const, const elem&>
  // {

  // };

  // friend class elem_iterator;




  public:

  /// ---------  Typedefs and Enums ------------ //


  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
  ///  ------------------------------------------------------ ///

  /// --- Constructor --- //
  explicit RHHAdjacencyMatrix(SegmentAllocator& allocator)
  {
    //DEBUG("RHHAdjacencyMatrix constructor");
    alloc(allocator, kInitialCapacity, true);
  }

  RHHAdjacencyMatrix(SegmentAllocator& allocator, const uint64_t initial_capasity, const bool is_allocate_valueblock)
  {
    //DEBUG("RHHAdjacencyMatrix constructor");
    alloc(allocator, initial_capasity, is_allocate_valueblock);
  }

  /// --- Copy constructor --- //

  /// --- Move constructor --- //
  RHHAdjacencyMatrix(RHHAdjacencyMatrix &&old_obj)
  {
    //DEBUG("RHHAdjacencyMatrix move-constructor");
    m_pos_head_ = old_obj.m_pos_head_;
    old_obj.m_pos_head_ = nullptr;
  }

  /// --- Destructor --- //
  ~RHHAdjacencyMatrix()
  {
    // DEBUG("RHHAdjacencyMatrix destructor");
    if (m_pos_head_ != nullptr) {
      assert(false);
      // DEBUG("RHHAdjacencyMatrix destructor : need free_buffer");
      //free_buffer(m_pos_head_, 0); /// XXX: memsize = 0
    }
  }

  /// ---  Move assignment operator --- //
  RHHAdjacencyMatrix &operator=(RHHAdjacencyMatrix&& old_obj)
  {
    //DEBUG("RHHAdjacencyMatrix move-assignment");
    m_pos_head_ = old_obj.m_pos_head_;
    old_obj.m_pos_head_ = nullptr;
    return *this;
  }


  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///

  bool insert_uniquely(SegmentAllocator& allocator, ElemType key, ElemType val)
  {
    // std::cout << key << "," << val << std::endl; //D
    // disp_elements(); //D

    ValueWrapperType value;
    value.value = val;
    const uint64_t mask = cal_mask();
    ProbeDistanceType dist = 0;
    const int64_t pos_first = lookup_first_position(key, &dist, mask);
    if (pos_first == kInvaridIndex) {
      /// 1. new vertex
      insert_directly_with_growing(allocator, std::move(key), std::move(value));
      return true;
    }

    if (is_adj_list(property(pos_first))) {
      /// 2. non-low-degree
      RHHAdjacencyMatrix& rhh_adj_list = m_pos_head_->pos_value_block[pos_first].adj_list;
      return rhh_adj_list.insert_key_unique(allocator, std::move(value.value));
    }

    /// Confirm unique insertion
    if (get_pos_equal_element_from_here(key, value.value, pos_first, dist, mask) != kInvaridIndex) {
    	return false;
    }

    int64_t moved_pos_list[kDirectInsertionThreshold];
    for (int64_t pos = pos_first, count_key = 0;; pos = (pos+1) & mask, ++dist) {
      PropertyBlockType existing_elem_property = property(pos);
      if (existing_elem_property == kClearedValue) {
        /// 3. low-degree
        construct(pos, hash_key(key), std::move(key), std::move(value), false, mask);
        return true;
      } else if (dist > extract_probedistance(existing_elem_property)) {
        /// 3. low-degree
        /// TODO: we can use value of 'pos' and 'dist' in order to reduce serching steps
        insert_directly_with_growing(allocator, std::move(key), std::move(value));
        return true;
      }else if (!is_deleted(existing_elem_property) && m_pos_head_->pos_key_block[pos] == key) {
        moved_pos_list[count_key] = pos;
        if (++count_key >= kDirectInsertionThreshold) break;
      }
    }

    /// 4. convert data strucure of value from direct-insertion model to adj-list-insertion model
    //DEBUG("Convert start, current initial alloc size = 2");
    ValueWrapperType value_wrapper;
    new(&value_wrapper.adj_list) RHHAdjacencyMatrix(allocator, kInitialCapacityAdjlist, false);
    value_wrapper.adj_list.insert_directly_with_growing(allocator, std::move(value.value));
    for (uint64_t k = 0; k < kDirectInsertionThreshold; ++k) {
      value_wrapper.adj_list.insert_directly_with_growing(allocator, std::move(m_pos_head_->pos_value_block[moved_pos_list[k]].value));
      delete_elem(allocator, moved_pos_list[k]);
    }
    /// insert adjacency-list
    //DEBUG("Move adj_list");
    construct(pos_first, hash_key(key), std::move(key), std::move(value_wrapper), true, mask);
    //DEBUG("Convert done");

    return true;
  }

  bool erase(SegmentAllocator& allocator, ElemType key, ElemType value)
  {
    const uint64_t mask = cal_mask();
  	ProbeDistanceType dist = 0;
    const int64_t pos_first = lookup_first_position(key, &dist, mask);
    if (pos_first == kInvaridIndex) {
      /// we don't have the key
      return false;
    }

    /// target value is stored as adj-list
    if (is_adj_list(property(pos_first))) {
      RHHAdjacencyMatrix& rhh_adj_list = m_pos_head_->pos_value_block[pos_first].adj_list;
      if (!rhh_adj_list.erase_key_only(value)) return false;
      if (rhh_adj_list.m_pos_head_->num_elems == kDirectInsertionThreshold) {
      	ElemType moved_value_list[kDirectInsertionThreshold];
      	rhh_adj_list.get_all_elem_key_list(moved_value_list);
      	delete_elem(allocator, pos_first);
        --m_pos_head_->num_elems;
      	for (int k=0; k<kDirectInsertionThreshold; ++k) {
          ElemType ky = key;
          ValueWrapperType val;
          val.value = moved_value_list[k];
      	  insert_directly_with_growing(allocator, std::move(ky), std::move(val));
      	}
      } else if (rhh_adj_list.m_pos_head_->num_elems < kDirectInsertionThreshold) {
        assert(false);
      }
      return true;
    }

    /// target value is stored as simple value
    int64_t pos = get_pos_equal_element_from_here(key, value, pos_first, dist, mask);
    if (pos != kInvaridIndex) {
    	delete_elem(allocator, pos);
    	--m_pos_head_->num_elems;
    	return true;
    }

    return false;
  }

  void free(SegmentAllocator& allocator)
  {
    free_buffer(allocator, reinterpret_cast<void *>(m_pos_head_), 0);
  }

  void dump_probedistance (const std::string& fname)
  {
    std::ofstream fout;
    fout.open(fname, std::ios::out | std::ios::app);
    for (uint64_t i = 0; i < m_pos_head_->capacity; ++i) {
      PropertyBlockType prop = m_pos_head_->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
      if (is_adj_list(prop)) {
        RHHAdjacencyMatrix& rhh_adj_list = m_pos_head_->pos_value_block[i].adj_list;
        rhh_adj_list.dump_probedistance(fname);
      } else {
        fout << m_pos_head_->capacity << "\t" << extract_probedistance(prop) << std::endl;
      }
    }
    fout.close();
  }

  void disp_status()
  {
    DISP_VAR(kDirectInsertionThreshold);
    DISP_VAR(kInitialCapacity);
    DISP_VAR(kInitialCapacityAdjlist);
    DISP_VAR(kResizeFactor);
    DISP_VAR(kCapacityGrowingFactor);
    DISP_VAR(m_pos_head_->num_elems);
    DISP_VAR(m_pos_head_->capacity);
  }

  void dump_elements(const std::string& fname)
  {
    std::ofstream fout;
    fout.open(fname, std::ios::out | std::ios::app);

    for (uint64_t i = 0; i < m_pos_head_->capacity; ++i) {
      PropertyBlockType prop = m_pos_head_->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
      if (is_adj_list(prop)) {
        RHHAdjacencyMatrix& rhh_adj_list = m_pos_head_->pos_value_block[i].adj_list;
        rhh_adj_list.dump_keys_with_prefix(fout, i);
      } else {
        fout << m_pos_head_->pos_key_block[i] << "\t" << m_pos_head_->pos_value_block[i].value << std::endl;
      }
    }
    fout.close();
  }

 void dump_keys_with_prefix(std::ofstream& fout, uint64_t prefix)
 {
    for (uint64_t i = 0; i < m_pos_head_->capacity; ++i) {
      PropertyBlockType prop = m_pos_head_->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
      fout << prefix << "\t" << m_pos_head_->pos_key_block[i] << std::endl;
    }
 }

  void disp_elements(void)
  {

    for (uint64_t i = 0; i < m_pos_head_->capacity; ++i) {
      PropertyBlockType prop = m_pos_head_->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
      if (is_adj_list(prop)) {
        RHHAdjacencyMatrix& rhh_adj_list = m_pos_head_->pos_value_block[i].adj_list;
        rhh_adj_list.disp_keys_with_prefix(i);
      } else {
        std::cout << m_pos_head_->pos_key_block[i] << "\t" << m_pos_head_->pos_value_block[i].value << std::endl;
      }
    }
    std::cout << "\n\n"; //D
  }

 void disp_keys_with_prefix(uint64_t prefix)
 {
    for (uint64_t i = 0; i < m_pos_head_->capacity; ++i) {
      PropertyBlockType prop = m_pos_head_->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
      std::cout << prefix << "\t" << m_pos_head_->pos_key_block[i] << std::endl;
    }
 }


 private:
  /// ---------  Typedefs and Enums ------------ ///
  typedef uint64_t PropertyBlockType;
  typedef int64_t ProbeDistanceType;
  typedef uint64_t HashElemType;

  union ValueWrapperType {
    ElemType value;
    RHHAdjacencyMatrix adj_list;
    ValueWrapperType() {
      //DEBUG("union constructor");
    }
    ValueWrapperType(ValueWrapperType &&old_obj)
    {
      //DEBUG("union move-constructor");
      value = old_obj.value;
      old_obj.value = NULL;
    }
    ~ValueWrapperType()
    {
      //DEBUG("union destructor");
    }
    ValueWrapperType &operator=(ValueWrapperType &&old_obj)
    {
     // DEBUG("union move-assigment");
      value = old_obj.value;
      old_obj.value = NULL;
      return *this;
    }
  };

  struct HeadBlock {
    PropertyBlockType *pos_property_block;
    ElemType *pos_key_block;
    ValueWrapperType *pos_value_block;

    uint64_t num_elems;
    uint64_t capacity;
  };



  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///

#define USE_TOMBSTONE 1
#if USE_TOMBSTONE == 0
  #error Backshift model has a bug
#endif

  static const uint64_t kDirectInsertionThreshold = 1ULL; /// NOTE: We assume that this value is small (may be less than 10)
  static const uint64_t kInitialCapacity = 2ULL; /// must be 1 or more
  static const uint64_t kInitialCapacityAdjlist = 4ULL; /// must be 1 or more
  static constexpr double kResizeFactor = 0.9;
  static const uint64_t kCapacityGrowingFactor = 2ULL;

  static const PropertyBlockType kTombstoneMask     = 0x8000000000000000LL; /// mask value to mark as deleted
  static const PropertyBlockType kAdjacencylistMask = 0x4000000000000000LL; /// mask value represents whether value is adjacency-list(RHH) or just a value
  static const PropertyBlockType kProbedistanceMask = 0x0000000000000FFFLL; /// mask value to extract probe distance
  static const PropertyBlockType kCapacityMask      = 0x3FFFFFFFFFFFF000LL; /// mask value to extract probe distance
  static const PropertyBlockType kClearedValue      = 0x3FFFFFFFFFFFFFFFLL; /// value repsents cleared space
  static const int64_t kInvaridIndex = -1LL;


  /// ------ Private member functions: algorithm core ----- ///
  /// NOTE: we assume that key can use static_cast to HashElemType
  inline HashElemType hash_key(ElemType& key)
  {
    return static_cast<HashElemType>(key);
  }

  inline uint64_t cal_mask() const
  {
    return m_pos_head_->capacity - 1ULL;
  }

  inline int64_t cal_desired_pos(HashElemType hash, const uint64_t mask)
  {
    return hash & mask;
  }

  /// XXX: capacity = mask + 1LL
  inline ProbeDistanceType cal_probedistance(HashElemType hash, const int64_t slot_index, const uint64_t mask)
  {
    return ((slot_index + (mask + 1LL) - cal_desired_pos(hash, mask)) & mask);
  }

  /// XXX: tombstone only
  inline PropertyBlockType& property(const int64_t ix)
  {
    return m_pos_head_->pos_property_block[ix];
  }

  inline PropertyBlockType property(const int64_t ix) const
  {
    return const_cast<RHHAdjacencyMatrix*>(this)->property(ix);
  }

  inline PropertyBlockType cal_property(HashElemType hash, const int64_t slot_index, const bool is_adj_list, const uint64_t mask)
  {
    return cal_probedistance(hash, slot_index, mask) | (static_cast<uint64_t>(is_adj_list) * kAdjacencylistMask);
  }

  inline PropertyBlockType cal_property(const ProbeDistanceType probedistance, const bool is_adj_list)
  {
    return probedistance | (static_cast<uint64_t>(is_adj_list) * kAdjacencylistMask);
  }

  inline bool extract_adj_list_flag(const PropertyBlockType prop)
  {
    return (prop & kAdjacencylistMask) == kAdjacencylistMask;
  }

  /// XXX: tombstone only
  inline void delete_elem(SegmentAllocator& allocator, const int64_t positon)
  {
  	if (is_adj_list(property(positon))) {
      m_pos_head_->pos_value_block[positon].adj_list.free(allocator);
    	m_pos_head_->pos_value_block[positon].adj_list.~RHHAdjacencyMatrix();
    }
    property(positon) |= kTombstoneMask;
  }

  /// XXX: tombstone only
  inline void delete_key(const int64_t positon)
  {
    property(positon) |= kTombstoneMask;
  }

  /// XXX: tombstone only
  inline static bool is_deleted(const PropertyBlockType prop)
  {
    return (prop & kTombstoneMask) == kTombstoneMask;
  }

  /// XXX: tombstone only
  inline static bool is_adj_list(const PropertyBlockType prop)
  {
    return (prop & kAdjacencylistMask) == kAdjacencylistMask;
  }

  inline static ProbeDistanceType extract_probedistance(PropertyBlockType prop)
  {
  	return prop & kProbedistanceMask;
  }

  inline static ProbeDistanceType extract_capacity(PropertyBlockType prop)
  {
    return prop & kCapacityMask;
  }


  /// ------ Private member functions: memory management ----- ///

  void grow(SegmentAllocator& allocator)
  {
    HeadBlock* old_pos_head = m_pos_head_;

    const uint64_t new_capacity = m_pos_head_->capacity * kCapacityGrowingFactor;
    alloc(allocator, new_capacity, true);
    m_pos_head_->num_elems = old_pos_head->num_elems;

    /// now copy over old elems
    const uint64_t mask = cal_mask();
    uint64_t old_capacity = old_pos_head->capacity;
    for(uint64_t i = 0; i < old_capacity; ++i) {
      PropertyBlockType prop = old_pos_head->pos_property_block[i];
      if (is_deleted(prop) || prop == kClearedValue) continue;
	    insert_helper(std::move(old_pos_head->pos_key_block[i]), std::move(old_pos_head->pos_value_block[i]), is_adj_list(prop), mask);
    }

    free_buffer(allocator, reinterpret_cast<void*>(old_pos_head), 0);
  }

  void grow_without_valueblock(SegmentAllocator& allocator)
  {
    HeadBlock* old_pos_head = m_pos_head_;

    const uint64_t new_capacity = m_pos_head_->capacity * kCapacityGrowingFactor;
    alloc(allocator, new_capacity, false);
    m_pos_head_->num_elems = old_pos_head->num_elems;

    /// now copy over old elems
    const uint64_t mask = cal_mask();
    uint64_t old_capacity = old_pos_head->capacity;
    for(uint64_t i = 0; i < old_capacity; ++i) {
      ProbeDistanceType dist = old_pos_head->pos_property_block[i];
      if (is_deleted(dist) || dist == kClearedValue)
		    continue;
  	  insert_helper(std::move(old_pos_head->pos_key_block[i]), mask);
    }

     free_buffer(allocator, reinterpret_cast<void*>(old_pos_head), 0);
  }

  void alloc(SegmentAllocator& allocator, const uint64_t capacity, const bool is_allocate_valueblock)
  {
    uint64_t mem_size = sizeof(HeadBlock);
    mem_size += capacity*sizeof(PropertyBlockType);
	  mem_size += capacity*sizeof(ElemType);
	  mem_size += capacity*sizeof(ValueWrapperType)*is_allocate_valueblock;

    bip::offset_ptr<void> ptr = allocator.allocate(mem_size);
    m_pos_head_ = reinterpret_cast<HeadBlock*>(ptr.get());
    m_pos_head_->pos_property_block   = reinterpret_cast<PropertyBlockType*>(&m_pos_head_[1]);
    m_pos_head_->pos_key_block        = reinterpret_cast<ElemType*>(&m_pos_head_->pos_property_block[capacity]);
  	if (is_allocate_valueblock)
  		 	m_pos_head_->pos_value_block  = reinterpret_cast<ValueWrapperType*>(&m_pos_head_->pos_key_block[capacity]);
  	else
  			m_pos_head_->pos_value_block  = nullptr;
  	m_pos_head_->capacity             = capacity;

    for (uint64_t k = 0; k < capacity; ++k) {
      m_pos_head_->pos_property_block[k] = kClearedValue;
    }
  }

  inline void free_buffer(SegmentAllocator& allocator, void* ptr, const size_t memsize)
  {
  	//allocator.mp_mngr->deallocate(ptr);
    allocator.deallocate(bip::offset_ptr<void>(ptr), memsize);
  }

  /// ------ Private member functions: insertion ----- ///

  inline bool insert_key_unique(SegmentAllocator& allocator, ElemType&& key)
  {
 	  ProbeDistanceType dist = 0;
    const uint64_t mask = cal_mask(); // sinse this function is called by parent RHH, we should calvurate mask at here
    const int64_t pos_first = lookup_first_position(key, &dist, mask);
    if (pos_first != kInvaridIndex) {
      return false;
    }
    if (get_pos_equal_element_from_here(key, pos_first, dist, mask) != kInvaridIndex) {
      return false;
    }
    insert_directly_with_growing(allocator, std::move(key));
    return true;
  }

  inline void insert_directly_with_growing(SegmentAllocator& allocator, ElemType&& key, ValueWrapperType&& value)
  {
  	if (++m_pos_head_->num_elems >= m_pos_head_->capacity*kResizeFactor) {
  		grow(allocator);
  	}
  	insert_helper(std::move(key), std::move(value), false, cal_mask());
  }

  inline void insert_directly_with_growing(SegmentAllocator& allocator, ElemType&& key)
  {
  	if (++m_pos_head_->num_elems >= m_pos_head_->capacity*kResizeFactor) {
  		grow_without_valueblock(allocator);
  	}
  	insert_helper(std::move(key), cal_mask());
  }

  void insert_helper(ElemType &&key, ValueWrapperType &&value, bool is_adj_list, const uint64_t mask)
  {
    int64_t pos = cal_desired_pos(hash_key(key), mask);
    ProbeDistanceType dist = 0;

    while(true) {
      PropertyBlockType existing_elem_property = property(pos);

      if(existing_elem_property == kClearedValue)
      {
        construct(pos, hash_key(key), std::move(key), std::move(value), is_adj_list, mask);
        return;
      }

      /// If the existing elem has probed less than or "equal to" us, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if(is_deleted(existing_elem_property))
        {
          construct(pos, hash_key(key), std::move(key), std::move(value), is_adj_list, mask);
          return;
        }
        m_pos_head_->pos_property_block[pos] = cal_property(dist, is_adj_list);
        std::swap(key, m_pos_head_->pos_key_block[pos]);
        std::swap(value, m_pos_head_->pos_value_block[pos]);
        dist = extract_probedistance(existing_elem_property);
        is_adj_list = extract_adj_list_flag(existing_elem_property);
      }

      pos = (pos+1) & mask;
      ++dist;
    }
  }

  void insert_helper(ElemType &&key, const uint64_t mask)
  {
    int64_t pos = cal_desired_pos(hash_key(key), mask);
    ProbeDistanceType dist = 0;

    while(true) {
      PropertyBlockType existing_elem_property = property(pos);

      if(existing_elem_property == kClearedValue)
      {
        construct(pos, hash_key(key), std::move(key), mask);
        return;
      }

      /// If the existing elem has probed less than or "equal to" us, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (extract_probedistance(existing_elem_property) <= dist)
      {
        if(is_deleted(existing_elem_property))
        {
          construct(pos, hash_key(key), std::move(key), mask);
          return;
        }
        m_pos_head_->pos_property_block[pos] = dist;
        std::swap(key, m_pos_head_->pos_key_block[pos]);
        dist = extract_probedistance(existing_elem_property);
      }

      pos = (pos+1) & mask;
      ++dist;
    }
  }

  inline void construct(const int64_t ix, const HashElemType hash, ElemType&& key, ValueWrapperType&& val, const bool is_adj_list, const uint64_t mask)
  {
    m_pos_head_->pos_property_block[ix] = cal_property(hash, ix, is_adj_list, mask);
    m_pos_head_->pos_key_block[ix] = std::move(key);
    m_pos_head_->pos_value_block[ix] = std::move(val);
  }

  inline void construct(const int64_t ix, const HashElemType hash, ElemType&& key, const uint64_t mask)
  {
    m_pos_head_->pos_property_block[ix] = cal_probedistance(hash, ix, mask);
    m_pos_head_->pos_key_block[ix] = std::move(key);
  }


  /// ------ Private member functions: search ----- ///

  int64_t lookup_first_position(ElemType& key, ProbeDistanceType* dist, const uint64_t mask)
  {
    const HashElemType hash = hash_key(key);
    int64_t pos = cal_desired_pos(hash, mask);

    while(true) {
      ProbeDistanceType existing_elem_property = property(pos);
      if (existing_elem_property == kClearedValue) { /// free space is found
        return kInvaridIndex;
      } else if (*dist > extract_probedistance(existing_elem_property)) {
        return kInvaridIndex;
      } else if (!is_deleted(existing_elem_property) && m_pos_head_->pos_key_block[pos] == key) {
        /// found !
        return pos;
      }
      pos = (pos+1) & mask;
      *dist = *dist + 1;
    }
  }


  /// ------ Private member functions: erase ----- ///

  bool erase_key_only(ElemType &key)
  {
    ProbeDistanceType dist = 0;
    const int64_t pos = lookup_first_position(key, &dist, cal_mask());
    if (pos == kInvaridIndex) {
      return false;
    }

    delete_key(pos);
    --m_pos_head_->num_elems;
    return true;
  }


  /// ------ Private member functions: utility ----- ///
    /// NOTE: this function require buffer, i.e. "pos_list", in order to store keys
  inline void get_all_elem_key_list(ElemType pos_list[])
  {
    uint64_t capacity = m_pos_head_->capacity;
    uint64_t num_elems = m_pos_head_->num_elems;
    for (int64_t pos=0, cout=0; pos < capacity && cout < num_elems; ++pos) {
      PropertyBlockType prop = property(pos);
      if (is_deleted(prop) || prop == kClearedValue) continue;
      pos_list[cout] = m_pos_head_->pos_key_block[pos];
      ++cout;
    }
  }

  inline int64_t get_pos_equal_element_from_here(const ElemType& key, const ElemType& value, int64_t pos, ProbeDistanceType dist, const uint64_t mask) const
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

  inline int64_t get_pos_equal_element_from_here(const ElemType& key, int64_t pos, ProbeDistanceType dist, const uint64_t mask) const
  {
    /// Check whether we have a same key element.
    PropertyBlockType existing_elem_property = property(pos);
    while(true) {
      if (existing_elem_property==kClearedValue || dist > extract_probedistance(existing_elem_property))
        return kInvaridIndex;
      if (is_equal_element(key, pos, existing_elem_property))
        return pos;
      pos=(pos+1)&mask;
      ++dist;
      existing_elem_property = property(pos);
    }
  }

  inline bool is_equal_element(const ElemType& key, const ElemType& value, const int64_t pos, const PropertyBlockType prop) const
  {
  	if (is_deleted(prop) || m_pos_head_->pos_key_block[pos] != key || m_pos_head_->pos_value_block[pos].value != value)
  		return false;
  	else
  		return true;
  }

  inline bool is_equal_element(const ElemType& key, const int64_t pos, const PropertyBlockType prop) const
  {
    if (is_deleted(prop) || m_pos_head_->pos_key_block[pos] != key)
      return false;
    else
      return true;
  }


  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
  HeadBlock *m_pos_head_;

};

} /// namespace mpi
} /// namespace havoqgt

#endif
