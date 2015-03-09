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
 
#ifndef HAVOQGT_MPI_ROBIN_HOOD_HASHING_HPP_INCLUDED
#define HAVOQGT_MPI_ROBIN_HOOD_HASHING_HPP_INCLUDED

#include <utility>
#include <boost/interprocess/allocators/allocator.hpp>


#define USE_SEPARATE_HASH_ARRAY 0

 namespace havoqgt {
namespace mpi {

namespace bip = boost::interprocess;

/// Note: Since we use 0 to indicate that the elem has never been used at all,
///       key 0 and 1 use same hash value (1).
///       This class is supporting duplicated-key, 
///       however not supporting duplicated-element.
template <typename Key, typename Value, typename SegementManager>


class robin_hood_hash {

  public:

  template<typename T>
  using SegmentAllocator = bip::allocator<T, SegementManager>;

  struct elem
  {
    Key key;
    Value value;
    elem(Key&& k, Value&& v) : key(std::move(k)), value(std::move(v)) {}
#if !USE_SEPARATE_HASH_ARRAY
    uint64_t hash;        
#endif
  };

  ///
  /// Iterator
  ///
  class elem_iterator : public std::iterator<std::input_iterator_tag, elem, ptrdiff_t, const elem* const, const elem&>
  {
  public:
    ///-----  Copy constructor -----///
    elem_iterator(const elem_iterator &obj)
    : hash_table_(obj.hash_table_)
    , key_(obj.key_) 
    {
      current_index_ = obj.current_index_;
      dist_ = obj.dist_;
    }


    ///-----  Operators -----///
    elem& operator*() const
    {
      return hash_table_->buffer_[current_index_];
    }

    elem_iterator& operator++() 
    {
      get_next();
      return *this;
    }

    elem_iterator operator++(int)
    {
      elem_iterator __tmp = *this;
      get_next();
      return __tmp;
    }

    elem *operator->() {
      return &hash_table_->buffer_[current_index_];
    }

    inline bool is_equal(const elem_iterator& x) const
    { return (x.current_index_ == current_index_); }

    ///  Return true if x and y are both end or not end, or x and y are the same.
    friend bool
    operator==(const elem_iterator& x, const elem_iterator& y)
    { return x.is_equal(y); }

    ///  Return false if x and y are both end or not end, or x and y are the same.
    friend bool
    operator!=(const elem_iterator& x, const elem_iterator& y)
    { return !x.is_equal(y); }


    ///-----  Public function -----///
    inline bool is_valid_index() const
    {
      return (current_index_ != kInvaridIndex);
    }

  private:
    friend class robin_hood_hash;

    ///----- Constructor -----///
    elem_iterator()
    : hash_table_(nullptr)
    , key_(0) 
    { 
      current_index_ = kInvaridIndex;
      dist_ = 0;
      has_key = false;
    }

    explicit elem_iterator(robin_hood_hash *_hash_table)
    : hash_table_(_hash_table)
    , key_(0)
    {
      dist_ = 0;
      current_index_ = -1;
      has_key = false;
      hash_table_->get_next_index(current_index_);
    }

    elem_iterator(robin_hood_hash *_hash_table, const Key& _key)
    : hash_table_(_hash_table)
    , key_(_key)
    {
      dist_ = 0;
      has_key = true;
      current_index_ = hash_table_->lookup_index_first(key_, dist_);
    }


    ///----- Private menber functions -----///
    void get_next()
    {
      if (has_key)
        hash_table_->get_next_index(key_, current_index_, dist_);
      else
        hash_table_->get_next_index(current_index_);
    }

    ///----- Private member variables -----///
    robin_hood_hash* hash_table_;
    const Key key_;
    int64_t dist_;
    int64_t current_index_;
    bool has_key;
  };

  friend class elem_iterator;


  /// ----------------------------------------------------------- ///
  ///                  robin_hood_hash class
  /// ----------------------------------------------------------- ///

  ///----- Constructor -----///
  explicit robin_hood_hash(SegmentAllocator<void>& _allocator)
  : allocator_(_allocator)
  , ptr_()
  {
    buffer_ = nullptr;
    num_elems_ = 0;
    capacity_ = kInitialCapacity;
    resize_threshold_ = (capacity_ * kLoadFactorPercent) / 100LL;
    mask_ = capacity_ - 1;
    alloc();
  }

  robin_hood_hash(SegmentAllocator<void>& _allocator, size_t initial_capasity) 
  : allocator_(_allocator)
  , ptr_()
  {
    buffer_ = nullptr;
    num_elems_ = 0;
    capacity_ = initial_capasity;
    resize_threshold_ = (capacity_ * kLoadFactorPercent) / 100LL;
    mask_ = capacity_ - 1;
    alloc();
  }

  ///-----  Copy constructor -----///
  robin_hood_hash(const robin_hood_hash &obj)
  : allocator_(obj.allocator_)
  , ptr_(obj.ptr_)
  {
    num_elems_ = obj.num_elems_;
    capacity_ = obj.capacity_;
    resize_threshold_ = obj.resize_threshold_;
    mask_ = obj.mask_;

    alloc();

    std::memcpy(buffer_, obj.buffer_, capacity_*sizeof(elem));
#if USE_SEPARATE_HASH_ARRAY
    std::memcpy(hashes_, obj.hashes_, capacity_*sizeof(uint64_t));
#endif
  }

  ///-----  Move constructor -----///
  robin_hood_hash(robin_hood_hash &&obj) 
  : allocator_(obj.allocator_)
  , ptr_(obj.ptr_)
  {
    buffer_ = obj.buffer_;
    obj.buffer_ = nullptr;
    num_elems_ = obj.num_elems_;
    obj.num_elems_ = 0;
    capacity_ = obj.capacity_;
    obj.capacity_ = 0;
    resize_threshold_ = obj.resize_threshold_;
    obj.resize_threshold_ = 0;
    mask_ = obj.mask_;
    obj.mask_ = 0;

#if USE_SEPARATE_HASH_ARRAY
    hashes_ = obj.hashes_;
    obj.hashes_ = nullptr;
#endif
  }

  ///-----  Copy assignment operator -----///
  robin_hood_hash& operator=( robin_hood_hash const& rhs )
  {
    return *this = robin_hood_hash(rhs);
  }

  ///-----  Move assignment operator -----///
  robin_hood_hash &operator=(robin_hood_hash&& obj)
  {
    ptr_ = obj.ptr_;
    buffer_ = obj.buffer_;
    obj.buffer_ = nullptr;
    num_elems_ = obj.num_elems_;
    obj.num_elems_ = 0;
    capacity_ = obj.capacity_;
    obj.capacity_ = 0;
    resize_threshold_ = obj.resize_threshold_;
    obj.resize_threshold_ = 0;
    mask_ = obj.mask_;
    obj.mask_ = 0;

#if USE_SEPARATE_HASH_ARRAY
    hashes_ = obj.hashes_;
    obj.hashes_ = nullptr;
#endif
    return *this;
  }

  void swap(robin_hood_hash& other) noexcept
  {
    using std::swap;
    swap(ptr_, other.ptr_);
    swap(num_elems_, other.num_elems_);
    swap(capacity_, other.capacity_);
    swap(resize_threshold_, other.resize_threshold_);
    swap(mask_, other.mask_);
    swap(buffer_, other.buffer_);
#if USE_SEPARATE_HASH_ARRAY
    swap(hashes_, other.hashes_);
#endif
  }
  void swap(robin_hood_hash& l, robin_hood_hash& r) noexcept {
    l.swap(r);
  }

  ///-----  deconstructors -----///
  ~robin_hood_hash()
  {
    for( int64_t i = 0; i < capacity_; ++i)
    {
      if (elem_hash(i) != 0) buffer_[i].~elem();
    }
    if (buffer_ != nullptr) free_buffer(ptr_, capacity_);
#if USE_SEPARATE_HASH_ARRAY
    delete [] hashes_;
#endif
  }

  ///----- Public memober functions ----- ///

  // Returns an iterator referring to the first element
  inline elem_iterator begin()
  {
    return elem_iterator(this);
  }

  inline elem_iterator end()
  {
    return elem_iterator();
  }


  inline void insert(Key key, Value val)
  {
    if (++num_elems_ >= resize_threshold_)
    {
      grow();
    }
    insert_helper(hash_key(key), std::move(key), std::move(val));
  }

  /// insert a element without duplicated key-value pair
  inline bool insert_unique(Key key, Value val)
  {
    if (has_data(key, val)) return false;

    if (++num_elems_ >= resize_threshold_)
    {
      grow();
    }
    insert_helper(hash_key(key), std::move(key), std::move(val));
    return true;
  }

  inline bool has_data(const Key& key, const Value& val) const
  {
    return (lookup_index(key, val) != kInvaridIndex);
  }

  // Searches for a first element with a key equivalent to 'key' and returns an iterator to it
  inline elem_iterator find(const Key& key)
  {
    return(elem_iterator(this, key));
  }

  inline void erase(elem_iterator& itr)
  {
    /// FIXME:
    //assert(itr.current_index_ >= 0 && itr.current_index_ < capacity_); //D

    erase_element(itr.current_index_);
#if USE_TOMBSTONE == 0
    if (!is_deleted(itr.current_index_)) {
      itr.current_index_ = (itr.current_index_ + capacity_ - 1) % capacity_;
      --(itr.dist_);
    }
#endif
    --num_elems_;
  }

  // inline size_t erase(const Key& key)
  // {
  //   int64_t dist = 0;
  //   int64_t pos = lookup_index_first(key, dist);
    
  //   size_t erased_count = 0;
  //   while(pos != kInvaridIndex) {
  //     get_next_index(key, pos, dist);
  //     erase_element(pos);
  //     ++erased_count;
  //   }
  //   num_elems_ -= erased_count;
  //   return erased_count;
  // }

  /// XXX: this function dose not supporting duplicated edges model
  inline bool erase(const Key& key, const Value& val)
  {
    const int64_t pos = lookup_index(key, val);
    if (pos == kInvaridIndex) return false;
    erase_element(pos);
    --num_elems_;
    return true;
  }

  inline size_t size() const
  {
    return num_elems_;
  }

  inline size_t allocated_size() const
  {
    return capacity_ * sizeof(elem);
  }

  // Count elements with a key equaivalent to 'key'
  inline size_t count(const Key& key)
  {
    int64_t dist = 0;
    int64_t pos = lookup_index_first(key, dist);
    
    size_t key_count = 0;
    while (pos != kInvaridIndex) { 
      get_next_index(key, pos, dist);
      ++key_count;
    }

    return key_count;
  }


  /// ----- Public member functions for debug ----- ///
  float average_probe_count() const
  {
    if (size() == 0) return 0;

    float probe_total = 0;
    for(int64_t i = 0; i < capacity_; ++i) {
      uint64_t hash = elem_hash(i);
      if (hash != 0 && !is_deleted(hash)) {
        probe_total += probe_distance(hash, i);
      }
    }
    return probe_total / size() + 1.0f;
  }

  void disp_elements()
  {
    const size_t length = capacity_;

    for (uint64_t i = 0; i < length; ++i) {
      if (elem_hash(i) == 0 || is_deleted(elem_hash(i))) continue;
      std::cout << buffer_[i].key << "\t" << buffer_[i].value << std::endl;
    }
  }

  void disp_hashes()
  {
    const size_t length = capacity_;

    for (uint64_t i = 0; i < length; ++i) {
      std::cout << buffer_[i].key << "\t" << elem_hash(i) << "\t" << desired_pos(elem_hash(i)) << "\t" << probe_distance(elem_hash(i), i) << std::endl;
    }
  }


  void dump_elements(const std::string& fname)
  {
    std::ofstream fout;
    fout.open(fname, std::ios::out | std::ios::app);
    const size_t length = capacity_;

    for (uint64_t i = 0; i < length; ++i) {
      if (elem_hash(i) == 0 || is_deleted(elem_hash(i))) continue;
      fout << buffer_[i].key << "\t" << buffer_[i].value << std::endl; 
    }
    fout.close();
  }

private:

#define USE_TOMBSTONE 1
#if USE_TOMBSTONE == 0
  #error Backshift model has a bug
#endif

  static const int64_t kInitialCapacity = 1; // must be 1 or more
  static const int64_t kLoadFactorPercent = 90LL;
  static const uint64_t kTombstoneMask = 0x8000000000000000ULL; // mask value to mark as deleted
  static const uint64_t kClearTombstoneMask = 0x7FFFFFFFFFFFFFFFULL;  // mask value to clear deleted flag
  static const int64_t  kCapacityGrowingFactor = 2LL;
  static const int64_t  kInvaridIndex = -1LL;

  /// ------ Private member functions: algorithm core ----- ///
  inline int64_t desired_pos(uint64_t hash) const
  {
    return hash & mask_;
  }

  inline int64_t probe_distance(const uint64_t hash, const uint64_t slot_index) const
  { 
    return (slot_index + capacity_ - desired_pos(hash)) & mask_;
  }

  inline uint64_t& elem_hash(int64_t ix)
  {
#if USE_SEPARATE_HASH_ARRAY
    return hashes_[ix];
#else
    return buffer_[ix].hash;
#endif
  }

  inline uint64_t elem_hash(int64_t ix) const
  {
    return const_cast<robin_hood_hash*>(this)->elem_hash(ix);
  }

  inline static uint64_t hash_key(const Key& key)
  {
    const std::hash<Key> hasher;
    uint64_t h = static_cast<uint64_t>(hasher(key));

#if USE_TOMBSTONE == 1
    // MSB is used to indicate a deleted elem, so
    // clear it
    h &= kClearTombstoneMask;
#endif

    // Ensure that we never return 0 as a hash,
    // since we use 0 to indicate that the elem has never
    // been used at all.
    h |= h==0ULL;
    return h;
  }

  inline void erase_element(const int64_t positon) {
    buffer_[positon].~elem();

// mark as deleted
#if USE_TOMBSTONE == 1
    elem_hash(positon) |= kTombstoneMask;
#else
    elem_hash(positon) = 0ULL;
    backshift_element(positon);
#endif

  }

  void backshift_element(const int64_t init_previous_pos) {
    int64_t previous_pos = init_previous_pos;
    int64_t swapped_pos = (init_previous_pos+1) & mask_;

    for(;;)
    {             
      if (elem_hash(swapped_pos) == 0) {// free space is found
        break ;
      } 
      if ( probe_distance(elem_hash(swapped_pos), swapped_pos) == 0) {
        break ;
      }

      std::swap(elem_hash(previous_pos), elem_hash(swapped_pos));
      std::swap(buffer_[previous_pos].key, buffer_[swapped_pos].key);
      std::swap(buffer_[previous_pos].value, buffer_[swapped_pos].value);
      previous_pos = swapped_pos;
      swapped_pos = (swapped_pos+1) & mask_;
    }
    
  }


  /// ----- Private funtions: memroy management ----- ///
  // alloc buffer according to currently set capacity
  void alloc()
  {
    /// This is temprorary codes to remove compile warning
    /// Todo: rewrite fllowing process
    // typedef bip::managed_mapped_file::segment_manager segment_manager_t;
    // segment_manager_t* segment_manager = allocator_.get_segment_manager();
    // bip::allocator<elem,segment_manager_t> alloc_inst(segment_manager);

    ptr_ = allocator_.allocate(capacity_ * sizeof(elem));
    buffer_ = reinterpret_cast<elem*>(ptr_.get());

#if USE_SEPARATE_HASH_ARRAY == 1
    hashes_ = new uint64_t[capacity_];
#endif

    // flag all elems as free
    for( uint64_t i = 0; i < capacity_; ++i)
    {
      elem_hash(i) = 0;
    }

    resize_threshold_ = (capacity_ * kLoadFactorPercent) / 100LL; 
    mask_ = capacity_ - 1;
  }

  inline void free_buffer(bip::offset_ptr<void>& _ptr, const size_t _capacity)
  {
    allocator_.deallocate(_ptr, _capacity * sizeof(elem));
  }

  void grow()
  {
    elem* old_elems = buffer_;
    bip::offset_ptr<void> old_ptr(ptr_);

    const uint64_t old_capacity = capacity_;
#if USE_SEPARATE_HASH_ARRAY
    auto old_hashes = hashes_;
#endif

    capacity_ *= kCapacityGrowingFactor;
    capacity_ |= capacity_== 0;
    alloc();

    // now copy over old elems
    for(uint64_t i = 0; i < old_capacity; ++i)
    {
      auto& e = old_elems[i];
#if USE_SEPARATE_HASH_ARRAY
      uint64_t hash = old_hashes[i];
#else
      uint64_t hash = e.hash;
#endif
      if (hash != 0 && !is_deleted(hash))
      {
        insert_helper(hash, std::move(e.key), std::move(e.value));
        e.~elem();
      }
    }
    free_buffer(old_ptr, old_capacity); 
#if USE_SEPARATE_HASH_ARRAY
    delete [] old_hashes;
#endif
  }


  /// ----- Private functions: search, delete, inset etc. ----- ///
  int64_t lookup_index_first(const Key& key, int64_t& dist) const
  {
    const uint64_t hash = hash_key(key);
    int64_t pos = desired_pos(hash);
    for(;;)
    {             
      if (elem_hash(pos) == 0) {// free space is found
        return kInvaridIndex;
      } else if (dist > probe_distance(elem_hash(pos), pos)) {
        return kInvaridIndex;
      } else if (elem_hash(pos) == hash && buffer_[pos].key == key) {
        return pos;
      }
      pos = (pos+1) & mask_;
      ++dist;
    }
    return kInvaridIndex;
  }

  int64_t lookup_index(const Key& key, const Value& val) const
  {
    const uint64_t hash = hash_key(key);
    int64_t dist = 0;
    int64_t pos = lookup_index_first(key, dist);
    if (pos == kInvaridIndex) return kInvaridIndex;
    for(;;) {
      if (elem_hash(pos) == 0) {// free space is found
        return kInvaridIndex;
      } else if (dist > probe_distance(elem_hash(pos), pos)) {
        return kInvaridIndex;
      } else if (elem_hash(pos) == hash && buffer_[pos].key == key && buffer_[pos].value == val) {
        return pos;
      }

      pos = (pos+1) & mask_;
      ++dist;
    }

    return kInvaridIndex;
  }


  /// ----- Private functions: search, delete, inset etc. ----- ///
  void get_next_index(int64_t& pos) const
  {
    ++pos;
    for(; pos < capacity_; ++pos)
    {
      if (elem_hash(pos) != 0 && !is_deleted(elem_hash(pos)))
        return;
    }
    pos = kInvaridIndex;
  }

  void get_next_index(const Key& key, int64_t& pos, int64_t& dist)
  {

    const uint64_t hash = hash_key(key);
    for(;;)
    { 
      pos = (pos+1) & mask_; // In order to prevent finding a start point,e increment positon at here.
      ++dist;
      if (elem_hash(pos) == 0) {// free space is found
        pos = kInvaridIndex;
        break;
      } else if (dist > probe_distance(elem_hash(pos), pos)) {
        pos = kInvaridIndex;
        break;
      } else if (elem_hash(pos) == hash && buffer_[pos].key == key) {
        break;
      }
    }
  }

  inline static bool is_deleted(const uint64_t hash)
  {
#if USE_TOMBSTONE == 1
    // MSB set indicates that this hash is a "tombstone"
    return (hash >> 63ULL) != 0;
#else
    return (hash == 0);
#endif
  }

  inline void construct(const int64_t ix, const uint64_t hash, Key&& key, Value&& val)
  {
    new (&buffer_[ix]) elem(std::move(key), std::move(val));  
    elem_hash(ix) = hash;
  }

  void insert_helper(uint64_t hash, Key&& key, Value&& val)
  {

    int64_t pos = desired_pos(hash);
    int64_t dist = 0;

    for(;;)
    {

      if(elem_hash(pos) == 0)
      {
        construct(pos, hash, std::move(key), std::move(val));
        return;
      }

      // If the existing elem has probed less than us, then swap places with existing
      // elem, and keep going to find another slot for that elem.
      int64_t existing_elem_probe_dist = probe_distance(elem_hash(pos), pos);
      if (existing_elem_probe_dist < dist)
      { 
        if(is_deleted(elem_hash(pos)))
        {
          construct(pos, hash, std::move(key), std::move(val));
          return;
        }
        std::swap(hash, elem_hash(pos));
        std::swap(key, buffer_[pos].key);
        std::swap(val, buffer_[pos].value);
        dist = existing_elem_probe_dist;        
      } else if (existing_elem_probe_dist == dist && is_deleted(elem_hash(pos))) {
        construct(pos, hash, std::move(key), std::move(val));
        return;        
      }

      pos = (pos+1) & mask_;
      ++dist;
    }
  }


  /// ---------- private menber variavles ---------- ///
  SegmentAllocator<void> allocator_;
  bip::offset_ptr<void> ptr_;

  elem* __restrict buffer_;

#if USE_SEPARATE_HASH_ARRAY
  uint64_t* __restrict hashes_;
#endif

  uint64_t num_elems_;
  uint64_t capacity_;
  uint64_t resize_threshold_;
  uint64_t mask_;

};

} // namespace mpi
} // namespace havoqgt

#endif
