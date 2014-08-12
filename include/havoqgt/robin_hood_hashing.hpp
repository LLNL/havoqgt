
#ifndef HAVOQGT_MPI_ROBIN_HOOD_HASHING_HPP_INCLUDED
#define HAVOQGT_MPI_ROBIN_HOOD_HASHING_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>


#define USE_SEPARATE_HASH_ARRAY 1

namespace havoqgt {
namespace mpi {

namespace bip = boost::interprocess;

template <typename Key, typename Value, typename SegementManager>
class hash_table {

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


  ///-----  Constructors -----///
  explicit hash_table(SegmentAllocator<void>& seg_allocator) :
  buffer_(nullptr),
  allocator_(seg_allocator)
  {
    init();
    alloc();
  }

  ///-----  deconstructors -----///
  ~hash_table()
  {   
    for( int64_t i = 0; i < capacity_; ++i)
    {
      if (elem_hash(i) != 0)
      {
        buffer_[i].~elem();
      }
    }
    free_buffer(buffer_, capacity_);
#if USE_SEPARATE_HASH_ARRAY
    delete [] hashes_;
#endif
  }


  ///----- Public memober functions ----- ///
  void resize(size_t new_size)
  {
    capacity_ = new_size;
    alloc();
  }

  inline void insert(Key key, Value val)
  {
    // std::cout << "Key = " << key << "\tval = " << val << std::endl;
    if (++num_elems_ >= resize_threshold_)
    {
      grow();
    }
    insert_helper(hash_key(key), std::move(key), std::move(val));   
  }

  inline void insert_unique(Key key, Value val)
  {
    if (is_duplicated(key, val)) return;
    insert(key, val);
  }

  inline bool is_duplicated(const Key& key, const Value& value)
  {
    const int64_t ix = lookup_index(key);
    return (ix != -1) && (buffer_[ix].value == value);
  }

  inline Value* find(const Key& key)
  {
    const int64_t ix = lookup_index(key);
    return ix != -1 ? &buffer_[ix].value : nullptr;
  }

  inline const Value* find(const Key& key) const
  {
    return const_cast<hash_table*>(this)->lookup(key);
  }

  bool erase(const Key& key)
  {
    const int64_t ix = lookup_index(key);

    if (ix == -1) return false;

    buffer_[ix].~elem();
    elem_hash(ix) |= 0x80000000; // mark as deleted
    --num_elems_;
    return true;
  }

  inline size_t size() const
  {
    return num_elems_;
  }


  /// ----- Public member functions: debug ----- ///
  float average_probe_count() const
  {
    if (size() == 0) {
      return 0;
    }
    float probe_total = 0;
    for(int64_t i = 0; i < capacity_; ++i)
    {
      uint64_t hash = elem_hash(i);
      if (hash != 0 && !is_deleted(hash))
      {
        probe_total += probe_distance(hash, i);
      }
    }
    return probe_total / size() + 1.0f;
  }

  void disp_elements()
  {
    const size_t length = capacity_;

    for (uint64_t i = 0; i < length; i++) {
      if (elem_hash(i) == 0 || is_deleted(elem_hash(i))) continue;
      std::cout << buffer_[i].key << "\t" << buffer_[i].value << std::endl;
    }
  }

  void dump_elements(const std::string& fname)
  {
    std::ofstream fout;
    fout.open(fname);
    const size_t length = capacity_;

    for (uint64_t i = 0; i < length; i++) {
      if (elem_hash(i) == 0 || is_deleted(elem_hash(i))) continue;
      fout << buffer_[i].key << "\t" << buffer_[i].value << std::endl; 
    }
    fout.close();
  }



private:
  static const int64_t INITIAL_SIZE = 1; // must be 1 or more
  static const int64_t LOAD_FACTOR_PERCENT = 90LL;

  void init()
  {
    num_elems_ = 0;
    capacity_ = INITIAL_SIZE;
    resize_threshold_ = (capacity_ * LOAD_FACTOR_PERCENT) / 100LL;
    mask_ = capacity_ - 1;
  }


  /// ------ Private member functions: algorithm core ----- ///
  inline int64_t desired_pos(uint64_t hash) const
  {
    return hash & mask_;
  }

  inline int64_t probe_distance(uint64_t hash, uint64_t slot_index) const
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
    return const_cast<hash_table*>(this)->elem_hash(ix);
  }

  inline static uint64_t hash_key(const Key& key)
  {
    const std::hash<Key> hasher;
    auto h = static_cast<uint64_t>(hasher(key));

    // MSB is used to indicate a deleted elem, so
    // clear it
    h &= 0x7fffffff;

    // Ensure that we never return 0 as a hash,
    // since we use 0 to indicate that the elem has never
    // been used at all.
    h |= h==0LL;
    return h; 
  }


  /// ----- Private funtions: memroy management ----- ///
  // alloc buffer according to currently set capacity
  void alloc()
  { 
    /// This is temprorary codes to remove compile warning
    /// Todo: rewrite fllowing process
    typedef bip::managed_mapped_file mapped_t;
    typedef mapped_t::segment_manager segment_manager_t;
    segment_manager_t* segment_manager = allocator_.get_segment_manager();
    bip::allocator<elem,segment_manager_t> alloc_inst(segment_manager);


    bip::offset_ptr<elem> ptr = alloc_inst.allocate(capacity_ * sizeof(elem));
    buffer_ = reinterpret_cast<elem*>(ptr.get());
#if USE_SEPARATE_HASH_ARRAY
    // hashes_ = reinterpret_cast<uint64_t*>( allocator_.allocate(capacity_ * sizeof(uint64_t)) );
    hashes_ = new uint64_t[capacity_];
#endif

    // flag all elems as free
    for( int64_t i = 0; i < capacity_; ++i)
    {
      elem_hash(i) = 0;
    }

    resize_threshold_ = (capacity_ * LOAD_FACTOR_PERCENT) / 100LL; 
    mask_ = capacity_ - 1;
  }

  void free_buffer(elem* buf, size_t capacity)
  {
    allocator_.deallocate(buf, capacity * sizeof(elem));
  }

  void grow()
  {
    elem* old_elems = buffer_;
    int64_t old_capacity = capacity_;
#if USE_SEPARATE_HASH_ARRAY
    auto old_hashes = hashes_;
#endif
    capacity_ *= 2LL;
    capacity_ |= capacity_==0;
    alloc();
    // now copy over old elems
    for(int64_t i = 0; i < old_capacity; ++i)
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

    free_buffer(old_elems, old_capacity);   
#if USE_SEPARATE_HASH_ARRAY
    delete [] old_hashes;
#endif
  }


  /// ----- Private functions: search, delete, inset etc. ----- ///
  int64_t lookup_index(const Key& key) const
  {
    const uint64_t hash = hash_key(key);
    int64_t pos = desired_pos(hash);
    int64_t dist = 0;
    for(;;)
    {             
      if (elem_hash(pos) == 0) 
        return -1;
      else if (dist > probe_distance(elem_hash(pos), pos)) 
        return -1;
      else if (elem_hash(pos) == hash && buffer_[pos].key == key) 
        return pos;       

      pos = (pos+1) & mask_;
      ++dist;
    }
  }

  inline static bool is_deleted(uint64_t hash)
  {
    // MSB set indicates that this hash is a "tombstone"
    return (hash >> 63LL) != 0;
  }

  inline void construct(int64_t ix, uint64_t hash, Key&& key, Value&& val)
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
      }

      pos = (pos+1) & mask_;
      ++dist;     
    }
  }


  /// ---------- private menber variavles ---------- ///
  SegmentAllocator<void> &allocator_;
  elem* __restrict buffer_;

#if USE_SEPARATE_HASH_ARRAY
  uint64_t* __restrict hashes_;
#endif  
  int64_t num_elems_;
  int64_t capacity_;
  int64_t resize_threshold_;
  uint64_t mask_;

};

} // namespace mpi
} // namespace havoqgt

#endif