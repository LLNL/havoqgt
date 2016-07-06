/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHHDA_RHHALLOCATOR_HOLDER_HPP_INCLUDED
#define RHHDA_RHHALLOCATOR_HOLDER_HPP_INCLUDED

#include <cstdio>
#include <memory>
#include <type_traits> // static_assert

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <havoqgt/graphstore/rhh/rhh_defs.hpp>

namespace graphstore {
namespace rhh {

///
/// \brief The static_node_allocator_wrp class
///   Wrapper class for static_node_allocator class to management the classes using an array data structure
class static_node_allocator_wrp
{
public:
  virtual void* allocate() =0;
  virtual void deallocate(void *) =0;
  virtual void deallocate_free_blocks() =0;
};

////
/// \brief The static_node_allocator class
///   actual class which holds a node allocater
template<typename segment_manager_type, size_t allocate_size>
class static_node_allocator : public static_node_allocator_wrp
{

public:

  explicit static_node_allocator(segment_manager_type* segment_manager)
    : m_node_allocator(segment_manager)
  {
    static_assert(kNodeAllocatorChunkSize >= allocate_size, "sizeof(rhhda_static) is larger than kNodeAllocatorChunkSize");
  }

  void* allocate()
  {
    return reinterpret_cast<void*>(m_node_allocator.allocate(1).get());
  }

  segment_manager_type* get_segment_manager() const
  {
    return m_node_allocator.get_segment_manager();
  }

  void deallocate(void *ptr)
  {
    m_node_allocator.deallocate(boost::interprocess::offset_ptr<dummy_node_type>(reinterpret_cast<dummy_node_type*>(ptr)),
                                1);
  }

  ///
  /// \brief deallocate_free_blocks
  ///   Actually deallocate node_allocater's free chunks.
  ///   The cost of this operation is O(n^2)
  void deallocate_free_blocks()
  {
    m_node_allocator.deallocate_free_blocks();
  }

private:
  struct dummy_node_type {
    char dummy[allocate_size];
  };
  enum NodesPerChunk : size_t {
    kNodesPerChunk  = kNodeAllocatorChunkSize / allocate_size,
  };
  using node_allocator_type  = boost::interprocess::node_allocator<dummy_node_type,  segment_manager_type, kNodesPerChunk>;

  node_allocator_type m_node_allocator;
};


///
/// \brief The static_node_allocators_holder class
///   array of node allocators (static_node_allocator_wrp class).
///   depends on the required size, uses proper sized node allocator
template<typename segment_manager_type, size_t element_size, size_t extra_size>
class static_node_allocators_holder
{
public:
  enum : size_t{
    kMaxCapacity = (1ULL << 18)
  };
  ///
  /// \brief rhhda_allocator_holder
  /// \param segment_manager
  ///
  explicit static_node_allocators_holder(segment_manager_type* segment_manager) :
    m_static_node_allocators{{
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 0] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 1] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 2] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 3] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 4] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 5] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 6] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 7] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 8] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[ 9] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[10] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[11] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[12] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[13] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[14] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[15] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[16] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[17] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[18] + extra_size>(segment_manager),
      new static_node_allocator<segment_manager_type, element_size * num_elements[19] + extra_size>(segment_manager)
    }}
  {}

  ~static_node_allocators_holder()
  {
    for (int i = 0; i < 20; ++i) {
      delete m_static_node_allocators[i];
    }
  }

  /// --- allocate & deallocate functions for rhhda_static allocator --- ///
  void* allocate(const size_t capacity)
  {
    const unsigned char idx = cal_allocator_index(capacity);
    return m_static_node_allocators[idx]->allocate();
  }

  void deallocate(void *ptr, const size_t capacity)
  {
    const unsigned char idx = cal_allocator_index(capacity);
    m_static_node_allocators[idx]->deallocate(ptr);
  }

  static inline unsigned char cal_allocator_index(const size_t size)
  {
    return cal_pos_msb_64bit(size) - 1;
  }


private:
  static constexpr size_t num_elements[20] = {
    1ULL << 0,
    1ULL << 1,
    1ULL << 2,
    1ULL << 3,
    1ULL << 4,
    1ULL << 5,
    1ULL << 6,
    1ULL << 7,
    1ULL << 8,
    1ULL << 9,
    1ULL << 10,
    1ULL << 11,
    1ULL << 12,
    1ULL << 13,
    1ULL << 14,
    1ULL << 15,
    1ULL << 16,
    1ULL << 17,
    1ULL << 18
  };

  ///
  /// \brief cal_pos_msb
  /// \param x
  /// \return postion of msb of x. 1 <= position <= 64;
  ///
  static inline unsigned int cal_pos_msb_64bit(const size_t x)
  {
    return 64 - __builtin_clzll(x-1) + 1;
  }

  std::array<static_node_allocator_wrp*, 20> m_static_node_allocators;
};


///
/// \brief The allocator_holder_sglt class
///  singleton pattern
template<typename _segment_manager_type, size_t element_size, size_t extra_size>
class allocator_holder_sglt
{

 public:

  using segment_manager_type = _segment_manager_type;
  using self_type = allocator_holder_sglt<segment_manager_type, element_size, extra_size>;

  static self_type& instance() {
    static self_type _instance;
    return _instance;
  }

  template <typename U1, size_t U2, size_t U3>
  struct rebind
  {
      using other = allocator_holder_sglt<U1, U2, U3>;
  };

  /// --- allocate & deallocate functions --- ///
  /// \brief allocate
  ///   allocate memory spaces equal to the size of required length of elements.
  ///   if required length is less than the number of elements per a chunk, use node_allocater.
  ///   otherwise, use raw allocator
  /// \param length
  ///   the number of elements want to allocate
  /// \return
  ///   void* pointer for allocated memory spaces
  void* allocate(const size_t capacity)
  {
    if (capacity <= node_allocators_type::kMaxCapacity) {
      return reinterpret_cast<void*>(m_node_allocators->allocate(capacity));
    } else {
      return reinterpret_cast<void*>(m_raw_allocator->allocate(capacity * element_size + extra_size).get());
    }
  }

  ///
  /// \brief deallocate
  ///   Deallocate memory spaces.
  ///   Note: In node_allocater, it dosen't actually deallocate the memory spaces.
  ///         To actually deallocate the memory spaceses, call deallocate_free_blocks().
  /// \param length
  /// \param ptr
  ///
  void deallocate(void *ptr, const size_t capacity)
  {
    if (capacity <= node_allocators_type::kMaxCapacity) {
      m_node_allocators->deallocate(ptr, capacity);
    } else {
      m_raw_allocator->deallocate(boost::interprocess::offset_ptr<char>(reinterpret_cast<char*>(ptr)),
                                 capacity * element_size + extra_size);
    }
  }

  ///
  /// \brief deallocate_free_blocks
  ///   Actually deallocate node_allocater's free chunks.
  ///   The cost of this operation is O(n^2)
  void deallocate_free_blocks()
  {
    m_node_allocators->deallocate_free_blocks();
  }

  void init(segment_manager_type* segment_manager)
  {
    m_raw_allocator = new raw_allocator_type(segment_manager);
    m_node_allocators = new node_allocators_type(segment_manager);
  }

  void destory()
  {
    delete m_raw_allocator;
    delete m_node_allocators;
  }

 private:
  /// raw allocator
  using raw_allocator_type    = boost::interprocess::allocator<char, segment_manager_type>;
  /// node allocators
  using node_allocators_type  = static_node_allocators_holder<segment_manager_type, element_size, extra_size>;

  allocator_holder_sglt() {}
  allocator_holder_sglt(const self_type &)  = delete;
  allocator_holder_sglt(const self_type &&) = delete;
  self_type &operator=(const self_type &)   = delete;
  self_type &operator=(const self_type &&)  = delete;

  /// ---- Private Variables ------ ///
  raw_allocator_type*  m_raw_allocator;
  node_allocators_type* m_node_allocators;
};


struct ALLOC_IN_CORE {};
///
/// \brief allocator_holder_sglt class
///  singleton pattern
///  specialized class for in core allocation
template<size_t element_size, size_t extra_size>
class allocator_holder_sglt <ALLOC_IN_CORE, element_size, extra_size>
{

 public:

  using segment_manager_type = ALLOC_IN_CORE;
  using self_type = allocator_holder_sglt<ALLOC_IN_CORE, element_size, extra_size>;

  static self_type& instance() {
    static self_type _instance;
    return _instance;
  }

  template <typename U1, size_t U2, size_t U3>
  struct rebind
  {
      using other = allocator_holder_sglt<U1, U2, U3>;
  };

  /// --- allocate & deallocate functions --- ///
  /// \brief allocate
  ///   allocate memory spaces equal to the size of required length of elements.
  ///   if required length is less than the number of elements per a chunk, use node_allocater.
  ///   otherwise, use raw allocator
  /// \param length
  ///   the number of elements want to allocate
  /// \return
  ///   void* pointer for allocated memory spaces
  void* allocate(const size_t capacity)
  {
    return malloc(capacity * element_size + extra_size);
  }

  ///
  /// \brief deallocate
  ///   Deallocate memory spaces.
  ///   Note: In node_allocater, it dosen't actually deallocate the memory spaces.
  ///         To actually deallocate the memory spaceses, call deallocate_free_blocks().
  /// \param length
  /// \param ptr
  ///
  void deallocate(void *ptr, const size_t dummy = 0)
  {
    free(ptr);
  }

//  void destory()
//  { }

 private:
  allocator_holder_sglt() {}
  allocator_holder_sglt(const self_type &)  = delete;
  allocator_holder_sglt(const self_type &&) = delete;
  self_type &operator=(const self_type &)   = delete;
  self_type &operator=(const self_type &&)  = delete;

};


///
/// \brief The allocator_holder class
///  --- not using now ---
template<typename segment_manager_type, size_t element_size, size_t extra_size>
class allocator_holder
{
public:
  ///
  /// \brief Constructor
  /// \param segment_manager
  ///
  explicit allocator_holder(segment_manager_type* segment_manager)
    : m_raw_allocator(segment_manager),
      m_node_allocators(segment_manager)
  {
    static_assert(kNodeAllocatorChunkSize >= element_size, "element_size is larger than kNodeAllocatorChunkSize");
  }

  /// --- allocate & deallocate functions --- ///
  /// \brief allocate
  ///   allocate memory spaces equal to the size of required length of elements.
  ///   if required length is less than the number of elements per a chunk, use node_allocater.
  ///   otherwise, use raw allocator
  /// \param length
  ///   the number of elements want to allocate
  /// \return
  ///   void* pointer for allocated memory spaces
  void* allocate(const size_t capacity)
  {
    if (capacity <= node_allocators_type::kMaxCapacity) {
      return reinterpret_cast<void*>(m_node_allocators.allocate(capacity));
    } else {
      return reinterpret_cast<void*>(m_raw_allocator.allocate(capacity * element_size + extra_size).get());
    }
  }

  ///
  /// \brief deallocate
  ///   Deallocate memory spaces.
  ///   Note: In node_allocater, it dosen't actually deallocate the memory spaces.
  ///         To actually deallocate the memory spaceses, call deallocate_free_blocks().
  /// \param length
  /// \param ptr
  ///
  void deallocate(void *ptr, const size_t capacity)
  {
    if (capacity <= node_allocators_type::kMaxCapacity) {
      m_node_allocators.deallocate(ptr, capacity);
    } else {
      m_raw_allocator.deallocate(boost::interprocess::offset_ptr<char>(reinterpret_cast<char*>(ptr)),
                                 capacity * element_size + extra_size);
    }
  }

  ///
  /// \brief deallocate_free_blocks
  ///   Actually deallocate node_allocater's free chunks.
  ///   The cost of this operation is O(n^2)
  void deallocate_free_blocks()
  {
    m_node_allocators.deallocate_free_blocks();
  }


private:
  /// raw allocator
  using raw_allocator_type    = boost::interprocess::allocator<char, segment_manager_type>;
  /// node allocators
  using node_allocators_type  = static_node_allocators_holder<segment_manager_type, element_size, extra_size>;


  /// ---- Private Variables ------ ///
  raw_allocator_type  m_raw_allocator;
  node_allocators_type m_node_allocators;
};



}}

#endif
