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
#include <havoqgt/graphstore/rhhda/rhhda_defs.hpp>
#include <havoqgt/graphstore/rhhda/rhhda_common.hpp>

namespace graphstore {
namespace rhhda {

/// Wrapper class for static_node_allocator to management the class using array data structure
class static_node_allocator_wrp
{
public:
  virtual void* allocate() =0;
  virtual void deallocate(void *) =0;
  virtual void deallocate_free_blocks() =0;
};

template<size_t allocate_size>
class static_node_allocator : public static_node_allocator_wrp
{

public:

  static_node_allocator(segment_manager_t* segment_manager)
    : m_node_allocator(segment_manager)
  {
    static_assert(kNodeAllocatorChunkSize >= allocate_size, "sizeof(rhhda_static) is larger than kNodeAllocatorChunkSize");
  }

  void* allocate()
  {
    return reinterpret_cast<void*>(m_node_allocator.allocate(1).get());
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
  using node_allocator_type  = boost::interprocess::node_allocator<dummy_node_type,  segment_manager_t, kNodesPerChunk>;

  node_allocator_type m_node_allocator;
};


template<size_t element_size, size_t extra_size>
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
  explicit static_node_allocators_holder(segment_manager_t* segment_manager) :
    m_static_node_allocators{{
      new static_node_allocator<element_size * num_elements[ 0] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 1] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 2] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 3] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 4] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 5] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 6] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 7] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 8] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[ 9] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[10] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[11] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[12] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[13] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[14] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[15] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[16] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[17] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[18] + extra_size>(segment_manager),
      new static_node_allocator<element_size * num_elements[19] + extra_size>(segment_manager)
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
    1ULL,
    2ULL,
    4ULL,
    8ULL,
    16ULL,
    32ULL,
    64ULL,
    128ULL,
    256ULL,
    512ULL,
    1024ULL,
    2048ULL,
    4096ULL,
    8192ULL,
    16384ULL,
    32768ULL,
    65536ULL,
    131072ULL,
    262144ULL
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



template<size_t element_size, size_t extra_size>
class allocator_holder_sglt
{

 public:

  using self_type = allocator_holder_sglt<element_size, extra_size>;

  static self_type& instance() {
    static self_type _instance;
    return _instance;
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

  void init(segment_manager_t* segment_manager)
  {
    m_raw_allocator = new raw_allocator_type(segment_manager);
    m_node_allocators = new node_allocators_type(segment_manager);
//    instance();
  }

  void destory()
  {
    delete m_raw_allocator;
    delete m_node_allocators;
  }

private:
  /// raw allocator
  using raw_allocator_type    = boost::interprocess::allocator<char, segment_manager_t>;
  /// node allocators
  using node_allocators_type  = static_node_allocators_holder<element_size, extra_size>;

  allocator_holder_sglt(){}
  allocator_holder_sglt(const self_type &other){}
  self_type &operator=(const self_type &other){}


  /// ---- Private Variables ------ ///
  raw_allocator_type*  m_raw_allocator;
  node_allocators_type* m_node_allocators;
};


template<typename allocator, typename segment_manager_type>
void init_allocator(segment_manager_type* segment_manager)
{
  allocator::instance().init(segment_manager);
}

template<typename allocator>
void destroy_allocator()
{
  allocator::instance().destory();
}


template<size_t element_size, size_t extra_size>
class allocator_holder
{
public:
  ///
  /// \brief Constructor
  /// \param segment_manager
  ///
  explicit allocator_holder(segment_manager_t* segment_manager)
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
    if (capacity >= node_allocators_type::kMaxCapacity) {
      return reinterpret_cast<void*>(m_raw_allocator.allocate(capacity * element_size + extra_size).get());
    } else {
      return reinterpret_cast<void*>(m_node_allocators.allocate(capacity));
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
    if (capacity >= node_allocators_type::kMaxCapacity) {
      m_raw_allocator.deallocate(boost::interprocess::offset_ptr<char>(reinterpret_cast<char*>(ptr)),
                                 capacity * element_size + extra_size);
    } else {
      m_node_allocators.deallocate(ptr, capacity);
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
  using raw_allocator_type    = boost::interprocess::allocator<char, segment_manager_t>;
  /// node allocators
  using node_allocators_type  = static_node_allocators_holder<element_size, extra_size>;


  /// ---- Private Variables ------ ///
  raw_allocator_type  m_raw_allocator;
  node_allocators_type m_node_allocators;
};



}}

#endif
















//  template<size_t element_size, size_t extra_size>
//  class allocator_holder
//  {
//  public:
//    ///
//    /// \brief Constructor
//    /// \param segment_manager
//    ///
//    explicit allocator_holder(segment_manager_t* segment_manager)
//      : m_raw_allocator(segment_manager),
//        m_node_allocator(segment_manager)
//    {
//      static_assert(kNodeAllocatorChunkSize >= element_size, "element_size is larger than kNodeAllocatorChunkSize");
//    }

//    /// --- allocate & deallocate functions --- ///
//    /// \brief allocate
//    ///   allocate memory spaces equal to the size of required length of elements.
//    ///   if required length is less than the number of elements per a chunk, use node_allocater.
//    ///   otherwise, use raw allocator
//    /// \param length
//    ///   the number of elements want to allocate
//    /// \return
//    ///   void* pointer for allocated memory spaces
//    void* allocate(const size_t length)
//    {
//      if (length >= kNodesPerChunk) {
//        return reinterpret_cast<void*>(m_raw_allocator.allocate(length * element_size).get());
//      } else {
//        return reinterpret_cast<void*>(m_node_allocator.allocate(length).get());
//      }
//    }

//    ///
//    /// \brief deallocate
//    ///   Deallocate memory spaces.
//    ///   Note: In node_allocater, it dosen't actually deallocate the memory spaces.
//    ///         To actually deallocate the memory spaceses, call deallocate_free_blocks().
//    /// \param length
//    /// \param ptr
//    ///
//    void deallocate(void *ptr, const size_t length)
//    {
//      if (length >= kNodesPerChunk) {
//        m_raw_allocator.deallocate(boost::interprocess::offset_ptr<char>(reinterpret_cast<char*>(ptr)),
//                                   length * element_size);
//      } else {
//        m_node_allocator.deallocate(boost::interprocess::offset_ptr<dummy_node_type>(reinterpret_cast<dummy_node_type*>(ptr)),
//                                    length);
//      }
//    }

//    ///
//    /// \brief deallocate_free_blocks
//    ///   Actually deallocate node_allocater's free chunks.
//    ///   The cost of this operation is O(n^2)
//    void deallocate_free_blocks()
//    {
//      m_node_allocator.deallocate_free_blocks();
//    }


//  private:
//    /// raw allocator
//    using raw_allocator_type = boost::interprocess::allocator<char, segment_manager_t>;

//    /// node allocator
//    struct dummy_node_type {
//      char dummy[element_size];
//    };
//    enum NodesPerChunk {
//      kNodesPerChunk  = kNodeAllocatorChunkSize / element_size,
//    };
//    using node_allocator_type  = boost::interprocess::node_allocator<dummy_node_type,  segment_manager_t, kNodesPerChunk>;


//    /// ---- Private Variables ------ ///
//    raw_allocator_type  m_raw_allocator;
//    node_allocator_type m_node_allocator;
//  };


//template<typename direct_value_type, size_t element_size>
//class rhhda_allocator_holder
//{
//public:

//  using segment_manager_t =  boost::interprocess::managed_mapped_file::segment_manager;

//  explicit rhhda_allocator_holder(segment_manager_t* segment_manager)
//    : m_raw_allocator(segment_manager),
//      m_key_array_allocator(segment_manager),
//      m_rhhda_satic_allocator1(segment_manager),
//      m_rhhda_satic_allocator2(segment_manager),
//      m_rhhda_satic_allocator3(segment_manager),
//      m_rhhda_satic_allocator4(segment_manager),
//      m_rhhda_satic_allocator5(segment_manager),
//      m_rhhda_satic_allocator6(segment_manager),
//      m_rhhda_satic_allocator7(segment_manager),
//      m_rhhda_satic_allocator8(segment_manager),
//      m_rhhda_satic_allocator9(segment_manager),
//      m_rhhda_satic_allocator10(segment_manager),
//      m_rhhda_satic_allocator11(segment_manager),
//      m_rhhda_satic_allocator12(segment_manager),
//      m_rhhda_satic_allocator13(segment_manager),
//      m_rhhda_satic_allocator14(segment_manager),
//      m_rhhda_satic_allocator15(segment_manager),
//      m_rhhda_satic_allocator16(segment_manager),
//      m_rhhda_satic_allocator17(segment_manager)
//  {
//    static_assert(kNodeAllocatorChunkSize >= sizeof(KeyType), "sizeof(KeyType) is larger than kNodeAllocatorChunkSize");
//    static_assert(kNodeAllocatorChunkSize >= sizeof(kRhhdaStaticCapacity_17), "sizeof(rhhda_static) is larger than kNodeAllocatorChunkSize");
//  }

//  /// --- allocate & deallocate functions for raw allocator --- ///
//  void* allocate_raw(size_t size)
//  {
//    return reinterpret_cast<void*>(m_raw_allocator.allocate(size).get());
//  }

//  void deallocate_raw(size_t size, void *ptr)
//  {
//    m_raw_allocator.deallocate(boost::interprocess::offset_ptr<unsigned char>(reinterpret_cast<unsigned char*>(ptr)), size);
//  }

//  /// --- allocate & deallocate functions for key array allocator --- ///
//  void* allocate_key_array(size_t length)
//  {
//    const size_t size = cal_next_highest_power_of_2_8bit(length);
//    return reinterpret_cast<void*>(m_key_array_allocator.allocate(size).get());
//  }

//  void deallocate_key_array(size_t length, void *ptr)
//  {
//    const size_t size = cal_next_highest_power_of_2_8bit(length);
//    m_key_array_allocator.deallocate(boost::interprocess::offset_ptr<direct_value_type>(reinterpret_cast<direct_value_type*>(ptr)), size);
//  }

//  /// --- allocate & deallocate functions for rhhda_satic allocator --- ///
//  void* allocate_rhhda_static_1()
//  {
//    reinterpret_cast<void*>(m_rhhda_satic_allocator1.allocate(1).get());
//  }
//  void deallocate_rhhda_static_1(void *ptr) {
//    m_rhhda_satic_allocator1.deallocate(boost::interprocess::offset_ptr<RhhdaStaticType>(reinterpret_cast<RhhdaStaticType*>(ptr)), 1);
//  }

//  /// =========    Private  ============ ///
//private:
//  /// raw allocator
//  using RawAllocatorType = boost::interprocess::allocator<unsigned char, segment_manager_t>;

//  /// key array's allocator
//  enum KeyArrayCapacity : size_t {
//    kKeyArrayCapacity1 = (1ULL) << 1,
//    kKeyArrayCapacity2 = (1ULL) << 2,
//    kKeyArrayCapacity3 = (1ULL) << 3
//  };
//  using KeyArrayAllocatorType = boost::interprocess::node_allocator<direct_value_type, segment_manager_t, kNodeAllocatorChunkSize / sizeof(direct_value_type)>;


//  /// rhhda_static's allocators
//  enum RhhdaStaticCapacity : size_t {
//    kRhhdaStaticCapacity_1  = (1ULL) << 3,
//    kRhhdaStaticCapacity_2  = (1ULL) << 4,
//    kRhhdaStaticCapacity_3  = (1ULL) << 5,
//    kRhhdaStaticCapacity_4  = (1ULL) << 6,
//    kRhhdaStaticCapacity_5  = (1ULL) << 7,
//    kRhhdaStaticCapacity_6  = (1ULL) << 8,
//    kRhhdaStaticCapacity_7  = (1ULL) << 9,
//    kRhhdaStaticCapacity_8  = (1ULL) << 10,
//    kRhhdaStaticCapacity_9  = (1ULL) << 11,
//    kRhhdaStaticCapacity_10 = (1ULL) << 12,
//    kRhhdaStaticCapacity_11 = (1ULL) << 13,
//    kRhhdaStaticCapacity_12 = (1ULL) << 14,
//    kRhhdaStaticCapacity_13 = (1ULL) << 15,
//    kRhhdaStaticCapacity_14 = (1ULL) << 16,
//    kRhhdaStaticCapacity_15 = (1ULL) << 17,
//    kRhhdaStaticCapacity_16 = (1ULL) << 18,
//    kRhhdaStaticCapacity_17 = (1ULL) << 19
//  };

//  struct dummy_node_type {
//    char dummy[element_size];
//  };

//  enum RhhdaStaticNodesPerChunk {
//    kRhhStaticNodesPerChunk1  = kNodeAllocatorChunkSize / sizeof(dummy_node_type<unsigned char, uint64_t, unsigned char, kRhhdaStaticCapacity_1>),
//  };
//  using RhhdaSaticAllocatorType1  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_1,  segment_manager_t, kRhhStaticNodesPerChunk1>;
//  using RhhdaSaticAllocatorType2  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_2,  segment_manager_t, kRhhStaticNodesPerChunk2>;
//  using RhhdaSaticAllocatorType3  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_3,  segment_manager_t, kRhhStaticNodesPerChunk3>;
//  using RhhdaSaticAllocatorType4  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_4,  segment_manager_t, kRhhStaticNodesPerChunk4>;
//  using RhhdaSaticAllocatorType5  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_5,  segment_manager_t, kRhhStaticNodesPerChunk5>;
//  using RhhdaSaticAllocatorType6  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_6,  segment_manager_t, kRhhStaticNodesPerChunk6>;
//  using RhhdaSaticAllocatorType7  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_7,  segment_manager_t, kRhhStaticNodesPerChunk7>;
//  using RhhdaSaticAllocatorType8  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_8,  segment_manager_t, kRhhStaticNodesPerChunk8>;
//  using RhhdaSaticAllocatorType9  = boost::interprocess::node_allocator<kRhhdaStaticCapacity_9,  segment_manager_t, kRhhStaticNodesPerChunk9>;
//  using RhhdaSaticAllocatorType10 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_10, segment_manager_t, kRhhStaticNodesPerChunk10>;
//  using RhhdaSaticAllocatorType11 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_11, segment_manager_t, kRhhStaticNodesPerChunk11>;
//  using RhhdaSaticAllocatorType12 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_12, segment_manager_t, kRhhStaticNodesPerChunk12>;
//  using RhhdaSaticAllocatorType13 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_13, segment_manager_t, kRhhStaticNodesPerChunk13>;
//  using RhhdaSaticAllocatorType14 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_14, segment_manager_t, kRhhStaticNodesPerChunk14>;
//  using RhhdaSaticAllocatorType15 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_15, segment_manager_t, kRhhStaticNodesPerChunk15>;
//  using RhhdaSaticAllocatorType16 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_16, segment_manager_t, kRhhStaticNodesPerChunk16>;
//  using RhhdaSaticAllocatorType17 = boost::interprocess::node_allocator<kRhhdaStaticCapacity_17, segment_manager_t, kRhhStaticNodesPerChunk17>;


//  /// ---- Private Functions ----- ///
//  inline unsigned char cal_rhhdastatic_allocator_index(size_t size)
//  {
//    if (node_capacity_17 >= size) return 16;
//    assert(size > node_capacity_1);
//    const int x = cal_pos_msb(size);
//    return x2 - 4;
//  }

//  static inline unsigned char cal_next_highest_power_of_2_8bit(unsigned char x)
//  {
//    --x;
//    x |= x >> 1;
//    x |= x >> 2;
//    x |= x >> 4;
//    ++x;
//    return x;
//  }

//  /// \return postion of msb of x. 1 <= position <= 64;
//  static inline unsigned int cal_pos_msb(size_t x)
//  {
//    return 64 - __builtin_clzll(size-1) + 1;
//  }


//  /// ---- Private Variables ------ ///
//  RawAllocatorType m_raw_allocator;
//  KeyArrayAllocatorType m_key_array_allocator;
//  RhhdaSaticAllocatorType1 m_rhhda_satic_allocator1;
//  RhhdaSaticAllocatorType2 m_rhhda_satic_allocator2;
//  RhhdaSaticAllocatorType3 m_rhhda_satic_allocator3;
//  RhhdaSaticAllocatorType4 m_rhhda_satic_allocator4;
//  RhhdaSaticAllocatorType5 m_rhhda_satic_allocator5;
//  RhhdaSaticAllocatorType6 m_rhhda_satic_allocator6;
//  RhhdaSaticAllocatorType7 m_rhhda_satic_allocator7;
//  RhhdaSaticAllocatorType8 m_rhhda_satic_allocator8;
//  RhhdaSaticAllocatorType9 m_rhhda_satic_allocator9;
//  RhhdaSaticAllocatorType10 m_rhhda_satic_allocator10;
//  RhhdaSaticAllocatorType11 m_rhhda_satic_allocator11;
//  RhhdaSaticAllocatorType12 m_rhhda_satic_allocator12;
//  RhhdaSaticAllocatorType13 m_rhhda_satic_allocator13;
//  RhhdaSaticAllocatorType14 m_rhhda_satic_allocator14;
//  RhhdaSaticAllocatorType15 m_rhhda_satic_allocator15;
//  RhhdaSaticAllocatorType16 m_rhhda_satic_allocator16;
//  RhhdaSaticAllocatorType17 m_rhhda_satic_allocator17;

//};

///// Wrapper class for rhhda_allocator to management the class using array data structure
//class rhhda_allocator_wrp
//{
//public:
//  virtual void* allocate() = 0;
//  virtual void deallocate(void *) = 0;
//};


//template<typename KeyType, typename  ValueType, size_t Capacity>
//class rhhda_static_node_allocator : rhhda_allocator_wrp
//{

//public:
//  using segment_manager_t =  boost::interprocess::managed_mapped_file::segment_manager;
//  using RhhdaStaticType  = rhhda_static<KeyType, ValueType, Capacity>;
//  enum { NodesPerChunk  = kNodeAllocatorChunkSize / sizeof(RhhdaStaticType) };
//  using RhhSaticAllocatorType = boost::interprocess::node_allocator<RhhdaStaticType,  segment_manager_t,  NodesPerChunk>;

//  rhhda_static_node_allocator(segment_manager_t* segment_manager)
//    : m_alloocator(segment_manager) {
//    static_assert(kNodeAllocatorChunkSize >= sizeof(RhhdaStaticType), "sizeof(rhhda_static) is larger than kNodeAllocatorChunkSize");
//  }

//  void* allocate() {
//    return reinterpret_cast<void*>(m_alloocator.allocate(1).get());
//  }

//  void deallocate(void *ptr) {
//    m_alloocator.deallocate(boost::interprocess::offset_ptr<RhhdaStaticType>(reinterpret_cast<RhhdaStaticType*>(ptr)), 1);
//  }

//private:
//  RhhSaticAllocatorType m_alloocator;
//};

/// Holds 3 type of allocators:
/// 1. raw allocator: each element is usingned char. not using node_allocator
/// 2. key-array's allocator: allocate array of KeyType element, where KeyType is a template parameter.
/// 3. rhhda_static's allocator: allocate rhhda_static.
///
/// The size of chunks in node_allocator is kNodeAllocatorChunkSize which is defined in rhhda_def.hpp
//template<typename KeyType, typename ValueType>
//class rhhda_allocator_holder
//{
//public:

//  using segment_manager_t =  boost::interprocess::managed_mapped_file::segment_manager;

//  /// raw allocator
//  using raw_allocator_t = boost::interprocess::allocator<unsigned char, segment_manager_t>;


//  /// key array's allocator
//  enum CapacityKeyArray {
//    kCapacityNormalArray1 = 2ULL,
//    kCapacityNormalArray2 = 4ULL,
//    kCapacityNormalArray3 = 8ULL
//  };
//  static constexpr size_t kCapacityKeyArray[3] = {kCapacityNormalArray1, kCapacityNormalArray2, kCapacityNormalArray3};
//  using key_array_allocator_t = boost::interprocess::node_allocator<KeyType, segment_manager_t, kNodeAllocatorChunkSize / sizeof(KeyType)>;

//  ///
//  /// \brief The CapacityRhhdaStatic enum
//  /// !!!: we assume that the number of threshhold values of capacity is less than 256
//  enum CapacityRhhdaStatic : size_t {
//    node_capacity_1 = 8ULL,
//    node_capacity_2 = 16ULL,
//    node_capacity_3 = 32ULL,
//    node_capacity_4 = 64ULL,
//    node_capacity_5 = 128ULL,
//    node_capacity_6 = 256ULL,
//    node_capacity_7 = 512ULL,
//    node_capacity_8 = 1024ULL,
//    node_capacity_9 = 2048ULL,
//    node_capacity_10 = 4096ULL,
//    node_capacity_11 = 8192ULL,
//    node_capacity_12 = 16384ULL,
//    node_capacity_13 = 32768ULL,
//    node_capacity_14 = 65536ULL,
//    node_capacity_15 = 131072ULL,
//    node_capacity_16 = 262144ULL,
//    node_capacity_17 = 524288ULL
//  };
//  static constexpr size_t node_capacity[17] = {node_capacity_1,
//                                               node_capacity_2,
//                                               node_capacity_3,
//                                               node_capacity_4,
//                                               node_capacity_5,
//                                               node_capacity_6,
//                                               node_capacity_7,
//                                               node_capacity_8,
//                                               node_capacity_9,
//                                               node_capacity_10,
//                                               node_capacity_11,
//                                               node_capacity_12,
//                                               node_capacity_13,
//                                               node_capacity_14,
//                                               node_capacity_15,
//                                               node_capacity_16,
//                                               node_capacity_17
//                                              };

//  ///
//  /// \brief rhhda_allocator_holder
//  /// \param segment_manager
//  ///
//  explicit rhhda_allocator_holder(segment_manager_t* segment_manager)
//    : m_raw_allocator(segment_manager),
//      m_key_array_allocator(segment_manager)
//  {
//    static_assert(kNodeAllocatorChunkSize >= sizeof(KeyType), "sizeof(KeyType) is larger than kNodeAllocatorChunkSize");

//    m_rhhda_static_node_allocators[0]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_1>(segment_manager);
//    m_rhhda_static_node_allocators[1]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_2>(segment_manager);
//    m_rhhda_static_node_allocators[2]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_3>(segment_manager);
//    m_rhhda_static_node_allocators[3]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_4>(segment_manager);
//    m_rhhda_static_node_allocators[4]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_5>(segment_manager);
//    m_rhhda_static_node_allocators[5]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_6>(segment_manager);
//    m_rhhda_static_node_allocators[6]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_7>(segment_manager);
//    m_rhhda_static_node_allocators[7]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_8>(segment_manager);
//    m_rhhda_static_node_allocators[8]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_9>(segment_manager);
//    m_rhhda_static_node_allocators[9]  = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_10>(segment_manager);
//    m_rhhda_static_node_allocators[10] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_11>(segment_manager);
//    m_rhhda_static_node_allocators[11] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_12>(segment_manager);
//    m_rhhda_static_node_allocators[12] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_13>(segment_manager);
//    m_rhhda_static_node_allocators[13] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_14>(segment_manager);
//    m_rhhda_static_node_allocators[14] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_15>(segment_manager);
//    m_rhhda_static_node_allocators[15] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_16>(segment_manager);
//    m_rhhda_static_node_allocators[16] = new rhhda_static_node_allocator<KeyType, ValueType, node_capacity_17>(segment_manager);

//  }


//  /// --- allocate & deallocate functions for raw allocator --- ///
//  void* allocate_raw(size_t length)
//  {
//    return reinterpret_cast<void*>(m_raw_allocator.allocate(length).get());
//  }

//  void deallocate_raw(size_t length, void *ptr)
//  {
//    m_raw_allocator.deallocate(boost::interprocess::offset_ptr<unsigned char>(reinterpret_cast<unsigned char*>(ptr)), length);
//  }

//  /// --- allocate & deallocate functions for key array allocator --- ///
//  void* allocate_key_array(int size_no)
//  {
//    return reinterpret_cast<void*>(m_key_array_allocator.allocate(kCapacityKeyArray[size_no]).get());
//  }

//  void deallocate_key_array(int size_no, void *ptr)
//  {
//    m_key_array_allocator.deallocate(boost::interprocess::offset_ptr<KeyType>(reinterpret_cast<KeyType*>(ptr)), kCapacityKeyArray[size_no]);
//  }

//  /// --- allocate & deallocate functions for rhhda_static allocator --- ///
//  void* allocate_rhhdasatic(const unsigned char allocator_no)
//  {
//    if (allocator_no < m_rhhda_static_node_allocators.size())
//      m_rhhda_static_node_allocators[allocator_no]->allocate();
//    else
//      return nullptr;
//  }

//  void deallocate_rhhdasatic(const unsigned char allocator_no, void *ptr)
//  {
//    m_rhhda_static_node_allocators[allocator_no]->deallocate(ptr);
//  }

//  inline unsigned char cal_allocator_index(size_t size)
//  {
//    if (node_capacity_17 >= size) return 16;
//    assert(size > node_capacity_1);
//    const int x = cal_pos_msb(size);
//    return x2 - 4;
//  }

//  ///
//  /// \brief cal_pos_msb
//  /// \param x
//  /// \return postion of msb of x. 1 <= position <= 64;
//  ///
//  static inline unsigned int cal_pos_msb(size_t x)
//  {
//    return 64 - __builtin_clzll(size-1) + 1;
//  }

//private:

//  raw_allocator_t m_raw_allocator;
//  key_array_allocator_t m_key_array_allocator;
//  std::array<std::unique_ptr<rhhda_allocator_wrp>, 17> m_rhhda_static_node_allocators;
//};



////template<typename ArrayType>
////class rhhda_normalarray_allocator : rhhda_allocator_wrp
////{

////public:
////  enum { NodesPerChunk  = kNodeAllocatorChunkSize / sizeof(ArrayType) };
////  using segment_manager_t =  boost::interprocess::managed_mapped_file::segment_manager;
////  using RhhdaNormalarrayAllocatorType = boost::interprocess::node_allocator<ArrayType, segment_manager_t, NodesPerChunk>;

////  rhhda_normalarray_allocator(segment_manager_t* segment_manager)
////    : m_alloocator(segment_manager) { }

////   void* allocate() {
////     return reinterpret_cast<void*>(m_alloocator.allocate(1).get());
////   }

////   void deallocate(void *ptr) {
////     m_alloocator.deallocate(boost::interprocess::offset_ptr<ArrayType>(reinterpret_cast<ArrayType*>(ptr)), 1);
////   }

////private:
////  RhhdaNormalarrayAllocatorType m_alloocator;
////};



///  =========================================================================== ///
///                             rhhda_allocator_holder
///  =========================================================================== ///
//template<typename KeyType, typename  ValueType>
//class rhhda_allocator_holder
//{
//public:

//  enum {
//    capacityNormalArray1 = 2ULL,
//    capacityNormalArray2 = 4ULL,
//    capacityNormalArray3 = 8ULL
//  };
//  using segment_manager_t =  boost::interprocess::managed_mapped_file::segment_manager;

//  enum {
//    node_capacity_1 = 8ULL,
//    node_capacity_2 = 16ULL,
//    node_capacity_3 = 32ULL,
//    node_capacity_4 = 64ULL,
//    node_capacity_5 = 128ULL,
//    node_capacity_6 = 256ULL,
//    node_capacity_7 = 512ULL,
//    node_capacity_8 = 1024ULL,
//    node_capacity_9 = 2048ULL,
//    node_capacity_10 = 4096ULL,
//    node_capacity_11 = 8192ULL,
//    node_capacity_12 = 16384ULL,
//    node_capacity_13 = 32768ULL,
//    node_capacity_14 = 65536ULL,
//    node_capacity_15 = 131072ULL,
//    node_capacity_16 = 262144ULL,
//    node_capacity_17 = 524288ULL
//  };

//  using RhhdaStaticType_1  = rhhda_static<KeyType, ValueType, node_capacity_1>;
//  using RhhdaStaticType_2  = rhhda_static<KeyType, ValueType, node_capacity_2>;
//  using RhhdaStaticType_3  = rhhda_static<KeyType, ValueType, node_capacity_3>;
//  using RhhdaStaticType_4  = rhhda_static<KeyType, ValueType, node_capacity_4>;
//  using RhhdaStaticType_5  = rhhda_static<KeyType, ValueType, node_capacity_5>;
//  using RhhdaStaticType_6  = rhhda_static<KeyType, ValueType, node_capacity_6>;
//  using RhhdaStaticType_7  = rhhda_static<KeyType, ValueType, node_capacity_7>;
//  using RhhdaStaticType_8  = rhhda_static<KeyType, ValueType, node_capacity_8>;
//  using RhhdaStaticType_9  = rhhda_static<KeyType, ValueType, node_capacity_9>;
//  using RhhdaStaticType_10 = rhhda_static<KeyType, ValueType, node_capacity_10>;
//  using RhhdaStaticType_11 = rhhda_static<KeyType, ValueType, node_capacity_11>;
//  using RhhdaStaticType_12 = rhhda_static<KeyType, ValueType, node_capacity_12>;
//  using RhhdaStaticType_13 = rhhda_static<KeyType, ValueType, node_capacity_13>;
//  using RhhdaStaticType_14 = rhhda_static<KeyType, ValueType, node_capacity_14>;
//  using RhhdaStaticType_15 = rhhda_static<KeyType, ValueType, node_capacity_15>;
//  using RhhdaStaticType_16 = rhhda_static<KeyType, ValueType, node_capacity_16>;
//  using RhhdaStaticType_17 = rhhda_static<KeyType, ValueType, node_capacity_17>;

//  enum {
//    NodesPerChunk_1  = kNodeAllocatorSize / sizeof(RhhdaStaticType_1),
//    NodesPerChunk_2  = kNodeAllocatorSize / sizeof(RhhdaStaticType_2),
//    NodesPerChunk_3  = kNodeAllocatorSize / sizeof(RhhdaStaticType_3),
//    NodesPerChunk_4  = kNodeAllocatorSize / sizeof(RhhdaStaticType_4),
//    NodesPerChunk_5  = kNodeAllocatorSize / sizeof(RhhdaStaticType_5),
//    NodesPerChunk_6  = kNodeAllocatorSize / sizeof(RhhdaStaticType_6),
//    NodesPerChunk_7  = kNodeAllocatorSize / sizeof(RhhdaStaticType_7),
//    NodesPerChunk_8  = kNodeAllocatorSize / sizeof(RhhdaStaticType_8),
//    NodesPerChunk_9  = kNodeAllocatorSize / sizeof(RhhdaStaticType_9),
//    NodesPerChunk_10 = kNodeAllocatorSize / sizeof(RhhdaStaticType_10),
//    NodesPerChunk_11 = kNodeAllocatorSize / sizeof(RhhdaStaticType_11),
//    NodesPerChunk_12 = kNodeAllocatorSize / sizeof(RhhdaStaticType_12),
//    NodesPerChunk_13 = kNodeAllocatorSize / sizeof(RhhdaStaticType_13),
//    NodesPerChunk_14 = kNodeAllocatorSize / sizeof(RhhdaStaticType_14),
//    NodesPerChunk_15 = kNodeAllocatorSize / sizeof(RhhdaStaticType_15),
//    NodesPerChunk_16 = kNodeAllocatorSize / sizeof(RhhdaStaticType_16),
//    NodesPerChunk_17 = kNodeAllocatorSize / sizeof(RhhdaStaticType_17)
//  };

//  typedef boost::interprocess::node_allocator<RhhdaStaticType_1,  segment_manager_t,  NodesPerChunk_1>   RhhSaticAllocatorType_1;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_2,  segment_manager_t,  NodesPerChunk_2>   RhhSaticAllocatorType_2;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_3,  segment_manager_t,  NodesPerChunk_3>   RhhSaticAllocatorType_3;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_4,  segment_manager_t,  NodesPerChunk_4>   RhhSaticAllocatorType_4;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_5,  segment_manager_t,  NodesPerChunk_5>   RhhSaticAllocatorType_5;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_6,  segment_manager_t,  NodesPerChunk_6>   RhhSaticAllocatorType_6;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_7,  segment_manager_t,  NodesPerChunk_7>   RhhSaticAllocatorType_7;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_8,  segment_manager_t,  NodesPerChunk_8>   RhhSaticAllocatorType_8;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_9,  segment_manager_t,  NodesPerChunk_9>   RhhSaticAllocatorType_9;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_10, segment_manager_t,  NodesPerChunk_10>  RhhSaticAllocatorType_10;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_11, segment_manager_t,  NodesPerChunk_11>  RhhSaticAllocatorType_11;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_12, segment_manager_t,  NodesPerChunk_12>  RhhSaticAllocatorType_12;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_13, segment_manager_t,  NodesPerChunk_13>  RhhSaticAllocatorType_13;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_14, segment_manager_t,  NodesPerChunk_14>  RhhSaticAllocatorType_14;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_15, segment_manager_t,  NodesPerChunk_15>  RhhSaticAllocatorType_15;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_16, segment_manager_t,  NodesPerChunk_16>  RhhSaticAllocatorType_16;
//  typedef boost::interprocess::node_allocator<RhhdaStaticType_17, segment_manager_t,  NodesPerChunk_17>  RhhSaticAllocatorType_17;
//  typedef boost::interprocess::node_allocator<uint64_t, segment_manager_t, 512> normalarray_allocator_t;
//  typedef boost::interprocess::allocator<unsigned char, segment_manager_t> raw_allocator_t;

//  explicit rhhda_allocator_holder(segment_manager_t* segment_manager)
//    : rhh_static_allocator_1(segment_manager)
//    , rhh_static_allocator_2(segment_manager)
//    , rhh_static_allocator_3(segment_manager)
//    , rhh_static_allocator_4(segment_manager)
//    , rhh_static_allocator_5(segment_manager)
//    , rhh_static_allocator_6(segment_manager)
//    , rhh_static_allocator_7(segment_manager)
//    , rhh_static_allocator_8(segment_manager)
//    , rhh_static_allocator_9(segment_manager)
//    , rhh_static_allocator_10(segment_manager)
//    , rhh_static_allocator_11(segment_manager)
//    , rhh_static_allocator_12(segment_manager)
//    , rhh_static_allocator_13(segment_manager)
//    , rhh_static_allocator_14(segment_manager)
//    , rhh_static_allocator_15(segment_manager)
//    , rhh_static_allocator_16(segment_manager)
//    , rhh_static_allocator_17(segment_manager)
//    , allocator_normalarray(segment_manager)
//    , allocator_raw(segment_manager)
//  { }


//  /// size = capacity * (1 + 8 + 1) + 8
//  /// probedistance = 1 byte, key = 8 byte, value block = 1 byte
//  RhhSaticAllocatorType_1 rhh_static_allocator_1;
//  RhhSaticAllocatorType_2 rhh_static_allocator_2;
//  RhhSaticAllocatorType_3 rhh_static_allocator_3;
//  RhhSaticAllocatorType_4 rhh_static_allocator_4;
//  RhhSaticAllocatorType_5 rhh_static_allocator_5;
//  RhhSaticAllocatorType_6 rhh_static_allocator_6;
//  RhhSaticAllocatorType_7 rhh_static_allocator_7;
//  RhhSaticAllocatorType_8 rhh_static_allocator_8;
//  RhhSaticAllocatorType_9 rhh_static_allocator_9;
//  RhhSaticAllocatorType_10 rhh_static_allocator_10;
//  RhhSaticAllocatorType_11 rhh_static_allocator_11;
//  RhhSaticAllocatorType_12 rhh_static_allocator_12;
//  RhhSaticAllocatorType_13 rhh_static_allocator_13;
//  RhhSaticAllocatorType_14 rhh_static_allocator_14;
//  RhhSaticAllocatorType_15 rhh_static_allocator_15;
//  RhhSaticAllocatorType_16 rhh_static_allocator_16;
//  RhhSaticAllocatorType_17 rhh_static_allocator_17;
//  normalarray_allocator_t allocator_normalarray;
//  raw_allocator_t allocator_raw;

//private:
//  // static inline uint64_t cal_next_highest_power_of_2(uint64_t x)
//  // {
//  //   --x;
//  //   x != x >> 1ULL;
//  //   x != x >> 2ULL;
//  //   x != x >> 4ULL;
//  //   x != x >> 8ULL;
//  //   x != x >> 16ULL;
//  //   x != x >> 32ULL;
//  //   return ++x;
//  // }

//};
