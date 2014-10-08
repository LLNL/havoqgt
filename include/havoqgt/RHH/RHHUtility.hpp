#ifndef HAVOQGT_MPI_RHHUTILITY_HPP_INCLUDED
#define HAVOQGT_MPI_RHHUTILITY_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>

namespace RHH {
  
  enum UpdateErrors {
    kSucceed,
    kDuplicated,
    kReachingFUllCapacity,
    kLongProbedistance
  };
  
  static const uint64_t kCapacityGrowingFactor = 2ULL;
  typedef void* NoValueType;
  
  ///  =========================================================================== ///
  ///                             Allocator Holder
  ///  =========================================================================== ///
  class AllocatorsHolder
  {
    namespace bip = boost::interprocess;
    typedef bip::managed_mapped_file  mapped_t;
    typedef mapped_t::segment_manager segment_manager_t;
    
  public:
    enum {
      capacityNormalArray1 = 2,
      capacityNormalArray2 = 4,
      capacityNormalArray3 = 8
      capacityRHHStatic1 = 8,
      capacityRHHStatic2 = 16,
      capacityRHHStatic3 = 32,
      capacityRHHStatic4 = 64,
      capacityRHHStatic5 = 128,
      capacityRHHStatic6 = 256
    }
    
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity1> RHHStaticNoVal1;
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity2> RHHStaticNoVal2;
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity3> RHHStaticNoVal3;
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity4> RHHStaticNoVal4;
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity5> RHHStaticNoVal5;
    typename RHHStatic<uint64_t, RHHStatic::NoValueType, capacity6> RHHStaticNoVal6;
    
    /// size = capacity * (1 + 8 + 1) + 8
    /// probedistance = 1 byte, key = 8 byte, value block = 1 byte
    bip::node_allocator<RHHStaticNoVal1, segment_manager_t,  47>  allocator_rhh_noval_1;
    bip::node_allocator<RHHStaticNoVal2, segment_manager_t,  24>  allocator_rhh_noval_2;
    bip::node_allocator<RHHStaticNoVal3, segment_manager_t,  12>  allocator_rhh_noval_3;
    bip::node_allocator<RHHStaticNoVal4, segment_manager_t,   6>  allocator_rhh_noval_4;
    bip::node_allocator<RHHStaticNoVal5, segment_manager_t,   3>  allocator_rhh_noval_5;
    bip::node_allocator<RHHStaticNoVal6, segment_manager_t,   1>  allocator_rhh_noval_6;
    bip::node_allocator<uint64_t, segment_manager_t, 2048> allocator_normalarray_1;
    bip::node_allocator<uint64_t, segment_manager_t, 1024> allocator_normalarray_2;
    bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_normalarray_3;
    bip::allocator<unsigned char, segment_manager_t> allocator_raw;
    
    explicit AllocatorsHolder(segment_manager_t segment_manager)
    : allocator_rhh_noval_1(segment_manager)
    , allocator_rhh_noval_2(segment_manager)
    , allocator_rhh_noval_3(segment_manager)
    , allocator_rhh_noval_4(segment_manager)
    , allocator_rhh_noval_5(segment_manager)
    , allocator_rhh_noval_6(segment_manager)
    , allocator_normalarray_1(segment_manager)
    , allocator_normalarray_2(segment_manager)
    , allocator_normalarray_3(segment_manager)
    , allocator_raw(segment_manager)
    { }
    
    void* allocate_normal_array(const uint64_t required_capacity)
    {
      if (required_capacity <= capacityNormalArray1) {
        return reinterpret_cast<void*>(allocator_normalarray_1.allocate(1).get());
      } else if (required_capacity <= capacityNormalArray2) {
        return reinterpret_cast<void*>(allocator_normalarray_2.allocate(1).get());
      } else if (required_capacity <= capacityNormalArray3) {
        return reinterpret_cast<void*>(allocator_normalarray_3.allocate(1).get());
      } else {
        assert(false);
      }
    }
    
    // void* allocate_rhh_static(const uint64_t required_capacity, uint64_t* allocated_capacity)
    // {
    //   if (required_capacity <= capacityRHHStatic1) {
    //     *allocated_capacity = capacityRHHStatic1;
    //     return reinterpret_cast<void*>(allocator_rhh_1.allocate(1).get());
    //   } else if (required_capacity <= capacityRHHStatic2){
    //     *allocated_capacity = capacityRHHStatic2;
    //     return reinterpret_cast<void*>(allocator_rhh_2.allocate(1).get());
    //   } else if (required_capacity <= capacityRHHStatic3){
    //     *allocated_capacity = capacityRHHStatic3;
    //     return reinterpret_cast<void*>(allocator_rhh_3.allocate(1).get());
    //   } else if (required_capacity <= capacityRHHStatic4){
    //     *allocated_capacity = capacityRHHStatic4;
    //     return reinterpret_cast<void*>(allocator_rhh_4.allocate(1).get());
    //   } else if (required_capacity <= capacityRHHStatic5){
    //     *allocated_capacity = capacityRHHStatic5;
    //     return reinterpret_cast<void*>(allocator_rhh_5.allocate(1).get());
    //   } else if (required_capacity <= capacityRHHStatic6){
    //     *allocated_capacity = capacityRHHStatic6;
    //     return reinterpret_cast<void*>(allocator_rhh_6.allocate(1).get());
    //   } else {
    //    assert(false);
    //   }
    // }
    
    void* allocate_rhh_main(const uint64_t required_capacity)
    {
      const uint64_t memsize = RHHMain<KeyType, ValueType>.cal_memsize(required_capacity);
      return reinterpret_cast<void*>(allocator_raw.allocate(memsize).get());
    }
    
    void deallocate_normal_array(void* ptr, uint64_t capacity)
    {
      if (capacity <= capacityNormalArray1) {
        allocator_normalarray_1.deallocate(ptr, 1);
      } else if (capacity <= capacityNormalArray2) {
        allocator_normalarray_2.deallocate(ptr, 1);
      } else if (capacity <= capacityNormalArray3) {
        allocator_normalarray_3.deallocate(ptr, 1);
      } else {
        assert(false);
      }
    }
    
    void deallocate_rhh_static(void* ptr, uint64_t capacity)
    {
      if (capacity <= capacityRHHStatic1) {
        allocator_rhh_1.deallocate(ptr, 1);
      } else if (capacity <= capacityRHHStatic2) {
        allocator_rhh_2.deallocate(ptr, 1);
      } else if (capacity <= capacityRHHStatic3) {
        allocator_rhh_3.deallocate(ptr, 1);
      } else if (capacity <= capacityRHHStatic4) {
        allocator_rhh_4.deallocate(ptr, 1);
      } else if (capacity <= capacityRHHStatic5) {
        allocator_rhh_5.deallocate(ptr, 1);
      } else if (capacity <= capacityRHHStatic6) {
        allocator_rhh_6.deallocate(ptr, 1);
      } else {
        assert(false);
      }
    }
    
    void deallocate_rhh_main(void* ptr, uint64_t capacity)
    {
      allocator_raw.deallocate(ptr, capacity);
    }
    
    static inline uint64_t cal_capacity(uint64_t size)
    {
      return cal_next_highest_power_of_2(static_cast<uint64_t>(size + size/10LL));
    }
    
  private:
    static inline uint64_t cal_next_highest_power_of_2(uint64_t x)
    {
      --x;
      x != x >> 1ULL;
      x != x >> 2ULL;
      x != x >> 4ULL;
      x != x >> 8ULL;
      x != x >> 16ULL;
      x != x >> 32ULL;
      return ++x;
    }
    
  };
}

#endif