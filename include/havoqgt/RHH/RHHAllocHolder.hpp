#ifndef RHHH_RHHALLOCHOLDER_HPP_INCLUDED
#define RHHH_RHHALLOCHOLDER_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include "RHHCommon.hpp"
#include "RHHStatic.hpp"

namespace RHH {

  namespace bip = boost::interprocess;


  enum {
    capacityNormalArray1 = 2ULL,
    capacityNormalArray2 = 4ULL,
    capacityNormalArray3 = 8ULL,
    capacityRHHStatic1 = 8ULL,
    capacityRHHStatic2 = 16ULL,
    capacityRHHStatic3 = 32ULL,
    capacityRHHStatic4 = 64ULL,
    capacityRHHStatic5 = 128ULL,
    capacityRHHStatic6 = 256ULL
  };

  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic1> RHHStaticNoVal1;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic2> RHHStaticNoVal2;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic3> RHHStaticNoVal3;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic4> RHHStaticNoVal4;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic5> RHHStaticNoVal5;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic6> RHHStaticNoVal6;


  ///  =========================================================================== ///
  ///                             Allocator Holder
  ///  =========================================================================== ///
  class AllocatorsHolder
  {
  public:
    typedef bip::managed_mapped_file  mapped_t;
    typedef mapped_t::segment_manager segment_manager_t;
    explicit AllocatorsHolder(segment_manager_t* segment_manager)
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

    // void* allocate_rhh_main(const uint64_t required_capacity)
    // {
    //   const uint64_t memsize = RHHMain<KeyType, ValueType>.cal_memsize(required_capacity);
    //   return reinterpret_cast<void*>(allocator_raw.allocate(memsize).get());
    // }

    void deallocate_normal_array(void* ptr, uint64_t capacity)
    {
      // if (capacity <= capacityNormalArray1) {
      //   allocator_normalarray_1.deallocate(bip::offset_ptr<void>(ptr), 1);
      // } else if (capacity <= capacityNormalArray2) {
      //   allocator_normalarray_2.deallocate(bip::offset_ptr<void>(ptr), 1);
      // } else if (capacity <= capacityNormalArray3) {
      //   allocator_normalarray_3.deallocate(bip::offset_ptr<void>(ptr), 1);
      // } else {
      //   assert(false);
      // }
    }

    void deallocate_rhh_static(void* ptr, uint64_t capacity)
    {
      if (capacity == capacityRHHStatic1) {
        allocator_rhh_noval_1.deallocate(bip::offset_ptr<RHHStaticNoVal1>(reinterpret_cast<RHHStaticNoVal1*>(ptr)), 1);
      } else if (capacity == capacityRHHStatic2) {
        allocator_rhh_noval_2.deallocate(bip::offset_ptr<RHHStaticNoVal2>(reinterpret_cast<RHHStaticNoVal2*>(ptr)), 1);
      } else if (capacity == capacityRHHStatic3) {
        allocator_rhh_noval_3.deallocate(bip::offset_ptr<RHHStaticNoVal3>(reinterpret_cast<RHHStaticNoVal3*>(ptr)), 1);
      } else if (capacity == capacityRHHStatic4) {
        allocator_rhh_noval_4.deallocate(bip::offset_ptr<RHHStaticNoVal4>(reinterpret_cast<RHHStaticNoVal4*>(ptr)), 1);
      } else if (capacity == capacityRHHStatic5) {
        allocator_rhh_noval_5.deallocate(bip::offset_ptr<RHHStaticNoVal5>(reinterpret_cast<RHHStaticNoVal5*>(ptr)), 1);
      } else if (capacity == capacityRHHStatic6) {
        allocator_rhh_noval_6.deallocate(bip::offset_ptr<RHHStaticNoVal6>(reinterpret_cast<RHHStaticNoVal6*>(ptr)), 1);
      } else {
        assert(false);
      }
    }

    void deallocate_rhh_main(void* ptr, uint64_t capacity)
    {
      //allocator_raw.deallocate(bip::offset_ptr<void>(ptr), capacity);
    }

    static inline uint64_t cal_capacity(uint64_t size)
    {
      return cal_next_highest_power_of_2(static_cast<uint64_t>(size + size/10LL));
    }


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