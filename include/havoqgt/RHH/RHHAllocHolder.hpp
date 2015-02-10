#ifndef RHHH_RHHALLOCHOLDER_HPP_INCLUDED
#define RHHH_RHHALLOCHOLDER_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include "RHHCommon.hpp"
#include "RHHStatic.hpp"

namespace RHH {

  namespace bip = boost::interprocess;
  typedef bip::managed_mapped_file  mapped_t;
  typedef mapped_t::segment_manager segment_manager_t;

  /// Memory size (Byte) = capacity * (  1  /// property block
  ///                                  + 8  /// key block
  ///                                  + 1) /// value block
  ///                                  + 8   /// next pointer
  enum {
    capacityNormalArray1 = 2ULL,
    capacityNormalArray2 = 4ULL,
    capacityNormalArray3 = 8ULL,
    capacityRHHStatic_1 = 8ULL,
    capacityRHHStatic_2 = 16ULL,
    capacityRHHStatic_3 = 32ULL,
    capacityRHHStatic_4 = 64ULL,
    capacityRHHStatic_5 = 128ULL,
    capacityRHHStatic_6 = 256ULL,
    capacityRHHStatic_7 = 512ULL,
    capacityRHHStatic_8 = 1024ULL,
    capacityRHHStatic_9 = 2048ULL,
    capacityRHHStatic_10 = 4096ULL,
    capacityRHHStatic_11 = 8192ULL,
    capacityRHHStatic_12 = 16384ULL,
    capacityRHHStatic_13 = 32768ULL,
    capacityRHHStatic_14 = 65536ULL,
    capacityRHHStatic_15 = 131072ULL,
    capacityRHHStatic_16 = 262144ULL,
    capacityRHHStatic_17 = 524288ULL
  };

  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_1> RHHStaticNoVal_1;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_2> RHHStaticNoVal_2;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_3> RHHStaticNoVal_3;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_4> RHHStaticNoVal_4;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_5> RHHStaticNoVal_5;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_6> RHHStaticNoVal_6;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_7> RHHStaticNoVal_7;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_8> RHHStaticNoVal_8;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_9> RHHStaticNoVal_9;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_10> RHHStaticNoVal_10;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_11> RHHStaticNoVal_11;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_12> RHHStaticNoVal_12;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_13> RHHStaticNoVal_13;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_14> RHHStaticNoVal_14;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_15> RHHStaticNoVal_15;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_16> RHHStaticNoVal_16;
  typedef RHHStatic<uint64_t, NoValueType, capacityRHHStatic_17> RHHStaticNoVal_17;

#if 0
  /// chuk size = 4MB
  typedef bip::node_allocator<RHHStaticNoVal_1, segment_manager_t,  47662>  Allocator_RHH_NoVal_1;
  typedef bip::node_allocator<RHHStaticNoVal_2, segment_manager_t,  24966>  Allocator_RHH_NoVal_2;
  typedef bip::node_allocator<RHHStaticNoVal_3, segment_manager_t,  12787>  Allocator_RHH_NoVal_3;
  typedef bip::node_allocator<RHHStaticNoVal_4, segment_manager_t,   6472>  Allocator_RHH_NoVal_4;
  typedef bip::node_allocator<RHHStaticNoVal_5, segment_manager_t,   3256>  Allocator_RHH_NoVal_5;
  typedef bip::node_allocator<RHHStaticNoVal_6, segment_manager_t,   1633>  Allocator_RHH_NoVal_6;
  typedef bip::node_allocator<RHHStaticNoVal_7, segment_manager_t,   817>  Allocator_RHH_NoVal_7;
  typedef bip::node_allocator<RHHStaticNoVal_8, segment_manager_t,   409>  Allocator_RHH_NoVal_8;
  typedef bip::node_allocator<RHHStaticNoVal_9, segment_manager_t,   204>  Allocator_RHH_NoVal_9;
  typedef bip::node_allocator<RHHStaticNoVal_10, segment_manager_t,  102>  Allocator_RHH_NoVal_10;
  typedef bip::node_allocator<RHHStaticNoVal_11, segment_manager_t,  51>  Allocator_RHH_NoVal_11;
  typedef bip::node_allocator<RHHStaticNoVal_12, segment_manager_t,  25>  Allocator_RHH_NoVal_12;
  typedef bip::node_allocator<RHHStaticNoVal_13, segment_manager_t,  12>  Allocator_RHH_NoVal_13;
  typedef bip::node_allocator<RHHStaticNoVal_14, segment_manager_t,  6>  Allocator_RHH_NoVal_14;
  typedef bip::node_allocator<RHHStaticNoVal_15, segment_manager_t,  3>  Allocator_RHH_NoVal_15;
  typedef bip::node_allocator<RHHStaticNoVal_16, segment_manager_t,  2>  Allocator_RHH_NoVal_16;
  typedef bip::node_allocator<RHHStaticNoVal_17, segment_manager_t,  1>  Allocator_RHH_NoVal_17;
  typedef bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_normalarray_t;
  typedef bip::allocator<unsigned char, segment_manager_t> allocator_raw_t;
#else
  /// chuk size = 4KB
  typedef bip::node_allocator<RHHStaticNoVal_1, segment_manager_t,  47>  Allocator_RHH_NoVal_1;
  typedef bip::node_allocator<RHHStaticNoVal_2, segment_manager_t,  24>  Allocator_RHH_NoVal_2;
  typedef bip::node_allocator<RHHStaticNoVal_3, segment_manager_t,  12>  Allocator_RHH_NoVal_3;
  typedef bip::node_allocator<RHHStaticNoVal_4, segment_manager_t,   6>  Allocator_RHH_NoVal_4;
  typedef bip::node_allocator<RHHStaticNoVal_5, segment_manager_t,   3>  Allocator_RHH_NoVal_5;
  typedef bip::node_allocator<RHHStaticNoVal_6, segment_manager_t,   1>  Allocator_RHH_NoVal_6;
  typedef bip::node_allocator<RHHStaticNoVal_7, segment_manager_t,   1>  Allocator_RHH_NoVal_7;
  typedef bip::node_allocator<RHHStaticNoVal_8, segment_manager_t,   1>  Allocator_RHH_NoVal_8;
  typedef bip::node_allocator<RHHStaticNoVal_9, segment_manager_t,   1>  Allocator_RHH_NoVal_9;
  typedef bip::node_allocator<RHHStaticNoVal_10, segment_manager_t,  1>  Allocator_RHH_NoVal_10;
  typedef bip::node_allocator<RHHStaticNoVal_11, segment_manager_t,  1>  Allocator_RHH_NoVal_11;
  typedef bip::node_allocator<RHHStaticNoVal_12, segment_manager_t,  1>  Allocator_RHH_NoVal_12;
  typedef bip::node_allocator<RHHStaticNoVal_13, segment_manager_t,  1>  Allocator_RHH_NoVal_13;
  typedef bip::node_allocator<RHHStaticNoVal_14, segment_manager_t,  1>  Allocator_RHH_NoVal_14;
  typedef bip::node_allocator<RHHStaticNoVal_15, segment_manager_t,  1>  Allocator_RHH_NoVal_15;
  typedef bip::node_allocator<RHHStaticNoVal_16, segment_manager_t,  1>  Allocator_RHH_NoVal_16;
  typedef bip::node_allocator<RHHStaticNoVal_17, segment_manager_t,  1>  Allocator_RHH_NoVal_17;
  typedef bip::node_allocator<uint64_t, segment_manager_t, 512> allocator_normalarray_t;
  typedef bip::allocator<unsigned char, segment_manager_t> allocator_raw_t;
#endif
  ///  =========================================================================== ///
  ///                             Allocator Holder
  ///  =========================================================================== ///
  class AllocatorsHolder
  {
  public:
    explicit AllocatorsHolder(segment_manager_t* segment_manager)
    : allocator_rhh_noval_1(segment_manager)
    , allocator_rhh_noval_2(segment_manager)
    , allocator_rhh_noval_3(segment_manager)
    , allocator_rhh_noval_4(segment_manager)
    , allocator_rhh_noval_5(segment_manager)
    , allocator_rhh_noval_6(segment_manager)
    , allocator_rhh_noval_7(segment_manager)
    , allocator_rhh_noval_8(segment_manager)
    , allocator_rhh_noval_9(segment_manager)
    , allocator_rhh_noval_10(segment_manager)
    , allocator_rhh_noval_11(segment_manager)
    , allocator_rhh_noval_12(segment_manager)
    , allocator_rhh_noval_13(segment_manager)
    , allocator_rhh_noval_14(segment_manager)
    , allocator_rhh_noval_15(segment_manager)
    , allocator_rhh_noval_16(segment_manager)
    , allocator_rhh_noval_17(segment_manager)
    , allocator_normalarray(segment_manager)
    , allocator_raw(segment_manager)
    { }


  /// size = capacity * (1 + 8 + 1) + 8
  /// probedistance = 1 byte, key = 8 byte, value block = 1 byte
  Allocator_RHH_NoVal_1 allocator_rhh_noval_1;
  Allocator_RHH_NoVal_2 allocator_rhh_noval_2;
  Allocator_RHH_NoVal_3 allocator_rhh_noval_3;
  Allocator_RHH_NoVal_4 allocator_rhh_noval_4;
  Allocator_RHH_NoVal_5 allocator_rhh_noval_5;
  Allocator_RHH_NoVal_6 allocator_rhh_noval_6;
  Allocator_RHH_NoVal_7 allocator_rhh_noval_7;
  Allocator_RHH_NoVal_8 allocator_rhh_noval_8;
  Allocator_RHH_NoVal_9 allocator_rhh_noval_9;
  Allocator_RHH_NoVal_10 allocator_rhh_noval_10;
  Allocator_RHH_NoVal_11 allocator_rhh_noval_11;
  Allocator_RHH_NoVal_12 allocator_rhh_noval_12;
  Allocator_RHH_NoVal_13 allocator_rhh_noval_13;
  Allocator_RHH_NoVal_14 allocator_rhh_noval_14;
  Allocator_RHH_NoVal_15 allocator_rhh_noval_15;
  Allocator_RHH_NoVal_16 allocator_rhh_noval_16;
  Allocator_RHH_NoVal_17 allocator_rhh_noval_17;
  allocator_normalarray_t allocator_normalarray;
  allocator_raw_t allocator_raw;

  private:
    // static inline uint64_t cal_next_highest_power_of_2(uint64_t x)
    // {
    //   --x;
    //   x != x >> 1ULL;
    //   x != x >> 2ULL;
    //   x != x >> 4ULL;
    //   x != x >> 8ULL;
    //   x != x >> 16ULL;
    //   x != x >> 32ULL;
    //   return ++x;
    // }

  };
}

#endif