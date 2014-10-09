#ifndef RHH_RHHCOMMON_HPP_INCLUDED
#define RHH_RHHCOMMON_HPP_INCLUDED

namespace RHH {

#define DEBUG(msg) do { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << msg << std::endl; } while (0)
#define DEBUG2(x) do  { std::cerr << "DEG: " << __FILE__ << "(" << __LINE__ << ") " << #x << " =\t" << static_cast<uint64_t>(x) << std::endl; } while (0)
#define DISP_VAR(x) do  { std::cout << #x << " =\t" << x << std::endl; } while (0)

  enum UpdateErrors {
    kSucceed,
    kDuplicated,
    kReachingFUllCapacity,
    kLongProbedistance
  };

  static const uint64_t kCapacityGrowingFactor = 2ULL;
  typedef unsigned char NoValueType;
}

#endif