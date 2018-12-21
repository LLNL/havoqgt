/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHHDA_DEFS_HPP
#define RHHDA_DEFS_HPP

#include <dynamic_graph_store/graphstore_utilities.hpp>

namespace graphstore {
namespace rhh {

#define RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE 0
#define RHH_CHAIN_AT_LARGE_TABLE_SIZE 0
#define RHH_USE_NUMA_ALLOC 0 /// use numa alloc instead of malloc
#define RHH_DETAILED_ANALYSYS 0 /// !!! will be invalid when the table grow !!!

static constexpr double kFullCapacitFactor  = 0.9;

enum {
    kNodeAllocatorChunkSize = (1 << 24),
    kCapacityGrowingFactor = 2
};

static void disp_configuration(void)
{
  DISP_LOG_VAR(RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE);
  DISP_LOG_VAR(RHH_CHAIN_AT_LARGE_TABLE_SIZE);
  DISP_LOG_VAR(RHH_USE_NUMA_ALLOC);
  //  DISP_LOG_VAR(RHH_DETAILED_ANALYSYS);
  DISP_LOG_VAR(kFullCapacitFactor);
  DISP_LOG_VAR(kNodeAllocatorChunkSize);
}

}}

#endif // RHHDA_DEFS_HPP

