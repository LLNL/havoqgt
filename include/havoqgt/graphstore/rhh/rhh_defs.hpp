/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHHDA_DEFS_HPP
#define RHHDA_DEFS_HPP

#include <havoqgt/graphstore/graphstore_utilities.hpp>

namespace graphstore {
namespace rhh {

#define RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE 1
#define RHH_CHAIN_AT_LARGE_TABLE_SIZE 0
#define RHH_DETAILED_ANALYSYS 0 /// !!! not supporting growing

static constexpr double kFullCapacitFactor  = 0.9;

enum {
    kNodeAllocatorChunkSize = (1 << 24),
    kCapacityGrowingFactor = 2
};

void disp_configuration(void)
{
  DISP_LOG_VAR(RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE);
  DISP_LOG_VAR(RHH_CHAIN_AT_LARGE_TABLE_SIZE);
  DISP_LOG_VAR(RHH_DETAILED_ANALYSYS);
  DISP_LOG_VAR(kFullCapacitFactor);
  DISP_LOG_VAR(kNodeAllocatorChunkSize);
}

}}

#endif // RHHDA_DEFS_HPP

