/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHHDA_DEFS_HPP
#define RHHDA_DEFS_HPP

#include <havoqgt/graphstore/graphstore_utilities.hpp>

namespace graphstore {
namespace rhh {

static constexpr double kFullCapacitFactor  = 0.9;

enum {
    kNodeAllocatorChunkSize = (1 << 24),
    kCapacityGrowingFactor = 2
};

void disp_configuration(void)
{
  DISP_LOG_VAR(kNodeAllocatorChunkSize);
}

}}

#endif // RHHDA_DEFS_HPP

