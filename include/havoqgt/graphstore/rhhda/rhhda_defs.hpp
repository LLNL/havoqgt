/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHHDA_DEFS_HPP
#define RHHDA_DEFS_HPP

#include <havoqgt/graphstore/graphstore_utilities.hpp>

namespace graphstore {
namespace rhhda {

/// If hash vertex ID
#define HASH_VERTEX_ID 0
static constexpr double kFullCapacitFactor  = 0.9;

enum {
    kNodeAllocatorChunkSize = (1 << 24),
    kCapacityGrwoingFactor = 2
};

void disp_configuration(void)
{
  DISP_LOG_VAR(HASH_VERTEX_ID);
  DISP_LOG_VAR(kNodeAllocatorChunkSize);
}

}}

#endif // RHHDA_DEFS_HPP

