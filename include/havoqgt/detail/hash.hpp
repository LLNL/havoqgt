// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_DETAIL_HASH_HPP_INCLUDED
#define HAVOQGT_DETAIL_HASH_HPP_INCLUDED

#include <assert.h>
#include <stdint.h>

namespace havoqgt {
namespace detail {
///
/// Hash functions
///
/// \todo requires documentation!
/// \todo requires testing!

inline uint32_t hash32(uint32_t a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

inline uint16_t hash16(uint16_t a) {
  a = (a + 0x5d16) + (a << 6);
  a = (a ^ 0xc23c) ^ (a >> 9);
  a = (a + 0x67b1) + (a << 5);
  a = (a + 0x646c) ^ (a << 7);
  a = (a + 0x46c5) + (a << 3);
  a = (a ^ 0x4f09) ^ (a >> 8);
  return a;
}

inline uint64_t shifted_n_hash32(uint64_t input, int n) {
  uint64_t to_hash = input >> n;
  uint64_t mask    = 0xFFFFFFFF;
  to_hash &= mask;
  to_hash = hash32(to_hash);

  to_hash <<= n;
  mask <<= n;
  // clear bits
  input &= ~mask;
  input |= to_hash;
  return input;
}

inline uint64_t shifted_n_hash16(uint64_t input, int n) {
  uint64_t to_hash = input >> n;
  uint64_t mask    = 0xFFFF;
  to_hash &= mask;
  to_hash = hash16(to_hash);

  to_hash <<= n;
  mask <<= n;
  // clear bits
  input &= ~mask;
  input |= to_hash;
  return input;
}

inline uint64_t hash_nbits(uint64_t input, int n) {
  // std::cout << "hash_nbits(" << input << ", " << n << ") = ";
  if (n == 32) {
    input = hash32(input);
  } else if (n > 32) {
    assert(n > 32);
    n -= 32;
    for (int i = 0; i <= n; ++i) {
      input = shifted_n_hash32(input, i);
    }
    for (int i = n; i >= 0; --i) {
      input = shifted_n_hash32(input, i);
    }
  } else if (n < 32) {
    assert(n < 32);
    assert(n > 16 && "Hashing less than 16bits is not supported");
    n -= 16;
    for (int i = 0; i <= n; ++i) {
      input = shifted_n_hash16(input, i);
    }
    for (int i = n; i >= 0; --i) {
      input = shifted_n_hash16(input, i);
    }
  }
  // std::cout << input << std::endl;
  return input;
}
}  // namespace detail
}  // namespace havoqgt

#endif
