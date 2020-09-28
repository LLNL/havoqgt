// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_RHH_DETAIL_UTILITY_HASH_HPP
#define HAVOQGT_RHH_DETAIL_UTILITY_HASH_HPP

namespace rhh {
namespace detail {
namespace utility {

// Forward declaration
uint64_t MurmurHash64A(const void *key, int len, uint64_t seed);

}  // namespace utility
}  // namespace detail
} // namespace rhh

namespace rhh {
template <class Key>
class hash {
 public:
  uint64_t operator()(Key const &key, uint64_t const &seed = 123) const {
    return detail::utility::MurmurHash64A(&key, sizeof(Key), seed);
  }
};
} // namespace rhh

namespace rhh {
namespace detail {
namespace utility {

// -----------------------------------------------------------------------------
// This file contains public domain code from MurmurHash2.
// From the MurmurHash2 header:
// -----------------------------------------------------------------------------
//  MurmurHash2, 64-bit versions, by Austin Appleby
//  MurmurHash2 was written by Austin Appleby, and is placed in the public
//  domain. The author hereby disclaims copyright to this source code.
// -----------------------------------------------------------------------------

uint64_t MurmurHash64A(const void *key, int len, uint64_t seed) {
  const uint64_t m = 0xc6a4a7935bd1e995ULL;
  const int      r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t *data = (const uint64_t *)key;
  const uint64_t *end  = data + (len / 8);

  while (data != end) {
    uint64_t k = *data++;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char *data2 = (const unsigned char *)data;

  switch (len & 7) {
    case 7:
      h ^= uint64_t(data2[6]) << 48;
    case 6:
      h ^= uint64_t(data2[5]) << 40;
    case 5:
      h ^= uint64_t(data2[4]) << 32;
    case 4:
      h ^= uint64_t(data2[3]) << 24;
    case 3:
      h ^= uint64_t(data2[2]) << 16;
    case 2:
      h ^= uint64_t(data2[1]) << 8;
    case 1:
      h ^= uint64_t(data2[0]);
      h *= m;
  };

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}

}  // namespace utility
}  // namespace detail
}  // namespace rhh

#endif  // HAVOQGT_RHH_DETAIL_UTILITY_HASH_HPP
