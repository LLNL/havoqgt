// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_INCLUDE_RHH_RHH_MAP_HPP
#define HAVOQGT_INCLUDE_RHH_RHH_MAP_HPP

#include <rhh/detail/basic_rhh.hpp>
#include <rhh/hash.hpp>

namespace rhh {
/// \brief Robin Hood Hashing class that takes a key and a value
/// \tparam key_type Type of the key values
/// \tparam mapped_value_type Type of the mapped values
/// \tparam hash Hash function
/// \tparam EqualTo Equal to function
/// \tparam AllocatorType Type of the allocator
template <typename key_type, typename mapped_value_type,
          typename hash      = hash<key_type>,
          typename key_equal = std::equal_to<key_type>,
          typename allocator = std::allocator<std::byte>>
using rhh_map = detail::basic_rhh<
    detail::utility::pair_value_type<key_type, mapped_value_type>, hash,
    key_equal, allocator>;
}  // namespace rhh
#endif  // HAVOQGT_INCLUDE_RHH_RHH_MAP_HPP
