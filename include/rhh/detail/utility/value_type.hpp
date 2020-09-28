// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_RHH_DYNAMI_VALUE_TYPE_HPP
#define HAVOQGT_RHH_DYNAMI_VALUE_TYPE_HPP

#include <memory>
#include <tuple>
#include <utility>

namespace rhh {
namespace detail {
namespace utility {

/// \brief Generalized class for the cases mapped_value_type is not
/// void_mapped_value_tag value_type is std::pair<key_type, mapped_value_type>
/// \tparam key_type Type of the key
/// \tparam mapped_value_type Type of the value
template <typename _key_type, typename _mapped_value_type>
class pair_value_type {
 public:
  using key_type             = _key_type;
  using value_type           = std::tuple<key_type, _mapped_value_type>;
  using const_key_value_type = std::tuple<const key_type, _mapped_value_type>;

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   const value_type &   value) {
    return value_type(std::allocator_arg, allocator, value);
  }

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   value_type &&        value) {
    return value_type(std::allocator_arg, allocator, std::move(value));
  }

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   const key_type &     key) {
    value_type value(std::allocator_arg, allocator);
    std::get<0>(value) = key;
    return value;
  }

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   key_type &&          key) {
    value_type value(std::allocator_arg, allocator);
    std::get<0>(value) = std::move(key);
    return value;
  }

  constexpr static std::size_t size() { return 2; }

  static key_type &key(value_type &value) { return std::get<0>(value); }
  static key_type &key(key_type &key) { return key; }

  static const key_type &key(const value_type &value) {
    return std::get<0>(value);
  }
  static const key_type &key(const key_type &key) { return key; }
};

/// \brief Specialized class whose second template type is void_mapped_value_tag
/// therefore, value_type is same to key_type
/// \tparam key_type Type of the key
template <typename _key_type>
class key_only_value_type {
 public:
  using key_type             = _key_type;
  using value_type           = key_type;
  using const_key_value_type = const key_type;

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   const key_type &     key) {
    return value_type(key);
  }

  template <typename allocator_type>
  static value_type allocate_value(const allocator_type allocator,
                                   key_type &&          key) {
    return value_type(std::move(key));
  }

  constexpr static std::size_t size() { return 1; }

  static key_type &      key(value_type &value) { return value; }
  static const key_type &key(const value_type &value) { return value; }
};

}  // namespace utility
}  // namespace detail
}  // namespace rhh
#endif  // HAVOQGT_RHH_DYNAMI_VALUE_TYPE_HPP
