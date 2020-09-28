// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_RHH_DETAIL_UTILITY_PACKED_TUPLE_HPP
#define HAVOQGT_RHH_DETAIL_UTILITY_PACKED_TUPLE_HPP

#include <iostream>

namespace rhh {
namespace detail {
namespace utility {

#pragma pack(1)
template <typename T1, typename T2>
struct packed_pair {
  T1 first;
  T2 second;

  packed_pair() = default;

  packed_pair(const T1 &t1, const T2 &t2) :
      first(t1),
      second(t2) {}

  packed_pair(T1 &&t1, T2 &&t2) :
      first(std::move(t1)),
      second(std::move(t2)) {}

  packed_pair(const packed_pair &) = default;
  packed_pair(packed_pair &&) noexcept = default;
  packed_pair &operator=(const packed_pair &) = default;
  packed_pair &operator=(packed_pair &&) noexcept = default;

  void swap(packed_pair<T1, T2> &other) {
    using std::swap;
    swap(first, other.first);
    swap(second, other.second);
  }

};

template <typename T1, typename T2>
inline bool operator==(const packed_pair<T1, T2> &lhs, const packed_pair<T1, T2> &rhs) {
  return lhs.first == rhs.first && lhs.second == rhs.second;
}
template <typename T1, typename T2>
inline bool operator!=(const packed_pair<T1, T2> &lhs, const packed_pair<T1, T2> &rhs) {
  return !(lhs == rhs);
}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &stream, const packed_pair<T1, T2> &obj) {
  stream << obj.first << ":" << obj.second;
  return stream;
}

#pragma pack(1)
template <typename T1, typename T2, typename T3, typename T4>
struct packed_tuple {
  T1 first;
  T2 second;
  T3 third;
  T4 fourth;

  using self_type = packed_tuple<T1, T2, T3, T4>;

  /// ---- Constructors ----
  packed_tuple() = default;

  packed_tuple(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4) :
      first(t1),
      second(t2),
      third(t3),
      fourth(t4) {}

  packed_tuple(T1 &&t1, T2 &&t2, T3 &&t3, T4 &&t4) :
      first(std::move(t1)),
      second(std::move(t2)),
      third(std::move(t3)),
      fourth(std::move(t4)) {}

  packed_tuple(const packed_tuple &) = default;
  packed_tuple(packed_tuple &&) noexcept = default;
  packed_tuple &operator=(const packed_tuple &) = default;
  packed_tuple &operator=(packed_tuple &&) noexcept = default;

  void swap(self_type &other) {
    using std::swap;
    swap(first, other.first);
    swap(second, other.second);
    swap(third, other.third);
    swap(fourth, other.fourth);
  }

};

template <typename T1, typename T2, typename T3, typename T4>
inline bool operator==(const packed_tuple<T1, T2, T3, T4> &lhs, const packed_tuple<T1, T2, T3, T4> &rhs) {
  return lhs.first == rhs.first && lhs.second == rhs.second &&
      lhs.third == rhs.third && lhs.fourth == rhs.fourth;
}

template <typename T1, typename T2, typename T3, typename T4>
inline bool operator!=(const packed_tuple<T1, T2, T3, T4> &lhs, const packed_tuple<T1, T2, T3, T4> &rhs) {
  return !(lhs == rhs);
}

template <typename T1, typename T2, typename T3, typename T4>
std::ostream &operator<<(std::ostream &stream, const packed_tuple<T1, T2, T3, T4> &obj) {
  stream << obj.first << ":" << obj.second << ":" << obj.third << ":" << obj.fourth;
  return stream;
}

#pragma pack(1)
template <typename T1, typename T2, typename T3>
struct packed_tuple<T1, T2, T3, void> {
  T1 first;
  T2 second;
  T3 third;

  using self_type = packed_tuple<T1, T2, T3, void>;

  /// ---- Constructors ----
  packed_tuple() {}

  packed_tuple(const T1 &t1, const T2 &t2, const T3 &t3) :
      first(t1),
      second(t2),
      third(t3) {}

  packed_tuple(T1 &&t1, T2 &&t2, T3 &&t3) :
      first(std::move(t1)),
      second(std::move(t2)),
      third(std::move(t3)) {}

  packed_tuple(const packed_tuple &) = default;
  packed_tuple(packed_tuple &&) noexcept = default;
  packed_tuple &operator=(const packed_tuple &) = default;
  packed_tuple &operator=(packed_tuple &&) noexcept = default;

  void swap(self_type &other) {
    using std::swap;
    swap(first, other.first);
    swap(second, other.second);
    swap(third, other.third);
  }

};

template <typename T1, typename T2, typename T3>
inline bool operator==(const packed_tuple<T1, T2, T3, void> &lhs, const packed_tuple<T1, T2, T3, void> &rhs) {
  return lhs.first == rhs.first && lhs.second == rhs.second &&
      lhs.third == rhs.third;
}

template <typename T1, typename T2, typename T3>
inline bool operator!=(const packed_tuple<T1, T2, T3, void> &lhs, const packed_tuple<T1, T2, T3, void> &rhs) {
  return !(lhs == rhs);
}

template <typename T1, typename T2, typename T3>
std::ostream &operator<<(std::ostream &stream, const packed_tuple<T1, T2, T3, void> &obj) {
  stream << obj.first << ":" << obj.second << ":" << obj.third;
  return stream;
}

} // namespace utility
} // namespace detail
} // namespace rhh
#endif /// HAVOQGT_RHH_DETAIL_UTILITY_PACKED_TUPLE_HPP
