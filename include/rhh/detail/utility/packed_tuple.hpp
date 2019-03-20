/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

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
