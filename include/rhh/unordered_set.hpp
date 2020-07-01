/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License. http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by the United States Government or Lawrence Livermore National Security, LLC.
 * The views and opinions of authors expressed herein do not necessarily state
 * or reflect those of the United States Government or Lawrence Livermore
 * National Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

#ifndef HAVOQGT_RHH_UNORDERED_SET_HPP
#define HAVOQGT_RHH_UNORDERED_SET_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <rhh/hash.hpp>
#include <rhh/rhh_set.hpp>

namespace rhh {

/// \brief Robin Hood hashing class
/// Following the interface of std::map
/// \tparam key_type Type of the key values
/// \tparam hash hash function
/// \tparam key_equal Equal to function
/// \tparam allocator Type of the allocator
template <typename _key_type, typename _hash = rhh::hash<_key_type>,
          typename _key_equal = std::equal_to<_key_type>,
          typename _allocator = std::allocator<_key_type>>
class unordered_set {
 private:
  using self_type = unordered_set<_key_type, _hash, _key_equal, _allocator>;
  using rhh_table_type = rhh_set<_key_type, _hash, _key_equal, _allocator>;

  template <bool is_const>
  class iterator_impl;

 public:
  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Public type
  // ----------------------------------------------------------------------------------------------------
  // //
  using key_type        = _key_type;
  using value_type      = key_type;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using hasher          = _hash;
  using key_equal       = _key_equal;
  using allocator_type  = _allocator;
  using reference       = value_type &;
  using const_reference = const value_type &;
  using pointer         = typename std::allocator_traits<_allocator>::pointer;
  using const_pointer =
      typename std::allocator_traits<_allocator>::const_pointer;
  using iterator       = iterator_impl<false>;
  using const_iterator = iterator_impl<true>;

  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Public function
  // ----------------------------------------------------------------------------------------------------
  // //

  // -------------------- Constructors -------------------- //
  explicit unordered_set(const allocator_type &allocator = allocator_type())
      : m_table(allocator) {}
  ~unordered_set() = default;

  unordered_set(const unordered_set &)     = default;
  unordered_set(unordered_set &&) noexcept = default;

  unordered_set &operator=(const unordered_set &) = default;
  unordered_set &operator=(unordered_set &&) noexcept = default;

  // -------------------- Capacity -------------------- //
  size_type size() const { return m_table.size(); }

  bool empty() const { return (m_table.size() == 0); }

  size_type max_size() const {
    return std::numeric_limits<difference_type>::max();
  }

  // -------------------- Iterator -------------------- //
  iterator begin() { return iterator(&m_table, 0); }

  const_iterator begin() const { return const_iterator(&m_table, 0); }

  iterator end() { return iterator(&m_table, m_table.capacity()); }

  const_iterator end() const {
    return const_iterator(&m_table, m_table.capacity());
  }

  const_iterator cbegin() const { return const_iterator(&m_table, 0); }

  const_iterator cend() const {
    return const_iterator(&m_table, m_table.capacity());
  }

  // -------------------- Modifiers -------------------- //
  void clear() { m_table.clear(); }

  // TODO: implement other versions
  std::pair<iterator, bool> insert(const value_type &value) {
    const size_type pos = m_table.find(value);
    if (pos != m_table.capacity()) {
      return std::make_pair(iterator(&m_table, pos), false);
    }

    return std::make_pair(iterator(&m_table, m_table.insert(value)), true);
  }

  std::pair<iterator, bool> insert(value_type &&value) {
    const size_type pos = m_table.find(value);
    if (pos != m_table.capacity()) {
      return std::make_pair(iterator(&m_table, pos), false);
    }

    return std::make_pair(iterator(&m_table, m_table.insert(std::move(value))),
                          true);
  }

  // TODO: implement emplace, emplace_hint

  iterator erase(const_iterator position) {
    assert(position.m_table == &m_table);
    m_table.erase_at(position.m_current_position);
    return iterator(&m_table, position.m_current_position);
  }

  iterator erase(const_iterator first, const_iterator last) {
    assert(first.m_table == &m_table);
    assert(last.m_table == &m_table);
    assert(first != last);
    for (; first != last; ++first) {
      m_table.erase_at(first.m_current_position);
    }
    return iterator(&m_table, last.m_current_position);
  }

  size_type erase(const key_type &key) {
    const size_type num_erased = m_table.erase(key);
    assert(num_erased <= 1);
    return num_erased;
  }

  void swap(self_type &other) {
    using std::swap;
    swap(m_table, other.m_table);
  }

  // -------------------- Element lookup -------------------- //
  size_type count(const key_type &key) const {
    const size_type num_elements =
        (m_table.find(key) == m_table.capacity()) ? 0 : 1;
    return num_elements;
  }

  iterator find(const key_type &key) {
    return iterator(&m_table, m_table.find(key));
  }

  const_iterator find(const key_type &key) const {
    return const_iterator(&m_table, m_table.find(key));
  }

  // -------------------- Statistic -------------------- //
  double load_factor() const { return m_table.load_factor(); }

 private:
  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Private variable
  // ----------------------------------------------------------------------------------------------------
  // //
  rhh_table_type m_table;
};

template <typename _key_type, typename _hash, typename _key_equal,
          typename _allocator>
template <bool is_const>
class unordered_set<_key_type, _hash, _key_equal, _allocator>::iterator_impl {
 public:
  using value_type =
      typename std::conditional<is_const, const _key_type, _key_type>::type;
  using pointer           = value_type *;
  using reference         = value_type &;
  using iterator_category = std::forward_iterator_tag;

  friend unordered_set<_key_type, _hash, _key_equal, _allocator>;

  iterator_impl() = default;

  iterator_impl(rhh_table_type *const table, const size_type position)
      : m_table(table), m_current_position(position) {
    m_current_position = m_table->find_next(m_current_position);
  }

  iterator_impl(const iterator_impl &) = default;
  iterator_impl(iterator_impl &&)      = default;
  iterator_impl &operator=(const iterator_impl &) = default;
  iterator_impl &operator=(iterator_impl &&) = default;

  iterator_impl(const iterator_impl<!is_const> &other)
      : m_table(other.m_table), m_current_position(other.m_current_position) {}

  iterator_impl(iterator_impl<!is_const> &&other)
      : m_table(other.m_table), m_current_position(other.m_current_position) {}

  iterator_impl &operator=(const iterator_impl<!is_const> &other) {
    m_table            = other.m_table;
    m_current_position = other.m_current_position;
    return *this;
  }

  iterator_impl &operator=(iterator_impl<!is_const> &&other) {
    m_table = other.m_table;

    m_current_position       = other.m_current_position;
    other.m_current_position = other.m_table->capacity();

    return *this;
  }

  iterator_impl operator++() {
    next_valid_value();
    return *this;
  }

  const iterator_impl operator++(int) {
    iterator_impl tmp = *this;
    next_valid_value();
    return tmp;
  }

  value_type &operator*() { return m_table->at(m_current_position); }

  const value_type &operator*() const {
    return m_table->at(m_current_position);
  }

  value_type *operator->() { return &(m_table->at(m_current_position)); }

  const value_type *operator->() const {
    return &(m_table->at(m_current_position));
  }

  bool operator==(const iterator_impl &rhs) const { return equal(rhs); }
  bool operator!=(const iterator_impl &rhs) const { return !equal(rhs); }

  bool operator==(const iterator_impl<!is_const> &rhs) const {
    return equal(rhs);
  }
  bool operator!=(const iterator_impl<!is_const> &rhs) const {
    return !equal(rhs);
  }

 private:
  void next_valid_value() {
    m_current_position = m_table->find_next(m_current_position + 1);
  }

  bool equal(const iterator_impl &rhs) const {
    return m_table == rhs.m_table &&
           m_current_position == rhs.m_current_position;
  }

  bool equal(const iterator_impl<!is_const> &rhs) const {
    return m_table == rhs.m_table &&
           m_current_position == rhs.m_current_position;
  }

  rhh_table_type *m_table;
  size_type       m_current_position;
};

// template <typename key_type, typename hash, typename key_equal,
//          typename allocator>
// inline bool operator==(
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::iterator &lhs,
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::iterator &rhs) {
//  return lhs.equal(rhs);
//}
//
// template <typename key_type, typename hash, typename key_equal,
//          typename allocator>
// inline bool operator!=(
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::iterator &lhs,
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::iterator &rhs) {
//  return !(lhs == rhs);
//}
//
// template <typename key_type, typename hash, typename key_equal,
//          typename allocator>
// inline bool operator==(
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::const_iterator &lhs,
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::const_iterator &rhs) {
//  return lhs.equal(rhs);
//}
//
// template <typename key_type, typename hash, typename key_equal,
//          typename allocator>
// inline bool operator!=(
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::const_iterator &lhs,
//    const typename unordered_set<key_type, hash, key_equal,
//                                 allocator>::const_iterator &rhs) {
//  return !(lhs == rhs);
//}

}  // namespace rhh
#endif  // HAVOQGT_RHH_UNORDERED_SET_HPP