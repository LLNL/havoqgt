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

#ifndef HAVOQGT_RHH_BASIC_RHH_HPP
#define HAVOQGT_RHH_BASIC_RHH_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <utility>

#include <rhh/detail/utility/dynamic_array.hpp>
#include <rhh/detail/utility/packed_tuple.hpp>
#include <rhh/detail/utility/value_type.hpp>

namespace rhh {
namespace detail {
//namespace utility {
//
//struct void_mapped_value_tag;
//
//// ----------------------------------------------------------------------------------------------------
//// //
////                                 Class value_type_selector
//// utility class to construct a rhh class which doesn't have mapped_value
//// ----------------------------------------------------------------------------------------------------
//// //
///// \brief Generalized class for the cases mapped_value_type is not
///// void_mapped_value_tag value_type is std::pair<key_type, mapped_value_type>
///// \tparam key_type Type of the key
///// \tparam mapped_value_type Type of the value
//template <typename key_type, typename mapped_value_type>
//class value_accessor {
// public:
//  using value_type = std::pair<key_type, mapped_value_type>;
//
//  static const key_type &get_key(value_type &value) { return value.first; }
//
//  static const key_type &get_key(const value_type &value) {
//    return value.first;
//  }
//
//  static void destroy_mapped_value(value_type *value) {
//    value->second.~mapped_value_type();
//  }
//
//  // ----- Explicitly delete constructors ----- //
//  using self_type  = value_accessor<key_type, mapped_value_type>;
//  value_accessor() = delete;
//  value_accessor(const self_type &) = delete;
//  value_accessor(self_type &&)      = delete;
//  ~value_accessor()                 = delete;
//};
//
///// \brief Specialized class whose second template type is void_mapped_value_tag
///// therefore, value_type is same to key_type
///// \tparam key_type Type of the key
//template <typename key_type>
//class value_accessor<key_type, void_mapped_value_tag> {
// public:
//  using value_type = key_type;
//
//  static key_type &get_key(value_type &value) { return value; }
//
//  static const key_type &get_key(const value_type &value) { return value; }
//
//  static void destroy_mapped_value(value_type *) {
//    // Does nothing
//  }
//
//  // ----- Explicitly delete constructors ----- //
//  using self_type  = value_accessor<key_type, void_mapped_value_tag>;
//  value_accessor() = delete;
//  value_accessor(const self_type &) = delete;
//  value_accessor(self_type &&)      = delete;
//  ~value_accessor()                 = delete;
//};
//
//}  // namespace utility

/// \brief Robin Hood Hashing class
/// \tparam key_type Type of the key values
/// \tparam mapped_value_type Type of the mapped values
/// \tparam hash Hash function
/// \tparam EqualTo Equal to function
/// \tparam AllocatorType Type of the allocator
template <typename _value_wrapper_type, typename _hash, typename _key_equal,
          typename _allocator>
class basic_rhh {
 public:
  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Public type
  // ----------------------------------------------------------------------------------------------------
  // //
  using key_type       = typename _value_wrapper_type::key_type;
  using value_type     = typename _value_wrapper_type::value_type;
  using hash           = _hash;
  using key_equal      = _key_equal;
  using allocator_type = _allocator;
  using size_type = typename std::allocator_traits<allocator_type>::size_type;

 private:
  using self_type =
      basic_rhh<_value_wrapper_type, _hash, _key_equal, _allocator>;

  struct element_type {
    class header_type;

    element_type()                         = default;
    ~element_type()                        = default;
    element_type(const element_type &)     = default;
    element_type(element_type &&) noexcept = default;
    element_type &operator=(const element_type &) = default;
    element_type &operator=(element_type &&) noexcept = default;

    header_type         header;
    _value_wrapper_type value_wrapper;
  };

  using internal_table_allocator_type = typename std::allocator_traits<
      _allocator>::template rebind_alloc<element_type>;

  using internal_table_type =
      utility::dynamic_array<element_type,
                             utility::scp_alloc<internal_table_allocator_type>>;

  // Internal data structure
  //  element {
  //    header {
  //      probe_distance
  //      tomb_stone
  //    };
  //    value {
  //      key
  //      mapped_value // or no mapped_value
  //    };
  //  };

 public:
  // --------------------------------------------------------------------------------
  // // Constructor & assign operator
  // --------------------------------------------------------------------------------
  // //
  explicit basic_rhh(const _allocator &allocator = _allocator())
      : m_num_elements(0), m_table(1, element_type(), allocator) {}

  basic_rhh(const size_type   initial_capacity,
            const _allocator &allocator = _allocator())
      : m_num_elements(0),
        m_table(initial_capacity, element_type(), allocator) {}

  ~basic_rhh() = default;

  // Copy constructor
  basic_rhh(const basic_rhh &) = default;

  // Move constructor
  basic_rhh(basic_rhh &&other) noexcept
      : m_num_elements(other.m_num_elements),
        m_table(std::move(other.m_table)) {
    other.m_num_elements = 0;
  }

  // Copy assignments
  basic_rhh &operator=(const basic_rhh &) = default;

  // Move assignments
  basic_rhh &operator=(basic_rhh &&other) noexcept {
    m_num_elements       = other.m_num_elements;
    other.m_num_elements = 0;
    m_table              = std::move(other.m_table);
    return *this;
  }

  // --------------------------------------------------------------------------------
  // // Public methods
  // --------------------------------------------------------------------------------
  // //
  // -------------------- Capacity -------------------- //
  /// \brief Returns the number of values in the container.
  /// \return the number of values in the container.
  size_type size() const { return m_num_elements; }

  /// \brief Returns the capacity of the container.
  /// \return The current capacity of the container.
  size_type capacity() const { return m_table.size(); }

  // -------------------- Element access -------------------- //
  /// \brief Accesses the element at 'position'.
  /// \param position The position of a value to access.
  /// \return A reference to the value at 'position'.
  /// Specifically, returns std::pair<key, mapped_value>&.
  /// Note that 'key' is not const;
  /// however, if key is modified, the container will be invalid.
  // TODO: need to return std::pair<const key, mapped_value> instead of
  // std::pair<key, mapped_value>
  //  when mapped_value is not utility::void_mapped_value_tag
  value_type &at(const size_type position) {
    assert(position < capacity());
    return m_table[position].value_wrapper.value();
  }

  /// \brief Accesses the element at 'position'.
  /// \param position The position of a value to access.
  /// \return A const reference to the element at 'position'.
  const value_type &at(const size_type position) const {
    assert(position < capacity());
    return m_table[position].value_wrapper;
  }

  // -------------------- Element lookup -------------------- //
  /// \brief Finds a value with key 'key'.
  /// \param key A key to search.
  /// \return The position of a value found.
  /// If not found, returns capacity().
  size_type find(const key_type &key) const {
    const auto ret = priv_locate_key(key);
    if (ret.second) {
      return ret.first;  // Found the key
    }
    return capacity();  // Didn't find the key
  }

  /// \brief Finds a value with key 'key'.
  /// \param key A key to search.
  /// \param hint_position A hint to find the value
  /// \return The position of a value found.
  size_type find(const key_type &key, const size_type hint_position) const {
    const auto ret = priv_locate_key(key, hint_position);
    if (ret.second) {
      return ret.first;  // Found the key
    } else {
      return capacity();  // Didn't find the key
    }
  }

  /// \brief Finds the next value with any key.
  /// \param start_position The current position.
  /// \return The position of the next value.
  size_type find_next(const size_type start_position) const {
    return priv_find_valid_element(start_position);
  }

  // -------------------- Modifiers -------------------- //
  /// \brief Inserts a value. Does not check duplicate element.
  /// \param value A value to insert.
  /// \return Returns the position the value was inserted.
  size_type insert(const value_type &value) {
    return priv_check_capacity_and_insert(value);
  }

  /// \brief Inserts a value. Does not check duplicate element.
  /// \param value A value to insert.
  /// \return Returns the position the value was inserted.
  size_type insert(value_type &&value) {
    return priv_check_capacity_and_insert(std::move(value));
  }

  /// \brief Erases values with the key.
  /// \param key A key to erase.
  /// \return The number of values erased.
  size_type erase(const key_type &key) {
    const auto old_num_elements = size();
    priv_erase_multiple(key);
    return old_num_elements - size();
  }

  /// \brief Erases the element at 'position'
  /// \param position The position of a value to erase.
  void erase_at(const size_type position) { priv_erase_at(position); }

  /// \brief Clears all values. Does not change the capacity.
  void clear() {
    m_num_elements = 0;
    for (size_type i = 0; i < capacity(); ++i) {
      if (!m_table[i].header.empty()) m_table[i].header.reset();
    }
  }

  void swap(basic_rhh &other) noexcept {
    using std::swap;
    swap(m_num_elements, other.m_num_elements);
    swap(m_table, other.m_table);
  }

  // -------------------- Etc -------------------- //
  allocator_type get_allocator() const noexcept {
    return m_table.get_allocator();
  }

  size_type probe_distance(const size_type &position) const {
    return priv_probe_distance(position);
  }

  // -------------------- Statistic -------------------- //
  /// \brief Return an average probe distance of valid (non empty) elements
  /// \return an average probe distance of valid (non empty) elements
  double load_factor() const {
    size_type sum = 0;
    for (size_type i = 0; i < capacity(); ++i) {
      if (!m_table[i].header.empty()) {
        sum += priv_probe_distance(i);
      }
    }

    return static_cast<double>(sum) / m_num_elements;
  }

 private:
  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Private static constant variables
  // ----------------------------------------------------------------------------------------------------
  // //
  static constexpr const size_type k_table_lenght_growing_factor = 2;

  // ----------------------------------------------------------------------------------------------------
  // //
  //                               Private functions
  // ----------------------------------------------------------------------------------------------------
  // //

  // -------------------- Hash -------------------- //
  size_type priv_hash_key(const key_type &key) const { return hash()(key); }

  size_type priv_ideal_position(const key_type &key) const {
    // std::cout << key << " -> " << priv_hash_key(key) << " : " <<
    // (priv_hash_key(key) & (table_length() - 1)) << std::endl; std::cout <<
    // key.first << " : " <<  key.second << " -> " << priv_hash_key(key) << " :
    // " << (priv_hash_key(key) & (capacity() - 1)) << std::endl;
    return priv_hash_key(key) & (capacity() - 1);
  }

  // -------------------- Probe distance -------------------- //
  size_type priv_probe_distance(const size_type current_position) const {
    const key_type &key = m_table[current_position].value_wrapper.key();
    return priv_probe_distance(key, current_position);
  }

  size_type priv_probe_distance(const key_type &key,
                                const size_type current_position) const {
    return (current_position + capacity() - priv_ideal_position(key)) &
           (capacity() - 1);
  }

  // -------------------- Capacity -------------------- //
  bool priv_enough_capacity() { return (m_num_elements < capacity() * 0.9); }

  void priv_grow_table() {
    const size_type new_capacity = capacity() * k_table_lenght_growing_factor;
    self_type       new_table(new_capacity, get_allocator());

    for (size_type i = 0; i < capacity(); ++i) {
      if (!m_table[i].header.empty() && !m_table[i].header.get_tomb_stone()) {
        new_table.priv_insert(std::move(m_table[i].value_wrapper.value()));
      }
    }

    (*this).swap(new_table);
  }

  // -------------------- Look up -------------------- //
  std::pair<size_type, bool> priv_locate_key(const key_type &key) const {
    return priv_locate_key(key, priv_ideal_position(key));
  }

  std::pair<size_type, bool> priv_locate_key(
      const key_type &key, const size_type hint_position) const {
    size_type current_position       = hint_position & (capacity() - 1);
    size_type current_probe_distance = 0;

    bool is_found_key = false;

    while (true) {
      if (m_table[current_position].header.empty()) {
        break;
      } else if (current_probe_distance >
                 priv_probe_distance(current_position)) {
        break;
      } else if (key_equal()(m_table[current_position].value_wrapper.key(),
                             key) &&
                 !m_table[current_position].header.get_tomb_stone()) {
        is_found_key = true;
        break;
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }

    return std::make_pair(current_position, is_found_key);
  }

  size_type priv_find_valid_element(const size_type start_position) const {
    if (start_position == capacity()) return capacity();

    size_type position = start_position;
    while (true) {
      if (!(m_table[position].header.empty() ||
            m_table[position].header.get_tomb_stone()))
        return position;

      position = (position + 1) & (capacity() - 1);

      if (position == start_position) {
        return capacity();
      }
    }
  }

  // -------------------- Insert -------------------- //
  template <typename T>
  size_type priv_check_capacity_and_insert(T &&value) {
    if (!priv_enough_capacity()) {
      priv_grow_table();
    }

    return priv_insert(std::forward<T>(value));
  }

  template <typename T>
  size_type priv_insert(T &&value) {
    // Find the position to insert the value
    auto insert_position = priv_locate_key(_value_wrapper_type::key(value));
    while (insert_position.second) {  // skip elements with the same key
      insert_position = priv_locate_key(_value_wrapper_type::key(value),
                                        insert_position.first + 1);
    }

    priv_insert_core(std::forward<T>(value), insert_position.first);

    return insert_position.first;
  }

  template <typename T>
  void priv_insert_core(T &&value, const size_type first_insert_position) {
    size_type current_position = first_insert_position;
    size_type current_probe_distance =
        priv_probe_distance(_value_wrapper_type::key(value), current_position);

    _value_wrapper_type wk_value(std::forward<T>(value));
    ++m_num_elements;

    while (true) {
      if (m_table[current_position].header.empty()) {
        priv_emplace_at(current_position, std::move(wk_value));
        return;
      }

      if (current_probe_distance > priv_probe_distance(current_position)) {
        if (m_table[current_position].header.get_tomb_stone()) {
          priv_emplace_at(current_position, std::move(wk_value));
          return;
        }

        using std::swap;
        swap(m_table[current_position].value_wrapper, wk_value);
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }
  }

  void priv_emplace_at(const size_type       position,
                       _value_wrapper_type &&value_wrapper) {
    m_table[position].value_wrapper = std::move(value_wrapper);
    m_table[position].header.reset();
    m_table[position].header.set();
  }

  // -------------------- Erase -------------------- //
  void priv_erase_multiple(const key_type &key) {
    const auto ret = priv_locate_key(key);
    if (!ret.second) {
      return;
    }

    priv_erase_multiple(key, ret.first);
  }

  void priv_erase_multiple(const key_type &key,
                           const size_type first_position) {
    size_type current_position       = first_position;
    size_type current_probe_distance = 0;

    while (true) {
      if (m_table[current_position].header.empty()) {
        break;
      } else if (current_probe_distance >
                 priv_probe_distance(current_position)) {
        break;
      } else if (key_equal()(m_table[current_position].value_wrapper.key(),
                             key) &&
                 !m_table[current_position].header.get_tomb_stone()) {
        priv_erase_at(current_position);
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }
  }

  void priv_erase_at(const size_type position) {
    assert(!m_table[position].header.get_tomb_stone());
    m_table[position].header.set_tomb_stone();
    m_table[position].value_wrapper.~_value_wrapper_type();

    assert(m_num_elements > 0);
    --m_num_elements;
  }

  size_type           m_num_elements;
  internal_table_type m_table;
};

template <typename _value_type, typename _hash, typename _key_equal,
          typename _allocator>
class basic_rhh<_value_type, _hash, _key_equal,
                _allocator>::element_type::header_type {
 public:
  using underlying_type = uint8_t;

  header_type()                        = default;
  ~header_type()                       = default;
  header_type(const header_type &)     = default;
  header_type(header_type &&) noexcept = default;
  header_type &operator=(const header_type &) = default;
  header_type &operator=(header_type &&) noexcept = default;

  bool get_tomb_stone() const { return static_cast<bool>(m_data & 0x1); }

  void set_tomb_stone() { m_data |= 0x1; }

  void reset_tomb_stone() {
    const underlying_type mask =
        ~static_cast<underlying_type>(0) - static_cast<underlying_type>(1);
    m_data &= mask;
  }

  void set() { m_data |= 0x10; }

  void reset() { m_data = static_cast<underlying_type>(0); }

  bool empty() const { return m_data == static_cast<underlying_type>(0); }

 private:
  underlying_type m_data{0};
};

}  // namespace detail
}  // namespace rhh

#endif  // HAVOQGT_RHH_BASIC_RHH_HPP
