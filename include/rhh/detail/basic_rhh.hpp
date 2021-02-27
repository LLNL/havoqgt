// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_RHH_BASIC_RHH_HPP
#define HAVOQGT_RHH_BASIC_RHH_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
#include <utility>

#include <boost/container/scoped_allocator.hpp>
#include <boost/container/vector.hpp>
#include <rhh/detail/utility/value_type.hpp>

namespace rhh {
namespace detail {

// --------------------------------------------------------------------------
//                        Forward Declaration
// --------------------------------------------------------------------------
template <typename value_wrapper, typename hash, typename key_equal,
          typename allocator>
class basic_rhh;

class rhh_header {
 private:
  static constexpr uint8_t k_tomb_stone_flag = 0x01;
  static constexpr uint8_t k_hold_value_flag = 0x10;

 public:
  using raw_type = uint8_t;

  rhh_header()                       = delete;
  ~rhh_header()                      = delete;
  rhh_header(const rhh_header &)     = delete;
  rhh_header(rhh_header &&) noexcept = delete;
  rhh_header &operator=(const rhh_header &) = delete;
  rhh_header &operator=(rhh_header &&) noexcept = delete;

  static void clear(raw_type &data) { data = static_cast<raw_type>(0); }

  static bool get_tomb_stone(const raw_type &data) {
    return static_cast<bool>(data & k_tomb_stone_flag);
  }

  static void set_tomb_stone(raw_type &data) { data |= k_tomb_stone_flag; }

  static void init_for_new_value(raw_type &data) {
    clear(data);
    data = k_hold_value_flag;
  }

  static bool empty(const raw_type &data) {
    return data == static_cast<raw_type>(0);
  }
};

/// \brief Robin Hood Hashing class
/// \tparam _value_wrapper Type of a value (key and mapped value)
/// \tparam hash Hash function
/// \tparam EqualTo Equal to function
/// \tparam AllocatorType Type of the allocator
template <typename _value_wrapper, typename _hash, typename _key_equal,
          typename _allocator_type>
class basic_rhh {
 private:
  // --------------------------------------------------------------------------
  // Private types & private static variables
  // --------------------------------------------------------------------------
  using self_type =
      basic_rhh<_value_wrapper, _hash, _key_equal, _allocator_type>;

  // Use tuple instead of pair because nested pair causes a problem with
  // scoped allocator adaptor.
  // C++ Standards Committee Library Working Group (LWG) issue# is 2975
  using block_type =
      std::tuple<rhh_header::raw_type, typename _value_wrapper::value_type>;
  using internal_table_allocator_type = typename std::allocator_traits<
      _allocator_type>::template rebind_alloc<block_type>;
  using internal_table_type =
      boost::container::vector<block_type,
                               boost::container::scoped_allocator_adaptor<
                                   internal_table_allocator_type>>;

  static constexpr double k_max_load_factor = 0.9;

 public:
  // --------------------------------------------------------------------------
  // Public type & public static variables
  // --------------------------------------------------------------------------
  using key_type       = typename _value_wrapper::key_type;
  using value_type     = typename _value_wrapper::value_type;
  using hash           = _hash;
  using key_equal      = _key_equal;
  using allocator_type = _allocator_type;
  using size_type = typename std::allocator_traits<allocator_type>::size_type;
  // friend key_iterator;
  // '- 1' for preventing an overflow with npos + 1.
  static constexpr size_type npos = std::numeric_limits<size_type>::max() - 1;

 public:
  // --------------------------------------------------------------------------
  // Constructor & assign operator
  // --------------------------------------------------------------------------

  // Note: always has at least one capacity to simplify the implementation

  basic_rhh() : m_num_blocks(0), m_table(1) {}

  explicit basic_rhh(const allocator_type &allocator)
      : m_num_blocks(0),
        m_table(1, block_type(std::allocator_arg, allocator), allocator) {}

  explicit basic_rhh(const size_type initial_capacity)
      : m_num_blocks(0), m_table(std::max(initial_capacity, (size_type)1)) {}

  basic_rhh(const size_type initial_capacity, const allocator_type &allocator)
      : m_num_blocks(0),
        m_table(std::max(initial_capacity, (size_type)1),
                block_type(std::allocator_arg, allocator), allocator) {}

  ~basic_rhh() = default;

  // Copy constructor
  basic_rhh(const basic_rhh &) = default;

  // Allocator-extended copy constructor
  basic_rhh(const basic_rhh &other, const allocator_type &allocator)
      : m_num_blocks(other.m_num_blocks),
        m_table(other.m_table, internal_table_allocator_type(allocator)) {}

  // Move constructor
  basic_rhh(basic_rhh &&other) noexcept
      : m_num_blocks(other.m_num_blocks), m_table(std::move(other.m_table)) {
    other.m_num_blocks = 0;
  }

  // Allocator-extended move constructor
  basic_rhh(basic_rhh &&other, const allocator_type &allocator)
      : m_num_blocks(other.m_num_blocks),
        m_table(std::move(other.m_table),
                internal_table_allocator_type(allocator)) {
    other.m_num_blocks = 0;
  }

  // Copy assignments
  basic_rhh &operator=(const basic_rhh &) = default;

  // Move assignments
  basic_rhh &operator=(basic_rhh &&other) noexcept {
    m_num_blocks       = other.m_num_blocks;
    other.m_num_blocks = 0;
    m_table            = std::move(other.m_table);
    return *this;
  }

  // --------------------------------------------------------------------------
  // Public methods
  // --------------------------------------------------------------------------

  // -------------------- Capacity -------------------- //
  /// \brief Returns the number of values in the container.
  /// \return the number of values in the container.
  size_type size() const { return m_num_blocks; }

  /// \brief Returns the capacity of the container.
  /// \return The current capacity of the container.
  size_type capacity() const { return m_table.size(); }

  /// \brief Returns the theoretical maximum capacity.
  /// \return The maximum capacity.
  size_type max_capacity() { return std::min(npos - 1, m_table.max_size()); }

  /// \brief Returns the maximum number of values can hold theoretically.
  /// \return The maximum number of values.
  size_type max_size() { return priv_max_num_blocks(max_capacity()); }

  /// \brief Sets the number of buckets to the number needed to accommodate at least count elements
  /// without exceeding maximum load factor and rehashes the container,
  /// i.e. puts the elements into appropriate buckets considering that total number of buckets has changed.
  /// \param count new size of the container
  void reserve(const size_type count) {
    if (count <= priv_max_num_blocks(capacity())) {
      return; // Enough capacity
    }

    size_type new_capacity = capacity();
    while (priv_max_num_blocks(new_capacity) < count) {
      new_capacity *= k_table_lenght_growing_factor;
    }

    priv_grow_table(capacity() * k_table_lenght_growing_factor);
  }

  // -------------------- Element access -------------------- //
  /// \brief Accesses the block at 'position'.
  /// \param position An position of a value to access.
  /// \return A reference to the value at 'position'.
  /// Specifically, returns std::pair<key, mapped_value>&.
  /// Note that 'key' is not const;
  /// however, if key is modified, the container will be invalid.
  value_type &at(const size_type position) {
    assert(position < capacity());
    return priv_value_at(position);
  }

  /// \brief Accesses the block at 'position'.
  /// \param position An position of a value to access.
  /// \return A const reference to the block at 'position'.
  const value_type &at(const size_type position) const {
    assert(position < capacity());
    return priv_value_at(position);
  }

  // -------------------- Element lookup -------------------- //
  /// \brief Finds a value with key 'key'.
  /// \param key A key to search.
  /// \return The position of a value found.
  /// If not found, returns npos.
  size_type find(const key_type &key) const {
    const auto ret = priv_locate_key(key);
    if (ret.second) {
      return ret.first;  // Found the key
    }
    return npos;  // Didn't find the key
  }

  /// \brief Finds a value with key 'key', starting from start_position.
  /// This function does not circulate the table, i.e., a returned
  /// value is always equal to or larger than a given start position.
  /// \param key A key to search.
  /// \param start_position A start position to find the value.
  /// \return The position of a value found.
  /// If there is no value found, returns the result of npos.
  size_type find(const key_type &key, const size_type start_position) const {
    const auto ret = priv_locate_key(key, start_position);
    if (ret.second && start_position <= ret.first) {
      return ret.first;  // Found the key
    }
    return npos;  // Didn't find the key
  }

  /// \brief Finds the first valid value with any key from a given start
  /// position. This function does not circulate the table, i.e., a returned
  /// value is always equal to or larger than a given start position.
  /// \param start_position A start position.
  /// \return The position of the first valid value.
  /// If there is no valid element, returns the result of npos.
  size_type find_any(const size_type start_position) const {
    return priv_find_any_valid_value(start_position);
  }

  /// \brief Returns the number of values that have the key.
  /// \param key A key to count.
  /// \return The number of values that have the key.
  size_type count(const key_type &key) const {
    return priv_count_same_keys(key);
  }

  // -------------------- Modifiers -------------------- //
  /// \brief Inserts a value. Does not check duplicate block.
  /// \param value A value to insert.
  /// \return Returns the position the value was inserted.
  size_type insert(const value_type &value) {
    return priv_check_capacity_and_insert(value);
  }

  /// \brief Inserts a value. Does not check duplicate block.
  /// \param value A value to insert.
  /// \return Returns the position the value was inserted.
  size_type insert(value_type &&value) {
    return priv_check_capacity_and_insert(std::move(value));
  }

  /// \brief Inserts a key. Does not check duplicate block.
  /// mapped_value will be uninitialized.
  /// \param key A key to insert.
  /// \return Returns the position the key was inserted.
  size_type insert_key(const key_type &key) {
    return priv_check_capacity_and_insert(key);
  }

  /// \brief Inserts a key. Does not check duplicate block.
  /// mapped_value will be uninitialized.
  /// \param key A key to insert.
  /// \return Returns the position the key was inserted.
  size_type insert_key(key_type &&key) {
    return priv_check_capacity_and_insert(std::move(key));
  }

  /// \brief Erases values with the key.
  /// \param key A key to erase.
  /// \return The number of values erased.
  size_type erase(const key_type &key) {
    const auto old_num_blocks = size();
    priv_erase_multiple(key);
    return old_num_blocks - size();
  }

  /// \brief Erases the block at 'position'
  /// \param position The position of a value to erase.
  void erase_at(const size_type position) { priv_erase_at(position); }

  /// \brief Clears all values. Does not change the capacity.
  void clear() {
    m_num_blocks = 0;
    for (size_type i = 0; i < capacity(); ++i) {
      if (!rhh_header::empty(priv_header_at(i)))
        rhh_header::clear(priv_header_at(i));
    }
  }

  void swap(basic_rhh &other) noexcept {
    using std::swap;
    swap(m_num_blocks, other.m_num_blocks);
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
  /// \brief Return an average probe distance of valid (non empty) blocks
  /// \return an average probe distance of valid (non empty) blocks
  auto load_factor() const {
    size_type sum = 0;
    for (size_type i = 0; i < capacity(); ++i) {
      if (!rhh_header::empty(priv_header_at(i))) {
        sum += priv_probe_distance(i);
      }
    }

    return static_cast<double>(sum) / m_num_blocks;
  }

 private:
  // --------------------------------------------------------------------------
  // Private static constant variables
  // --------------------------------------------------------------------------
  static constexpr const size_type k_table_lenght_growing_factor = 2;

  // --------------------------------------------------------------------------
  // Private functions
  // --------------------------------------------------------------------------
  auto &priv_header_at(const size_type position) {
    return std::get<0>(m_table[position]);
  }

  const auto &priv_header_at(const size_type position) const {
    return std::get<0>(m_table[position]);
  }

  value_type &priv_value_at(const size_type position) {
    return std::get<1>(m_table[position]);
  }

  const value_type &priv_value_at(const size_type position) const {
    return std::get<1>(m_table[position]);
  }

  // -------------------- Hash -------------------- //
  size_type priv_hash_key(const key_type &key) const { return hash()(key); }

  size_type priv_ideal_position(const key_type &key) const {
    return priv_hash_key(key) & (capacity() - 1);
  }

  // -------------------- Probe distance -------------------- //
  size_type priv_probe_distance(const size_type position) const {
    return priv_probe_distance(_value_wrapper::key(priv_value_at(position)),
                               position);
  }

  size_type priv_probe_distance(const key_type &key,
                                const size_type position) const {
    return (position + capacity() - priv_ideal_position(key)) &
           (capacity() - 1);
  }

  // -------------------- Capacity -------------------- //
  static constexpr size_type priv_max_num_blocks(const size_type capacity) {
    return capacity * k_max_load_factor;
  }

  void priv_grow_table(const size_type new_capacity) {
    self_type       new_table(new_capacity, get_allocator());

    for (size_type i = 0; i < capacity(); ++i) {
      if (!rhh_header::empty(priv_header_at(i)) &&
          !rhh_header::get_tomb_stone(priv_header_at(i))) {
        new_table.priv_insert(std::move(priv_value_at(i)));
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
      const auto &header = priv_header_at(current_position);
      const auto &value  = priv_value_at(current_position);
      if (rhh_header::empty(header)) {
        break;
      } else if (current_probe_distance >
                 priv_probe_distance(current_position)) {
        break;
      } else if (key_equal()(_value_wrapper::key(value), key) &&
                 !rhh_header::get_tomb_stone(header)) {
        is_found_key = true;
        break;
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }

    return std::make_pair(current_position, is_found_key);
  }

  size_type priv_find_any_valid_value(const size_type start_position) const {
    for (auto position = start_position; position < capacity(); ++position) {
      if (!rhh_header::empty(priv_header_at(position)) &&
          !rhh_header::get_tomb_stone(priv_header_at(position)))
        return position;
    }
    return npos;
  }

  size_type priv_count_same_keys(const key_type &key) const {
    auto      position = priv_locate_key(key);
    size_type count    = 0;
    while (position.second) {
      ++count;
      position = priv_locate_key(key, position.first + 1);
    }

    return count;
  }

  // -------------------- Insert -------------------- //
  template <typename T>
  size_type priv_check_capacity_and_insert(T &&value) {
    if (m_num_blocks >= priv_max_num_blocks(capacity())) {
      priv_grow_table(capacity() * k_table_lenght_growing_factor);
    }

    return priv_insert(std::forward<T>(value));
  }

  template <typename T>
  size_type priv_insert(T &&value) {
    // Find the position to insert the value
    auto insert_position = priv_locate_key(_value_wrapper::key(value));
    while (insert_position.second) {  // skip blocks with the same key
      insert_position = priv_locate_key(_value_wrapper::key(value),
                                        insert_position.first + 1);
    }

    priv_insert_core(std::forward<T>(value), insert_position.first);

    return insert_position.first;
  }

  template <typename T>
  void priv_insert_core(T &&value, const size_type first_insert_position) {
    size_type current_position = first_insert_position;
    size_type current_probe_distance =
        priv_probe_distance(_value_wrapper::key(value), current_position);

    auto wk_value =
        _value_wrapper::allocate_value(get_allocator(), std::forward<T>(value));

    ++m_num_blocks;

    while (true) {
      if (rhh_header::empty(priv_header_at(current_position))) {
        priv_set_value_at(current_position, std::move(wk_value));
        return;
      }

      if (current_probe_distance > priv_probe_distance(current_position)) {
        if (rhh_header::get_tomb_stone(priv_header_at(current_position))) {
          priv_set_value_at(current_position, std::move(wk_value));
          return;
        }

        using std::swap;
        swap(priv_value_at(current_position), wk_value);
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }
  }

  void priv_set_value_at(const size_type position, value_type &&value) {
    priv_value_at(position) = std::move(value);
    rhh_header::init_for_new_value(priv_header_at(position));
  }

  /// \brief Insert only a key.
  /// This function is enabled if value_wrapper has both key and value.
  template <typename T, std::enable_if_t<std::is_same<T, key_type>::value &&
                                         _value_wrapper::size() == 2> = 0>
  void priv_set_value_at(const size_type position, T &&key) {
    _value_wrapper::key(priv_value_at(position)) = std::move(key);
    rhh_header::init_for_new_value(priv_header_at(position));
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
      const auto &header = priv_header_at(current_position);
      const auto &value  = priv_value_at(current_position);

      if (rhh_header::empty(header)) {
        break;
      } else if (current_probe_distance >
                 priv_probe_distance(current_position)) {
        break;
      } else if (key_equal()(_value_wrapper::key(value), key) &&
                 !rhh_header::get_tomb_stone(header)) {
        priv_erase_at(current_position);
      }

      current_position = (current_position + 1) & (capacity() - 1);
      ++current_probe_distance;
    }
  }

  void priv_erase_at(const size_type position) {
    assert(!rhh_header::get_tomb_stone(priv_header_at(position)));
    rhh_header::set_tomb_stone(priv_header_at(position));
    priv_value_at(position).~value_type();

    assert(m_num_blocks > 0);
    --m_num_blocks;
  }

  size_type           m_num_blocks{0};
  internal_table_type m_table;
};

template <typename value_type, typename hash, typename key_equal,
          typename allocator>
inline void swap(
    basic_rhh<value_type, hash, key_equal, allocator> &lhs,
    basic_rhh<value_type, hash, key_equal, allocator> &rhs) noexcept {
  lhs.swap(rhs);
}

}  // namespace detail
}  // namespace rhh

#endif  // HAVOQGT_RHH_BASIC_RHH_HPP
