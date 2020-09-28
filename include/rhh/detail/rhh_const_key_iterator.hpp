// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_INCLUDE_RHH_DETAIL_RHH_CONST_KEY_ITERATOR_HPP_
#define HAVOQGT_INCLUDE_RHH_DETAIL_RHH_CONST_KEY_ITERATOR_HPP_

#include <iterator>

namespace rhh {
namespace detail {

/// \brief Const key iterator.
/// Any update to the container makes this iterator invalid.
template <typename rhh_type>
class rhh_key_iterator {
 public:
  using value_type        = const typename rhh_type::key_type;
  using pointer           = value_type *;
  using reference         = value_type &;
  using size_type         = typename rhh_type::size_type;
  using iterator_category = std::forward_iterator_tag;

  rhh_key_iterator() : m_rhh_table(nullptr), m_current_position(0) {}

  explicit rhh_key_iterator(const rhh_type *const table,
                            const size_type       position = 0)
      : m_rhh_table(table), m_current_position(position) {
    m_current_position = m_rhh_table->find_any(m_current_position);
  }

  rhh_key_iterator(const rhh_key_iterator &)     = default;
  rhh_key_iterator(rhh_key_iterator &&) noexcept = default;
  rhh_key_iterator &operator=(const rhh_key_iterator &) = default;
  rhh_key_iterator &operator=(rhh_key_iterator &&) = default;

  rhh_key_iterator operator++() {
    next_valid_value();
    return *this;
  }

  rhh_key_iterator operator++(int) {
    rhh_key_iterator tmp(*this);
    next_valid_value();
    return tmp;
  }

  value_type &operator*() const {
    return std::get<0>(m_rhh_table->at(m_current_position));
  }

  const value_type *operator->() const {
    return &(std::get<0>(m_rhh_table->at(m_current_position)));
  }

  size_type position() const { return m_current_position; }

  void move_to_end() {
    m_current_position = rhh_type::npos;
  }

  bool operator==(const rhh_key_iterator &rhs) const { return equal(rhs); }

  bool operator!=(const rhh_key_iterator &rhs) const { return !equal(rhs); }

 private:
  void next_valid_value() {
    m_current_position = m_rhh_table->find_any(m_current_position + 1);
  }

  bool equal(const rhh_key_iterator &other) const {
    return m_rhh_table == other.m_rhh_table &&
           m_current_position == other.m_current_position;
  }

  const rhh_type *m_rhh_table{nullptr};
  size_type       m_current_position{0};
};

}  // namespace detail
}  // namespace rhh
#endif  // HAVOQGT_INCLUDE_RHH_DETAIL_RHH_CONST_KEY_ITERATOR_HPP_
