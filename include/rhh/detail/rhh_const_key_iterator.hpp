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
