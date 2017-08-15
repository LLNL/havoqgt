//
// Created by Iwabuchi, Keita on 8/14/17.
//
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

#ifndef HAVOQGT_MPI_FIXED_SIZE_UNORDERED_MAP_HPP
#define HAVOQGT_MPI_FIXED_SIZE_UNORDERED_MAP_HPP

#include <tuple>
#include <limits>
#include <vector>

namespace havoqgt
{

template <typename key_t,
  typename value_t,
  typename hash = std::hash<key_t>,
  typename key_equal = std::equal_to<key_t>,
  typename allocator = std::allocator<std::tuple<bool, key_t, value_t>>>
class fixed_size_unordered_map
{
 public:
  using size_type = typename allocator::size_type;

 private:
  using self_type = fixed_size_unordered_map<key_t, value_t, hash, key_equal, allocator>;

  /// enum class table_element { empty_flag = 0, key = 1, value = 2 };

  using table_element_t = std::tuple<bool, key_t, value_t>; /// 0:is_empty flag, 1:key, 2:value
  using table_t = std::vector<table_element_t, allocator>;

  static const bool k_empty = true;
  static const bool k_not_empty = false;

 public:

  /// ---------- Constructors ---------- ///
  /// Default constructor
  fixed_size_unordered_map() :
    m_size(0),
    m_table()
  { }

  /// constructor
  explicit fixed_size_unordered_map(const size_type size) :
    m_size(size),
    m_table(m_size, std::make_tuple(k_empty, key_t(), value_t()))
  { }

  /// Copy constructor
  fixed_size_unordered_map(const self_type& other) :
    m_size(other.m_size),
    m_table(other.m_table)
  { }

  /// Move constructor
  fixed_size_unordered_map(self_type&& other) :
    m_size(std::move(other.m_size)),
    m_table(std::move(other.m_table))
  {
    other.m_size = 0;
  }

  /// Copy and move assignments
  self_type& operator=(self_type rhs)
  {
    rhs.swap(*this);
    return (*this);
  }

  /// ---------- Operators ---------- ///
  value_t& operator[](const key_t& key)
  {
    const size_type position = search_key_or_empty_position(key);

    if (is_empty(position))
      set_key_at(key, position);

    return std::get<2>(m_table[position]);
  }

  /// ---------- Capacity ---------- ///
  inline size_type size() const
  {
    return m_size;
  }

  inline void resize(const size_type size)
  {
    m_size = size;
    m_table.resize(m_size, std::make_tuple(k_empty, key_t(), value_t()));
  }


  /// ---------- Iterators ---------- ///
  /// This implementation is not correct because it access every element in the table even its empty flag is true
  inline typename table_t::iterator begin()
  {
    return m_table.begin();
  }

  inline typename table_t::iterator end()
  {
    return m_table.end();
  }

  inline typename table_t::const_iterator cbegin()
  {
    return m_table.cbegin();
  }

  inline typename table_t::const_iterator cend()
  {
    return m_table.cend();
  }


 private:

  /// TODO: make sure this function is called
  inline void swap(self_type& other)
  {
    using std::swap;
    swap(m_size, other.m_size);
    swap(m_table, other.m_table);
  }

  inline bool is_empty(const size_type position) const
  {
    return std::get<0>(m_table[position]);
  }

  inline size_type search_key_or_empty_position(const key_t& key) const
  {
    const size_type ideal_position = static_cast<size_type>(hash()(key)) % m_size;

    for (size_type i = 0; i < m_size; ++i) {
      const size_type position = (ideal_position + i) % m_size;

      if (std::get<1>(m_table[position]) == key || is_empty(position))
        return position;
    }

    exit(1);
  }

  inline void set_key_at(const key_t& key, const size_type position)
  {
    m_table[position] = std::make_tuple(k_not_empty, key, value_t());
  }

  size_type m_size;
  table_t m_table;
};

} /// namespace havoqgt

#endif //HAVOQGT_MPI_FIXED_SIZE_UNORDERED_MAP_HPP
