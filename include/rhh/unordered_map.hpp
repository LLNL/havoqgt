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

#ifndef UNORDERED_MAP_HPP
#define UNORDERED_MAP_HPP

#include "rhh_no_delete.hpp"

namespace rhh {

/// \brief Robin Hood Hashing class, however not supporting deletes
/// Following the interface of std::map
/// \tparam KeyType Type of the key values
/// \tparam MappedValueType Type of the mapped values
/// \tparam Hasher Hash function
/// \tparam EqualTo Equal to function
/// \tparam AllocatorType Type of the allocator
template <typename KeyType,
  typename MappedValueType,
  typename Hash = std::hash<typename std::remove_const<KeyType>::type>,
  typename EqualTo = std::equal_to<KeyType>,
  typename AllocatorType = std::allocator<void>>
class
unordered_map
{
 private:
  using self_type = unordered_map<KeyType, MappedValueType, Hash, EqualTo, AllocatorType>;

  using table_element_type = typename rhh::detail::rhh_no_delete_internal_type<KeyType, MappedValueType>::table_element_type;
  using allocator_type = typename std::allocator_traits<AllocatorType>::template rebind_alloc<table_element_type>;
  template <typename T>
  using underling_table_type = std::vector<T, allocator_type>;
  using rhh_table_type = rhh::detail::rhh_no_delete<KeyType, MappedValueType, underling_table_type, Hash, EqualTo>;

 public:
  // ---------------------------------------------------------------------------------------------------- //
  //                               Public type
  // ---------------------------------------------------------------------------------------------------- //
  using key_type = typename rhh_table_type::key_type;
  using mapped_value_type = typename rhh_table_type::mapped_value_type ;
  using value_type = typename rhh_table_type::value_type;
  using size_type = typename rhh_table_type::size_type;

  // TODO: implement const iterator and move to an another file
  class iterator
  {
   public:
    // iterator() {} // TODO: implement

    iterator(rhh_table_type* const table, const size_type hint_position) :
      m_table(table),
      m_current_position(m_table->find_valid_value_with_hint(hint_position))
    {  }

    iterator operator++()
    {
      iterator tmp = *this;
      next_valid_value();
      return tmp;
    }

    iterator operator++(int)
    {
      next_valid_value();
      return *this;
    }

    value_type& operator*()
    {
      return m_table->value_at_by_position(m_current_position);
    }

    value_type* operator->()
    {
      return &(m_table->value_at_by_position(m_current_position));
    }

    bool operator==(const iterator& rhs)
    {
      return equal(rhs);
    }

    bool operator!=(const iterator& rhs)
    {
      return !equal(rhs);;
    }

   private:

    bool equal(const iterator& rhs)
    {
      return m_table == rhs.m_table && m_current_position == rhs.m_current_position;
    }

    void next_valid_value()
    {
      m_current_position = m_table->find_valid_value_with_hint(m_current_position + 1);
    }

    rhh_table_type* m_table;
    size_type m_current_position;
  };

  // ---------------------------------------------------------------------------------------------------- //
  //                               Public function
  // ---------------------------------------------------------------------------------------------------- //

  // -------------------- Constructors -------------------- //
  unordered_map() = default;
  unordered_map(const self_type&) = default;
  unordered_map(self_type&&) = default;
  self_type& operator=(const self_type&) = default;
  self_type& operator=(self_type&&) = default;

//  unordered_map() :
//    m_table() {}
//
//  // Copy constructor
//  unordered_map(const self_type& other) :
//    m_table(other.m_table) { }
//
//  // Move constructor
//  unordered_map(self_type&& other) :
//    m_table(std::move(other.m_table))
//  { }
//
//  /// Copy and move assignments
//  /// When this function work as a move assignment,
//  /// a move constructor will be called when recieve "self_type rhs" instead of a copy constructor
//  /// due to the optimization of C++11 compilers
//  self_type& operator=(self_type rhs)
//  {
//    rhs.swap(*this);
//    return (*this);
//  }


  // -------------------- Capacity -------------------- //
  size_type size() const
  {
    return m_table.size();
  }

  size_type table_length() const
  {
    return m_table.table_length();
  }


  // -------------------- Iterator -------------------- //
  // TODO: implement cbegin and cend
  iterator begin()
  {
    return iterator(&m_table, 0);
  }

  iterator end()
  {
    return iterator(&m_table, table_length());
  }


  // -------------------- Element access -------------------- //
  mapped_value_type& operator[](const key_type& key)
  {
    size_type position = m_table.find(key);
    if (position >= m_table.table_length()) {
      const auto ret = m_table.insert(value_type(key, mapped_value_type()));
      position = ret.first;
    }
    return m_table.value_at_by_position(position).second;
  }

  mapped_value_type& operator[](key_type&& key)
  {
    return (*this)[key];
  }

  mapped_value_type& at(const key_type& key)
  {
    return m_table.value_at_by_key(key).second;
  }

  // -------------------- Element lookup -------------------- //
  iterator find(const key_type& key)
  {
    return iterator(&m_table, m_table.find(key));
  }

  // TODO: how to avoid duplicated implementations?
  const iterator find(const key_type& key) const
  {
    return iterator(&m_table, m_table.find(key));
  }

  size_type count(const key_type& key) const
  {
    size_type position = m_table.find(key);
    if (position < m_table.table_length()) {
      return 1;
    }
    return 0;
  }


  // -------------------- Modifiers -------------------- //
  std::pair<iterator, bool> insert(value_type value)
  {
    auto ret = m_table.insert(std::move(value));
    return std::make_pair(iterator(&m_table, ret.first), ret.second);
  }

  void clear()
  {
    m_table.clear();
  }

  void swap(self_type& other)
  {
    using std::swap;
    swap(m_table, other.m_table);
  }

  // -------------------- Statistic -------------------- //
  double load_factor() const
  {
    return m_table.load_factor();
  }


 private:

  // ---------------------------------------------------------------------------------------------------- //
  //                               Private variable
  // ---------------------------------------------------------------------------------------------------- //
  rhh_table_type m_table;
};

} // namespace rhh
#endif //UNORDERED_MAP_HPP
