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
