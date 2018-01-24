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

#ifndef HAVOQGT_BITMAP_HPP
#define HAVOQGT_BITMAP_HPP

#include <type_traits>

namespace havoqgt {

namespace detail {
/// \brief calculate log2
/// n must be larger than 0
/// \param n
/// \return log2 of n
inline constexpr size_t cal_log2(const size_t n)
{
  return (n < 2) ? 0 : 1 + cal_log2(n / 2);
}

/// examples
/// input 0 ~ 63 -> return 0; input 64 ~ 127 -> return 1;
template <typename bitmap_base_type>
inline constexpr size_t bitmap_global_pos(const size_t pos)
{
  return (pos >> cal_log2(sizeof(bitmap_base_type) * 8ULL));
}

template <typename bitmap_base_type>
inline constexpr size_t bitmap_local_pos(const size_t pos)
{
  return pos & (sizeof(bitmap_base_type) * 8 - 1);
}


/// \brief define a proper bitmap underling type based on the number of bits needed
/// \example
/// 0 ~ 8   bits -> uint8_t
/// 9 ~ 16  bits -> uint16_t
/// 17 ~ 32 bits -> uint32_t
/// 33~     bits -> uint64_t
template <size_t num_bits>
using bitmap_base_type =
typename std::conditional<num_bits <= 8,
                          uint8_t,
                          typename std::conditional<num_bits <= 16,
                                                   uint16_t,
                                                   typename std::conditional<num_bits <= 32,
                                                                            uint32_t,
                                                                            uint64_t
                                                   >::type
                          >::type
>::type;

/// exapmles: bitmap_base_type = uint64_t
/// input 1 ~ 64 -> return 1;  input 65 ~ 128 -> return 2
template <typename bitmap_base_type>
inline constexpr size_t bitmap_size(const size_t size)
{
  return (size == 0) ? 0 : (size - 1ULL) / (sizeof(bitmap_base_type) * 8ULL) + 1ULL;
}

template <typename bitmap_base_type>
inline constexpr bool get_bit(const bitmap_base_type* const bitmap, const size_t pos)
{
  return bitmap[bitmap_global_pos<bitmap_base_type>(pos)]
         & (0x1ULL << bitmap_local_pos<bitmap_base_type>(pos));
}

template <typename bitmap_base_type>
inline void set_bit(bitmap_base_type* const bitmap, const size_t pos)
{
  bitmap[bitmap_global_pos<bitmap_base_type>(pos)] |= (0x1ULL << bitmap_local_pos<bitmap_base_type>(pos));
}

template <size_t _num_bit>
class static_bitmap
{
 public:
  static constexpr size_t num_bit = _num_bit;
  using bitmap_base_t = havoqgt::detail::bitmap_base_type<num_bit>;
  using const_bitmap_base_t = const havoqgt::detail::bitmap_base_type<num_bit>;
  static constexpr size_t size = havoqgt::detail::bitmap_size<bitmap_base_t>(num_bit);

  static_bitmap<_num_bit>() = default;
  static_bitmap<_num_bit>(const static_bitmap<_num_bit>&) = default;
  static_bitmap<_num_bit>(static_bitmap<_num_bit>&&) = default;

  static_bitmap<_num_bit> &operator=(const static_bitmap &rhs)
  {
    for (size_t i = 0; i < size; ++i) {
      m_map[i] = rhs.m_map[i];
    }
    return (*this);
  }

  const static_bitmap<_num_bit> &operator|=(const static_bitmap &rhs)
  {
    for (size_t i = 0; i < size; ++i) {
      m_map[i] |= rhs.m_map[i];
    }
    return (*this);
  }

  bool get(const size_t pos) const
  {
    return havoqgt::detail::get_bit(m_map, pos);
  }

  void set(const size_t pos)
  {
    havoqgt::detail::set_bit(m_map, pos);
  }

  bitmap_base_t* get_map()
  {
    return m_map;
  }

  const_bitmap_base_t* get_map() const
  {
    return m_map;
  }

  bool has_bit() const
  {
    for (size_t i = 0; i < size; ++i) {
      if (m_map[i]) return true;
    }
    return false;
  }

 private:
  bitmap_base_t m_map[size] = {0};
} __attribute__ ((packed));

template <size_t num_bit>
inline bool operator==(const static_bitmap<num_bit> &lhs, const static_bitmap<num_bit> &rhs)
{
  typename static_bitmap<num_bit>::const_bitmap_base_t *const lmap = lhs.get_map();
  typename static_bitmap<num_bit>::const_bitmap_base_t *const rmap = rhs.get_map();
  for (size_t i = 0; i < static_bitmap<num_bit>::size; ++i) {
    if (lmap[i] != rmap[i]) return false;
  }
  return true;
}

template <size_t num_bit>
inline bool operator!=(const static_bitmap<num_bit> &lhs, const static_bitmap<num_bit> &rhs)
{
  return !(lhs == rhs);
}

} // namespace detail

} // namespace havoqgt
#endif //HAVOQGT_BITMAP_HPP
