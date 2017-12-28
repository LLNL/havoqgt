//
// Created by Iwabuchi, Keita on 10/27/17.
//

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
  bitmap[bitmap_global_pos<bitmap_base_type>(pos)]
    |= (0x1ULL << bitmap_local_pos<bitmap_base_type>(pos));
}

template <size_t _num_bit>
class static_bitmap
{
 public:
  static constexpr size_t num_bit = _num_bit;
  using bitmap_base_t = havoqgt::detail::bitmap_base_type<num_bit>;
  using const_bitmap_base_t = const havoqgt::detail::bitmap_base_type<num_bit>;
  static constexpr size_t size = havoqgt::detail::bitmap_size<bitmap_base_t>(num_bit);

  const static_bitmap &operator|=(const static_bitmap &rhs)
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
