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
  return bitmap[detail::bitmap_global_pos<bitmap_base_type>(pos)]
         & (0x1ULL << detail::bitmap_local_pos<bitmap_base_type>(pos));
}

template <typename bitmap_base_type>
inline void set_bit(bitmap_base_type* const bitmap, const size_t pos)
{
  bitmap[detail::bitmap_global_pos<bitmap_base_type>(pos)]
    |= (0x1ULL << detail::bitmap_local_pos<bitmap_base_type>(pos));
}

} // end namespace havoqgt
#endif //HAVOQGT_BITMAP_HPP
