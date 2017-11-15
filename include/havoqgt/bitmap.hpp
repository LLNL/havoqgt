//
// Created by Iwabuchi, Keita on 10/27/17.
//

#ifndef HAVOQGT_BITMAP_HPP
#define HAVOQGT_BITMAP_HPP

namespace havoqgt {

using bitmap_base_type = uint64_t;

/// exapmles
/// input 1 ~ 64 -> return 1;  input 65 ~ 128 -> return 2
inline constexpr size_t bitmap_size (const size_t size) noexcept
{
  return (size == 0) ? 0 :
         (size - 1ULL) / (sizeof(bitmap_base_type) * 8ULL) + 1ULL;
}

/// examples
/// input 0 ~ 63 -> return 0; input 64 ~ 127 -> return 1;
inline constexpr size_t bitmap_global_pos(const size_t pos) noexcept
{
  return (pos >> 6ULL);
}

inline constexpr size_t bitmap_local_pos(const size_t pos) noexcept
{
  return pos & 0x3FULL;
}

inline constexpr bool get_bit(const bitmap_base_type* const bitmap, const size_t pos) noexcept
{
  return bitmap[bitmap_global_pos(pos)] & (0x1ULL << bitmap_local_pos(pos));
}

inline void set_bit(bitmap_base_type* const bitmap, const size_t pos)
{
  bitmap[bitmap_global_pos(pos)] |= 0x1ULL << bitmap_local_pos(pos);
}

} // end namespace havoqgt
#endif //HAVOQGT_BITMAP_HPP
