//
// Created by Iwabuchi, Keita on 12/14/17.
//

#ifndef HAVOQGT_K_VISIT_BITMAP_HPP
#define HAVOQGT_K_VISIT_BITMAP_HPP

#include <havoqgt/detail/bitmap.hpp>
#include <cstdio>


namespace havoqgt
{

template <size_t _num_sources>
struct k_visit_bitmap
{

  static constexpr size_t num_sources = _num_sources;
  using bitmap_base_t = havoqgt::detail::bitmap_base_type<num_sources>;
  static constexpr size_t size = havoqgt::detail::bitmap_size<bitmap_base_t>(num_sources);
  bitmap_base_t bitmap[size] = {0};

  inline bool equal(const k_visit_bitmap &rhs) const
  {
    for (size_t i = 0; i < size; ++i) {
      if (bitmap[i] != rhs.bitmap[i]) return false;
    }
    return true;
  }

  const k_visit_bitmap &operator|=(const k_visit_bitmap &rhs)
  {
    for (size_t i = 0; i < size; ++i) {
      bitmap[i] |= rhs.bitmap[i];
    }
    return (*this);
  }

  bool get(const size_t pos) const
  {
    return havoqgt::detail::get_bit(bitmap, pos);
  }

  void set(const size_t pos)
  {
    havoqgt::detail::set_bit(bitmap, pos);
  }

} __attribute__ ((packed));

template <size_t num_sources>
inline bool operator==(const k_visit_bitmap<num_sources> &lhs, const k_visit_bitmap<num_sources> &rhs)
{
  return lhs.equal(rhs);
}

template <size_t num_sources>
inline bool operator!=(const k_visit_bitmap<num_sources> &lhs, const k_visit_bitmap<num_sources> &rhs)
{
  return !(lhs.equal(rhs));
}

template <size_t num_sources>
bool is_contain(const k_visit_bitmap<num_sources> &lhs, const k_visit_bitmap<num_sources> &rhs)
{
  for (size_t i = 0; i < k_visit_bitmap<num_sources>::size; ++i) {
    if ((lhs.bitmap[i] | rhs.bitmap[i]) > lhs.bitmap[i]) return false;
  }
  return true;
}

} // namespace havoqgt
#endif //HAVOQGT_K_VISIT_BITMAP_HPP
