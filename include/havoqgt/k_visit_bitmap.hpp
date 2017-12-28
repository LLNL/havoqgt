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
using k_visit_bitmap = detail::static_bitmap<_num_sources>;

} // namespace havoqgt
#endif //HAVOQGT_K_VISIT_BITMAP_HPP
