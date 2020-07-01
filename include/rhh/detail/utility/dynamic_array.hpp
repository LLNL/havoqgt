//
// Created by Iwabuchi, Keita on 2019-02-12.
//

#ifndef HAVOQGT_RHH_DETAIL_UTILITY_DYNAMIC_ARRAY_HPP
#define HAVOQGT_RHH_DETAIL_UTILITY_DYNAMIC_ARRAY_HPP

#include <memory>
#include <boost/container/vector.hpp>
#include <boost/container/scoped_allocator.hpp>

namespace rhh {
namespace detail {
namespace utility {

template <typename T, typename allocator>
using dynamic_array = boost::container::vector<T, typename std::allocator_traits<allocator>::template rebind_alloc<T>>;

template <typename allocator>
using scp_alloc = boost::container::scoped_allocator_adaptor<allocator>;

} // namespace rhh
} // namespace detail
} // namespace utility
#endif //HAVOQGT_RHH_DETAIL_UTILITY_DYNAMIC_ARRAY_HPP
