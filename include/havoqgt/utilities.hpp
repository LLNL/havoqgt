// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_UTILIES_HPP
#define HAVOQGT_UTILIES_HPP

#include <memory>

namespace havoqgt {

/// Frees the container of edges
template <typename Container>
void free_edge_container(Container &edges) {};

template<>
void free_edge_container<std::vector<std::pair<uint64_t, uint64_t> > >(std::vector<std::pair<uint64_t, uint64_t> > &edges){
  std::vector< std::pair<uint64_t, uint64_t> >empty(0);
  edges.swap(empty);
};

template <typename Allocator, typename T2>
using other_allocator = typename std::allocator_traits<Allocator>::template rebind_alloc<T2>;

#endif // HAVOQGT_UTILIES_HPP

} // namespace havoqgt
