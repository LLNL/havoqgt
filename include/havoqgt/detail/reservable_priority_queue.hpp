// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_DETAIL_RESERVABLE_PRIORITY_QUEUE_HPP_INCLUDED
#define HAVOQGT_DETAIL_RESERVABLE_PRIORITY_QUEUE_HPP_INCLUDED


#include <vector>
#include <queue>

namespace havoqgt { namespace detail {

template <typename  T, typename CONT=std::vector<T>, typename COMP=std::less<T> >
class reservable_priority_queue: public std::priority_queue<T, CONT, COMP>
{
public:
    typedef typename std::priority_queue<T>::size_type size_type;
    reservable_priority_queue(size_type capacity = 0) { reserve(capacity); };
    void reserve(size_type capacity) { this->c.reserve(capacity); } 
    size_type capacity() const { return this->c.capacity(); } 
};

} } //end namespace havoqgt { namespace detail {

#endif //end HAVOQGT_DETAIL_RESERVABLE_PRIORITY_QUEUE_HPP_INCLUDED
