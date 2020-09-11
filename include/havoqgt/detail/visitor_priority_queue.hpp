// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_DETAIL_VISITOR_PRIORITY_QUEUE_HPP_INCLUDED
#define HAVOQGT_DETAIL_VISITOR_PRIORITY_QUEUE_HPP_INCLUDED


#include <deque>
#include <queue>

namespace havoqgt { namespace detail {

  template <typename Visitor>
  class visitor_priority_queue
  {

  protected:
    std::priority_queue< Visitor, std::deque<Visitor>, 
                                 std::greater<Visitor> > m_data;
  public:
    visitor_priority_queue() { }

    bool push(Visitor const & task)
    {
      m_data.push(task);
      return true;
    }

    void pop()
    {
      m_data.pop();
    }

    Visitor const & top() //const
    {
      return m_data.top();
    }

    size_t size() const
    {
      return m_data.size();;
    }

    bool empty() const
    {
      return m_data.empty();
    }

    void clear()
    {
      m_data.clear();
    }
  };


} } //end namespace havoqgt { namespace detail {

#endif //end HAVOQGT_DETAIL_VISITOR_PRIORITY_QUEUE_HPP_INCLUDED
