// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_UPPER_TRIANGLE_EDGE_GENERATOR_INCLUDED
#define HAVOQGT_UPPER_TRIANGLE_EDGE_GENERATOR_INCLUDED

#include <boost/random.hpp>
#include <havoqgt/detail/hash.hpp>

#include <utility>
#include <stdint.h>

namespace havoqgt {

class upper_triangle_edge_generator {
 public:
  // typedef uint64_t                      vertex_descriptor;
  typedef std::pair<uint64_t, uint64_t> edge_type;


  ///
  /// InputIterator class for upper_triangle_edge_generator

  class input_iterator_type : public std::iterator<std::input_iterator_tag,
      edge_type, ptrdiff_t, const edge_type*, const edge_type&> {
   private:
    void get_next() {
      if (m_count == m_num_edges) {
        m_current = std::make_pair(0, 0);
        assert(false);
        return;
      }

      if (m_undirected && m_make_undirected) {
        std::swap(m_current.first, m_current.second);
        m_make_undirected = false;
      } else {
        m_current = generate_edge();
        m_make_undirected = true;
        ++m_count;
      }
    };  // get_next

   public:
    input_iterator_type(uint64_t max_vertex, uint64_t num_edges,
          bool undirected)
      : m_max_vertex(max_vertex)
      , m_num_edges(num_edges)
      , m_undirected(undirected)
      , m_count(num_edges) {
      m_current = std::make_pair(0, 0);
    }

    input_iterator_type(uint64_t max_vertex, uint64_t num_edges,
          bool undirected, uint64_t i, uint64_t j)
      : m_max_vertex(max_vertex)
      , m_num_edges(num_edges)
      , m_undirected(undirected)
      , m_i(i)
      , m_j(j) {
      // Initilize the first current value
      get_next();
    }

    const edge_type& operator*() const { return m_current; }
    // const uint64_t* operator->() const { return &(operator*()); }
    input_iterator_type& operator++() {
      get_next();
      return *this;
    }

    input_iterator_type operator++(int) {
      input_iterator_type __tmp = *this;
      get_next();
      return __tmp;
    }

    edge_type *operator->() {
      return &m_current;
    }

    bool is_equal(const input_iterator_type& _x) const {
      return m_count == (_x.m_count);
    }

    // Return true if x and y are both end or not end, or x and y are the same.
    friend bool
    operator==(const input_iterator_type& x, const input_iterator_type& y)
    { return x.is_equal(y); }

    // Return false if x and y are both end or not end, or x and y are the same.
    friend bool
    operator!=(const input_iterator_type& x, const input_iterator_type& y)
    { return !x.is_equal(y); }

   private:
    input_iterator_type();

    edge_type generate_edge() {
      if (m_j < m_max_vertex) {
        m_j++;
      } else {
        m_i++;
        m_j = m_i;
      }
      assert(m_i <= m_max_vertex);
      assert(m_j <= m_max_vertex);
      return std::make_pair(m_j + 1, m_i + 1);
    }  // generate_edge

    const uint64_t m_max_vertex, m_num_edges;
    const bool m_undirected;
    bool m_make_undirected {false};
    uint64_t m_count {0}, m_i, m_j;
    edge_type m_current;
  };  // class input_iterator_type

  uint64_t find_max_vertex(uint64_t num_edges) {
    // (1/2)(x)(x+1) <= num_edges.
    uint64_t x = 0, y = 0;
    num_edges *= 2;
    while (y < num_edges) {
      x++;
      y = x * (x + 1);
    }

    uint64_t count = 0;
    for (uint64_t i = 0; i < x; i++) {
      for (uint64_t j = i; j < x; j++) {
        count++;
      }
    }
    assert(count * 2 >= num_edges);

    return x + 1; // +1 because we skip 0 vertex
  };

  upper_triangle_edge_generator(uint64_t total_edges, int rank, int size,
        bool undirected)
      : m_total_edges(total_edges)
      , m_undirected(undirected)
      , m_max_vertex(find_max_vertex(total_edges)) {
    // Determine the max_vertex id
    m_lower_edge_id = (total_edges / size) * rank;
    m_upper_edge_id = (total_edges / size) * (rank + 1);
    m_upper_edge_id = std::min(total_edges, m_upper_edge_id);

    uint64_t start_count = (total_edges / size) * rank;

    uint64_t count = 0;
    for (m_i = 0; count < start_count && m_i < m_max_vertex; m_i++) {
      for (m_j = m_i; count < start_count && m_j < m_max_vertex; m_j++) {
        count++;
      }
    }
    sanity_check();
  }

  void sanity_check() {
    auto itr1 = begin();
    auto itr2 = begin();
    auto itr_end = end();
    while (itr1 != itr_end) {
      assert(itr1 == itr2);
      assert(*itr1 == *itr2);
      itr1++;
      itr2++;
    }

  }


  /// Returns the begin of the input iterator
  input_iterator_type begin() {
    const uint64_t num_edges = m_upper_edge_id - m_lower_edge_id;
    return input_iterator_type(m_max_vertex, num_edges, m_undirected, m_i, m_j);
  }

  /// Returns the end of the input iterator
  input_iterator_type end() {
    const uint64_t num_edges = m_upper_edge_id - m_lower_edge_id;
    return input_iterator_type(m_max_vertex, num_edges, m_undirected);
  }

  uint64_t max_vertex_id() {
    return m_max_vertex;
  }

  size_t size() {
    const uint64_t num_edges = m_upper_edge_id - m_lower_edge_id;
    return num_edges;
  }

  uint64_t m_upper_edge_id, m_lower_edge_id;
  const uint64_t m_total_edges;
  const bool m_undirected;

  const uint64_t m_max_vertex;
  uint64_t m_i, m_j;

};



}  // name space havoqgt

#endif  // HAVOQGT_UPPER_TRIANGLE_EDGE_GENERATOR_INCLUDED

