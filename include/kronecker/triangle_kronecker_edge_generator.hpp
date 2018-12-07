#pragma once

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>
#include <ygm/mpi.hpp>

template <typename C1, typename C2, typename Small_Vertex = uint32_t>
class triangle_kronecker_edge_generator {
  typedef uint64_t                      vertex_descriptor;
  typedef std::pair<uint64_t, uint64_t> value_type;
  typedef value_type                    edge_type;

  class input_iterator_type
      : public std::iterator<std::input_iterator_tag, edge_type, ptrdiff_t,
                             const edge_type*, const edge_type&> {
   public:
    input_iterator_type(triangle_kronecker_edge_generator* ptr_kron,
                        uint64_t                           count)
        : m_ptr_kron(ptr_kron), m_count(count), m_make_undirected(false) {
      if (m_count == 0) {
        get_next();
        m_count = 0;
      }
    }

    const edge_type&     operator*() const { return m_current; }
    input_iterator_type& operator++() {
      get_next();
      return *this;
    }

    input_iterator_type operator++(int) {
      input_iterator_type __tmp = *this;
      get_next();
      return __tmp;
    }

    inline edge_type* operator->() { return &m_current; }

    bool is_equal(const input_iterator_type& _x) const {
      return m_count == (_x.m_count);
    }

    /// Return true if x and y are both end or not end, or x and y are the same.
    friend bool operator==(const input_iterator_type& x,
                           const input_iterator_type& y) {
      return x.is_equal(y);
    }

    /// Return false if x and y are both end or not end, or x and y are the
    /// same.
    friend bool operator!=(const input_iterator_type& x,
                           const input_iterator_type& y) {
      return !x.is_equal(y);
    }

   private:
    input_iterator_type();

    void get_next() {
      if (m_ptr_kron->m_undirected && m_make_undirected) {
        std::swap(std::get<0>(m_current), std::get<1>(m_current));
        m_make_undirected = false;
      } else {
        do {
          ++m_count;
          if (*this == m_ptr_kron->end())
            break;  /// Ends even if last edge of Kronecker is skipped
          m_current = m_ptr_kron->generate_edge();
        } while (!valid_edge());
        m_make_undirected = true;
      }
    }

    bool valid_edge() {
      return true;  /// Can filter edges here
    }

   protected:
    triangle_kronecker_edge_generator* m_ptr_kron;
    uint64_t                           m_count;
    edge_type                          m_current;
    bool                               m_make_undirected;
  };

 public:
  triangle_kronecker_edge_generator(C1 graph1, C2 graph2,
                                    Small_Vertex num_vertices_graph1,
                                    Small_Vertex num_vertices_graph2,
                                    bool         undirected = false)
      : m_graph1(graph1),
        m_graph2(graph2),
        m_num_vertices_graph1(num_vertices_graph1),
        m_num_vertices_graph2(num_vertices_graph2),
        m_undirected(undirected),
        m_local_edge_count(
            graph1.size() * graph2.size() / ygm::comm_world().size() +
            (ygm::comm_world().rank() <
             (graph1.size() % ygm::comm_world().size()) * graph2.size())),
        m_graph1_itr(m_graph1.begin()),
        m_graph2_itr(m_graph2.begin()),
        m_graph1_pos(ygm::comm_world().rank()),
        m_graph1_adj_lists(num_vertices_graph1),
        m_graph2_adj_lists(num_vertices_graph2) {
    m_graph1_itr += ygm::comm_world().rank();

    Small_Vertex row, col;
    int64_t      val;
    for (auto& edge : m_graph1) {
      row = std::get<0>(edge);
      col = std::get<1>(edge);
      val = std::get<2>(edge);

      m_graph1_adj_lists[row][col] = val;
    }

    for (auto& edge : m_graph2) {
      row = std::get<0>(edge);
      col = std::get<1>(edge);
      val = std::get<2>(edge);

      m_graph2_adj_lists[row][col] = val;
    }
  }

  input_iterator_type begin() { return input_iterator_type(this, 0); }

  input_iterator_type end() {
    return input_iterator_type(this, m_local_edge_count);
  }

  uint64_t max_vertex_id() {
    return (m_num_vertices_graph1 * m_num_vertices_graph2) - 1;
  }

  size_t size() { return m_local_edge_count; }

  bool undirected() { return m_undirected; }

  uint64_t query(uint64_t row, uint64_t col) {
    Small_Vertex row1, row2, col1, col2;

    row1 = row / m_num_vertices_graph2;
    row2 = row % m_num_vertices_graph2;

    col1 = col / m_num_vertices_graph2;
    col2 = col % m_num_vertices_graph2;

    // std::cout << "Searching for (" << row1 << ", " << col1 << ")\n\t(" <<
    // row2 << ", " << col2 << ")" <<  std::endl;

    auto it1 = m_graph1_adj_lists[row1].find(col1);
    auto it2 = m_graph2_adj_lists[row2].find(col2);
    if ((it1 == m_graph1_adj_lists[row1].end()) ||
        (it2 == m_graph2_adj_lists[row2].end())) {
      return 0;
    } else
      return (*it1).second * (*it2).second;
  }

 private:
  edge_type generate_edge() {
    uint64_t     row, col;
    Small_Vertex row1, col1, row2, col2;

    row1 = std::get<0>(*m_graph1_itr);
    col1 = std::get<1>(*m_graph1_itr);
    row2 = std::get<0>(*m_graph2_itr);
    col2 = std::get<1>(*m_graph2_itr);

    row = row1 * m_num_vertices_graph2 + row2;
    col = col1 * m_num_vertices_graph2 + col2;

    m_graph2_itr++;
    if (m_graph2_itr == m_graph2.end()) {
      m_graph2_itr = m_graph2.begin();

      int increments = std::min(Small_Vertex(ygm::comm_world().size()),
                                m_num_vertices_graph1 - m_graph1_pos);
      m_graph1_itr += increments;
    }

    return std::make_pair(row, col);
  }

  C1           m_graph1;
  C2           m_graph2;
  uint64_t     m_local_edge_count;  /// Local edge count
  bool         m_undirected;
  Small_Vertex m_num_vertices_graph1;
  Small_Vertex m_num_vertices_graph2;
  Small_Vertex m_graph1_pos;

  std::vector<std::map<Small_Vertex, int32_t> > m_graph1_adj_lists;
  std::vector<std::map<Small_Vertex, int32_t> > m_graph2_adj_lists;

  typename C1::iterator m_graph1_itr;
  typename C2::iterator m_graph2_itr;
};
