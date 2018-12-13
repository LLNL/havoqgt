/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License. http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by the United States Government or Lawrence Livermore National Security, LLC.
 * The views and opinions of authors expressed herein do not necessarily state
 * or reflect those of the United States Government or Lawrence Livermore
 * National Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

#pragma once

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <havoqgt/detail/hash.hpp>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <tuple>
#include <vector>
#include <ygm/mpi.hpp>

namespace kronecker {

template <typename T, typename W>
T read_graph_file(std::string                       filename,
                  std::vector<std::tuple<T, T, W>>& edge_list) {
  std::ifstream filestream(filename);
  T             num_vertices;

  if (filestream.is_open()) {
    std::string line;
    if (!std::getline(filestream, line)) {
      std::cerr << "Empty file\n";
      exit(-1);
    }
    std::istringstream iss(line);
    iss >> num_vertices;
    if (line.find(" ") != line.npos) {
      std::cout << iss.str() << std::endl;
      std::cerr << "First line of input has too many values\n";
      exit(-1);
    }
    while (std::getline(filestream, line)) {
      std::istringstream iss2(line);
      T src, dest;
      W wgt;
      if (!(iss2 >> src >> dest >> wgt)) {
        std::cerr << "Malformed line in input\n";
        exit(-1);
      } else {
        edge_list.push_back(std::make_tuple(src, dest, wgt));
      }
    }
    filestream.close();
  } else {
    std::cerr << "Unable to open file " << filename << std::endl;
    exit(-1);
  }

  return num_vertices;
}

}  // end namespace kronecker

namespace havoqgt {

template <typename edge_data_type = uint8_t, typename Small_Index = uint32_t,
          typename C1 =
              std::vector<std::tuple<Small_Index, Small_Index, edge_data_type>>,
          typename C2 =
              std::vector<std::tuple<Small_Index, Small_Index, edge_data_type>>>
class kronecker_edge_generator {
 public:
  typedef uint64_t                                       vertex_descriptor;
  typedef std::tuple<uint64_t, uint64_t, edge_data_type> value_type;
  typedef value_type                                     edge_type;
  typedef edge_data_type                                 edge_data_value_type;

  class input_iterator_type
      : public std::iterator<std::input_iterator_tag, edge_type, ptrdiff_t,
                             const edge_type*, const edge_type&> {
   public:
    input_iterator_type(kronecker_edge_generator* ptr_kron, uint64_t count)
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

    kronecker_edge_generator* m_ptr_kron;
    uint64_t                  m_count;
    edge_type                 m_current;
    bool                      m_make_undirected;
  };

 public:
  kronecker_edge_generator(C1 graph1, C2 graph2,
                           Small_Index num_vertices_graph1,
                           Small_Index num_vertices_graph2,
                           bool scramble = false, bool undirected = false)
      : m_graph1(graph1),
        m_graph2(graph2),
        m_num_vertices_graph1(num_vertices_graph1),
        m_num_vertices_graph2(num_vertices_graph2),
        m_scramble(scramble),
        m_undirected(undirected),
        m_local_edge_count(
            graph1.size() * graph2.size() / ygm::comm_world().size() +
            (ygm::comm_world().rank() <
             (graph1.size() % ygm::comm_world().size()) * graph2.size())),
        m_graph1_itr(m_graph1.begin()),
        m_graph2_itr(m_graph2.begin()),
        m_graph1_pos(ygm::comm_world().rank()),
        m_has_edge_data(true) {
    m_vertex_scale =
        (uint64_t)ceil(log2(m_num_vertices_graph1 * m_num_vertices_graph2));
    m_graph1_itr += ygm::comm_world().rank();
  }

  kronecker_edge_generator(std::string filename1, std::string filename2,
                           bool scramble = false, bool undirected = false)
      : m_scramble(scramble),
        m_undirected(undirected),
        m_graph1_pos(ygm::comm_world().rank()),
        m_has_edge_data(true) {
    m_num_vertices_graph1 = kronecker::read_graph_file(filename1, m_graph1);
    m_num_vertices_graph2 = kronecker::read_graph_file(filename2, m_graph2);
    m_local_edge_count =
        m_graph1.size() * m_graph2.size() / ygm::comm_world().size() +
        (ygm::comm_world().rank() <
         (m_graph1.size() % ygm::comm_world().size()) * m_graph2.size());
    m_graph1_itr = m_graph1.begin();
    m_graph2_itr = m_graph2.begin();
    m_graph1_itr += ygm::comm_world().rank();

    m_vertex_scale =
        (uint64_t)ceil(log2(m_num_vertices_graph1 * m_num_vertices_graph2));
    std::cout << "Vertex Scale: " << m_vertex_scale << std::endl;
  }

  input_iterator_type begin() {
    // Reset iterators to prepare to read
    m_graph1_itr = m_graph1.begin();
    m_graph2_itr = m_graph2.begin();
    m_graph1_pos = ygm::comm_world().rank();
    m_graph1_itr += ygm::comm_world().rank();
    return input_iterator_type(this, 0);
  }

  input_iterator_type end() {
    return input_iterator_type(this, m_local_edge_count);
  }

  uint64_t max_vertex_id() {
    return (m_num_vertices_graph1 * m_num_vertices_graph2) - 1;
  }

  size_t size() { return m_local_edge_count; }

  bool undirected() { return m_undirected; }

  bool has_edge_data() { return m_has_edge_data; }

 private:
  edge_type generate_edge() {
    uint64_t       row, col;
    Small_Index    row1, col1, row2, col2;
    edge_data_type val1, val2, val;

    row1 = std::get<0>(*m_graph1_itr);
    col1 = std::get<1>(*m_graph1_itr);
    val1 = std::get<2>(*m_graph1_itr);
    row2 = std::get<0>(*m_graph2_itr);
    col2 = std::get<1>(*m_graph2_itr);
    val2 = std::get<2>(*m_graph2_itr);

    row = row1 * m_num_vertices_graph2 + row2;
    col = col1 * m_num_vertices_graph2 + col2;
    val = val1 * val2;

    m_graph2_itr++;
    if (m_graph2_itr == m_graph2.end()) {
      m_graph2_itr = m_graph2.begin();

      int increments = std::min(Small_Index(ygm::comm_world().size()),
                                m_num_vertices_graph1 - m_graph1_pos);
      m_graph1_itr += increments;
    }

    if (m_scramble) {
      row = havoqgt::detail::hash_nbits(row, m_vertex_scale);
      col = havoqgt::detail::hash_nbits(col, m_vertex_scale);
    }

    return std::make_tuple(row, col, val);
  }

  C1          m_graph1;
  C2          m_graph2;
  uint64_t    m_local_edge_count;  /// Local edge count
  bool        m_scramble;
  bool        m_undirected;
  bool        m_has_edge_data;
  Small_Index m_num_vertices_graph1;
  Small_Index m_num_vertices_graph2;
  Small_Index m_graph1_pos;
  uint64_t    m_vertex_scale;

  typename C1::iterator m_graph1_itr;
  typename C2::iterator m_graph2_itr;
};

}  // end namespace havoqgt
