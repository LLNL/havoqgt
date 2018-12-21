/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef BASELINE_MAP_VERTEX_ITERATOR_HPP
#define BASELINE_MAP_VERTEX_ITERATOR_HPP

#include <dynamic_graph_store/baseline/baseline.hpp>

namespace graphstore {

template <typename vertex_type,
          typename vertex_property_data_type,
          typename edge_property_data_type,
          typename segment_manager_type>
class graphstore_baseline_map<vertex_type,
                          vertex_property_data_type,
                          edge_property_data_type,
                          segment_manager_type>::vertex_iterator
{
 private:
  using graphstore_type       = graphstore_baseline_map<vertex_type,
                                                    vertex_property_data_type,
                                                    edge_property_data_type,
                                                    segment_manager_type>;
  using self_type             = graphstore_type::vertex_iterator;
  using table_vertex_iterator = typename graphstore_type::vertex_map_table_type::iterator;

 public:

  explicit vertex_iterator (table_vertex_iterator&& iterator) :
    m_iterator(iterator)
  { }


  void swap(self_type &other) noexcept
  {
    using std::swap;
    swap(m_iterator, other.m_iterator);
  }

  self_type &operator++ () // Pre-increment
  {
    find_next_value();
    return *this;
  }

  self_type operator++ (int) // Post-increment
  {
    self_type tmp(*this);
    find_next_value();
    return tmp;
  }

  // two-way comparison: v.begin() == v.cbegin() and vice versa
  bool operator == (const self_type &rhs) const
  {
    return is_equal(rhs);
  }

  bool operator != (const self_type &rhs) const
  {
    return !is_equal(rhs);
  }

  const vertex_type& source_vertex()
  {
    return m_iterator->first;
  }

  vertex_property_data_type& property_data()
  {
    return m_iterator->second.first;
  }


private:

  inline bool is_equal(const self_type &rhs) const
  {
    return (m_iterator == rhs.m_iterator);
  }

  inline void find_next_value()
  {
    ++m_iterator;
  }

  table_vertex_iterator m_iterator;
};

}
#endif // BASELINE_MAP_VERTEX_ITERATOR_HPP
