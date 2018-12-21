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

#ifndef BASELINE_MAP_HPP
#define BASELINE_MAP_HPP

#include <boost/unordered_map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <dynamic_graph_store/graphstore_utilities.hpp>


namespace graphstore {

template<typename map_t>
void bucket_size_histogram(const map_t& map,
                           std::vector<size_t>& histogram)
{
  histogram.clear();
  for (size_t i = 0; i < map.bucket_count(); ++i) {
    size_t s = map.bucket_size(i);
    if (histogram.size() < s + 1) histogram.resize(s+1, 0);
    ++histogram[s];
  }
}

}


namespace graphstore {
template <typename _vertex_type,
          typename _vertex_property_data_type,
          typename _edge_property_data_type,
          typename _segment_manager_type>
class graphstore_baseline_map
{
public:
  using vertex_type                 = _vertex_type;
  using vertex_property_data_type   = _vertex_property_data_type;
  using edge_property_data_type     = _edge_property_data_type;
  using segment_manager_type        = _segment_manager_type;

  /// iterator
  class vertex_iterator;
  class adjacent_edge_iterator;

private:
  using size_type             = size_t;
  template <typename t1, typename t2>
  using pair_type             = graphstore::utility::packed_pair<t1, t2>;
  using allocator_type        = boost::interprocess::allocator<void, segment_manager_type>;


  /// adjacency list
  using edge_map_value_type     = edge_property_data_type;
  using edge_map_element_type   = std::pair<const vertex_type, edge_map_value_type>;
  using edge_map_allocator_type = boost::interprocess::allocator<edge_map_element_type, segment_manager_type>;
  using edge_map_table_type     = boost::unordered_multimap<vertex_type,
                                                            edge_map_value_type,
                                                            boost::hash<vertex_type>,
                                                            std::equal_to<vertex_type>,
                                                            edge_map_allocator_type>;

  /// vertex table
  using vertex_map_value_type        = std::pair<vertex_property_data_type, edge_map_table_type>;
  using vertex_map_element_type      = std::pair<const vertex_type, vertex_map_value_type>;
  using vertex_map_allocator_type    = boost::interprocess::allocator<vertex_map_element_type, segment_manager_type>;
  using vertex_map_table_type        = boost::unordered_map<vertex_type,
                                                            vertex_map_value_type,
                                                            boost::hash<vertex_type>,
                                                            std::equal_to<vertex_type>,
                                                            vertex_map_allocator_type>;

public:

  explicit graphstore_baseline_map(segment_manager_type* segment_manager) :
    m_allocator(segment_manager),
    m_map_table(m_allocator),
    m_num_edges(0),
    m_segment_manager(segment_manager)
  {}


  /// -------- Lookup -------- ///
  std::pair<vertex_iterator, vertex_iterator> vertices()
  {
    return std::make_pair(vertices_begin(), vertices_end());
  }

  vertex_iterator vertices_begin()
  {
    return vertex_iterator(m_map_table.begin());
  }

  vertex_iterator vertices_end()
  {
    return vertex_iterator(m_map_table.end());
  }

  adjacent_edge_iterator adjacent_edge_begin(const vertex_type& src_vrt)
  {
    const auto itr = m_map_table.find(src_vrt);
    edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*itr));
    return adjacent_edge_iterator(edge_map.begin());
  }

  adjacent_edge_iterator adjacent_edge_end(const vertex_type& src_vrt)
  {
    const auto itr = m_map_table.find(src_vrt);
    edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*itr));
    return adjacent_edge_iterator(edge_map.end());
  }

  // Returns the degree of a vertex.
  size_type degree(const vertex_type& vertex)
  {
    const auto itr = m_map_table.find(vertex);
    if (itr == m_map_table.end()) {
      return 0;
    } else {
      edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*itr));;
      return edge_map.size();
    }
  }

  /// -------- Modifiers ------- ////
  bool insert_edge(const vertex_type& src, const vertex_type& trg, const edge_property_data_type& edge_property)
  {

    auto vmp_itr = m_map_table.find(src);
    if (vmp_itr == m_map_table.end()) { // new vertex
      insert_new_vertex(src, trg, vertex_property_data_type(), edge_property);
    } else {
      edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));

      /// --- find duplicated edge --- ///
      auto emp_itr = edge_map.find(trg);
      if (emp_itr != edge_map.end()) {
        return false; /// found duplicated edge
      }

      /// --- since there is no empty spase, add it at end --- ///
      edge_map.insert(edge_map_element_type(trg, edge_property));
    }

    ++m_num_edges;
    return true;
  }

  /// allows duplicated insertion
  void insert_edge_dup(const vertex_type& src, const vertex_type& trg, const edge_property_data_type& edge_property)
  {

    auto vmp_itr = m_map_table.find(src);
    if (vmp_itr == m_map_table.end()) { // new vertex
      insert_new_vertex(src, trg, vertex_property_data_type(), edge_property);
    } else {
      edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
      edge_map.insert(edge_map_element_type(trg, edge_property));
    }

    ++m_num_edges;
  }

  bool erase_edge(const vertex_type& src, const vertex_type& trg)
  {
    auto vmp_itr = m_map_table.find(src);
    if (vmp_itr == m_map_table.end()) {
      return false;
    }

    edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
    auto edg_itr = edge_map.find(trg);
    if (edg_itr != edge_map.end()) {
      edge_map.erase(edg_itr);
      --m_num_edges;

      if (edge_map.size() == 0) {
        m_map_table.erase(vmp_itr);
      }
      return true;
    }

    return false;
  }

  /// this fuction deletes all duplicated edges
  size_type erase_edge_dup(const vertex_type& src, const vertex_type& trg)
  {
    auto vmp_itr = m_map_table.find(src);
    if (vmp_itr == m_map_table.end()) {
      return false;
    }

    /// erase edges
    edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
    const size_t num_erased = edge_map.erase(trg);
    m_num_edges -= num_erased;

    if (edge_map.size() == 0) {
      m_map_table.erase(vmp_itr);
    }

    return num_erased;
  }

  inline vertex_property_data_type& vertex_property_data(const vertex_type& vertex)
  {
    auto vmp_itr = m_map_table.find(vertex);
    return std::get<vtx_prp>(std::get<vmp_val>(*vmp_itr));
  }

  edge_property_data_type& edge_property_data(const vertex_type& src, const vertex_type& trg) {
    auto vmp_itr = m_map_table.find(src);
//    if (vmp_itr == m_map_table.end()) {
//      std::cout << src << ":" << trg;
//      std::cout << " Bad edge: vertex lookup failed. ";
//      ::exit(1);
//    }

    edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
    auto edg_itr = edge_map.find(trg);
//    if (edg_itr == edge_map.end()) {
//      std::cout << src << ":" << trg;
//      std::cout << " Bad edge: could not find destination. ";
//      ::exit(1);
//    }

    return std::get<edg_prp>(*edg_itr);
  }


  void opt()
  { }

  inline size_type num_edges() const
  {
    return m_num_edges;
  }

  segment_manager_type* get_segment_manager() const
  {
      return m_segment_manager;
  }

  void clear()
  { }

  void print_status(const int level = 0) const
  {
    size_t cnt = 0;
    double ave_num_buckets = 0;
    for (auto vmp_itr = m_map_table.begin(), end = m_map_table.end(); vmp_itr != end; ++vmp_itr) {
      const edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
      ave_num_buckets += edge_map.bucket_count();
      ++cnt;
    }
    std::cout << "edge-maps' average num buckets: " << ave_num_buckets / cnt << std::endl;

    std::cout << "vertex-map's num buckets: " << m_map_table.bucket_count() << std::endl;
    std::cout << "vertex-map's size: " << m_map_table.size() << std::endl;

    if (level < 1) return;

    /// --- print vertex table's detailed status --- ///
    {
      std::vector<size_t> cnt;
      bucket_size_histogram(m_map_table, cnt);
      for (size_t i = 0; i < cnt.size(); ++i) {
        std::cout << cnt[i] << "[" << i << "] ";
      }
      std::cout << std::endl;
    }

    /// --- print edge table's detailed status --- ///
    {
      std::vector<size_t> cnt;
      for (auto vmp_itr = m_map_table.begin(), end = m_map_table.end(); vmp_itr != end; ++vmp_itr) {
        const edge_map_table_type& edge_map = std::get<edg_map>(std::get<vmp_val>(*vmp_itr));
        std::vector<size_t> wk;
        bucket_size_histogram(edge_map, wk);
        if (cnt.size() < wk.size()) cnt.resize(wk.size(), 0);
        for (size_t i = 0; i < cnt.size(); ++i) {
          cnt[i] += wk[i];
        }
      }
      for (size_t i = 0; i < cnt.size(); ++i) {
        std::cout << cnt[i] << "[" << i << "] ";
      }
      std::cout << std::endl;
    }

  }

  void fprint_all_elements(std::ofstream& of)
  {
    for (auto vitr = vertices_begin(), end = vertices_end();
         vitr != end;
         ++vitr) {
      for (auto eitr = adjacent_edge_begin(vitr.source_vertex()), end = adjacent_edge_end(vitr.source_vertex());
           eitr != end;
           ++eitr) {
        of << vitr.source_vertex() << " " << eitr.target_vertex() << " " << eitr.property_data() << "\n";
      }
    }
  }


 private:
  /// vert_map key, val (prop, edge_map)
  /// edge_map key, val
  enum {
    vtx_src = 0,
    vmp_val = 1,
    vtx_prp = 0,
    edg_map = 1,
    edg_vrt = 0,
    edg_prp = 1
  };

  void insert_new_vertex(const vertex_type& src,
                         const vertex_type& trg,
                         const vertex_property_data_type& vertex_property,
                         const edge_property_data_type& edge_property)
  {
    edge_map_table_type edge_map(m_allocator);
    edge_map.insert(edge_map_element_type(trg, edge_property));
    vertex_map_value_type vmp_value(vertex_property, edge_map);
    m_map_table.insert(vertex_map_element_type(src, vmp_value));
  }

  allocator_type m_allocator;
  vertex_map_table_type m_map_table;
  size_type m_num_edges;
  segment_manager_type* m_segment_manager;
};

}

#include <dynamic_graph_store/baseline/baseline_map_vertex_iterator.hpp>
#include <dynamic_graph_store/baseline/baseline_map_adjacent_edge_iterator.hpp>

#endif // BASELINE_MAP_HPP
