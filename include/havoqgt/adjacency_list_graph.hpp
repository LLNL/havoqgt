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


#ifndef HAVOQGT_OMP_ADJACENCY_LIST_GRAPH_HPP_INCLUDED
#define HAVOQGT_OMP_ADJACENCY_LIST_GRAPH_HPP_INCLUDED

#include <boost/interprocess/containers/vector.hpp>
#include <stdint.h>
#include <utility>
#include <sstream>
#include <omp.hpp>
#include <vector>

namespace havoqgt { namespace omp {

namespace bip = boost::interprocess;

template<typename MemoryArena>
class adjacency_list_graph 
{
public:
  struct vertex_descriptor
  {
    uint64_t thread_id : 16;
    uint64_t local_id  : 48;
  };
  typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;

  template <typename T>
  class vertex_data
  {
  public:
    vertex_data(adjacency_list_graph<MemoryArena>* graph, T init) 
      : m_data(num_threads())
    {
      assert_sequential();
      #pragma omp parallel
      {
        assert_parallel();
        m_data[thread_num()].resize(graph->m_thread_data[thread_num()]->num_vertices(), init);
      }
    }

    T& operator[](vertex_descriptor vert) 
    {
      assert(vert.thread_id == thread_num());
      return m_data[thread_num()][vert.local_id];
    }

    uint64_t count_equal(T x) 
    {
      assert_sequential();
      uint64_t to_return = 0;
      #pragma omp parallel reduction(+: to_return)
      {
        for(size_t i=0; i<m_data[thread_num()].size(); ++i) 
        {
          if(m_data[thread_num()][i] == x) ++to_return;
        }
      }
      return to_return;
    }
  private:
    std::vector< std::vector< T > > m_data;
  };
  template <typename T> friend class vertex_data;

private:
  class thread_data {
    typedef typename MemoryArena::template allocator<vertex_descriptor>::type allocator_vd;
    typedef bip::vector< vertex_descriptor, allocator_vd >     adjacency_type;
    typedef typename MemoryArena::template allocator<adjacency_type>::type    allocator_adj_list;
    typedef bip::vector< adjacency_type, allocator_adj_list > adj_list_type;
  public:
    typedef typename adjacency_type::iterator edge_iterator;

    thread_data(typename MemoryArena::segment_manager* ptr_sm, size_t num_vertices)
      : m_adj_list(num_vertices,adjacency_type(ptr_sm),ptr_sm)
    { }
    
    void add_edge(edge_descriptor edge)
    {
      m_adj_list[edge.first.local_id].push_back(edge.second);
    }

    size_t num_edges(vertex_descriptor vert) const
    {
      return m_adj_list[vert.local_id].size();
    }

    edge_iterator edges_begin(vertex_descriptor vert)
    {
      return m_adj_list[vert.local_id].begin();
    }

    edge_iterator edges_end(vertex_descriptor vert)
    {
      return m_adj_list[vert.local_id].end(); 
    }

    size_t num_vertices() const { return m_adj_list.size(); }
  private:
    adj_list_type m_adj_list;
    char __pad[128 - sizeof(adj_list_type)];
  };

public:
  typedef typename thread_data::edge_iterator edge_iterator;

  adjacency_list_graph(MemoryArena& arena, const char* gname, size_t num_vertices)
  {
    assert_sequential();
    m_thread_data.resize(num_threads());
    size_t tl_num_vertices = num_vertices / num_threads();
    if(num_vertices % num_threads() > 0) ++tl_num_vertices;
    #pragma omp parallel
    {
      assert_parallel();
      std::stringstream ssname;
      ssname << gname << thread_num();
      m_thread_data[thread_num()] = 
        arena.get_sm()->template construct<thread_data>(ssname.str().c_str())(arena.get_sm(), tl_num_vertices);
    }
  }

  adjacency_list_graph(MemoryArena& arena, const char* gname)
  {
    assert_sequential();
    m_thread_data.resize(num_threads());
    #pragma omp parallel
    {
      assert_parallel();
      std::stringstream ssname;
      ssname << gname << thread_num();
      m_thread_data[thread_num()] = 
        arena.get_sm()->template find<thread_data>(ssname.str().c_str()).first;
    }
  }

  void add_edge(edge_descriptor edge)
  {
    assert(thread_num() == edge.first.thread_id);
    m_thread_data[thread_num()]->add_edge(edge);
  }

  size_t num_edges(vertex_descriptor vert) const
  {
    assert(thread_num() == vert.thread_id);
    return m_thread_data[thread_num()]->num_edges(vert);
  }

  edge_iterator edges_begin(vertex_descriptor vert)
  {
    assert(thread_num() == vert.thread_id);
    return m_thread_data[thread_num()]->edges_begin(vert);
  }

  edge_iterator edges_end(vertex_descriptor vert)
  {
    assert(thread_num() == vert.thread_id);
    return m_thread_data[thread_num()]->edges_end(vert);
  }

  template <typename T>
  vertex_data<T> init_vertex_data(T init = T()) 
  {
    assert_sequential();
    return vertex_data<T>(this, init);
  }

protected:
  std::vector<thread_data*> m_thread_data;
};


} } //end havoqgt::omp

#endif //HAVOQGT_OMP_ADJACENCY_LIST_GRAPH_HPP_INCLUDED
