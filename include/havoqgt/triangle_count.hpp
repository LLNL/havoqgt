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

#ifndef HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
#define HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>


namespace havoqgt { namespace mpi {


template <typename Visitor>
class triangle_priority_queue
{

protected:
  std::priority_queue< Visitor, std::deque<Visitor>, 
                               std::greater<Visitor> > m_data;
public:
  triangle_priority_queue() { }

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


template<typename Graph>
class triangle_count_visitor {
public:
  typedef triangle_count_visitor<Graph>    my_type;
  typedef typename Graph::vertex_locator                 vertex_locator;

  triangle_count_visitor():
      vertex(),
      first(),
      second(),
      last_degree(0) { }
  
  triangle_count_visitor(vertex_locator v):
        vertex(v),
        first(),
        second(),
        last_degree(0) { }

  triangle_count_visitor(vertex_locator v, vertex_locator f, uint32_t deg):
        vertex(v), first(f),
        second(),
        last_degree(deg) { }

  triangle_count_visitor(vertex_locator v, vertex_locator f, vertex_locator s, uint32_t deg):
        vertex(v), first(f), second(s), last_degree(deg) { }

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    //return true;//(second == std::numeric_limits<vertex_descriptor_type>::max());
    uint32_t my_degree = std::get<0>(alg_data).degree(vertex); 
    if( last_degree < my_degree ) {
      return true;
    } else if( last_degree > my_degree) {
      return false;
    }
    return first < vertex;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {

    if(!first.is_valid()) {
      typedef typename Graph::edge_iterator eitr_type;
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        my_type new_visitor( neighbor, vertex, g.degree(vertex) );
        vis_queue->queue_visitor(new_visitor);
      }
    } else if(!second.is_valid()) {
      typedef typename Graph::edge_iterator eitr_type;
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        my_type new_visitor( neighbor, vertex, first, g.degree(vertex) );
        vis_queue->queue_visitor(new_visitor);
      }
    } else {
      //Make a binary search here!
      typedef typename Graph::edge_iterator eitr_type;
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        if(neighbor == second) {
          std::get<1>(alg_data)[vertex] = std::get<1>(alg_data)[vertex] + 1;
        }
      }
    }
    return true;
  }

  int get_state() const {
    if(!first.is_valid()) {
      return 0;
    }
    if(!second.is_valid()) {
      return 1;
    }
    return 2;
  }
  
  friend inline bool operator>(const triangle_count_visitor& v1, const triangle_count_visitor& v2) {
    return v1.get_state() < v2.get_state();
  }

  friend inline bool operator<(const triangle_count_visitor& v1, const triangle_count_visitor& v2) {
    return v1.get_state() > v2.get_state();
  }

  vertex_locator vertex, first, second;
  uint32_t last_degree;
};

template <typename TGraph>
uint64_t triangle_count(TGraph& g, typename TGraph::vertex_locator s) {
  typedef TGraph                                             graph_type;
  typedef typename mpi::triangle_count_visitor<TGraph>             visitor_type;
  typename graph_type::template vertex_data<uint64_t, std::allocator<uint64_t> >   tc_data(g);
  auto alg_data = std::forward_as_tuple(g, tc_data);
  auto vq = create_visitor_queue<visitor_type, triangle_priority_queue>(&g, alg_data);
  vq.init_visitor_traversal(s);
  return tc_data.global_accumulate();  
}

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
