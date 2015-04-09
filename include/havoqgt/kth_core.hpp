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


#ifndef HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED
#define HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <vector>

namespace havoqgt { namespace mpi {


template <typename Visitor>
class kcore_queue
{

protected:
  std::vector<Visitor> m_data;
public:
  kcore_queue() { }

  bool push(Visitor const & task)
  {
    m_data.push_back(task);
    return true;
  }

  void pop()
  {
    m_data.pop_back();
  }

  Visitor const & top() //const
  {
    return m_data.back();
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


/// This is the kth-core data stored per vertex in the Graph
class kth_core_data {
public:  
  kth_core_data() : core_bound(std::numeric_limits<uint32_t>::max()),
               alive(false) { }
  
  uint32_t get_core_bound() const { return core_bound; }
  void set_core_bound(uint32_t bound) { core_bound = bound; }

  bool get_alive() const { return alive; }
  void set_alive( bool b) { alive = b; }

private:
  uint32_t core_bound;
  bool alive;
};

/// This is the kth-core visitor
template<typename Graph, typename KCoreData>
class kth_core_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  typedef kth_core_visitor<Graph, KCoreData>    my_type;

  kth_core_visitor() {}
  kth_core_visitor(vertex_locator _vertex):vertex(_vertex) { }

  bool pre_visit() const {
    if((*k_core_data)[vertex].get_alive()) {
      (*k_core_data)[vertex].set_core_bound( (*k_core_data)[vertex].get_core_bound() - 1);
      if((*k_core_data)[vertex].get_core_bound() < kth_core) {
        (*k_core_data)[vertex].set_alive(false);
        return true;
      }
    }
    return false;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    //(*k_core_data)[vertex].set_alive(false);
    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vis_queue->queue_visitor( my_type(eitr.target()));
    }
    return true;
  }

  vertex_locator vertex;
  static uint32_t kth_core;
  static KCoreData* k_core_data;
};

template<typename Graph, typename KCoreData>
uint32_t kth_core_visitor<Graph, KCoreData>::kth_core;

template<typename Graph, typename KCoreData>
KCoreData* kth_core_visitor<Graph,KCoreData>::k_core_data;


template <typename TGraph, typename KCoreData>
uint64_t kth_core(TGraph& graph, KCoreData& k_core_data, uint32_t kthcore) {
  
  uint64_t to_return(0);
  typedef kth_core_visitor<TGraph, KCoreData>           visitor_type;
  visitor_type::kth_core = kthcore;
  visitor_type::k_core_data = &k_core_data;

  //reset graph data
  for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    k_core_data[*vitr].set_alive(true);
    k_core_data[*vitr].set_core_bound(graph.degree(*vitr) + 1);
  }
  for(auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    k_core_data[*citr].set_alive(true);
    k_core_data[*citr].set_core_bound(graph.degree(*citr) + 1);
  } 

  typedef visitor_queue< visitor_type, kcore_queue, TGraph >    visitor_queue_type;

  visitor_queue_type vq(&graph);
  vq.init_visitor_traversal();
  
  for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    if(k_core_data[*vitr].get_alive()) ++to_return;
  }
  for(auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    if(k_core_data[*citr].get_alive()) ++to_return;
  }
  return mpi_all_reduce(to_return,std::plus<uint64_t>(), MPI_COMM_WORLD);
}

} } // end havoqgt::mpi


#endif //HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED