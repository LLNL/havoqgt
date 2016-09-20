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
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {


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
template<typename Graph>
class kth_core_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  typedef kth_core_visitor<Graph>    my_type;

  kth_core_visitor() {}
  kth_core_visitor(vertex_locator _vertex):vertex(_vertex) { }

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(std::get<0>(alg_data)[vertex].get_alive()) {
      std::get<0>(alg_data)[vertex].set_core_bound( std::get<0>(alg_data)[vertex].get_core_bound() - 1);
      if(std::get<0>(alg_data)[vertex].get_core_bound() < std::get<1>(alg_data)) {
        std::get<0>(alg_data)[vertex].set_alive(false);
        return true;
      }
    }
    return false;
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(std::get<0>(alg_data)[vertex].get_alive() && 
       std::get<0>(alg_data)[vertex].get_core_bound() < std::get<1>(alg_data)) {
      std::get<0>(alg_data)[vertex].set_alive(false);
      return visit(g, vis_queue, alg_data);
    }
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //std::get<0>(alg_data)[vertex].set_alive(false);
    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vis_queue->queue_visitor( my_type(eitr.target()));
    }
    return true;
  }
  
  friend inline bool operator>(const kth_core_visitor& v1, const kth_core_visitor& v2) {
    if(v1.vertex == v2.vertex) return false;
    return !(v1.vertex < v2.vertex);
  }

  vertex_locator vertex;
};

template <typename TGraph, typename KCoreData>
void kth_core(TGraph& graph, KCoreData& k_core_data) {
  
  uint64_t to_return(0);
  typedef kth_core_visitor<TGraph>           visitor_type;
  int32_t kth_core_count = 0;
  
  auto alg_data = std::forward_as_tuple(k_core_data, kth_core_count); 

  //reset graph data
  for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    k_core_data[*vitr].set_alive(true);
    k_core_data[*vitr].set_core_bound(graph.degree(*vitr));
  }
  for(auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    k_core_data[*citr].set_alive(true);
    k_core_data[*citr].set_core_bound(graph.degree(*citr));
  } 

  auto vq = create_visitor_queue<visitor_type, detail::visitor_priority_queue>(&graph, alg_data);
  uint64_t count_alive = 0;
  do {
    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
      vq.init_visitor_traversal_new();
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();
    uint64_t local_alive(0);
    for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
      if(k_core_data[*vitr].get_alive()) ++local_alive;
    }
    for(auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
      if(k_core_data[*citr].get_alive()) ++local_alive;
    }
    count_alive = mpi_all_reduce(local_alive,std::plus<uint64_t>(), MPI_COMM_WORLD);
    if(havoqgt_env()->world_comm().rank() == 0) {
      std::cout << "Core " << std::get<1>(alg_data) << ", size = " << count_alive << ", time = " << time_end-time_start << std::endl;
    }
  } while(count_alive && ++std::get<1>(alg_data) < 50);
  
  //
  // for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
  //   if(k_core_data[*vitr].get_alive()) ++to_return;
  // }
  // for(auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
  //   if(k_core_data[*citr].get_alive()) ++to_return;
  // }
  // return mpi_all_reduce(to_return,std::plus<uint64_t>(), MPI_COMM_WORLD);
}

} } // end havoqgt::mpi


#endif //HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED
