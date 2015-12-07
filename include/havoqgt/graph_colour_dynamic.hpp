/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce, Scott Sallinen <{rpearce, sallinen1}@llnl.gov>.
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

#ifndef HAVOQGT_MPI_GRAPH_COLOUR_DYN_HPP_INCLUDED
#define HAVOQGT_MPI_GRAPH_COLOUR_DYN_HPP_INCLUDED


#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>


namespace havoqgt { namespace mpi {


enum visit_t { BAD, INI, ADD, CHK, DEL };


template<typename Graph, typename prop_t>
class gc_dynamic {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor. Needs to be defined, but should not be used.
  gc_dynamic() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor. Needs to be defined, but should not be used.
  explicit gc_dynamic(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), colour(0), vis_type(INI) {  }

  // Who I am, who notified me, and what type of visit it is.
  gc_dynamic(vertex_locator _vertex, vertex_locator _caller, visit_t type) :
      vertex(_vertex), caller(_caller), colour(0), vis_type(type) {  }

  // Who I am, who notified me, the incoming colour, and what visit type it is.
  gc_dynamic(vertex_locator _vertex, vertex_locator _caller, prop_t _colour, visit_t type) :
      vertex(_vertex), caller(_caller), colour(_colour), vis_type(type) {  }


  // TODO
  bool pre_visit() const {
    // Perform an action dependant on the visit type.
    switch(vis_type) {
      case ADD:
        graph_ref()->insert_edge(vertex.id(), caller.id(), 0);

        if (vertex.id() <= 50) {
          std::cout << havoqgt::havoqgt_env()->world_comm().rank() << ":" << vertex.id() << "," << caller.id() << " ";
        }

        break;
      case CHK:
        break;
      default:
        std::cerr << "ERROR:  Bad visit type." << std::endl; exit(-1);
        break;
    }

    return false;
  }


  // A visit will send the colour of the vertex to its neighbours.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t mycolour = graph.vertex_property_data(vertex.id());
    // Send to all nbrs our current colour.
    for (auto nbr  = graph.adjacent_edge_begin(vertex.id());
              nbr != graph.adjacent_edge_end(vertex.id()); nbr++) {
      vertex_locator vl_nbr = vertex_locator(nbr.target_vertex());

      // Send the neighbour a visitor with our colour.
      gc_dynamic new_visitor(vl_nbr, vertex, mycolour, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
    //
    return false;
  }

  // TODO
  friend inline bool operator > (const gc_dynamic& v1,
                                 const gc_dynamic& v2) {
    return false;
  }


  template<typename VisitorQueueHandle>
  static void add_edge(havoqgt::parallel_edge_list_reader::edge_type edge,
                       VisitorQueueHandle vis_queue) {
    auto src = std::get<0>(edge);
    auto dst = std::get<1>(edge);

    vertex_locator vl_src(src);
    vertex_locator vl_dst(dst);

    // Make one src to dst, and another the opposite.
    gc_dynamic srcdst(vl_src, vl_dst, ADD);
    vis_queue->queue_visitor(srcdst);

    gc_dynamic dstsrc(vl_dst, vl_src, ADD);
    vis_queue->queue_visitor(dstsrc);
  }


  // Static graph reference.
  static void set_graph_ref(Graph* _graph_ref) {
    graph_ref() = _graph_ref;
  }
  static Graph*& graph_ref() {
    static Graph* _graph_ref;
    return _graph_ref;
  }


  // Instance variables.
  vertex_locator vertex;
  vertex_locator caller;
  prop_t         colour;
  visit_t        vis_type;
} __attribute__((packed));



// Launch point for graph colouring.
template <typename TGraph, typename prop_t>
void graph_colour_dynamic(TGraph* graph) {
  double time_start = MPI_Wtime();

  {
    typedef gc_dynamic<TGraph, prop_t> colourer_t;
    colourer_t::set_graph_ref(graph);

    typedef visitor_queue<colourer_t, havoqgt::detail::visitor_priority_queue,
                          TGraph> colourer_queue_t;
    colourer_queue_t colourer(graph);

    // Begin traversal.
    colourer.init_dynamic_traversal();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  if (mpi_rank == 0) {
    std::cout << "Colour time = " << time_end - time_start << std::endl;
  }
}



}  // namespace mpi
}  // namespace havoqgt



#endif  // HAVOQGT_MPI_GRAPH_COLOUR_DYN_HPP_INCLUDED
