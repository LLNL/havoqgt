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
#include <boost/dynamic_bitset.hpp>


namespace havoqgt { namespace mpi {


enum visit_t { BAD, ADD, REVERSEADD, CHK, DEL };



template <typename Visitor>
class fifo_queue {
 protected:
  boost::container::deque<Visitor> m_data;
 public:
  fifo_queue() { }

  bool push(Visitor const & task) {
    m_data.push_back(task);
    return true;
  }

  void pop() {
    m_data.pop_front();
  }

  Visitor const & top() {
    return m_data.front();
  }

  size_t size() const {
    return m_data.size();;
  }

  bool empty() const {
    return m_data.empty();
  }

  void clear() {
    m_data.clear();
  }
};



template<typename Graph, typename prop_t>
class gc_dynamic {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor. Needs to be defined, but should not be used.
  gc_dynamic() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor. Needs to be defined, but should not be used.
  explicit gc_dynamic(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_colour(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming colour, and what visit type it is.
  gc_dynamic(vertex_locator _vertex, vertex_locator _caller, prop_t _colour, visit_t type) :
      vertex(_vertex), caller(_caller), caller_colour(_colour), vis_type(type) {  }

  // Creates an ADD visitor, based on the priority of the given vertices.
  static gc_dynamic<Graph, prop_t> create_visitor_add_type(vertex_locator vl_dst, vertex_locator vl_src) {
    // Start the add -> reverse add process with the higher priority vertex.
    if (vertex_locator::lesser_hash_priority(vl_dst, vl_src)) {
      // Make one src to dst.
      return gc_dynamic(vl_src, vl_dst, 0, ADD);
    } else {
      // Make one dst to src.
      return gc_dynamic(vl_dst, vl_src, 0, ADD);
    }
  }

  // Creates a DELETE visitor. TODO(Scott): Confirm it does not need priority.
  static gc_dynamic<Graph, prop_t> create_visitor_del_type(vertex_locator vl_dst, vertex_locator vl_src) {
    return gc_dynamic(vl_dst, vl_src, 0, DEL);
  }

  // Recolours our vertex based on knowledge of our neighbours/edges colours.
  prop_t recolour() const {
    boost::dynamic_bitset<> bitmap;
    bitmap.resize(graph_ref()->degree(vertex) + 1);  // 0 Colour offset.

    // Collect neighbour colours.
    for (auto nbr  = graph_ref()->adjacent_edge_begin(vertex);
              nbr != graph_ref()->adjacent_edge_end(vertex); nbr++) {
      auto nbr_col = nbr.property_data();
      if (nbr_col < bitmap.size() && nbr_col != 0) {
        bitmap.set(nbr_col - 1);  // 0 Colour offset.
      }
    }

    // Find first unused out of neighbour colours.
    prop_t colour = 1;
    if (bitmap.size() != 0) {  // TODO(Scott): double check this logic.
      bitmap.flip();
      colour = (bitmap.find_first() + 1);  // 0 Colour offset.
    }

    // Set colour.
    graph_ref()->vertex_property_data(vertex) = colour;
    return colour;
  }


  bool pre_visit() {
    // TODO(Scott): Consider self edges instead of just dropping them.
    if (vertex.id() == caller.id()) return false;

    // Perform an action depending on the visit type.
    switch (vis_type) {
      case ADD: {
        return true;

      } case REVERSEADD: {
        return true;

      } case CHK: {
        return true;

      } case DEL: {
        return true;

      } default: {
        std::cerr << "ERROR: Bad visit type (DEFAULT IN PRE)." << std::endl; exit(-1);
      }
    }

    return false;
  }


  // A visit will send the colour of the vertex to one or all of its neighbours.
  // Performs a specific action depending on the visit type.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    switch (vis_type) {
      case ADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == false);
        graph_ref()->insert_edge(vertex, caller, 0);
        prop_t* our_colour = &(graph_ref()->vertex_property_data(vertex));

        // If we are uncoloured (new), colour us (the first colour).
        if (*our_colour == 0) {
          *our_colour = 1;
        }

        // Note: visit only needs to be to new vertex, doesn't need to be all.
        gc_dynamic new_visitor(caller, vertex, *our_colour, REVERSEADD);
        vis_queue->queue_visitor(new_visitor);
        return false;

      } case REVERSEADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == true);
        graph_ref()->insert_edge(vertex, caller, 0);
        prop_t* our_colour = &(graph_ref()->vertex_property_data(vertex));
        graph_ref()->edge_property_data(vertex, caller) = caller_colour;

        // If we are uncoloured (new), colour us.
        if (*our_colour == 0) {
          if (caller_colour == 1) {
            *our_colour = 2;
          } else {
            *our_colour = 1;
          }

          // Note: visit only needs to be to new vertex, doesn't need to be all.
          gc_dynamic new_visitor(caller, vertex, *our_colour, CHK);
          vis_queue->queue_visitor(new_visitor);
          return false;
        // Not a new vertex.
        } else {
          // Check whether or not we conflict.
          if (*our_colour == caller_colour) {
            // Reverse add is always lower priority: we need to recolour.
            recolour();
            // Need to send our new colour to all neighbours.
            visitAllNbrs(graph, vis_queue);
            return false;
          }
          // No conflict, we can remain as-is.
        }
        return false;

      } case CHK: {
        // Set associated edge (vertex -> caller) with the caller's colour.
        assert(caller_colour != 0);
        graph_ref()->edge_property_data(vertex, caller) = caller_colour;

        // Check whether or not we conflict.
        if (graph_ref()->vertex_property_data(vertex) == caller_colour) {
          // Conflict: check if their are a higher priority than us.
          if (vertex_locator::lesser_hash_priority(vertex, caller)) {
            // Yes: we need to recolour.
            recolour();
            visitAllNbrs(graph, vis_queue);  // Need to send colour to all neighbours.
            return false;
          } else {
            // TODO(Scott): This should not be necessary!?!?!??!?! wat
            // gc_dynamic new_visitor(caller, vertex, graph_ref()->vertex_property_data(vertex), CHK);
            // vis_queue->queue_visitor(new_visitor);
            return false;
          }
        }
        // No conflict, we can remain as-is.
        return false;

      } case DEL: {
        // TODO(Scott): Might want to check that the edge actually exists, first.
        //              In case the graph is buggy.
        // Find what colour we will no longer be connected to, first.
        auto deleted_colour = graph_ref()->edge_property_data(vertex, caller);
        // Remove the edge as instructed.
        graph_ref()->erase_edge(vertex, caller);
        // Short circuit if we wouldn't want to take the colour anyway (it's a higher colour).
        auto our_colour = &(graph_ref()->vertex_property_data(vertex));
        if (deleted_colour > *our_colour) {  return false;  }

        for (auto nbr  = graph_ref()->adjacent_edge_begin(vertex);
                  nbr != graph_ref()->adjacent_edge_end(vertex); nbr++) {
          auto nbr_col = nbr.property_data();
          if (nbr_col == deleted_colour) {  return false;  } // Can't do anything: a nbr has the colour.
        }  // Check passed: no neighbours have the colour we want to take.
        *our_colour = deleted_colour;
        visitAllNbrs(graph, vis_queue);  // Need to send check colour to all neighbours.
        return false;

      } default: {
        std::cerr << "ERROR: Bad visit type (DEFAULT IN VISIT)." << std::endl; exit(-1);
      }
    }
    return false;
  }


  // This "Checks" all neighbours, with our colour.
  template<typename VisitorQueueHandle>
  inline void visitAllNbrs(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t mycolour = graph.vertex_property_data(vertex);
    // Send to all nbrs our current colour.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our colour.
      gc_dynamic new_visitor(vl_nbr, vertex, mycolour, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
  }


  friend inline bool operator > (const gc_dynamic& v1,
                                 const gc_dynamic& v2) {
    // Have to do adds first, for some reason that is completely unknown.
    //if (v2.vis_type == ADD || v2.vis_type == REVERSEADD) {
    //  return true;
    //}
    return false;
  }

//  /// --- Keita commented out --- ///
//  template<typename VisitorQueueHandle>
//  static void add_edge(havoqgt::parallel_edge_list_reader::edge_type edge,
//                       VisitorQueueHandle vis_queue) {
//    auto src = std::get<0>(edge);
//    auto dst = std::get<1>(edge);

//    vertex_locator vl_src(src);
//    vertex_locator vl_dst(dst);

//    // Start the add -> reverse add process with the higher priority vertex.
//    if (vertex_locator::lesser_hash_priority(vl_dst, vl_src)) {
//      // Make one src to dst.
//      gc_dynamic srcdst(vl_src, vl_dst, 0, ADD);
//      vis_queue->queue_visitor(srcdst);
//    } else {
//      // Make one dst to src.
//      gc_dynamic dstsrc(vl_dst, vl_src, 0, ADD);
//      vis_queue->queue_visitor(dstsrc);
//    }
//  }


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
  prop_t         caller_colour;
  visit_t        vis_type;
} __attribute__((packed));



template<typename Graph, typename prop_t>
class gc_dynamic_tester {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor.
  gc_dynamic_tester() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor.
  explicit gc_dynamic_tester(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_colour(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming colour, and what visit type it is.
  gc_dynamic_tester(vertex_locator _vertex, vertex_locator _caller, prop_t _colour, visit_t type) :
      vertex(_vertex), caller(_caller), caller_colour(_colour), vis_type(type) {  }


  bool pre_visit() {
    if (vis_type == CHK) {
      if (vertex.id() == caller.id()) return false;
      if (caller_colour == graph_ref()->vertex_property_data(vertex)) {
        std::cout << "Bad colouring! Failure from: \n";
        std::cout << vertex.id() << ":" << caller.id() << " with colour " << caller_colour << "\n";
        exit(0);
      }
      if (graph_ref()->vertex_property_data(vertex) == 0) {
        std::cout << "Bad colouring! Uncoloured vertex: \n";
        std::cout << vertex.id() << "\n";
        exit(0);
      }
      return false;
    } else {
      return true;
    }
  }


  // A visit will send the colour of the vertex to its neighbours.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t mycolour = graph.vertex_property_data(vertex);
    // Send to all nbrs our current colour.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our colour.
      gc_dynamic_tester new_visitor(vl_nbr, vertex, mycolour, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
    return false;
  }


  friend inline bool operator > (const gc_dynamic_tester& v1,
                                 const gc_dynamic_tester& v2) {
    return false;
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
  prop_t         caller_colour;
  visit_t        vis_type;
} __attribute__((packed));



// Launch point for graph colouring.
template <typename TGraph, typename edgelist_t, typename prop_t>
void graph_colour_dynamic(TGraph* graph, edgelist_t* edgelist, bool test_result) {
  double time_start = MPI_Wtime();

  {
    typedef gc_dynamic<TGraph, prop_t> colourer_t;
    colourer_t::set_graph_ref(graph);

    typedef visitor_queue<colourer_t, fifo_queue,
                          TGraph> colourer_queue_t;
    colourer_queue_t colourer(graph);

    colourer.init_dynamic_traversal(edgelist);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  if (mpi_rank == 0) {
    std::cout << "Colour time = " << time_end - time_start << std::endl;
  }

  // Finished colouring: test.
  if (test_result) {
    if (mpi_rank == 0) {  std::cout << "Verifying...\n";  }
    typedef gc_dynamic_tester<TGraph, prop_t> tester_t;
    tester_t::set_graph_ref(graph);

    typedef visitor_queue<tester_t, havoqgt::detail::visitor_priority_queue,
                          TGraph> tester_queue_t;
    tester_queue_t tester(graph);

    // Begin traversal.
    tester.init_dynamic_test_traversal();
    if (mpi_rank == 0) {  std::cout << "Valid colouring.\n";  }
  }
}



}  // namespace mpi
}  // namespace havoqgt



#endif  // HAVOQGT_MPI_GRAPH_COLOUR_DYN_HPP_INCLUDED
