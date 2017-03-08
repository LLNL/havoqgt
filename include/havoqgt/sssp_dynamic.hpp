/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Scott Sallinen, Roger Pearce <rpearce@llnl.gov, scottsallinen@gmail.com>.
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

#ifndef HAVOQGT_MPI_SSSP_DYN_HPP_INCLUDED
#define HAVOQGT_MPI_SSSP_DYN_HPP_INCLUDED


#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <boost/dynamic_bitset.hpp>


namespace havoqgt { namespace mpi {


enum visit_t { BAD, INIT, ADD, REVERSEADD, CHK };



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
class sssp_dynamic {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor. Needs to be defined, but should not be used.
  sssp_dynamic() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor. Needs to be defined, but should not be used.
  explicit sssp_dynamic(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_cost(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming weight, and what visit type it is.
  sssp_dynamic(vertex_locator _vertex, vertex_locator _caller, prop_t _cost, visit_t type) :
      vertex(_vertex), caller(_caller), caller_cost(_cost), vis_type(type) {  }

  // Creates an ADD visitor, based on the priority of the given vertices.
  static sssp_dynamic<Graph, prop_t> create_visitor_add_type(vertex_locator vl_dst, vertex_locator vl_src) {
    // Start the add -> reverse add process with the higher priority vertex.
    if (vertex_locator::lesser_hash_priority(vl_dst, vl_src)) {
      // Make one src to dst.
      return sssp_dynamic(vl_src, vl_dst, 0, ADD);
    } else {
      // Make one dst to src.
      return sssp_dynamic(vl_dst, vl_src, 0, ADD);
    }
  }

  // Not implemented.
  static sssp_dynamic<Graph, prop_t> create_visitor_del_type(vertex_locator vl_dst, vertex_locator vl_src) {
    assert(false);
  }

  // For init
  static sssp_dynamic<Graph, prop_t> create_visitor_init_type(vertex_locator vl_src, vertex_locator vl_dst) {
    return sssp_dynamic(vl_src, vl_dst, 0, INIT);
  }

  bool pre_visit() {
    // Self edges
    if (vertex.id() == caller.id()) return false;

    // Perform an action depending on the visit type.
    switch (vis_type) {
      case INIT: {
        return true;

      }case ADD: {
        return true;

      } case REVERSEADD: {
        return true;

      } case CHK: {
        return true;

      } default: {
        std::cerr << "ERROR: Bad visit type (DEFAULT IN PRE)." << std::endl; exit(-1);
      }
    }

    return false;
  }


  // A visit will send the cost of the vertex to one or all of its neighbours.
  // Performs a specific action depending on the visit type.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    switch (vis_type) {
      case INIT: {
        prop_t* our_cost = &(graph_ref()->vertex_property_data(vertex)); // segfault if vertex doesn't exist yet..
        prop_t edge_weight = graph_ref()->edge_property_data(vertex, caller);
        *our_cost = edge_weight;
        // Need to send our new cost to all neighbours.
        visitAllNbrs(graph, vis_queue);
        return false;

      } case ADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == false);
        prop_t edge_weight = 2; // todo : add edge weight
        graph_ref()->insert_edge(vertex, caller, edge_weight);
        prop_t* our_cost = &(graph_ref()->vertex_property_data(vertex));
        // init cost
        if (*our_cost == 0) {  *our_cost = std::numeric_limits<prop_t>::max();  }

        // Do nothing but send our cost.

        // Note: visit only needs to be to new vertex, doesn't need to be all.
        sssp_dynamic new_visitor(caller, vertex, *our_cost, REVERSEADD);
        vis_queue->queue_visitor(new_visitor);
        return false;

      } case REVERSEADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == true);
        prop_t edge_weight = 2; // todo : add edge weight
        graph_ref()->insert_edge(vertex, caller, edge_weight);
        prop_t* our_cost = &(graph_ref()->vertex_property_data(vertex));
        // init cost
        if (*our_cost == 0) {  *our_cost = std::numeric_limits<prop_t>::max();  }

        if (*our_cost == caller_cost) {
          // do nothing

        // Check whether we have a lower cost. (hop away offset)
        } else if (*our_cost < caller_cost - edge_weight) {
          // Chk back the caller.
          sssp_dynamic new_visitor(caller, vertex, *our_cost, CHK);
          vis_queue->queue_visitor(new_visitor);

        // Check if they have a lower cost. (hop away offset)
        } else if (*our_cost > caller_cost + edge_weight) {
          *our_cost = caller_cost + edge_weight;
          // Need to send our new cost to all neighbours.
          visitAllNbrs(graph, vis_queue);
        }
        return false;

      } case CHK: {
        prop_t* our_cost = &(graph_ref()->vertex_property_data(vertex));
        prop_t edge_weight = graph_ref()->edge_property_data(vertex, caller);

        if (*our_cost == caller_cost) {
          // do nothing

        // Check whether we have a lower cost. (hop away offset)
        } if (*our_cost < caller_cost - edge_weight) {
          // Chk back the caller.
          sssp_dynamic new_visitor(caller, vertex, *our_cost, CHK);
          vis_queue->queue_visitor(new_visitor);

        // Check if they have a lower cost. (hop away offset)
        } else if (*our_cost > caller_cost + edge_weight) {
          *our_cost = caller_cost + edge_weight;
          // Need to send our new cost to all neighbours.
          visitAllNbrs(graph, vis_queue);
        }
        return false;

      } default: {
        std::cerr << "ERROR: Bad visit type (DEFAULT IN VISIT)." << std::endl; exit(-1);
      }
    }
    return false;
  }


  // This "Checks" all neighbours, with our cost.
  template<typename VisitorQueueHandle>
  inline void visitAllNbrs(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t my_cost = graph.vertex_property_data(vertex);
    // Send to all nbrs.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our cost.
      sssp_dynamic new_visitor(vl_nbr, vertex, my_cost, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
  }


  friend inline bool operator > (const sssp_dynamic& v1,
                                 const sssp_dynamic& v2) {
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
  prop_t         caller_cost;
  visit_t        vis_type;
} __attribute__((packed));



template<typename Graph, typename prop_t>
class sssp_dynamic_tester {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor.
  sssp_dynamic_tester() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor.
  explicit sssp_dynamic_tester(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_cost(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming cost, and what visit type it is.
  sssp_dynamic_tester(vertex_locator _vertex, vertex_locator _caller, prop_t _cost, visit_t type) :
      vertex(_vertex), caller(_caller), caller_cost(_cost), vis_type(type) {  }

  // TODO(Scott) : figure out a good verify test
  bool pre_visit() {
    return false;
  }


  // TODO(Scott) : figure out a good verify test
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t mycost = graph.vertex_property_data(vertex);
    // Send to all nbrs our cost.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our cost.
      sssp_dynamic_tester new_visitor(vl_nbr, vertex, mycost, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
    return false;
  }


  friend inline bool operator > (const sssp_dynamic_tester& v1,
                                 const sssp_dynamic_tester& v2) {
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
  prop_t         caller_cost;
  visit_t        vis_type;
} __attribute__((packed));



// Launch point for dynamic sssp.
template <typename TGraph, typename edgelist_t, typename prop_t>
void single_source_shortest_path_dynamic(TGraph* graph, edgelist_t* edgelist, typename TGraph::vertex_locator src, bool test_result) {
  double time_start = MPI_Wtime();

  {
    typedef sssp_dynamic<TGraph, prop_t> sssp_t;
    sssp_t::set_graph_ref(graph);

    typedef visitor_queue<sssp_t, fifo_queue, TGraph> sssp_queue_t;
    sssp_queue_t sssp(graph);

    sssp.init_dynamic_traversal(edgelist, src);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  if (mpi_rank == 0) {
    std::cout << "sssp time = " << time_end - time_start << std::endl;
  }

  // Finished computing, now test.
  if (test_result) {
    if (mpi_rank == 0) {  std::cout << "Verifying...\n";  }
    typedef sssp_dynamic_tester<TGraph, prop_t> tester_t;
    tester_t::set_graph_ref(graph);

    typedef visitor_queue<tester_t, havoqgt::detail::visitor_priority_queue,
                          TGraph> tester_queue_t;
    tester_queue_t tester(graph);

    // Begin traversal.
    tester.init_dynamic_test_traversal();
    if (mpi_rank == 0) {  std::cout << "Valid.\n";  }
  }
}



}  // namespace mpi
}  // namespace havoqgt



#endif  // HAVOQGT_MPI_SSSP_DYN_HPP_INCLUDED
