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

#ifndef HAVOQGT_MPI_BFS_DYN_HPP_INCLUDED
#define HAVOQGT_MPI_BFS_DYN_HPP_INCLUDED


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

template <typename Visitor>
class bfs_queue {
public:
  typedef uint32_t level_number_type;
  typedef typename boost::container::deque<Visitor>::size_type size_type;

protected:
  std::vector<boost::container::deque<Visitor> > m_vec_bfs_level_stack;
  level_number_type m_cur_min_level;
  size_type m_size;
public:
  bfs_queue() : m_vec_bfs_level_stack(20), m_cur_min_level(std::numeric_limits<level_number_type>::max()), m_size(0) { }

  bool push(Visitor const & task) {
    while(task.level() >= m_vec_bfs_level_stack.size()) {
      m_vec_bfs_level_stack.push_back(boost::container::deque<Visitor>());
    }
    m_vec_bfs_level_stack[task.level()].push_back(task);
    ++m_size;
    m_cur_min_level = std::min(m_cur_min_level, (uint32_t)task.level());
    return true;
  }

  void pop() {
    m_vec_bfs_level_stack[m_cur_min_level].pop_back();
    --m_size;
    if(m_vec_bfs_level_stack[m_cur_min_level].empty()) {
      //if now empty, find next level;
      for(;m_cur_min_level < m_vec_bfs_level_stack.size(); ++m_cur_min_level) {
        if(!m_vec_bfs_level_stack[m_cur_min_level].empty()) break;
      }
    }
  }

  Visitor const & top() {
    return m_vec_bfs_level_stack[m_cur_min_level].back();
  }

  size_type size() const {
    return m_size;
  }

  bool empty() const {
    return (m_size == 0);
  }

  void clear() {
     for(typename std::vector<boost::container::deque<Visitor> >::iterator itr = m_vec_bfs_level_stack.begin();
       itr != m_vec_bfs_level_stack.end(); ++itr) {
       itr->clear();
     }
     m_size = 0;
     m_cur_min_level = std::numeric_limits<level_number_type>::max();
   }
};



template<typename Graph, typename prop_t>
class bfs_dynamic {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor. Needs to be defined, but should not be used.
  bfs_dynamic() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor. Needs to be defined, but should not be used.
  explicit bfs_dynamic(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_level(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming level, and what visit type it is.
  bfs_dynamic(vertex_locator _vertex, vertex_locator _caller, prop_t _level, visit_t type) :
      vertex(_vertex), caller(_caller), caller_level(_level), vis_type(type) {  }

  // Creates an ADD visitor, based on the priority of the given vertices.
  static bfs_dynamic<Graph, prop_t> create_visitor_add_type(vertex_locator vl_dst, vertex_locator vl_src) {
    // Start the add -> reverse add process with the higher priority vertex.
    if (vertex_locator::lesser_hash_priority(vl_dst, vl_src)) {
      // Make one src to dst.
      return bfs_dynamic(vl_src, vl_dst, 0, ADD);
    } else {
      // Make one dst to src.
      return bfs_dynamic(vl_dst, vl_src, 0, ADD);
    }
  }

  // Not implemented.
  static bfs_dynamic<Graph, prop_t> create_visitor_del_type(vertex_locator vl_dst, vertex_locator vl_src) {
    assert(false);
  }

  // For init
  static bfs_dynamic<Graph, prop_t> create_visitor_init_type(vertex_locator vl_src, vertex_locator vl_dst) {
    return bfs_dynamic(vl_src, vl_dst, 0, INIT);
  }

  bool pre_visit() {
    // Self edges
    if (vertex.id() == caller.id()) return false;

    // Perform an action depending on the visit type.
    switch (vis_type) {
      case INIT: {
        return true;

      } case ADD: {
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


  // A visit will send the level of the vertex to one or all of its neighbours.
  // Performs a specific action depending on the visit type.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    switch (vis_type) {
      case INIT: {
        prop_t* our_level = &(graph_ref()->vertex_property_data(vertex)); // segfault if vertex doesn't exist yet..
        *our_level = 1;
        // Need to send our new level to all neighbours.
        visitAllNbrs(graph, vis_queue);
        return false;

      } case ADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == false);
        graph_ref()->insert_edge(vertex, caller, 0);
        prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));
        // init level
        if (*our_level == 0) {  *our_level = std::numeric_limits<prop_t>::max();  }

        // Do nothing but send our label.

        // Note: visit only needs to be to new vertex, doesn't need to be all.
        bfs_dynamic new_visitor(caller, vertex, *our_level, REVERSEADD);
        vis_queue->queue_visitor(new_visitor);
        return false;

      } case REVERSEADD: {
        assert(vertex_locator::lesser_hash_priority(vertex, caller) == true);
        graph_ref()->insert_edge(vertex, caller, 0);
        prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));
        // init level
        if (*our_level == 0) {  *our_level = std::numeric_limits<prop_t>::max();  }

        if (*our_level == caller_level) {
          // do nothing

        // Check whether we have a lower level. (hop away offset)
        } else if (*our_level < caller_level - 1) {
          // Chk back the caller.
          bfs_dynamic new_visitor(caller, vertex, *our_level, CHK);
          vis_queue->queue_visitor(new_visitor);

        // Check if they have a lower level. (hop away offset)
        } else if (*our_level > caller_level + 1) {
          *our_level = caller_level + 1;
          // Need to send our new level to all neighbours.
          visitAllNbrs(graph, vis_queue);
        }
        return false;

      } case CHK: {
        prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));

        if (*our_level == caller_level) {
          // do nothing

        // Check whether we have a lower level. (hop away offset)
        } if (*our_level < caller_level - 1) {
          // Chk back the caller.
          bfs_dynamic new_visitor(caller, vertex, *our_level, CHK);
          vis_queue->queue_visitor(new_visitor);

        // Check if they have a lower level. (hop away offset)
        } else if (*our_level > caller_level + 1) {
          *our_level = caller_level + 1;
          // Need to send our new level to all neighbours.
          visitAllNbrs(graph, vis_queue);
        }
        return false;

      } default: {
        std::cerr << "ERROR: Bad visit type (DEFAULT IN VISIT)." << std::endl; exit(-1);
      }
    }
    return false;
  }


  // This "Checks" all neighbours, with our level.
  template<typename VisitorQueueHandle>
  inline void visitAllNbrs(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t my_level = graph.vertex_property_data(vertex);
    // Send to all nbrs.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our level.
      bfs_dynamic new_visitor(vl_nbr, vertex, my_level, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
  }


  friend inline bool operator > (const bfs_dynamic& v1,
                                 const bfs_dynamic& v2) {
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
  prop_t         caller_level;
  visit_t        vis_type;
} __attribute__((packed));



template<typename Graph, typename prop_t>
class bfs_dynamic_tester {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor.
  bfs_dynamic_tester() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor.
  explicit bfs_dynamic_tester(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_level(0), vis_type(BAD) {  }

  // Who I am, who notified me, the incoming level, and what visit type it is.
  bfs_dynamic_tester(vertex_locator _vertex, vertex_locator _caller, prop_t _level, visit_t type) :
      vertex(_vertex), caller(_caller), caller_level(_level), vis_type(type) {  }


  bool pre_visit() {
    if (vis_type == CHK) {
      if (vertex.id() == caller.id()) return false;  // self edge
      auto our_level = graph_ref()->vertex_property_data(vertex);
      if (caller_level == 0) {
        return false;
      }
      if (caller_level == our_level) {
        return false;
      }
      if (our_level != caller_level + 1 && our_level != caller_level - 1) {
        std::cout << "Bad result! Failure from: \n";
        std::cout << vertex.id() << ":" << caller.id() << " with caller level " << caller_level << ", our level " << our_level << "\n";
        exit(0);
      }
      return false;
    } else {
      return true;
    }
  }


  // A visit will send the level of the vertex to its neighbours. They should all be +- 1, ==, or 0.
  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    const prop_t mylevel = graph.vertex_property_data(vertex);
    // Send to all nbrs our level.
    for (auto nbr  = graph.adjacent_edge_begin(vertex);
              nbr != graph.adjacent_edge_end(vertex); nbr++) {
      auto edge = nbr.target_vertex();
      vertex_locator vl_nbr = vertex_locator(edge);

      // Send the neighbour a visitor with our level.
      bfs_dynamic_tester new_visitor(vl_nbr, vertex, mylevel, CHK);
      vis_queue->queue_visitor(new_visitor);
    }
    return false;
  }


  friend inline bool operator > (const bfs_dynamic_tester& v1,
                                 const bfs_dynamic_tester& v2) {
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
  prop_t         caller_level;
  visit_t        vis_type;
} __attribute__((packed));



// Static BFS in the dynamic world. To be run on a graph that is already constructed / construction is paused.
template<typename Graph, typename prop_t>
class bfs_static {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor.
  bfs_static() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor -- INIT?
  explicit bfs_static(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_level(0), vis_type(INIT) {  }

  // Who I am, who notified me, the incoming level, and what visit type it is.
  bfs_static(vertex_locator _vertex, vertex_locator _caller, prop_t _level, visit_t type) :
      vertex(_vertex), caller(_caller), caller_level(_level), vis_type(type) {  }
  
  // For init
  static bfs_static<Graph, prop_t> create_visitor_init_type(vertex_locator vl_src, vertex_locator vl_dst) {
    return bfs_static(vl_src, vl_dst, 0, INIT);
  }


  bool pre_visit() {
    if (vis_type == INIT) {
      prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));
      *our_level = 1;
      return true;
    } else if (vis_type == CHK) {
      if (vertex.id() == caller.id()) {  return false;  }
      prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));

      // Check if they have a lower level. (hop away offset)
      if (*our_level > caller_level + 1) {
        *our_level = caller_level + 1;
        return true;
      }
    }
    return false;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    prop_t* our_level = &(graph.vertex_property_data(vertex));
    if (*our_level == caller_level + 1) {
      
      // Send to all nbrs our level.
      for (auto nbr  = graph.adjacent_edge_begin(vertex);
                nbr != graph.adjacent_edge_end(vertex); nbr++) {
        auto edge = nbr.target_vertex();
        vertex_locator vl_nbr = vertex_locator(edge);
      
        // Send the neighbour a visitor with our level.
        bfs_static new_visitor(vl_nbr, vertex, *our_level, CHK);
        vis_queue->queue_visitor(new_visitor);
      }
      return false;
    }
  }


  uint64_t level() const {  return caller_level; }

  friend inline bool operator > (const bfs_static& v1,
                                 const bfs_static& v2) {
    if(v1.level() > v2.level()) {
      return true;
    } else if (v1.level() < v2.level()) {
      return false;
    }
    return !(v1.vertex < v2.vertex);
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
  prop_t         caller_level;
  visit_t        vis_type;
} __attribute__((packed));


// Resets vertex_property_data to MAX
template<typename Graph, typename prop_t>
class reset_vertex_prop {
 public:
  typedef typename Graph::vertex_locator vertex_locator;
  // Default constructor.
  reset_vertex_prop() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor
  explicit reset_vertex_prop(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), caller_level(0), vis_type(INIT) {  }

  // Who I am, who notified me, the incoming level, and what visit type it is.
  reset_vertex_prop(vertex_locator _vertex, vertex_locator _caller, prop_t _level, visit_t type) :
      vertex(_vertex), caller(_caller), caller_level(_level), vis_type(type) {  }


  bool pre_visit() {
    return true;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& graph, VisitorQueueHandle vis_queue) const {
    prop_t* our_level = &(graph_ref()->vertex_property_data(vertex));
    *our_level = std::numeric_limits<prop_t>::max();
    return false;
  }


  friend inline bool operator > (const reset_vertex_prop& v1,
                                 const reset_vertex_prop& v2) {
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
  prop_t         caller_level;
  visit_t        vis_type;
} __attribute__((packed));


// Launch point for dynamic bfs.
template <typename TGraph, typename edgelist_t, typename prop_t>
void breadth_first_search_dynamic(TGraph* graph, edgelist_t* edgelist, typename TGraph::vertex_locator src, bool test_result) {
  double time_start = MPI_Wtime();

  {
    typedef bfs_dynamic<TGraph, prop_t> bfs_t;
    bfs_t::set_graph_ref(graph);

    typedef visitor_queue<bfs_t, fifo_queue, TGraph> bfs_queue_t;
    bfs_queue_t bfs(graph);

    bfs.init_dynamic_traversal(edgelist, src);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_end = MPI_Wtime();
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  if (mpi_rank == 0) {
    std::cout << "Dynamic BFS time = " << time_end - time_start << std::endl;
  }

  {
  // Finished computing, now test.
  if (test_result) {
    if (mpi_rank == 0) {  std::cout << "Verifying...\n";  }
    typedef bfs_dynamic_tester<TGraph, prop_t> tester_t;
    tester_t::set_graph_ref(graph);

    typedef visitor_queue<tester_t, havoqgt::detail::visitor_priority_queue,
                          TGraph> tester_queue_t;
    tester_queue_t tester(graph);

    // Begin traversal.
    tester.init_dynamic_test_traversal();
    if (mpi_rank == 0) {  std::cout << "Valid.\n";  }
  }
  }

  // reset
  {
    if (mpi_rank == 0) {  std::cout << "Resetting vertex properties...\n";  }
    typedef reset_vertex_prop<TGraph, prop_t> reset_t;
    reset_t::set_graph_ref(graph);
    typedef visitor_queue<reset_t, fifo_queue, TGraph> reset_queue_t;
    reset_queue_t reset(graph);
    reset.init_dynamic_test_traversal();
    if (mpi_rank == 0) {  std::cout << "Reset.\n";  }
  }
  
  // static bfs
  double time_start_s = MPI_Wtime();
  {
    if (mpi_rank == 0) {  std::cout << "Init static...\n";  }
    typedef bfs_static<TGraph, prop_t> bfs_static_t;
    bfs_static_t::set_graph_ref(graph);
    typedef visitor_queue<bfs_static_t, detail::visitor_priority_queue, TGraph> bfs_static_queue_t;
    bfs_static_queue_t bfs_s(graph);
    bfs_s.init_dynamic_test_traversal_src(src);
    if (mpi_rank == 0) {  std::cout << "Done static...\n";  }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  double time_end_s = MPI_Wtime();
  int mpi_rank_s(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_s));
  if (mpi_rank_s == 0) {
    std::cout << "Static BFS time = " << time_end_s - time_start_s << std::endl;
  }
}



}  // namespace mpi
}  // namespace havoqgt



#endif  // HAVOQGT_MPI_BFS_DYN_HPP_INCLUDED
