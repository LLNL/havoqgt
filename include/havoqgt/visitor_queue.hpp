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

#ifndef HAVOQGT_MPI_VISITOR_QUEUE_HPP_INCLUDED
#define HAVOQGT_MPI_VISITOR_QUEUE_HPP_INCLUDED

#include <havoqgt/new_mailbox.hpp>
#include <havoqgt/mailbox.hpp>
#include <havoqgt/termination_detection.hpp>
#include <havoqgt/detail/reservable_priority_queue.hpp>
#include <vector>
#include <iterator>
#include <sched.h>

namespace havoqgt { namespace mpi {

class oned_blocked_partitioned_t { };
class el_partitioned_t { };

template <typename TVisitor, template<typename T> class Queue, typename TGraph>
class visitor_queue {
  typedef TVisitor              visitor_type;

  typedef termination_detection<uint64_t> termination_detection_type;
  typedef TGraph                graph_type;
  typedef typename TGraph::vertex_locator vertex_locator;
  //typedef typename havoqgt::detail::reservable_priority_queue<visitor_type,
  //    std::vector<visitor_type>, std::greater<visitor_type> > local_queue_type;
  typedef  Queue<visitor_type> local_queue_type;

#ifdef __bgp__
  typedef mailbox_bgp_torus<visitor_type> mailbox_type;
#else
  //typedef mailbox_routed<visitor_type> mailbox_type;
  typedef mailbox<visitor_type> mailbox_type;
#endif


public:
  visitor_queue(TGraph* _graph)
    : m_mailbox(/*MPI_COMM_WORLD,*/ 0)
    , m_termination_detection(MPI_COMM_WORLD, 2, 2, 3, 4)
    , m_ptr_graph(_graph) {
    //m_localqueue_owned.reserve(_graph->num_local_vertices());
    //m_localqueue_delegates.reserve(_graph->num_delegates() * 4);
  }

  ~visitor_queue() {
    // if(m_mailbox.comm_rank() == 0) {
    //   std::cout << "***************  Visitor Queue Statistics *****************" << std::endl;
    //   std::cout << "m_localqueue_owned.reserve  = " << m_ptr_graph->num_local_vertices() << std::endl;
    //   std::cout << "m_localqueue_owned.capacity = " << m_localqueue_owned.capacity() << std::endl;
    //   std::cout << "m_localqueue_delegates.reserve  = " << 0/*m_ptr_graph->num_delegates()*/ << std::endl;
    //   std::cout << "m_localqueue_delegates.capacity = " << m_localqueue_delegates.capacity() << std::endl;
    //   std::cout << "***********************************************************" << std::endl;
    // }
  }

  class visitor_queue_inserter
        : public std::iterator<std::output_iterator_tag, void, void, void, void> {
  public:
    visitor_queue_inserter(visitor_queue* _vq):m_vq(_vq) { }
    visitor_queue_inserter& operator=(const visitor_type& __value) {
      m_vq->handle_mailbox_receive(__value);
      return *this;
    }

    bool intercept(visitor_type& __value) {
      // assert(m_vq->m_ptr_graph->master(__value.m_visitor.vertex) != uint32_t(m_vq->m_mailbox.comm_rank()));
      bool ret = __value.pre_visit();
      if(!ret) {
        m_vq->m_termination_detection.inc_completed();
      }
      return ret;
    }

    /// Simply returns *this.
    visitor_queue_inserter&
    operator*()
    { return *this; }

    /// Simply returns *this.  (This %iterator does not "move".)
    visitor_queue_inserter&
    operator++()
    { return *this; }

    /// Simply returns *this.  (This %iterator does not "move".)
    visitor_queue_inserter
    operator++(int)
    { return *this; }

  private:
    visitor_queue* m_vq;
  };

  void init_visitor_traversal(vertex_locator _source_v) {
    if(0 /*_source_v.owner()*/ == m_mailbox.comm_rank()) {
      queue_visitor(visitor_type(_source_v));
    }
    do {
      do {
      process_pending_controllers();
      while(!empty()) {
        process_pending_controllers();
        visitor_type this_visitor = pop_top();
        vertex_locator v = this_visitor.vertex;
        bool ret = this_visitor.visit(*m_ptr_graph, this);
        if(ret && v.is_delegate() && m_ptr_graph->master(v) == m_mailbox.comm_rank()) {
          m_mailbox.bcast(this_visitor, visitor_queue_inserter(this));
          m_termination_detection.inc_queued(m_mailbox.comm_size());
        }
        m_termination_detection.inc_completed();
      }
      m_mailbox.flush_buffers_if_idle();
      } while(!m_local_controller_queue.empty() || !m_mailbox.is_idle() );
      sched_yield();
    } while(!m_termination_detection.test_for_termination());
  }


  void do_visit(visitor_type& this_visitor) {
    vertex_locator v = this_visitor.vertex;
    bool ret = this_visitor.visit(*m_ptr_graph, this);
    if(ret && v.is_delegate() && m_ptr_graph->master(v) == m_mailbox.comm_rank()) {
      visitor_type v = this_visitor;
      m_mailbox.bcast(v, visitor_queue_inserter(this));
      m_termination_detection.inc_queued(m_mailbox.comm_size());
    }
  }

  void do_init_visit(visitor_type& this_visitor) {
    vertex_locator v = this_visitor.vertex;
    bool ret = this_visitor.init_visit(*m_ptr_graph, this);
    if(ret && v.is_delegate() && m_ptr_graph->master(v) == m_mailbox.comm_rank()) {
      visitor_type v = this_visitor;
      m_mailbox.bcast(v, visitor_queue_inserter(this));
      m_termination_detection.inc_queued(m_mailbox.comm_size());
    }
  }

  // For dynamic graph traversal: this currently does one edge at a time per
  // process.
  void init_dynamic_traversal() {
    do {
      do {
        process_pending_controllers();
        while(!empty()) {
          process_pending_controllers();
          visitor_type this_visitor = pop_top();
          do_visit(this_visitor);
          m_termination_detection.inc_completed();
        }
        m_mailbox.flush_buffers_if_idle();
      } while(!m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
    // std::cout << havoqgt::havoqgt_env()->world_comm().rank() << "TERM ";
  }


  // Note: similar to below, but uses graphstore iterator and no delegates.
  void init_dynamic_test_traversal() {
    for(auto vitr = m_ptr_graph->vertices_begin(); vitr != m_ptr_graph->vertices_end(); vitr++) {
      vertex_locator vl(vitr.source_vertex());
      visitor_type v(vl);
      if(v.pre_visit()) {
        do_visit( v );
        check_mailbox();
      }
    }
    do {
      do {
      process_pending_controllers();
      while(!empty()) {
        process_pending_controllers();
        visitor_type this_visitor = pop_top();
        do_visit(this_visitor);
        m_termination_detection.inc_completed();
      }
      m_mailbox.flush_buffers_if_idle();
      } while(!m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }


  void init_visitor_traversal() {
    typename TGraph::controller_iterator citr = m_ptr_graph->controller_begin();
    for(; citr != m_ptr_graph->controller_end(); ++citr) {
      visitor_type v(*citr);
      if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
        do_visit( v );
        check_mailbox();
      }
    }
    typename TGraph::vertex_iterator vitr = m_ptr_graph->vertices_begin();
    for(; vitr != m_ptr_graph->vertices_end(); ++vitr) {
      visitor_type v(*vitr);
      if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
        do_visit( v );
        check_mailbox();
      }
    }
    do {
      do {
      process_pending_controllers();
      while(!empty()) {
        process_pending_controllers();
        visitor_type this_visitor = pop_top();
        do_visit(this_visitor);
        m_termination_detection.inc_completed();
      }
      m_mailbox.flush_buffers_if_idle();
      } while(!m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }

  void init_visitor_traversal_new() {
    auto citr = m_ptr_graph->controller_begin();
    auto vitr = m_ptr_graph->vertices_begin();

    do {
      do {
        if(citr != m_ptr_graph->controller_end()) {
          visitor_type v(*citr);
          do_init_visit(v);
          ++citr;
        }
        if(vitr != m_ptr_graph->vertices_end()) {
          visitor_type v(*vitr);
          do_init_visit(v);
          ++vitr;
        }
        process_pending_controllers();
        while(!empty()) {
          process_pending_controllers();
          visitor_type this_visitor = pop_top();
          do_visit(this_visitor);
          m_termination_detection.inc_completed();
        }
        m_mailbox.flush_buffers_if_idle();
      } while(citr != m_ptr_graph->controller_end() || vitr != m_ptr_graph->vertices_end()
              || !empty() || !m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }

  void init_visitor_traversal_new_alt() {
    auto citr = m_ptr_graph->controller_begin();
    auto vitr = m_ptr_graph->vertices_begin();

    do {
      do {
        do {
          if(citr != m_ptr_graph->controller_end()) {
            visitor_type v(*citr);
            if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
              do_visit( v );
              check_mailbox();
            }
            ++citr;
          }
          if(vitr != m_ptr_graph->vertices_end()) {
            visitor_type v(*vitr);
            if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
              do_visit( v );
              check_mailbox();
            }
            ++vitr;
          }
          process_pending_controllers();
          while(!empty()) {
            process_pending_controllers();
            visitor_type this_visitor = pop_top();
            do_visit(this_visitor);
            m_termination_detection.inc_completed();
          }
        } while( citr != m_ptr_graph->controller_end() ||  vitr != m_ptr_graph->vertices_end() );
        m_mailbox.flush_buffers_if_idle();
      } while( !empty() || !m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }

  void init_visitor_traversal_local_first() {
    typename TGraph::controller_iterator citr = m_ptr_graph->controller_begin();
    for(; citr != m_ptr_graph->controller_end(); ++citr) {
      visitor_type v(*citr);
      if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
        do_visit( v );
        // check_mailbox();
      }
    }
    typename TGraph::vertex_iterator vitr = m_ptr_graph->vertices_begin();
    for(; vitr != m_ptr_graph->vertices_end(); ++vitr) {
      visitor_type v(*vitr);
      if(v.pre_visit()) {  //RECENTLY ADDED 2013.10.10
        do_visit( v );
        // check_mailbox();
      }
    }
    do {
      do {
      process_pending_controllers();
      while(!empty()) {
        process_pending_controllers();
        visitor_type this_visitor = pop_top();
        do_visit(this_visitor);
        m_termination_detection.inc_completed();
      }
      m_mailbox.flush_buffers_if_idle();
      } while(!m_local_controller_queue.empty() || !m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }

  void queue_visitor(visitor_type& v) {
    if(v.vertex.is_delegate()) {
      local_delegate_visit(v);
    } else {
      if(v.vertex.owner() == uint32_t(m_mailbox.comm_rank())) {
        if(v.pre_visit()) {
          push(v);
          m_termination_detection.inc_queued();
        }
      } else {
        visitor_type vw = v;
        //vw.m_dest = v.vertex.owner();
        m_mailbox.send(v.vertex.owner(), vw, visitor_queue_inserter(this), false);
        m_termination_detection.inc_queued();
      }
    }
  }

private:
  // This occurs when the local process first encounters a delegate
  void local_delegate_visit(visitor_type& v) {
    if(v.pre_visit()) {
      if(m_ptr_graph->master(v.vertex) == uint32_t(m_mailbox.comm_rank())) {
        //delegate_bcast(v);
        push(v);
        m_termination_detection.inc_queued();
        /*  THIS was working, but trying to change how delegates are bcast from the master
        visitor_wrapper vw;
        vw.m_visitor = v;
        vw.set_bcast(true);
        m_mailbox.bcast(vw, visitor_queue_inserter(this));
        m_termination_detection.inc_queued(m_mailbox.comm_size());*/
      } else { //send interceptable to parent
        uint32_t master_rank = m_ptr_graph->master(v.vertex);
        m_mailbox.send(master_rank, v, visitor_queue_inserter(this), true);
        //delegate_parent(v);
        m_termination_detection.inc_queued();
      }
    }
  }

  /*oid delegate_parent(const visitor_type& v) {
    uint64_t parent = offset_tree_parent(m_mailbox.comm_size(), 2, m_ptr_graph->master(v.vertex),
                                m_mailbox.comm_rank());
    visitor_wrapper vw;
    vw.m_visitor = v;
    vw.m_visitor.vertex.set_dest(parent);
    m_mailbox.send(parent, vw, visitor_queue_inserter(this));
  }*/


  /*void delegate_bcast(const visitor_type& v) {
    uint32_t root = m_ptr_graph->master(v.vertex);
    uint32_t num_bcast_children = offset_tree_num_children(m_mailbox.comm_size(),
                                                           2, root,
                                                           m_mailbox.comm_rank());
    if(num_bcast_children > 0) {
      uint32_t first_bcast_child = offset_tree_first_child(m_mailbox.comm_size(),
                                                           2, root,
                                                           m_mailbox.comm_rank());
      for(uint32_t i=0; i<num_bcast_children; ++i) {
        uint32_t child = (first_bcast_child + i)%m_mailbox.comm_size();
        visitor_wrapper vw;
        vw.m_visitor = v;
        vw.m_visitor.vertex.set_dest(child);
        vw.m_visitor.vertex.set_bcast(true);
        m_mailbox.send(child, vw, visitor_queue_inserter(this));
      }
    }
  }*/

  void process_pending_controllers()
  {
    while(!m_local_controller_queue.empty()) {
      TVisitor v = m_local_controller_queue.front();
      m_local_controller_queue.pop();
      v.visit(*m_ptr_graph, this);
      m_termination_detection.inc_completed();
    }
  }

  void handle_mailbox_receive(visitor_type v) {
    if(v.vertex.is_delegate()) {
      if(v.vertex.get_bcast()) {
        //delegate_bcast(vw.m_visitor);
        if(m_ptr_graph->master(v.vertex) == uint32_t(m_mailbox.comm_rank())) {
          //This is because the mailbox bcast returns to self -- this should be fixed!
          m_termination_detection.inc_completed();
        } else {
          //vw.m_visitor.pre_visit();
          //push(vw.m_visitor);
          /*  2013.10.11 -- this causes too much recursion in mailbox, trying something new....
          vw.m_visitor.visit(*m_ptr_graph, this);
          m_termination_detection.inc_completed();
          */
          m_local_controller_queue.push(v);
        }
      } else {
        assert(m_ptr_graph->master(v.vertex) == uint32_t(m_mailbox.comm_rank()));
        if(v.pre_visit()) {
          //if(m_ptr_graph->master(vw.m_visitor.vertex) == m_mailbox.comm_rank()) {
            //delegate_bcast(vw.m_visitor);
            push(v);
            /* This was working, trying new way for master bcast
            vw.set_bcast(true);
            vw.set_intercept(false);
            m_mailbox.bcast(vw, visitor_queue_inserter(this));
            m_termination_detection.inc_queued(m_mailbox.comm_size());*/
          //} else {
          //  delegate_parent(vw.m_visitor);
          //}
        } else {
          m_termination_detection.inc_completed();
        }
      }
    } else {
      assert(v.vertex.owner() == uint32_t(m_mailbox.comm_rank()));
      //
      // Now handle owned vertices
      if(v.pre_visit()) {
        push(v);
      } else {
        m_termination_detection.inc_completed();
      }
    }
  }

  void push(const visitor_type& v) {
    /*if(v.vertex.is_delegate()) {
      m_localqueue_delegates.push(v);
    } else {
      m_localqueue_owned.push(v);
    }*/
    m_localqueue_owned.push(v);
  }

  visitor_type pop_top() {
    check_mailbox();
    visitor_type to_return;
    assert(!(m_localqueue_delegates.empty() && m_localqueue_owned.empty()));

    /*if(m_localqueue_delegates.empty()) {
      to_return = m_localqueue_owned.top();
      m_localqueue_owned.pop();
    } else {
      if(m_localqueue_owned.empty()) {
        to_return = m_localqueue_delegates.top();
        m_localqueue_delegates.pop();
      } else {
        if(m_localqueue_delegates.top() < m_localqueue_owned.top()) {
          to_return = m_localqueue_delegates.top();
          m_localqueue_delegates.pop();
        } else {
          to_return = m_localqueue_owned.top();
          m_localqueue_owned.pop();
        }
      }
    }*/
    to_return = m_localqueue_owned.top(); m_localqueue_owned.pop();
    return to_return;
  }


  void check_mailbox() {
    m_mailbox.receive(visitor_queue_inserter(this));
  }

  bool empty() {
    if(m_localqueue_owned.empty() && m_localqueue_delegates.empty()) {
      check_mailbox();
    }
    return m_localqueue_owned.empty() && m_localqueue_delegates.empty();
  }

/*  void init_visitor_traversal() {
    typedef typename graph_type::vertex_iterator vitr_type;
    std::pair< vitr_type, vitr_type > vitr = vertices(*m_ptr_graph);

    do {
      do {
        do {
          if(vitr.first != vitr.second) {
            uint64_t sampled_vertex = *(vitr.first);
            uint64_t my_rank = m_mailbox.comm_rank();
            if(m_ptr_graph->is_local(vertex(sampled_vertex, *m_ptr_graph))) {
              sampled_vertex |= (my_rank << owner_shift);
              queue_visitor( visitor_type( sampled_vertex ) );
            }
            ++vitr.first;
          }
          while(!empty()) {
            visitor_type this_visitor = top(); pop();
            vertex_locator v = this_visitor.vertex;
            this_visitor.visit((*m_ptr_graph)[v], *m_ptr_graph, this);
            m_termination_detection.inc_completed();
          }
       } while( vitr.first != vitr.second);
        m_mailbox.flush_buffers_if_idle();
      } while(!m_mailbox.is_idle() );
    } while(!m_termination_detection.test_for_termination());
  }*/





  mailbox_type           m_mailbox;
  termination_detection_type m_termination_detection;
  local_queue_type       m_localqueue_owned;
  local_queue_type       m_localqueue_delegates;
  TGraph*                m_ptr_graph;
  std::queue<TVisitor> m_local_controller_queue;
};


}} //namespace havoqgt::mpi

#endif //HAVOQGT_MPI_VISITOR_QUEUE_HPP_INCLUDED
