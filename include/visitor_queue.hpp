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


#ifndef HAVOQGT_OMP_VISITOR_QUEUE_HPP_INCLUDED
#define HAVOQGT_OMP_VISITOR_QUEUE_HPP_INCLUDED

#include <mailbox.hpp>
#include <termination_detection.hpp>
#include <omp.hpp>
#include <deque>
#include <queue>

namespace havoqgt { namespace omp {

template <typename Graph, typename Visitor>
class visitor_queue 
{
  typedef typename Graph::vertex_descriptor vertex_descriptor;
private:
  class thread_data
  {
  public:
    std::priority_queue<Visitor, std::vector<Visitor>, std::greater<Visitor> > local_queue;
  private:
    char __pad[128-sizeof(local_queue)];
  };

public:
  visitor_queue(Graph* ptr_graph)
    : m_mailbox(this)
    , m_thread_data(num_threads())
    , m_ptr_graph(ptr_graph)
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      m_thread_data[thread_num()] = new thread_data();
    }
  }

  ~visitor_queue()
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      delete m_thread_data[thread_num()];
    }
  }

  void init_visitor_traversal(vertex_descriptor _source_v)
  {
    assert_sequential();
    #pragma omp parallel
    {
      assert_parallel();
      if(thread_num() == _source_v.thread_id) 
      { 
        queue_visitor( Visitor(_source_v) ); 
      }

      do {
        while(!m_thread_data[thread_num()]->local_queue.empty())
        {
          Visitor this_visitor = m_thread_data[thread_num()]->local_queue.top();
          m_thread_data[thread_num()]->local_queue.pop();
          this_visitor.visit(m_ptr_graph, this);
          m_td.inc_completed();
        }
        m_mailbox.idle_receive();
      } while(!m_mailbox.is_idle() || !m_td.test_for_termination());
    }
  }

  void queue_visitor(const Visitor& v) 
  {
    if(v.vertex.thread_id == thread_num()) {
      if(v.pre_visit()) {
        m_thread_data[thread_num()]->local_queue.push(v);
        m_td.inc_sent();
      }
    } else {
      m_mailbox.send(v.vertex.thread_id, v);
      m_td.inc_sent();
    }
  }

protected:
  template <typename T, typename U> friend class havoqgt::omp::mailbox;
  template <typename T>
  void mailbox_receive(const T& data)
  {
    // do something with data, previsit!
    if(data.pre_visit())
    {
      m_thread_data[thread_num()]->local_queue.push(data);
    } else
    {
      m_td.inc_completed();
    }
  }

private:
  termination_detection m_td;
  mailbox<Visitor, visitor_queue> m_mailbox;
  std::vector<thread_data*> m_thread_data;
  Graph* m_ptr_graph;
  
};

} } //end havoqgt::omp

#endif //HAVOQGT_OMP_VISITOR_QUEUE_HPP_INCLUDED
