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


#ifndef HAVOQGT_OMP_GENERATE_GRAPH500_RMAT_HPP_INCLUDED
#define HAVOQGT_OMP_GENERATE_GRAPH500_RMAT_HPP_INCLUDED

#include <mailbox.hpp>
#include <termination_detection.hpp>
#include <rmat_edge_generator.hpp>

namespace havoqgt { namespace omp {

template <typename AdjListGraph>
class generate_graph500_rmat {
  typedef AdjListGraph adj_list_graph_type;
  typedef typename adj_list_graph_type::edge_descriptor edge_descriptor;

public:
  generate_graph500_rmat(adj_list_graph_type* ptr_graph, uint64_t scale)
    : m_mailbox(this)
    , m_ptr_graph(ptr_graph)
    , m_scale(scale)
  {
    run();
  }
private:
  void run()
  {
    using namespace havoqgt;
    using namespace havoqgt::omp;
    assert_sequential();

    uint64_t num_vertices = uint64_t(1) << m_scale;
    uint64_t num_uedges = num_vertices * 16; 
    #pragma omp parallel
    {
      assert_parallel();
      rmat_edge_generator rmat(thread_num() * 13, m_scale, num_uedges,
                               0.57, 0.19, 0.19, 0.05, true, true);
      rmat_edge_generator::input_iterator_type itr = rmat.begin();

      #pragma omp for nowait
      for (uint64_t i = 0; i < num_uedges; ++i)
      {
        std::pair<uint64_t,uint64_t> raw_edge = *itr; ++itr;
        edge_descriptor edge;
        edge.first.thread_id  = raw_edge.first % num_threads();
        edge.first.local_id   = raw_edge.first / num_threads();
        edge.second.thread_id = raw_edge.second % num_threads();
        edge.second.local_id  = raw_edge.second / num_threads();

        m_mailbox.send(edge.first.thread_id, edge);
        m_td.inc_sent();

        // make edge undirected
        std::swap(edge.first, edge.second);
        m_mailbox.send(edge.first.thread_id, edge);
        m_td.inc_sent();
      }
      
      while(!m_mailbox.is_idle() || !m_td.test_for_termination()) {
        m_mailbox.idle_receive();
      }
    }

    assert(m_td.global_completed() == num_uedges*2);
  }

protected:
  template <typename T, typename U> friend class havoqgt::omp::mailbox;
  template <typename T>
  void mailbox_receive(const T& data)
  {
    assert(data.first.thread_id == thread_num());
    m_ptr_graph->add_edge(data);
    m_td.inc_completed(1);
  }
private:
  havoqgt::omp::termination_detection m_td;
  havoqgt::omp::mailbox<edge_descriptor, generate_graph500_rmat> m_mailbox;
  adj_list_graph_type* m_ptr_graph;
  uint64_t m_scale;
};

}} //end havoqgt::omp


#endif // HAVOQGT_OMP_GENERATE_GRAPH500_RMAT_HPP_INCLUDED
