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

#ifndef HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED
#define HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED


#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Graph>
class cc_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  cc_visitor() { }
  cc_visitor(vertex_locator _vertex, vertex_locator _cc)
    : vertex(_vertex)
    , m_cc(_cc) { }

  cc_visitor(vertex_locator _vertex)
    : vertex(_vertex)
    , m_cc(_vertex) { }

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(m_cc < (*alg_data)[vertex]) { 
      (*alg_data)[vertex] = m_cc; 
      return true;
    }
    return false;
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if((*alg_data)[vertex] == m_cc) {
      for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        auto neighbor = eitr.target();
        if(m_cc < neighbor) {
          cc_visitor new_visitor(neighbor, m_cc);
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return true;
    }
    return false;
  }


  friend inline bool operator>(const cc_visitor& v1, const cc_visitor& v2) {
    if(v2.m_cc < v1.m_cc)
    {
      return true;
    } else if(v1.m_cc < v2.m_cc)
    {
      return false;
    }
    if(v1.vertex == v2.vertex) return false;
    return !(v1.vertex < v2.vertex);
  }

  vertex_locator   vertex;
  vertex_locator  m_cc;
};

template <typename TGraph, typename CCData>
void connected_components(TGraph* g, CCData& cc_data) {

  typedef  cc_visitor<TGraph>    visitor_type;
 
  for(auto vitr = g->vertices_begin(); vitr != g->vertices_end(); ++vitr) {
    cc_data[*vitr] = *vitr;
  }
  for(auto citr = g->controller_begin(); citr != g->controller_end(); ++citr) {
    cc_data[*citr] = *citr;
  } 
  
  auto alg_data = &cc_data;
  auto vq = create_visitor_queue<visitor_type, detail::visitor_priority_queue>(g, alg_data); 
  vq.init_visitor_traversal_new();
}



}} //end namespace havoqgt::mpi




#endif //HAVOQGT_CONNECTED_COMPONENTS_HPP_INCLUDED
