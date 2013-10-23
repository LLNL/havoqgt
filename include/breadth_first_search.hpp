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


#ifndef HAVOQGT_OMP_BREADTH_FIRST_SEARCH_HPP_INCLUDED
#define HAVOQGT_OMP_BREADTH_FIRST_SEARCH_HPP_INCLUDED

#include <mailbox.hpp>
#include <termination_detection.hpp>
#include <visitor_queue.hpp>

namespace havoqgt { namespace omp {

template <typename Graph>
class bfs_visitor 
{
public:
  typedef typename Graph::vertex_descriptor    vertex_descriptor;
  typedef typename Graph::template vertex_data<uint8_t> level_data;

  bfs_visitor(vertex_descriptor vert)
    : vertex (vert)
    , m_level(0) { }

  bfs_visitor(vertex_descriptor vert, uint8_t level)
    : vertex (vert)
    , m_level (level) { }

  bool pre_visit() const 
  {
    bool do_visit = m_level < (*sp_level_data)[vertex];
    if(do_visit)
    {
      (*sp_level_data)[vertex] = m_level;
    }
    return do_visit;
  }

  template <typename VisitorQueue>
  void visit(Graph* ptr_graph, VisitorQueue* ptr_vq)
  {
    if(m_level == (*sp_level_data)[vertex])
    {
      for(typename Graph::edge_iterator itr = ptr_graph->edges_begin(vertex);
          itr != ptr_graph->edges_end(vertex); ++itr)
      {
        vertex_descriptor neighbor = *itr;
        bfs_visitor new_visitor(neighbor, m_level + 1);
        ptr_vq->queue_visitor(new_visitor);
      }
    }
  }

  friend inline bool operator >(const bfs_visitor& v1, const bfs_visitor& v2)
  {
    return v1.m_level > v2.m_level;
  }

  static void set_level_data(level_data* ld)
  {
    sp_level_data = ld;
  }
public:
  vertex_descriptor vertex;
  uint64_t          m_level : 8;
  static level_data* sp_level_data;
};

template <typename Graph>
typename bfs_visitor<Graph>::level_data* bfs_visitor<Graph>::sp_level_data;


template <typename Graph>
void breadth_first_search(Graph& graph, typename Graph::vertex_descriptor source,
                          typename Graph::template vertex_data<uint8_t>& level_data)
{
    bfs_visitor<Graph>::set_level_data(&level_data);
    visitor_queue<Graph, bfs_visitor<Graph> > vq(&graph);
    vq.init_visitor_traversal(source);
}


}}

#endif //HAVOQGT_OMP_BREADTH_FIRST_SEARCH_HPP_INCLUDED