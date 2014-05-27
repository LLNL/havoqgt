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

#ifndef HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED
#define HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED



#include <havoqgt/visitor_queue.hpp>
#include <queue>

namespace havoqgt { namespace mpi {


template <typename Visitor>
class sssp_queue
{

protected:
  std::priority_queue< Visitor, std::vector<Visitor>, 
                               std::greater<Visitor> > m_data;
public:
  sssp_queue() { }

  bool push(Visitor const & task)
  {
    m_data.push(task);
    return true;
  }

  void pop()
  {
    m_data.pop();
  }

  Visitor const & top() //const
  {
    return m_data.top();
  }

  size_t size() const
  {
    return m_data.size();;
  }

  bool empty() const
  {
    return m_data.empty();
  }

  void clear()
  {
    m_data.clear();
  }
};


template<typename Graph, typename PathData, typename EdgeWeight>
class sssp_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  typedef typename PathData::value_type                  path_type;
  sssp_visitor(): m_path(std::numeric_limits<path_type>::max())  { }
  sssp_visitor(vertex_locator _vertex, path_type _level):
              vertex(_vertex), m_path(_level) { }

  sssp_visitor(vertex_locator _vertex) :
              vertex(_vertex), m_path(0) { }            

  
  bool pre_visit() const {
    bool do_visit  = (*path_data())[vertex] > m_path;
    if(do_visit) {
      (*path_data())[vertex] = m_path;
    }
    return do_visit;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    if(m_path == (*path_data())[vertex]) 
    {
      typedef typename Graph::edge_iterator eitr_type;
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        path_type weight = (*edge_data())[eitr];
        //std::cout << "Visiting neighbor: " << g.locator_to_label(neighbor) << std::endl;
        sssp_visitor new_visitor( neighbor, weight + m_path);
        vis_queue->queue_visitor(new_visitor);
      }
      return true;
    }
    return false;
  }

  friend inline bool operator>(const sssp_visitor& v1, const sssp_visitor& v2) {
    return v1.m_path > v2.m_path;
  }

  friend inline bool operator<(const sssp_visitor& v1, const sssp_visitor& v2) {
    return v1.m_path < v2.m_path;
  }

  static void set_path_data(PathData* _data) { path_data() = _data; }

  static PathData*& path_data() {
    static PathData* data;
    return data;
  }

  static void set_edge_weight(EdgeWeight* _data) { edge_data() = _data; }

  static EdgeWeight*& edge_data() {
    static EdgeWeight* data;
    return data;
  }

  vertex_locator   vertex;
  path_type        m_path;
};

 
template <typename TGraph, typename PathData, typename EdgeWeight>
void single_source_shortest_path(TGraph& g, 
                                 PathData& path_data, 
                                 EdgeWeight& edge_data,
                                 typename TGraph::vertex_locator s) {
  
  typedef  sssp_visitor<TGraph, PathData, EdgeWeight>    visitor_type;
  visitor_type::set_path_data(&path_data);
  visitor_type::set_edge_weight(&edge_data);
  typedef visitor_queue< visitor_type, sssp_queue, TGraph >    visitor_queue_type;
   
  visitor_queue_type vq(&g);
  vq.init_visitor_traversal(s);
}



}} //end namespace havoqgt::mpi




#endif //HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED
