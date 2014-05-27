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
 

#ifndef HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED
#define HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED



#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>
#include <vector>

namespace havoqgt { namespace mpi {

template <typename Visitor>
class pr_queue
{

protected:
  std::vector< Visitor > m_data;
public:
  pr_queue() { }

  bool push(Visitor const & task)
  {
    m_data.push_back(task);
    return true;
  }

  void pop()
  {
    m_data.pop_back();
  }

  Visitor const & top() //const
  {
    return m_data.back();
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



template<typename Graph, typename PRData>
class pr_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  pr_visitor(): rank(0)  { }

  pr_visitor(vertex_locator _vertex, double _rank)
    : vertex(_vertex)
    , rank(_rank) { }

  pr_visitor(vertex_locator _vertex)
    : vertex(_vertex)
    , rank(0) { }      

  
  bool pre_visit() const {
    rank_data()[vertex] += rank;
    return false;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    double old_rank = 1;
    uint64_t degree = g.degree(vertex);
    double send_rank = 1;//old_rank / double(degree);


    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      pr_visitor new_visitor( neighbor, send_rank);
      vis_queue->queue_visitor(new_visitor);
    }
    return true;

  }


  friend inline bool operator>(const pr_visitor& v1, const pr_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const pr_visitor& v1, const pr_visitor& v2) {
    return false;
  }

  static PRData& rank_data(PRData* _data = NULL) {
    static PRData* data;
    if(_data) data = _data;
    return *data;
  }


  vertex_locator   vertex;
  double           rank;
};

 
template <typename TGraph, typename PRData>
void page_rank(TGraph& g, PRData& pr_data) {
  typedef  pr_visitor<TGraph, PRData>    visitor_type;
  visitor_type::rank_data(&pr_data);
  typedef visitor_queue< visitor_type, pr_queue, TGraph >    visitor_queue_type;
   
  visitor_queue_type vq(&g);
  vq.init_visitor_traversal();
  pr_data.all_reduce();
}



}} //end namespace havoqgt::mpi




#endif //HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED
