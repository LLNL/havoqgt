// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED
#define HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED



#include <havoqgt/visitor_queue.hpp>
#include <queue>

namespace havoqgt {


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

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    bool do_visit = std::get<0>(alg_data)[vertex] > m_path;
    if(do_visit) {
      std::get<0>(alg_data)[vertex] = m_path; 
    }
    return do_visit;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {  
    if(m_path <= std::get<0>(alg_data)[vertex]) 
    {
      std::stringstream line;
      //line << "Visiting vertex: " << g.locator_to_label(vertex);
      typedef typename Graph::edge_iterator eitr_type;
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        path_type weight = std::get<1>(alg_data)[eitr];
        //line << ", Sending to neighbor: " << g.locator_to_label(neighbor) << " dist = " << weight+m_path;
        sssp_visitor new_visitor( neighbor, weight + m_path);
        vis_queue->queue_visitor(new_visitor);
      }
      //std::cout << line.str() << std::endl;
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

  vertex_locator   vertex;
  path_type        m_path;
};

 
template <typename TGraph, typename PathData, typename EdgeWeight>
void single_source_shortest_path(TGraph& g, 
                                 PathData& path_data, 
                                 EdgeWeight& edge_data,
                                 typename TGraph::vertex_locator s) {
  
  typedef  sssp_visitor<TGraph, PathData, EdgeWeight>    visitor_type;
  auto alg_data = std::forward_as_tuple(path_data, edge_data);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(&g, alg_data);
  vq.init_visitor_traversal(s);
}



} //end namespace havoqgt




#endif //HAVOQGT_MPI_SINGLE_SOURCE_SHORTEST_PATH_HPP_INCLUDED
