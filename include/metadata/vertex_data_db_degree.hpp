// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { 

template<typename Visitor>
class vertex_data_degree_queue {

public:

  vertex_data_degree_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    data.pop_back();
  }
 
  Visitor const& top() {
    return data.back();
  } 
  
  size_t size() const {
    return data.size();;
  }

  bool empty() const {
    return data.empty();
  }

  void clear() {
    data.clear();
  }

protected:

  std::vector<Visitor> data;

};

template<typename Graph, typename VertexData>
class vertex_data_degree_visitor {

public:

  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  vertex_data_degree_visitor() : 
    vertex_data(0),
    do_update_vertex_data(false) {}

  vertex_data_degree_visitor(vertex_locator _vertex) : 
    vertex(_vertex), 
    vertex_data(0), 
    do_update_vertex_data(false) {}

  vertex_data_degree_visitor(vertex_locator _vertex, VertexData _vertex_data) :
    vertex(_vertex), 
    vertex_data(_vertex_data), 
    do_update_vertex_data(true) {}
 
  ~vertex_data_degree_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    return true;
  } 

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //if (!do_update_vertex_data) {
      //std::cout << "False" << std::endl; // Test 
    //  return false;
    //} 
  
    // log2  
    std::get<0>(alg_data)[vertex] = static_cast<VertexData>(ceil(log2(g.degree(vertex) + 1)));
    
    //log6  
    //auto base = 6;
    //auto log_base_of_x = log10(g.degree(vertex)) / log10(base); 
    //std::get<0>(alg_data)[vertex] = static_cast<VertexData>(ceil(log_base_of_x + 1));         

    //std::cout << "Visiting " << g.locator_to_label(vertex) << " " << std::get<0>(alg_data)[vertex] << std::endl; // Test
    return true;
  } 

  friend inline bool operator>(const vertex_data_degree_visitor& v1, const vertex_data_degree_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const vertex_data_degree_visitor& v1, const vertex_data_degree_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  VertexData vertex_data;
  bool do_update_vertex_data;
};

template <typename TGraph, typename VertexMetaData, typename Vertex, 
  typename VertexData>
void vertex_data_db_degree(TGraph* g, VertexMetaData& vertex_metadata) {

  int mpi_rank = comm_world().rank();

  if (mpi_rank == 0) {
    std::cout << "Building Distributed Vertex Metadata Store (Vertex Degree) ... " << std::endl;
  }

  typedef vertex_data_degree_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  vq.init_visitor_traversal();
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "Done Building Vertex Metadata Store." << std::endl;
  }   
}

} //end namespace havoqgt
