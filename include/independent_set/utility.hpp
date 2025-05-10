#pragma once

#include <math.h>

#include <cstdint>
#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_map>

#include <havoqgt/distributed_db.hpp>

namespace maximal_independent_set {

const uint64_t in = std::numeric_limits<uint64_t>::max();
const uint64_t out = std::numeric_limits<uint64_t>::min();
const double float_factor = 0.00000000023283064365386962890625l;

size_t average_degree;
double scaled_average_degree;

uint64_t hash(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return x;  
} // hash

template <typename Graph, typename VertexID, 
  typename VertexPriority>
void graph_statistic(Graph* graph) {

  typedef typename Graph::vertex_locator vloc_type;	
  typedef typename Graph::vertex_iterator vitr_type;
  typedef typename Graph::controller_iterator citr_type;
  typedef typename Graph::edge_iterator eitr_type;   

  int mpi_rank = havoqgt::comm_world().rank();

  // vertex and edge count
   
  size_t local_vertex_count(0);
  size_t local_edge_count(0);
  size_t local_max_degree(0);

  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    ++local_vertex_count;
    local_edge_count+=graph->degree(vertex);
    if (graph->degree(vertex) > local_max_degree) {
      local_max_degree = graph->degree(vertex);	    
    } 	    
  } // for	  

  for (citr_type citr = graph->controller_begin();
    citr != graph->controller_end(); ++citr) {
    vloc_type vertex = *citr;
    ++local_vertex_count;
    local_edge_count+=graph->degree(vertex);
    if (graph->degree(vertex) > local_max_degree) {
      local_max_degree = graph->degree(vertex);	    
    } 
  } // for

  const size_t vertex_count = havoqgt::mpi_all_reduce(local_vertex_count,
    std::plus<size_t>(), MPI_COMM_WORLD);
  const size_t edge_count = havoqgt::mpi_all_reduce(local_edge_count,
    std::plus<size_t>(), MPI_COMM_WORLD);
  const size_t max_degree = havoqgt::mpi_all_reduce(local_max_degree,
    std::greater<size_t>(), MPI_COMM_WORLD); 		  

  average_degree = static_cast<size_t>(ceill(edge_count / vertex_count));
  scaled_average_degree = ((maximal_independent_set::in / 2) - 1) * average_degree;

  if (mpi_rank == 0) {
    std::cout << "#Vertices: " << vertex_count << std::endl;
    std::cout << "#Edges: " << edge_count << std::endl;
    std::cout << "Max degree: " << max_degree << std::endl;
    std::cout << "Average degree: " << average_degree << std::endl;
    std::cout << "Scaled average degree: " << scaled_average_degree <<
    std::endl;
  }  
} // graph_statistic

template <typename VertexID, typename EdgeID, typename VertexPriority>
VertexPriority vertex_priority_ecl(VertexID vertex_ID, EdgeID vertex_degree) {
  if (vertex_degree > 0) {
    auto x = vertex_degree - 
      (maximal_independent_set::hash(vertex_ID) * 
      maximal_independent_set::float_factor);	    
    VertexPriority res = 
      static_cast<VertexPriority>
      (maximal_independent_set::scaled_average_degree / 
      (maximal_independent_set::average_degree + x));
    return (res + res) | 1;  
  } else { 	  
    return maximal_independent_set::in;
  }
} // vertex_priority_ecl	

template <typename VertexID, typename EdgeID, typename VertexPriority>
VertexPriority vertex_priority(VertexID vertex_ID, EdgeID vertex_degree) {
  return maximal_independent_set::vertex_priority_ecl
    <VertexID, EdgeID, VertexPriority>(vertex_ID, vertex_degree);
} // vertex_priority	

} // end namespace maximal_independent_set	
