// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT


#ifndef HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED
#define HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt {

template <typename TGraph, typename DegreeData>
void compute_degree_no_selfloops(TGraph& graph, DegreeData& degree) {
  degree.reset(0);
  //reset graph data
  for(auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    for(auto eitr = graph.edges_begin(*vitr); eitr != graph.edges_end(*vitr); ++eitr) {
      if(eitr.source() != eitr.target()) {
        degree[*vitr]++;
      }
    }
  }
  for(auto ditr = graph.delegate_vertices_begin(); ditr != graph.delegate_vertices_end(); ++ditr) {
    for(auto eitr = graph.edges_begin(*ditr); eitr != graph.edges_end(*ditr); ++eitr) {
      if(eitr.source() != eitr.target()) {
        degree[*ditr]++;
      }
    }
  }
  degree.all_reduce();
  /*
  {  // Computes & prints number of self loops found
  uint64_t local_count(0);
  for (auto vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    if(deg_no_self[*vitr] > graph.degree(*vitr)) {
      std::cerr << "LOGIC ERROR" << std::endl;  exit(-1);
    }
    local_count += graph.degree(*vitr) - deg_no_self[*vitr];
  }
  for (auto citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    if(deg_no_self[*citr] > graph.degree(*citr)) {
      std::cerr << "LOGIC ERROR" << std::endl;  exit(-1);
    }
    local_count += graph.degree(*citr) - deg_no_self[*citr];
  }
  
  uint64_t global_count = mpi_all_reduce(local_count,std::plus<uint64_t>(), MPI_COMM_WORLD);
  if(mpi_rank == 0) {
    std::cout << "Number of self loops = " << global_count << std::endl;
  }
  }
  */
  
}


  //kth_core(TGraph& graph, typename TGraph::vertex_data<bool>& , uint32_t k)
  
  //kth_core(TGraph& graph, typename TGraph::vertex_data<uint32_t>& , )
  
  //log2_kth_core(TGraph& graph, typename TGraph::vertex_data<uint8_t>& log2_core)


template <typename TGraph>
void kth_core(TGraph& graph) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  

  
  typename TGraph::template vertex_data<uint32_t, std::allocator<uint32_t>>  deg_no_self(graph);
  compute_degree_no_selfloops(graph, deg_no_self);
  


}


} // end havoqgt


#endif //HAVOQGT_MPI_KTH_CORE_HPP_INCLUDED
