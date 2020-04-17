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
