// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <include/havoqgt_setup.hpp>

graph_type* read_delegate_graph(std::string graph_input, const int mpi_rank) {
  distributed_db ddb(db_open_read_only(), graph_input.c_str());

  auto graph =  ddb.get_manager()->find<graph_type>("graph_obj").first;
  std::cout << "MPI Rank: " << mpi_rank << " Loading Graph: " << graph_input 
  << "." <<  std::endl; 
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "MPI Rank: " << mpi_rank << " Graph Loaded Ready." << std::endl;
  }
  //graph->print_graph_statistics(); 
  // causes MPI runtime exception on flash.llnl.gov
  MPI_Barrier(MPI_COMM_WORLD);

//  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
//       ++vitr) {
//    vloc_type vertex = *vitr;
//    std::cout << "MPI Rank -> " << mpi_rank << " local vertex " 
//    << graph->locator_to_label(vertex) << std::endl;
//  }
  return graph;
}
