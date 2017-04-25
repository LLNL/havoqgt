#include <include/havoqgt_setup.hpp>

graph_type* read_delegate_graph(std::string graph_input, const int mpi_rank) {
  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  graph_type* graph =  ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
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
