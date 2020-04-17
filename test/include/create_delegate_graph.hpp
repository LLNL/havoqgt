#include <fstream>
#include <include/havoqgt_setup.hpp>
#include <include/util.hpp>

bool create_edge_list_file(
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& input_graph,
  std::string output_filename) {

  std::ofstream output_file(output_filename, std::ofstream::out);   
  for (auto itr = input_graph.begin(); itr != input_graph.end(); ++itr) {
    output_file << std::get<0>(*itr) << " " 
                << std::get<1>(*itr) << " "
                << std::get<2>(*itr) << "\n";     
  }  
  output_file.flush();
  output_file.close();     

  std::ifstream input_file(output_filename);
  if (input_file.good()) {
    input_file.close();
    return SUCCESS;
  } 
  else return FAILURE;
}

//template <typename edge_data_type>
void create_delegate_graph(
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& input_graph, 
  std::string graph_filename, std::string graph_unique_instance_name, 
  std::string edge_data_unique_instance_name, int mpi_rank) {

  if (mpi_rank == 0) {
     //assert(create_edge_list_file(input_graph, graph_filename));
     create_edge_list_file(input_graph, graph_filename);
  } 
  
  MPI_Barrier(MPI_COMM_WORLD); 

  uint64_t                   delegate_threshold = 4;
  uint64_t                   partition_passes = 1;
  double                     gbyte_per_rank = 0.25;
  uint64_t                   chunk_size = 8*1024;
  bool                       undirected = false;   
 
  std::vector< std::string > input_filenames;
  input_filenames.clear();
  input_filenames.push_back(graph_filename);
  
  std::string                output_filename = graph_filename;

  if (mpi_rank == 0) {
      std::cout << "Ingesting graph from " << input_filenames.size() 
      << " files." << std::endl;
  }

  havoqgt::distributed_db ddb(havoqgt::db_create(), 
    output_filename.c_str(), gbyte_per_rank);  

  segment_manager_t* segment_manager = ddb.get_segment_manager();
  
  bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  graph_type::edge_data<edge_data_type, 
    bip::allocator<edge_data_type, segment_manager_t>> edge_data(alloc_inst); 

  //Setup edge list reader
  havoqgt::parallel_edge_list_reader<edge_data_type> 
    pelr(input_filenames, undirected);
  bool has_edge_data = pelr.has_edge_data(); 

  if (mpi_rank == 0) {
    std::cout << "Generating new graph." << std::endl;
  } 

  graph_type *graph = segment_manager->construct<graph_type>
    (graph_unique_instance_name.c_str())
    (alloc_inst, MPI_COMM_WORLD, pelr, pelr.max_vertex_id(), 
     delegate_threshold, partition_passes, chunk_size, edge_data);

  if (has_edge_data) {
      graph_type::edge_data<edge_data_type, 
      bip::allocator<edge_data_type, segment_manager_t>>* edge_data_ptr
      = segment_manager->construct<graph_type::edge_data<edge_data_type, 
      bip::allocator<edge_data_type, segment_manager_t>>>
        (edge_data_unique_instance_name.c_str())
        (edge_data);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "Graph Ready, Calculating Stats. " << std::endl;
  }

  //for (int i = 0; i < mpi_size; i++) {
  //  if (i == mpi_rank) {
  //    double percent = double(segment_manager->get_free_memory()) /
  //    double(segment_manager->get_size());
  //    std::cout << "[" << mpi_rank << "] " << segment_manager->get_free_memory()
  //              << "/" << segment_manager->get_size() << " = " << percent << std::endl;
  //  }
  //  MPI_Barrier(MPI_COMM_WORLD);
  //}
 
  MPI_Barrier(MPI_COMM_WORLD);

  if(mpi_rank == 0) {
    sync();
  }

  MPI_Barrier(MPI_COMM_WORLD);  
}

