#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/pattern_util.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/fuzzy_pattern_matching.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
using namespace havoqgt;

template <typename Vertex, typename VertexData, typename VertexDataList>
    void read_vertex_data_list(std::string vertex_data_input_filename,
      VertexDataList& vertex_data) {
      std::ifstream vertex_data_input_file(vertex_data_input_filename,
        std::ifstream::in);
      std::string line;
      while (std::getline(vertex_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        VertexData v_data(0);
        iss >> v_source >> v_data;
        vertex_data.push_back(v_data);
      }
      vertex_data_input_file.close();
}

int main(int argc, char** argv) {
  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
  havoqgt::havoqgt_init(&argc, &argv);
  {

  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);  

  std::string graph_input = argv[1];
  std::string backup_filename; 

  // for fuzzy pattern matching
  std::string vertex_data_input_filename = argv[2];
  std::string pattern_input_filename = argv[3];
  std::string vertex_rank_output_filename = argv[4];

  // parse commandline

  MPI_Barrier(MPI_COMM_WORLD);
  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }

  //graph->print_graph_statistics(); // causes MPI error
  MPI_Barrier(MPI_COMM_WORLD);  

  // fuzzy pattern matching
  {
  // types for reading from file
  typedef uint64_t Vertex;
  typedef uint16_t VertexData;   // assuming metadata is a 16-bit uint
  typedef uint64_t VertexRank;

  // types used the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  //typedef typename graph_type::edge_iterator eitr_type;

  typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetaData; 

  std::cout << "Distributed fuzzy pattern matching." << std::endl;

  // read vertex data from file
  std::cout << "Reading vertex data list ... " << std::endl;
  std::vector<VertexData> vertex_metadata_all; // for now
  read_vertex_data_list<Vertex, VertexData>(vertex_data_input_filename, vertex_metadata_all);  
  std::cout << "Size of vertex data list: " << vertex_metadata_all.size() << std::endl;

  // vertex metadata container
  VertexMetaData vertex_metadata(*graph);
  //graph_type::vertex_data<VertexRank, std::allocator<VertexRank> >  vertex_rank(*graph); 
  
  // iterate over the vertices and populate metadata
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();       
    ++vitr) {
    vloc_type vertex = *vitr;
    vertex_metadata[vertex] = vertex_metadata_all[graph->locator_to_label(vertex)];
  }

  for (vitr_type vitr = graph->delegate_vertices_begin(); 
    vitr != graph->delegate_vertices_end(); ++vitr) { 
    vloc_type vertex = *vitr;
    vertex_metadata[vertex] = vertex_metadata_all[graph->locator_to_label(vertex)];
  }

  // setup patterns
  std::cout << "Setting up patterns to search ... " << std::endl;
  pattern_util<VertexData> ptrn_util_two(pattern_input_filename, true);
  auto pattern = std::get<0>(ptrn_util_two.input_patterns[3]);
  auto pattern_indices = std::get<1>(ptrn_util_two.input_patterns[3]);

  MPI_Barrier(MPI_COMM_WORLD);

  double time_start = MPI_Wtime();
  
  // run application
  fuzzy_pattern_matching(graph, vertex_metadata, pattern, pattern_indices);

  MPI_Barrier(MPI_COMM_WORLD);

  double time_end = MPI_Wtime();

  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time = " << time_end - time_start << std::endl;
  }     

  } // fuzzy pattern matching  

  } // havoqgt_init
  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;  
}
