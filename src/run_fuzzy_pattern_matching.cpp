#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/fuzzy_pattern_matching.hpp>
#include <havoqgt/pattern_util.hpp>
#include <havoqgt/vertex_data_db.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
using namespace havoqgt;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
template<typename T>
    using SegmentAllocator = bip::allocator<T, segment_manager_t>;

template <typename Vertex, typename VertexData, typename VertexDataList, typename VertexEntryList>
    void read_vertex_data_list(std::string vertex_data_input_filename,
      VertexDataList& vertex_data, VertexEntryList& vertex_entry) {
      std::ifstream vertex_data_input_file(vertex_data_input_filename,
        std::ifstream::in);
      std::string line;
      while (std::getline(vertex_data_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        VertexData v_data(0);
        iss >> v_source >> v_data;
        vertex_data.push_back(v_data);
        vertex_entry.push_back(std::forward_as_tuple(v_source, v_data)); 
      }
      vertex_data_input_file.close();
}

int main(int argc, char** argv) {
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

  // TODO: parse commandline

  MPI_Barrier(MPI_COMM_WORLD);
  if(backup_filename.size() > 0) {
    distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  }

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  segment_manager_t* segment_manager = ddb.get_segment_manager();
    bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  graph_type *graph = segment_manager->
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

  // user defined types for the application
  typedef uint64_t Vertex;
  typedef uint16_t VertexData; // assuming metadata is a 16-bit uint
  typedef uint64_t VertexRankType;

  // types used by the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  //typedef typename graph_type::edge_iterator eitr_type;

  typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetaData;
  typedef graph_type::vertex_data<VertexRankType, SegmentAllocator<VertexRankType> > VertexRank; 
  typedef std::vector<std::tuple<Vertex, VertexData>> VertexEntry;

  std::cout << "Distributed fuzzy pattern matching." << std::endl;

  // vertex containers
  VertexMetaData vertex_metadata(*graph, alloc_inst);
  VertexRank vertex_rank(*graph, alloc_inst); 
 
  // build the distributed vertex data db
  vertex_data_db<graph_type, VertexMetaData, Vertex, VertexData>
    (graph, vertex_metadata, vertex_data_input_filename, 10000); // each rank reads 10K lines at a time
 
  // test print
/*  int set_mpi_rank = 4;
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    if (mpi_rank == set_mpi_rank)
      std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_metadata[vertex] << std::endl; 
  }
  
  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    if (mpi_rank == set_mpi_rank)
      std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_metadata[vertex] << std::endl; 
  }*/ 
  // test print
  
  // setup patterns
  std::cout << "Setting up patterns to search ... " << std::endl;
  pattern_util<VertexData> ptrn_util_two(pattern_input_filename, true);

  for(size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {

  auto pattern = std::get<0>(ptrn_util_two.input_patterns[pl]);
  auto pattern_indices = std::get<1>(ptrn_util_two.input_patterns[pl]);

  if(mpi_rank == 0) {
    std::cout << "[" << pl << "] Searching pattern: "; 
    pattern_util<VertexData>::output_pattern(pattern);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double time_start = MPI_Wtime();

  vertex_rank.reset(0);
  
  // run application
  fuzzy_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank);

  MPI_Barrier(MPI_COMM_WORLD);

  double time_end = MPI_Wtime();

  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time = " << time_end - time_start << std::endl;
  }     

  // test print
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_rank[vertex] > 0)
    std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl; 

  }
  
  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_rank[vertex] > 0)
    std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl; 
  }
  // test print 

  } // for

  } // fuzzy pattern matching  

  } // havoqgt_init
  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;  
}
