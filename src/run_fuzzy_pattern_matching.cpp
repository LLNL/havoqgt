#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/fuzzy_pattern_matching.hpp>
#include <havoqgt/graph.hpp>
//#include <havoqgt/label_propagation_pattern_matching.hpp>
#include <havoqgt/label_propagation_pattern_matching_bsp.hpp> 
#include <havoqgt/pattern_util.hpp>
#include <havoqgt/token_passing_pattern_matching.hpp>
#include <havoqgt/vertex_data_db.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
using namespace havoqgt;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;

template<typename T>
  using SegmentAllocator = bip::allocator<T, segment_manager_t>;

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

  // types used by the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  //typedef typename graph_type::edge_iterator eitr_type;

  // user defined types
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint16_t VertexData; // assuming metadata is a 16-bit uint
  typedef uint64_t VertexRankType;

  typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetaData;
  typedef graph_type::vertex_data<VertexRankType, SegmentAllocator<VertexRankType> > VertexRank;
  typedef graph_type::vertex_data<bool, SegmentAllocator<bool> > VertexActive;
  typedef graph_type::vertex_data<uint64_t, SegmentAllocator<uint64_t> > VertexIteration;

  typedef vertex_state<uint64_t> VertexState;
  typedef std::unordered_map<Vertex, VertexState> VertexStateMap;

  if(mpi_rank == 0) {
    std::cout << "Distributed fuzzy pattern matching." << std::endl;
  }

  // per-rank containers
  VertexStateMap vertex_state_map; 

  // vertex containers
  VertexMetaData vertex_metadata(*graph, alloc_inst);
  VertexRank vertex_rank(*graph, alloc_inst);
  VertexActive vertex_active(*graph, alloc_inst);
  VertexIteration vertex_iteration(*graph, alloc_inst);  

  MPI_Barrier(MPI_COMM_WORLD); 
 
  // build the distributed vertex data db
  // each rank reads 10K lines at a time
  vertex_data_db<graph_type, VertexMetaData, Vertex, VertexData>
    (graph, vertex_metadata, vertex_data_input_filename, 10000);

  // barrier 

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
  if(mpi_rank == 0) { 
    std::cout << "Setting up patterns to search ... " << std::endl;
  }
  pattern_util<VertexData> ptrn_util_two(pattern_input_filename, true);

  for(size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {

  auto pattern = std::get<0>(ptrn_util_two.input_patterns[pl]);
  auto pattern_indices = std::get<1>(ptrn_util_two.input_patterns[pl]);

  if(mpi_rank == 0) {
    std::cout << "[" << pl << "] Searching pattern: "; 
    pattern_util<VertexData>::output_pattern(pattern);
  } // for

  // setup patterns - version 2
  typedef ::graph<Vertex, Edge, VertexData> PatternGraph; // TODO: ::graph class name conflict!
  PatternGraph pattern_graph(pattern_input_filename + "_edge",
    pattern_input_filename + "_vertex",
    pattern_input_filename + "_vertex_data", true, true);

  // test print
  if(mpi_rank == 0) {
  std::cout << "Searching pattern (version 2): " << std::endl;
  for (Vertex v = 0; v < pattern_graph.vertex_count; v++) {
    std::cout << v << " : " << pattern_graph.vertices[v] << " " << pattern_graph.vertex_data[v] 
    << " " << pattern_graph.vertex_degree[v] << std::endl;
    for (auto e = pattern_graph.vertices[v]; e < pattern_graph.vertices[v + 1]; e++) {
    auto v_nbr = pattern_graph.edges[e];
      std::cout << v_nbr << ", " ;
    }
    std::cout << std::endl; 
  } 
  }
  // test print

  // initialize containers
  vertex_rank.reset(0);
  vertex_active.reset(true);
  vertex_iteration.reset(0); // TODO: -1 ?
  std::cout << "Vertex state map size (initially): " << vertex_state_map.size() << std::endl; // Test 

  MPI_Barrier(MPI_COMM_WORLD);

  double time_start = MPI_Wtime();
 
  // run application

  // clone pattern matchng
  //fuzzy_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank);

  // label propagation pattern matching 
//  label_propagation_pattern_matching<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    VertexRank, VertexActive, VertexIteration, VertexStateMap, PatternGraph>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph);

  // label propagation pattern matching bsp 
//  label_propagation_pattern_matching_bsp<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    VertexRank, VertexActive, VertexIteration, VertexStateMap, PatternGraph>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph);

  // toekn passing
  token_passing_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank, pattern_graph);  

  // barrier
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

  double time_end = MPI_Wtime();

  if(mpi_rank == 0) {
    std::cout << "Fuzzy Pattern Matching Time | label Propagation = " << time_end - time_start << std::endl;
  }    

  std::cout << "Vertex state map size (finally): " << vertex_state_map.size() << std::endl; // Test 

  // test print
  uint64_t vertex_active_count = 0;
  uint64_t vertex_inactive_count = 0;
  uint64_t vertex_iteration_valid_count = 0;  
  for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_iteration[vertex] >= 2*pattern_graph.vertex_data.size() + 1) { 
      //std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_iteration[vertex] << std::endl;
      vertex_iteration_valid_count++;
    }   
    if (vertex_rank[vertex] > 0) {
      std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl;
    }
    if (!vertex_active[vertex]) {
      //std::cout << mpi_rank << " l " << graph->locator_to_label(vertex) << " " << vertex_active[vertex] << std::endl;
      vertex_inactive_count++;
    } else {
      vertex_active_count++;
    }
  }
  
  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex_iteration[vertex] >= 2*pattern_graph.vertex_data.size() + 1) {
      //std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_iteration[vertex] << std::endl;
      vertex_iteration_valid_count++;
    }
    if (vertex_rank[vertex] > 0) {
      std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_rank[vertex] << std::endl; 
    }  
    if (!vertex_active[vertex]) {   
      //std::cout << mpi_rank << " d " << graph->locator_to_label(vertex) << " " << vertex_active[vertex] << std::endl;
      vertex_inactive_count++;
    } else {
      vertex_active_count++;
    }
  }
  std::cout << mpi_rank << " # active vertices " << vertex_active_count << std::endl;
  std::cout << mpi_rank << " # inactive vertices " << vertex_inactive_count << std::endl; 
  std::cout << mpi_rank << " # vertices reached max-iterations " << vertex_iteration_valid_count << std::endl;
  // test print
  
  // cleanup memeory
  //vertex_rank.clear(); // TODO: add clear() method to vertex_data.cpp   
  //vertex_active.clear(); // TODO: add clear() method to vertex_data.cpp
  vertex_state_map.clear();

  // toekn passing
  //token_passing_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank);
  //

  } // for - loop over query patterns

  } // fuzzy pattern matching  

  } // havoqgt_init
  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;  
}
