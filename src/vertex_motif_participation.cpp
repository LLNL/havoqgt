// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <bitset>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/container/vector.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
///#include <havoqgt/environment.hpp>

#include <metall/metall.hpp>
#include <metall_utility/metall_mpi_adaptor.hpp>

#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>

#include <prunejuice/template.hpp>
#include <prunejuice/nonlocal_constraint.hpp>
#include <prunejuice/algorithm_state.hpp>
#include <prunejuice/local_constraint_checking.hpp>
#include <prunejuice/nonlocal_constraint_checking_basic.hpp>
#include <prunejuice/nonlocal_constraint_checking_tds_batch.hpp>
#include <prunejuice/motif_container.hpp>

#define OUTPUT_RESULT
 
using namespace havoqgt;

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

template<typename T>
  using DelegateGraphVertexDataSTDAllocator = graph_type::vertex_data
  <T, std::allocator<T>>; 

template<typename T>
  using DelegateGraphEdgeDataSTDAllocator = graph_type::edge_data
  <T, std::allocator<T>>;  

void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -p <string> -o <string>\n"
      << " -i <string>   - input graph base filename (required)\n"
      << " -b <string>   - backup graph base filename. If set, \"input\" graph will be deleted if it exists\n"
      << " -v <string>   - vertex metadata base filename (optional, default is degree based metadata)\n"
      << " -e <string>   - edge metadata base filename (optional, N/A)\n" 
      << " -p <string>   - pattern base directory (required)\n"
      << " -o <string>   - output base directory (optional)\n"
      << " -x <int>      - batch count (optional, default/min batch count is 1, max batch count is "
        << comm_world().size() << "\n"   
      << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_input, 
  std::string& backup_graph_input, std::string& vertex_metadata_input, 
  std::string& edge_metadata_input, std::string& pattern_input, 
  std::string& result_output, uint64_t& tp_batch_size) {

  if(comm_world().rank() == 0) {
    std::cout << "CMD Line :";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool print_help = false;
  std::bitset<3> required_input;
  required_input.reset();
  required_input.set(2); // TODO: improve

  char c;
  while ((c = getopt(argc, argv, "i:b:v:e:p:o:x:h ")) != -1) {
    switch (c) {
      case 'h' :  
        print_help = true;
        break;  
      case 'i' :
        graph_input = optarg;
        required_input.set(0);
        break;
      case 'b' :
        backup_graph_input = optarg;
        break;			
      case 'v' :
        vertex_metadata_input = optarg;
        break;
      case 'e' :
        edge_metadata_input = optarg;
        break;					
      case 'p' :
        pattern_input = optarg;
        required_input.set(1);   
        break;
      case 'o' :         
        result_output = optarg;
        required_input.set(2);
        break;
      case 'x' :		 	 
        tp_batch_size = std::stoull(optarg);
        if (tp_batch_size < 1 || tp_batch_size > comm_world().size()) {
          print_help = true;
        } else if (tp_batch_size > 1) {
          tp_batch_size = comm_world().size() / tp_batch_size;
        } else {
          tp_batch_size = comm_world().size();
        } 		
        break;
      default:
        std::cerr << "Unrecognized Option : " << c << ", Ignore."<<std::endl;
        print_help = true;
        break;
    } 		
  }

  if (print_help || !required_input.all()) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {
 
  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
   
  havoqgt::init(&argc, &argv);
  {

  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  //havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI Initialized With " << mpi_size << " Ranks" << std::endl;
    //havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);  

  //////////////////////////////////////////////////////////////////////////////

  // parse command line

  std::string graph_input;
  std::string backup_graph_input;
  std::string vertex_metadata_input;
  std::string edge_metadata_input;
  std::string pattern_input;
  std::string result_output;  

  uint64_t tp_batch_size = comm_world().size();

  parse_cmd_line(argc, argv, graph_input, backup_graph_input, 
    vertex_metadata_input, edge_metadata_input, pattern_input, result_output, 
    tp_batch_size); 

  std::string pattern_dir = pattern_input; 
  std::string result_dir = "/todo/"; //result_output;
  std::string motif_result_dir = result_output;

  MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////

  // load graph

  double hgt_time_start = MPI_Wtime();
  double hgt_time_end = MPI_Wtime();

  if (mpi_rank == 0) {
    std::cout << "HavoqGT | Loading Graph ... " << std::endl;
  }

  if (backup_graph_input.size() > 0) {
    distributed_db::transfer(backup_graph_input.c_str(), graph_input.c_str());
  }

  distributed_db ddb(db_open_read_only(), graph_input.c_str());

  auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  // edge data
  //if (mpi_rank == 0) {
    //std::cout << "Loading / Initializing Edge Data ... " << std::endl;
  //}

  typedef uint8_t edge_data_type; 
  // TODO: figure out a way to get it from graph_type;
  // see edge_data_value_type in parallel_edge_list_reader.hpp

  typedef graph_type::edge_data<edge_data_type, 
    distributed_db::allocator<edge_data_type>> edge_data_t;

  auto edge_data_ptr = ddb.get_manager()->find<edge_data_t>
    ("graph_edge_data_obj").first;
  //assert(edge_data_ptr != nullptr); 

  MPI_Barrier(MPI_COMM_WORLD);
  hgt_time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "HavoqGT | Done Loading Graph" << std::endl;
    std::cout << "HavoqGT Time | Graph Loading : "
      << hgt_time_end - hgt_time_start << std::endl;
  }

  //graph->print_graph_statistics(); // causes MPI error
  //MPI_Barrier(MPI_COMM_WORLD);  

  //////////////////////////////////////////////////////////////////////////////

  // pattern matching
   
  {

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching ... " << std::endl;
  }

  // types used by the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  typedef typename graph_type::edge_iterator eitr_type;
 
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData; // for string hash
  //typedef uint8_t VertexData; // for log binning 
  //typedef edge_data_type EdgeData;
  typedef uint8_t EdgeData; 

  static constexpr size_t max_bit_vector_size = 16; // TODO:
  static constexpr size_t max_template_vertex_count = 16;  
  typedef std::bitset<max_bit_vector_size> BitSet; // TODO: rename to TemplateVertexSet
  typedef BitSet TemplateVertexBitSet; 
  typedef uint16_t TemplateVertexType; // TODO: rename to TemplateVertexBitSetToUint
  
  typedef uint8_t Boolean; // TODO: replace all bool with Boolean ?

  // TODO: mmap
  //typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetadata;
  //typedef graph_type::vertex_data<bool, SegmentAllocator<bool> > VertexActive;

  typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetadata;
  typedef graph_type::vertex_data<Boolean, std::allocator<Boolean> > VertexActive; 
  // TODO: solution_graph; you are mixing bool and uint, fix it
  typedef graph_type::vertex_data<TemplateVertexType, std::allocator<TemplateVertexType> > TemplateVertex; 
  // TODO: solution_graph, rename to VertexTemplateVertexBitSetToUint

  typedef prunejuice::vertex_state<Vertex, VertexData, BitSet> VertexState;
  typedef std::unordered_map<Vertex, VertexState> VertexStateMap; // TODO: solution_graph

  typedef std::unordered_set<Vertex> VertexSet;  
  typedef graph_type::vertex_data<VertexSet, std::allocator<VertexSet> > VertexSetCollection; 
  
  typedef std::unordered_map<Vertex, uint8_t> VertexUint8Map; 
  typedef graph_type::vertex_data<VertexUint8Map, std::allocator<VertexUint8Map> > VertexUint8MapCollection;    

  typedef graph_type::edge_data<EdgeData, std::allocator<EdgeData> > EdgeMetadata;
  typedef graph_type::edge_data<Boolean, std::allocator<Boolean> > EdgeActive; // TODO: solution_graph 

  typedef std::vector<Boolean> VectorBoolean;

  ////////////////////////////////////////////////////////////////////////////// 

  double time_start = MPI_Wtime();
  double time_end = MPI_Wtime();

  // per rank containers    
  VertexStateMap vertex_state_map; 

  // vertex containers    
  // TODO: need a new alloc_inst to use bip/mmap
  // VertexMetadata vertex_metadata(*graph, alloc_inst);
  // VertexActive vertex_active(*graph, alloc_inst);
 
  VertexMetadata vertex_metadata(*graph); 
  VertexActive vertex_active(*graph);
  TemplateVertex template_vertices(*graph);
  VertexUint8MapCollection vertex_active_edges_map(*graph);
  VertexSetCollection vertex_token_source_set(*graph); // per vertex set

  // edge containers
  //EdgeMetadata edge_metadata(*graph);
  //EdgeActive edge_active(*graph);

  uint8_t vertex_rank; // TODO: dummy, remove
  uint8_t vertex_iteration; // TODO: dummy, remove  

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching | Allocated Vertex and Edge Containers" 
    << std::endl;
  }

  ////////////////////////////////////////////////////////////////////////////// 

  // application parameters // TODO: command line input

  bool do_output_metadata = true; //false;

  MPI_Barrier(MPI_COMM_WORLD);
 
  ////////////////////////////////////////////////////////////////////////////// 

  // build the distributed vertex data store
  
  time_start = MPI_Wtime();

  if (vertex_metadata_input.size() > 0) {
    vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata, vertex_metadata_input, 10000);
      // TODO: each rank reads 10K lines from file at a time
  } else {
    vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata);
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Vertex Data Store : " 
      << time_end - time_start << std::endl;
  }

  if (do_output_metadata) {
    std::string vertex_data_filename = result_dir +  
      "/all_ranks_vertex_data/vertex_data_" + std::to_string(mpi_rank);
    std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

    for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
      ++vitr) {
      vloc_type vertex = *vitr;
      vertex_data_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
        << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  

      //std::cout << mpi_rank << ", l, " << graph->locator_to_label(vertex)
      //  << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n"; // Test
    } 	
    	
    for(vitr_type vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vloc_type vertex = *vitr;

      if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
        vertex_data_file << mpi_rank << ", c, " << graph->locator_to_label(vertex) 
          << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
      } else {	
        vertex_data_file << mpi_rank << ", d, " << graph->locator_to_label(vertex) 
          << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
      }
    }
  
    vertex_data_file.close();
  }
 
  ////////////////////////////////////////////////////////////////////////////// 
 
  // result
   
  std::string pattern_set_result_filename = result_dir + "/result_pattern_set";
  std::ofstream pattern_set_result_file;  
  if (mpi_rank == 0) {
    pattern_set_result_file = std::ofstream(pattern_set_result_filename, 
      std::ofstream::out);
  }

  //////////////////////////////////////////////////////////////////////////////

  // setup pattrens   

  std::string pattern_input_filename = pattern_dir + "/pattern";
   
  typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData,
      EdgeData> PatternGraph;

  size_t pattern_count = PatternGraph::read_vertex_data_list_2
    (pattern_input_filename + "_vertex_data_list");

//  size_t ps = 0;

  // loop over the pattern set
//  for (size_t ps = 0; ps < pattern_count; ps++) { 
   
    // setup pattern
    if(mpi_rank == 0) { 
//      std::cout << "Setting up Pattern [" << ps << "] ... " << std::endl;
    }
     
    // setup pattern for local constraint checking 
    //std::string pattern_input_filename = pattern_dir + "/pattern";
    
    //typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData, 
    //  EdgeData> PatternGraph;
    PatternGraph pattern_graph(
      pattern_input_filename + "_edge",
      pattern_input_filename + "_vertex",
      pattern_input_filename + "_vertex_data",
      pattern_input_filename + "_edge_data",
      pattern_input_filename + "_stat",
      false, false); // TODO: improve

//    pattern_graph.set_vertex_data(ps);

    // TODO: create a graphical representation of the pattern 
    // using a 2D graphics library
    // Test print
    /*if(mpi_rank == 0) {
      std::cout << "Pattern Matching | Searching Pattern [" << ps
	<< "] : " << std::endl;
      for (auto v = 0; v < pattern_graph.vertex_count; v++) {
	std::cout << v << " : offset " << pattern_graph.vertices[v] 
	<< ", vertex_data " << pattern_graph.vertex_data[v] 
	<< ", vertex_degree " << pattern_graph.vertex_degree[v] << std::endl;
	std::cout << " neighbors : "; 
	for (auto e = pattern_graph.vertices[v]; 
	  e < pattern_graph.vertices[v + 1]; e++) {
	  auto v_nbr = pattern_graph.edges[e];
	  std::cout << v_nbr << ", " ;
	}
	std::cout << std::endl;
	std::cout << " neighbor vertex data count : ";
	for (auto& nd : pattern_graph.vertex_neighbor_data_count_map[v]) {
	  std::cout << "(" << nd.first << ", " << nd.second << "), ";
	}
	std::cout << std::endl; 
      }
      //std::cout << "diameter : " << pattern_graph.diameter << std::endl; 
    }*/
    // Test print

    // setup pattern for token passing
    // TODO: remove from here ?    
   
    typedef prunejuice::pattern_nonlocal_constraint<Vertex, Edge, VertexData, 
      PatternGraph> PatternNonlocalConstraint;

    PatternNonlocalConstraint ptrn_util_two(pattern_graph,
      pattern_dir + "/pattern_nonlocal_constraint");

  //////////////////////////////////////////////////////////////////////////////
  
  // motif

  {
  //std::string motif_db_name = "reddit_motifs";

  using MotifContainer = prunejuice::container::motif_container<Vertex, VertexData,
    uint64_t, metall::manager::allocator_type<char>>;
   
  std::string motif_db_dir = "/p/lustre2/reza2/graph_learning/reddit/metall_motif_collection/motif_A_"; // TODO
  std::string motif_db_local_partiton_dir = motif_db_dir + std::to_string(mpi_rank);

  //metall::manager manager(metall::create_only, motif_db_local_partiton_dir.c_str());
  //MotifContainer *motif_db =
  //  manager.construct<MotifContainer>("motif_container")(manager.get_allocator()); 
 
  metall::manager manager(metall::open_only, motif_db_local_partiton_dir.c_str());
  MotifContainer *motif_db = manager.find<MotifContainer>("motif_container").first;

  MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////

  // loop over the pattern set
   
  for (size_t ps = 0; ps < pattern_count; ps++) { 

    pattern_graph.set_vertex_data(ps);
    //ptrn_util_two // TODO 

    if (ps >= 0) {   
    //if (ps == 39) {
    //if (ps == pattern_count - 1) {

    // Test print
    if(mpi_rank == 0) {
      std::cout << "Pattern Matching | Searching Pattern [" << ps
	<< "] : " << std::endl;
      for (auto v = 0; v < pattern_graph.vertex_count; v++) {
	std::cout << v << " : offset " << pattern_graph.vertices[v] 
	<< ", vertex_data " << pattern_graph.vertex_data[v] 
	<< ", vertex_degree " << pattern_graph.vertex_degree[v] << std::endl;
	std::cout << " neighbors : "; 
	for (auto e = pattern_graph.vertices[v]; 
	  e < pattern_graph.vertices[v + 1]; e++) {
	  auto v_nbr = pattern_graph.edges[e];
	  std::cout << v_nbr << ", " ;
	}
	std::cout << std::endl;
	std::cout << " neighbor vertex data count : ";
	for (auto& nd : pattern_graph.vertex_neighbor_data_count_map[v]) {
	  std::cout << "(" << nd.first << ", " << nd.second << "), ";
	}
	std::cout << std::endl; 
      }
      //std::cout << "diameter : " << pattern_graph.diameter << std::endl; 
    }
    // Test print

    //MPI_Barrier(MPI_COMM_WORLD);
    //continue; 
    //return 0; 

    } else {
      MPI_Barrier(MPI_COMM_WORLD); 
      continue;
    }
 
  //////////////////////////////////////////////////////////////////////////////

    // per pattern initialization

    // initialize containers
    vertex_state_map.clear(); // Important
    vertex_active.reset(true); // initially all vertices are active
    //edge_active.reset(true); //  initially all edges are active
    vertex_active_edges_map.clear(); // Important
    vertex_token_source_set.clear(); // clear all the sets on all the vertices

    // initialize application parameters  
    bool global_init_step = true; // TODO: Boolean 
    bool global_not_finished = false; // TODO: Boolean

    bool do_nonlocal_constraint_checking = true; // TODO: Boolean

    uint64_t global_itr_count = 0;

    uint64_t active_vertices_count = 0;
    uint64_t active_edges_count = 0;

    uint64_t message_count = 0;

  ////////////////////////////////////////////////////////////////////////////// 

    // result

    //std::string itr_result_filename = result_dir +  
    //  "/result_iteration"; // TODO: improve
    //std::ofstream itr_result_file(itr_result_filename, std::ofstream::out);

    std::string step_result_filename = result_dir +  
      "/result_step";
    std::ofstream step_result_file(step_result_filename, std::ofstream::out); 

    std::string superstep_result_filename = result_dir + 
      "/result_superstep";
    std::ofstream superstep_result_file(superstep_result_filename, std::ofstream::out);

    std::string active_vertices_count_result_filename = result_dir +  
      "/all_ranks_active_vertices_count/active_vertices_" + std::to_string(mpi_rank); 
    std::ofstream active_vertices_count_result_file(active_vertices_count_result_filename, 
      std::ofstream::out);

    std::string active_vertices_result_filename = result_dir + 
      "/all_ranks_active_vertices/active_vertices_" + std::to_string(mpi_rank);
    std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);

    std::string active_edges_count_result_filename = result_dir + 
      "/all_ranks_active_edges_count/active_edges_" + std::to_string(mpi_rank);
    std::ofstream active_edges_count_result_file(active_edges_count_result_filename, std::ofstream::out);

    std::string active_edges_result_filename = result_dir + 
      "/all_ranks_active_edges/active_edges_" + std::to_string(mpi_rank);
    std::ofstream active_edges_result_file(active_edges_result_filename, std::ofstream::out);

    std::string message_count_result_filename = result_dir + 
      "/all_ranks_messages/messages_" + std::to_string(mpi_rank); // TODO:message_count
    std::ofstream message_count_result_file(message_count_result_filename, std::ofstream::out);

    MPI_Barrier(MPI_COMM_WORLD);

  //////////////////////////////////////////////////////////////////////////////

    // run constraint checking
    
    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Running Constraint Checking ... " << std::endl;
    }

    double pattern_time_start = MPI_Wtime();
    //double itr_time_start = MPI_Wtime();

  //////////////////////////////////////////////////////////////////////////////
  
    // local constraint checking  
      
    double label_propagation_time_start = MPI_Wtime();
    double label_propagation_time_end = MPI_Wtime();

    prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData, 
      graph_type, VertexMetadata, VertexStateMap, VertexActive, 
      VertexUint8MapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>
      (graph, vertex_metadata, vertex_state_map, vertex_active, 
      vertex_active_edges_map, template_vertices, pattern_graph, global_init_step, 
      global_not_finished, global_itr_count, superstep_result_file, 
      active_vertices_count_result_file, active_edges_count_result_file,
      message_count_result_file);

    MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
    label_propagation_time_end = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "Pattern Matching Time | Local Constraint Checking : " 
	<< label_propagation_time_end - label_propagation_time_start << std::endl;
    }

    // result
    if(mpi_rank == 0) {
      step_result_file << global_itr_count << ", LP, " << "0"  << ", "
	<< (label_propagation_time_end - label_propagation_time_start) << "\n";
    }

  //////////////////////////////////////////////////////////////////////////////

    if (global_init_step) { // Important
      global_init_step = false;
    } 

    // global termination detection

    // global finished status
    // TODO: in lcc, if no vertices and edges are deleted then 
    // global_not_finished = false; change it to true ? 
    global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, 
      std::greater<uint8_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
    if(mpi_rank == 0) {
      std::cout << "Pattern Matching | Global Finished Status : "; 
      if (global_not_finished) { 
	std::cout << "Continue" << std::endl;
      } else {
	std::cout << "Stop" << std::endl;
	do_nonlocal_constraint_checking = false;
      } 
    }

    // global active vertex count
    size_t global_active_vertices_count =  
      havoqgt::mpi_all_reduce(vertex_state_map.size(),
      std::greater<size_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
    if (global_active_vertices_count < 1) {
      //break;
      global_not_finished = false;
      do_nonlocal_constraint_checking = false;
    } else {
      global_not_finished = true;
      do_nonlocal_constraint_checking = true;
    }
    //if(mpi_rank == 0) {
    //  std::cout << "Pattern Matching | Global Active Vertex Count : "
    //  << global_active_vertices_count << std::endl;
    //}

    global_itr_count++; 

    // end of local constraint chceking  

  //////////////////////////////////////////////////////////////////////////////

    // nonlocal constraint checking
     
    //double token_passing_time_start = MPI_Wtime(); 

    if (ptrn_util_two.input_patterns.size() < 1) {
      global_not_finished = false;  
      do_nonlocal_constraint_checking = false; 
    }   

    // Test
    // forced token passing
    //if (global_itr_count == 0) {
    //  global_not_finished = true;  
    //  do_nonlocal_constraint_checking = true; 
    //}
    // Test  

    if (do_nonlocal_constraint_checking && global_not_finished) { // TODO: do we need this? 

      VectorBoolean pattern_found(ptrn_util_two.input_patterns.size(), 0); // per rank state
      VectorBoolean pattern_token_source_found(ptrn_util_two.input_patterns.size(), 0); // per rank state

      VertexUint8Map token_source_map; // per vertex state

      global_not_finished = false;  

      // loop over the set of nonlocal constraints and 
      // run token passing and local constarint checking
      for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {

	// setup subpattern 
	// TODO: This is actually a bad design. 
	// It should be one object (ptrn_util_two) per entry.

        // TODO: imporove 	
        // replace with updated pattern_graph.vertex_data
        for (size_t i = 0; i < std::get<1>(ptrn_util_two.input_patterns[pl]).size(); i++) {
          auto v = std::get<1>(ptrn_util_two.input_patterns[pl])[i]; 
          std::get<0>(ptrn_util_two.input_patterns[pl])[i] = pattern_graph.vertex_data[v];
        } 

	auto pattern_tp = std::get<0>(ptrn_util_two.input_patterns[pl]);
	auto pattern_indices_tp = std::get<1>(ptrn_util_two.input_patterns[pl]);
	auto pattern_cycle_length_tp = std::get<2>(ptrn_util_two.input_patterns[pl]); // uint
	auto pattern_valid_cycle_tp = std::get<3>(ptrn_util_two.input_patterns[pl]); // boolean
	auto pattern_is_tds_tp = std::get<4>(ptrn_util_two.input_patterns[pl]); // boolean
	auto pattern_interleave_label_propagation_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean
	auto pattern_enumeration_tp = ptrn_util_two.enumeration_patterns[pl]; 
	auto pattern_aggregation_steps_tp = ptrn_util_two.aggregation_steps[pl]; 

	// TODO: read from file / remove
	auto pattern_selected_vertices_tp = 0; // TODO: remove 
	auto pattern_selected_edges_tp = false; // boolean
	auto pattern_mark_join_vertex_tp = false; // boolean
	auto pattern_ignore_join_vertex_tp = false; // boolean  
	auto pattern_join_vertex_tp = 0; // TODO: 
	//bool do_tds_tp = false;
    
	// Test print 
	if(mpi_rank == 0) {
	  std::cout << "Token Passing [" << pl << "] | Searching Subpattern : ";
	  PatternNonlocalConstraint::output_pattern(pattern_tp);
	  std::cout << "Token Passing [" << pl << "] | Vertices : ";
	  PatternNonlocalConstraint::output_pattern(pattern_indices_tp);
	  std::cout << "Token Passing [" << pl << "] | Arguments : " 
	    << pattern_cycle_length_tp << " " 
	    << pattern_valid_cycle_tp << " " 
	    << pattern_interleave_label_propagation_tp << std::endl; 
	  std::cout << "Token Passing [" << pl << "] | Enumeration Indices : ";
	  PatternNonlocalConstraint::output_pattern(pattern_enumeration_tp);
	  //std::cout << "Token Passing [" << pl << "] | Agreegation Steps : TODO" 
	  //  << std::endl;
	}
	// Test print

	// initialize containers  
	token_source_map.clear(); // Important
	vertex_token_source_set.clear(); // Important
     
	// initialize application parameters 
	bool token_source_deleted = false;
	message_count = 0;
     
	// result
	// TODO: only output subgraphs when doing enumeration  
	std::string paths_result_filename = result_dir + 
	  "/all_ranks_subgraphs/subgraphs_" +
	  std::to_string(pl) + "_" + std::to_string(mpi_rank);
	  std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);
      
	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
  
  ////////////////////////////////////////////////////////////////////////////// 
  
	// run token passing

	time_start = MPI_Wtime();

	// token passing - distributed processing

	if (pattern_is_tds_tp) {
	  prunejuice::token_passing_pattern_matching<graph_type, Vertex, Edge, VertexData, 
	    EdgeData, VertexMetadata, EdgeMetadata, VertexActive, 
	    VertexUint8MapCollection, TemplateVertex, VertexStateMap, PatternGraph, 
	    PatternNonlocalConstraint, VertexUint8Map, VertexSetCollection, 
	    DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>
	    (graph, vertex_metadata, vertex_active, vertex_active_edges_map, 
	    template_vertices, vertex_state_map, pattern_graph, ptrn_util_two, pl,
	    token_source_map, vertex_token_source_set, 
	    pattern_found[pl], tp_batch_size, paths_result_file, message_count);

	} else {     
	  prunejuice::token_passing_pattern_matching<graph_type, VertexMetadata, 
	    decltype(pattern_tp), decltype(pattern_indices_tp), uint8_t, PatternGraph,
	    VertexStateMap, VertexUint8Map, edge_data_t,
	    VertexSetCollection, VertexActive, TemplateVertex, VertexUint8MapCollection, BitSet>
	    (graph, vertex_metadata, pattern_tp,
	    pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map,
	    token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
	    pattern_found[pl], *edge_data_ptr, vertex_token_source_set, vertex_active, 
	    template_vertices, vertex_active_edges_map, pattern_selected_vertices_tp,
	    pattern_selected_edges_tp, pattern_mark_join_vertex_tp,
	    pattern_ignore_join_vertex_tp, pattern_join_vertex_tp, message_count);
	} 
       
	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?    
	time_end = MPI_Wtime();
	if(mpi_rank == 0) {
	  std::cout << "Pattern Matching Time | Token Passing (Traversal) [" 
	    << pl << "] : " << time_end - time_start << std::endl;
	}

	// result
	paths_result_file.close();

	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

  //////////////////////////////////////////////////////////////////////////////

	// token passing - node local processing

	// remove invalid (token source) vertices from the vertex_state_map
	// for delegates, set vertex_active to false

	// TODO: in the case, a vertex is on multiple cycles/chains (not as the token source)
	// only invalidate it as a token source, but do not remove it from the vertex_state_map

	size_t token_source_pattern_indices_tp_index = 0; // TODO: this is confusing, update 
	
	bool is_token_source_map_not_empty = havoqgt::mpi_all_reduce
	  (!token_source_map.empty(), std::greater<uint8_t>(), MPI_COMM_WORLD); // TODO: less did not work ?
	MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
	pattern_token_source_found[pl] = is_token_source_map_not_empty;

	uint64_t remove_count = 0;      

	for (auto& s : token_source_map) {
	  if (!s.second) {
	    auto v_locator = graph->label_to_locator(s.first);
	    BitSet v_template_vertices(template_vertices[v_locator]);

	    if (v_template_vertices.none()) {
	      continue;
	    } 

	    //pattern_indices_tp[0]; // token source template vertex ID

	    if (v_template_vertices.test(pattern_indices_tp[token_source_pattern_indices_tp_index])) {
	      assert(pattern_indices_tp[token_source_pattern_indices_tp_index] < max_template_vertex_count); // Test
	      v_template_vertices.reset(pattern_indices_tp[token_source_pattern_indices_tp_index]);
	      template_vertices[v_locator] = v_template_vertices.to_ulong();
	    }

	    if (v_template_vertices.none()) {
	      vertex_active[graph->label_to_locator(s.first)] = false;
	    }

	    if (!global_not_finished) {
	      global_not_finished = true;
	    }
       
	    if (!token_source_deleted) {  
	      token_source_deleted = true;
	    }

	    //if (global_itr_count > 0) { // Test
	    //  std::cout << "MPI Rank: " << mpi_rank << " template vertices " 
	    //    << v_template_vertices << " global_not_finished " 
	    //    << global_not_finished << std::endl;
	    //} 
	     
	  } // if
	} // for

	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

	vertex_active.all_min_reduce(); 
	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

	// TODO: this is a temporary patch, forcing all the delegates 
	// to have no identity
	for(vitr_type vitr = graph->delegate_vertices_begin();
	  vitr != graph->delegate_vertices_end(); ++vitr) {
	  auto vertex = *vitr;
	  if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
	    continue;  // skip the controller
	  }
	  else {
	    auto find_vertex = vertex_state_map.find(graph->locator_to_label(vertex));
	    if (find_vertex == vertex_state_map.end()) {
	      template_vertices[vertex] = 0;
	    }
	  }
	} // for
	MPI_Barrier(MPI_COMM_WORLD);

	template_vertices.all_max_reduce(); 
	// ensure all the delegates have the same value as the controller
	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

	// remove deleted token sources from vertex_state_map
	for (auto& s : token_source_map) {
	  auto v_locator = graph->label_to_locator(s.first);
	  if (!vertex_active[v_locator]) {
	    auto find_vertex = 
	      vertex_state_map.find(s.first);
       
	    if (find_vertex != vertex_state_map.end()) { 
	    
	      if (vertex_state_map.erase(s.first) < 1) { // s.first is the vertex
		std::cerr << "Error: failed to remove an element from the map."
		<< std::endl;
	      } else {
		remove_count++;
		//if (!global_not_finished) {
		//  global_not_finished = true;
		//}
	      }   
	    } // if
	  } // if 
	} // for

	//  if(mpi_rank == 0) {
	//    std::cout << "Token Passing [" << pl << "] | MPI Rank [" << mpi_rank 
	//      << "] | Removed " << remove_count << " vertices."<< std::endl; 
	//      TODO: not useful informtion
	//  }
 
  //////////////////////////////////////////////////////////////////////////////

        //MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
	time_end = MPI_Wtime();
	if(mpi_rank == 0) {
	  std::cout << "Pattern Matching Time | Token Passing [" << pl 
	    << "] : " << time_end - time_start << std::endl;
	}

#ifdef OUTPUT_RESULT
	// result
	if(mpi_rank == 0) {
	  superstep_result_file << global_itr_count << ", TP, "
	    << pl << ", "
	    << time_end - time_start << "\n";
	}

        if(mpi_rank == 0) {
          step_result_file << global_itr_count << ", TP, "
            << pl << ", "
            << time_end - time_start << "\n";
        } 

	// Important : this may slow down things - only for presenting results
	active_vertices_count = 0;
	active_edges_count = 0;   

	for (auto& v : vertex_state_map) {
	  auto v_locator = graph->label_to_locator(v.first);
	  if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
	    active_vertices_count++;

	    // edges
	    active_edges_count+=vertex_active_edges_map[v_locator].size();   
	  } else if (!v_locator.is_delegate()) {
	    active_vertices_count++;

	    // edges
	    active_edges_count+=vertex_active_edges_map[v_locator].size();   
	  }
	}

	// vertices
	active_vertices_count_result_file << global_itr_count << ", TP, "
	  << pl << ", "  
	  << active_vertices_count << "\n";  

	// edges   
	active_edges_count_result_file << global_itr_count << ", TP, "
	  << pl << ", "
	  << active_edges_count << "\n";

	// messages
	message_count_result_file << global_itr_count << ", TP, "
	  << pl << ", "
	  << message_count << "\n";
#endif
  
	MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

	// token source found or not; token source deleted or not
	if (is_token_source_map_not_empty) {

	  // subpattren (nonlocal constraint) was found or not 
	  // TODO: write output to file 
	  havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), 
	    MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
	  if(mpi_rank == 0) {   
	    std::string s = pattern_found[pl] == 1 ? "True" : "False";
	    std::cout << "Token Passing [" << pl << "] | Found Subpattern : " 
	      << s << std::endl;
	  }

	  //MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

	  // verify global token source deleted status
	  token_source_deleted = havoqgt::mpi_all_reduce(token_source_deleted,
	    std::greater<uint8_t>(), MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 
	  if(mpi_rank == 0) {
	    std::cout << "Token Passing [" << pl 
	      << "] | Token Source Deleted Status : ";
	    if (token_source_deleted) {
	      std::cout << "Deleted" << std::endl;
	    } else {
	      std::cout << "Not Deleted" << std::endl;
	    }
	  } 

	} else {
	  if(mpi_rank == 0) {
	    std::cout << "Token Passing [" << pl << "] | No Token Source Found"  
	      << std::endl; 
	  }
	} // if - is_token_source_map_not_empty
   
        global_itr_count++; 

        // end of token passing 

  //////////////////////////////////////////////////////////////////////////////

	// interleave nonlocal constraint checking with local constraint checking 
  
	if (token_source_deleted && pattern_interleave_label_propagation_tp) {

	  label_propagation_time_start = MPI_Wtime();

	  prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData, 
	    graph_type, VertexMetadata, VertexStateMap, VertexActive, 
	    VertexUint8MapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>
	    (graph, vertex_metadata, vertex_state_map, vertex_active, 
	    vertex_active_edges_map, template_vertices, pattern_graph, global_init_step, 
	    global_not_finished, global_itr_count, superstep_result_file, 
	    active_vertices_count_result_file, active_edges_count_result_file,
	    message_count_result_file);  

	  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
	  label_propagation_time_end = MPI_Wtime();
	  if(mpi_rank == 0) {
	    std::cout << "Pattern Matching Time | Local Constraint Checking (Interleaved) : " 
	      << label_propagation_time_end - label_propagation_time_start << std::endl;
	  }

	  // result
	  if(mpi_rank == 0) {
	    step_result_file << global_itr_count << ", LP, " << pl << ", "
	      << (label_propagation_time_end - label_propagation_time_start) << "\n";
	  }

	  // global termination detection
	   
	  // global active vertices count 
	  size_t global_active_vertices_count =
	    havoqgt::mpi_all_reduce(vertex_state_map.size(),
	    std::greater<size_t>(), MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
	  if (global_active_vertices_count < 1) {
	    global_not_finished = false;
	    do_nonlocal_constraint_checking = false;
	    break; // exits the current loop and print results
	  } else {
	    global_not_finished = true;
	    do_nonlocal_constraint_checking = true;
	  }

          global_itr_count++;  
	
	} else { 
	  if(mpi_rank == 0) {
	    std::cout << "Pattern Matching | Skipping Local Constraint Checking" 
	      << " (Interleaved)" << std::endl;
	  }        
	}
        
	// end of local constraint checking

  //////////////////////////////////////////////////////////////////////////////

	// end of the current nonlocal constraint

      } // for - loop over nonlocal constraints 

  ////////////////////////////////////////////////////////////////////////////// 

      // end of nonlocal constraint checking  

      // subpattrens (nonlocal constraints) were found or not // TODO: write output to file 
      havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
      if(mpi_rank == 0) {   
	for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
	  if (pattern_token_source_found[pl]) {
	    std::string s = pattern_found[pl] == 1 ? "True" : "False";
	    std::cout << "Token Passing [" << pl << "] | Found Subpattern : " << s << std::endl;
	  } else {
	    std::cout << "Token Passing [" << pl << "] | No Token Source Found" << std::endl;
	  }
	}
      }

      // result
      // Important : this may slow down things - only for presenting results
      //--  active_vertices_count = 0;
      //--  for (auto& v : vertex_state_map) {
      //--    auto v_locator = graph->label_to_locator(v.first);
      //--    if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
      //--      active_vertices_count++;
      //--    } else if (!v_locator.is_delegate()) {
      //--      active_vertices_count++;
      //--    }
      //--  }

      //--  active_vertices_count_result_file << global_itr_count << ", TP, "
      //--    << "0, "  
      //--    << active_vertices_count << "\n";    
   
    } else {
      if(mpi_rank == 0) { 
	std::cout << "Pattern Matching | Skipping Token Passing" << std::endl;
      }
      global_not_finished = false;

    } // if - do_nonlocal_constraint_checking 

    // done token passing // TODO: not useful anymore 
    //double token_passing_time_end = MPI_Wtime();
    //if(mpi_rank == 0) {
    //  std::cout << "Fuzzy Pattern Matching Time | Token Passing : " 
    //    << (token_passing_time_end - token_passing_time_start) << std::endl;
    //}

    // result // TODO: not useful anymore
    //if(mpi_rank == 0) {
    //  step_result_file << global_itr_count << ", TP, "  
    //    << (token_passing_time_end - token_passing_time_start) << "\n"; 
    //}

    global_not_finished = false;

    // end of constraint checking for the current pattern

  ///////////////////////////////////////////////////////////////////////////// 

    MPI_Barrier(MPI_COMM_WORLD);  
    double pattern_time_end = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "Pattern Matching Time | Pattern [" << ps << "] : " 
      << pattern_time_end - pattern_time_start << std::endl;
    }

    // global termination detection // TODO: remove ? 
    global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, 
      std::greater<uint8_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

    if(mpi_rank == 0) {
      std::cout << "Pattern Matching | Global Finished Status : ";
      if (global_not_finished) {
	std::cout << "Continue" << std::endl;
      } else {
	std::cout << "Stop" << std::endl;
      }
    }

    //double itr_time_end = MPI_Wtime(); 
    //if(mpi_rank == 0) { //TODO: sum of LCC and NLCC iterations
    //  std::cout << "Pattern Matching Time | Pattern [" << ps 
    //    << "] | Iteration [" << global_itr_count << "] : " 
    //    << itr_time_end - itr_time_start << std::endl;  
    //}

    // result 
    //if(mpi_rank == 0) {
      // iteration number, time
      //itr_result_file << global_itr_count << ", "
	//<< (itr_time_end - itr_time_start) << "\n";
    //}

    //global_itr_count++; //TODO: sum of LCC and NLCC iterations

    //if(mpi_rank == 0) { //TODO: sum of LCC and NLCC iterations
    //std::cout << "Pattern Matching | Pattern [" << ps 
    //  << "] | # Iterations : " << global_itr_count << std::endl;
    //}

    // TODO: mention whether the pattern was found or not  

  //////////////////////////////////////////////////////////////////////////////

    // result

    // Important : this may slow things down - only for presenting results

    size_t active_edge_count_local = 0;

    for (auto& v : vertex_state_map) {
      auto v_locator = graph->label_to_locator(v.first);
      if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
	BitSet v_template_vertices(template_vertices[v_locator]); 
	active_vertices_result_file << mpi_rank << ", "
	  << v.first << ", " 
	  //<< v.second.vertex_pattern_index << ", "
	  << "0, " 
	  << vertex_metadata[v_locator] << ", "
	  << v_template_vertices << "\n";
       
	// edges 
	for (auto& n : vertex_active_edges_map[v_locator]) { 	 
	  active_edges_result_file << mpi_rank << ", "
	    << v.first << ", "
	    << n.first //<< ", "
	    //<< (size_t)n.second << ", "
	    //<< vertex_active_edges_map[v_locator].size() 
	    << "\n";
	}
   
	active_edge_count_local += vertex_active_edges_map[v_locator].size(); 
    
      } else if (!v_locator.is_delegate()) {
	BitSet v_template_vertices(template_vertices[v_locator]);
	active_vertices_result_file << mpi_rank << ", "
	  << v.first << ", "
	  //<< v.second.vertex_pattern_index << ", "
	  << "0, "
	  << vertex_metadata[v_locator] << ", "
	  << v_template_vertices << "\n";

	// edges
	for (auto& n : vertex_active_edges_map[v_locator]) {  
	  active_edges_result_file << mpi_rank << ", "
	    << v.first << ", "
	    << n.first //<< ", "          
	    //<< (size_t)n.second << ", " 
	    //<< vertex_active_edges_map[v_locator].size()
	    << "\n"; 
	}

	active_edge_count_local += vertex_active_edges_map[v_locator].size(); 

      }

      // motif 
      motif_db->vertex_motif_map[v.first].insert(pattern_graph.pattern_ID);

    } // for

    MPI_Barrier(MPI_COMM_WORLD);

    size_t vertex_state_map_set_size_global =
    havoqgt::mpi_all_reduce(vertex_state_map.size(), std::plus<size_t>(), 
      MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    size_t active_edge_count_global =
    havoqgt::mpi_all_reduce(active_edge_count_local, std::plus<size_t>(), 
      MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Pattern ["
        << pattern_graph.pattern_ID << "] " << std::endl;
      std::cout << "Pattern Matching | Pattern Vertex Match Count | Pattern [" 
	<< ps << "] : "
	<< vertex_state_map_set_size_global << std::endl;
      std::cout << "Pattern Matching | Pattern Edge Match Count | Pattern [" 
	<< ps << "] : "
	<< active_edge_count_global << std::endl;
    }  

    // result
    if(mpi_rank == 0) {  
      // pattern set element ID, number of MPI ranks, 
      // total number of iterations (lcc + nlcc), total time, 
      // #vertices in the pattern, #edges in the pattern, 
      // #nonlocal constraints,
      // #vertices in the solution subgraph, #edges in the solution subgraph
      pattern_set_result_file << ps << ", " 
	<< mpi_size << ", "
	<< global_itr_count << ", " 
	<< (pattern_time_end - pattern_time_start) << ", "
	<< pattern_graph.vertex_count << ", "
	<< pattern_graph.edge_count << ", " 
	<< ptrn_util_two.input_patterns.size() << ", "
        << vertex_state_map_set_size_global << ", "
        << active_edge_count_global << ", "
        << "\n"; 
    }

    // motif  
    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Motif Database (Local) Item Count : "
        << motif_db->vertex_motif_map.size() << std::endl;
    }  

  //////////////////////////////////////////////////////////////////////////////

    // cleanup memeory
     
    //vertex_active.clear(); // TODO: add clear() method to vertex_data.cpp
    //vertex_state_map.clear();

    // close files
     
    //itr_result_file.close();
    step_result_file.close();
    superstep_result_file.close();
    active_vertices_count_result_file.close(); 
    active_vertices_result_file.close();
    active_edges_count_result_file.close();
    active_edges_result_file.close();
    message_count_result_file.close();

    // end of the current pattern

    MPI_Barrier(MPI_COMM_WORLD);

  } // for - loop over the pattern set

  } // motif

  //////////////////////////////////////////////////////////////////////////////

  // motif

  {
 
  if (mpi_rank == 0) {
    std::cout << "Pattern Matching | Writing Output to File ... " << std::endl; 
  }


  // TODO: same as all_ranks_, directories are created ahead of time
  std::string motif_db_dir = "/p/lustre2/reza2/graph_learning/reddit/metall_motif_collection/motif_A_"; // TODO
  //std::string motif_result_dir = "/p/lustre2/reza2/graph_learning/reddit/lbann_motif_collection/";

  using MotifContainer = prunejuice::container::motif_container<uint64_t, uint64_t,
    uint64_t, metall::manager::allocator_type<char>>;
 
  MPI_Comm shm_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);     

  int shm_rank, shm_size;
  MPI_Comm_rank(shm_comm, &shm_rank);
  MPI_Comm_size(shm_comm, &shm_size);
  //std::cout << mpi_rank << " " << shm_rank << " " << shm_size << std::endl;
  
  std::vector<int> recv_buffer;
  havoqgt::mpi_all_gather(mpi_rank, recv_buffer, shm_comm); // works
  MPI_Barrier(shm_comm); // TODO: may not need this here 
  //std::cout << recv_buffer.size() << std::endl;

  if (shm_rank == 0) {

    int mpi_node = 0;
    if (mpi_rank > 0) {
      mpi_node = static_cast<int>(floor(mpi_rank/shm_size));
    }

    std::string motif_result_filename = motif_result_dir + "/motifs_" +
      std::to_string(mpi_node);
    std::ofstream motif_result_file(motif_result_filename, std::ofstream::out);
   
    for (auto& r : recv_buffer) { // node local ranks
 
      std::string motif_db_partition_dir = motif_db_dir + std::to_string(r);
    
      metall::manager partition_manager(metall::open_read_only, 
        motif_db_partition_dir.c_str());
      MotifContainer *motif_db_partition = 
        partition_manager.find<MotifContainer>("motif_container").first;

      for (auto& v : motif_db_partition->vertex_motif_map) {
        if (motif_db_partition->vertex_motif_map[v.first].size() > 0) {
          motif_result_file << v.first << " ";
          for (auto& m : v.second) {
            motif_result_file << m << " ";
          }

          //size_t pad_count = global_max_motif_count - v.second.size();
          //for (size_t i = 0; i < pad_count; i++) {
          //  motif_result_file << "-1 ";
          //}

          motif_result_file << "\n";
        }
      } // for       
       
    } // for

    motif_result_file.close(); 

  } // if

  MPI_Comm_free(&shm_comm); 

  } // motif

  //////////////////////////////////////////////////////////////////////////////

  // end of all patterns 

  if (mpi_rank == 0) {
      pattern_set_result_file.close(); // close file
      std::cout << "Pattern Matching | Done" << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  } // pattern matching  

  } // havoqgt_init
  ;
  // END Main MPI
  ///havoqgt::havoqgt_finalize();

  return 0;  
} // main
