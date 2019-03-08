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
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
///#include <havoqgt/environment.hpp>

#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>

#include <prunejuice/template.hpp>
#include <prunejuice/non_local_constraint.hpp>
#include <prunejuice/algorithm_state.hpp>
#include <prunejuice/local_constraint_checking.hpp>
#include <prunejuice/non_local_constraint_checking_unique.hpp>
#include <prunejuice/non_local_constraint_checking_tds_batch.hpp> // TODO: optimize for batch_size = mpi_size

/*
//#include <havoqgt/graph.hpp>
#include <havoqgt/pattern_graph.hpp>
#include <havoqgt/pattern_util.hpp>
//#include <havoqgt/label_propagation_pattern_matching.hpp>
//#include <havoqgt/label_propagation_pattern_matching_bsp.hpp> 
//#include <havoqgt/label_propagation_pattern_matching_iterative.hpp>
#include <havoqgt/label_propagation_pattern_matching_nonunique_ee.hpp>
//#include <havoqgt/label_propagation_pattern_matching_nonunique_counting_ee.hpp>
//#include <havoqgt/token_passing_pattern_matching.hpp>
//--#include <havoqgt/token_passing_pattern_matching_new.hpp> // TP_ASYNC
//#include <havoqgt/token_passing_pattern_matching_batch.hpp> // TP_BATCH
//#include <havoqgt/token_passing_pattern_matching_iterative.hpp>
//#include <havoqgt/token_passing_pattern_matching_edge_aware.hpp> // TP_ASYNC
//#include <havoqgt/token_passing_pattern_matching_path_checking.hpp> // TP_ASYNC
//#include <havoqgt/token_passing_pattern_matching_nonunique_ee.hpp> // TP_ASYNC // nonunique path checking
//#include <havoqgt/token_passing_pattern_matching_nonunique_nem.hpp> // TP_ASYNC // nonunique path checking, edge monocyclic 
#include <havoqgt/token_passing_pattern_matching_nonunique_nem_1.hpp> // TP_ASYNC // nonunique path checking, edge monocyclic 
//#include <havoqgt/token_passing_pattern_matching_nonunique_tds_1.hpp> // TP_ASYNC // nonunique path checking, edge monocyclic, tds
///#include <havoqgt/token_passing_pattern_matching_nonunique_tds_batch_1.hpp>
#include <havoqgt/token_passing_pattern_matching_nonunique_tds_batch_4.hpp>
//#include <havoqgt/token_passing_pattern_matching_nonunique_iterative_tds_1.hpp> // TP_ASYNC // nonunique path checking, edge monocyclic, tds
//#include <havoqgt/token_passing_pattern_matching_nonunique_iterative_tds_batch_1.hpp> // TP_ASYNC // nonunique path checking, edge monocyclic, tds, batching
//#include <havoqgt/update_edge_state.hpp>
///#include <havoqgt/vertex_data_db.hpp>
#include <havoqgt/vertex_data_db_degree.hpp>*/

//#define OUTPUT_RESULT
//#define ENABLE_BLOCK

#define TP_ASYNC
//#define TP_BATCH

///namespace hmpi = havoqgt::mpi;
///using namespace havoqgt::mpi;
using namespace havoqgt;
//using namespace prunejuice;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;

template<typename T>
  using SegmentAllocator = bip::allocator<T, segment_manager_t>;

///typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
//typedef havoqgt::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef havoqgt::delegate_partitioned_graph
  <typename segment_manager_t::template allocator<void>::type> graph_type;

template<typename T>
  using DelegateGraphVertexDataSTDAllocator = graph_type::vertex_data
  <T, std::allocator<T>>; 

template<typename T>
  using DelegateGraphEdgeDataSTDAllocator = graph_type::edge_data
  <T, std::allocator<T>>;  

void usage()  {
  //if(havoqgt_env()->world_comm().rank() == 0) {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -p <string> -o <string>\n"
      << " -i <string>   - input graph base filename (required)\n"
      << " -b <string>   - backup graph base filename. If set, \"input\" graph will be deleted if it exists\n"
      << " -v <string>   - vertex metadata base filename (optional, Default is degree based metadata)\n"
      << " -e <string>   - edge metadata base filename (optional)\n" 
      << " -p <string>   - pattern base directory (required)\n"
      << " -o <string>   - output base directory (required)\n"
      //<< " -x <int>      - Token Passing batch size (optional, Default/max batch size is " 
      //  << havoqgt_env()->world_comm().size() << ", Min batch size is 1)\n"
      //  << comm_world().size() << " , min batch size is 1)\n" 
      << " -x <int>      - Token Passing batch count (optional, Default/min batch count is 1, max batch count is "
        << comm_world().size() << "\n"   
      << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_input, 
  std::string& backup_graph_input, std::string& vertex_metadata_input, 
  std::string& edge_metadata_input, std::string& pattern_input, 
  std::string& result_output, uint64_t& tp_vertex_batch_size) {

  //if(havoqgt_env()->world_comm().rank() == 0) {
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
        tp_vertex_batch_size = std::stoull(optarg);
        if (tp_vertex_batch_size < 1 || tp_vertex_batch_size > comm_world().size()) {
          print_help = true;
        } else if (tp_vertex_batch_size > 1) {
          tp_vertex_batch_size = comm_world().size() / tp_vertex_batch_size;
        } else {
          tp_vertex_batch_size = comm_world().size();
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
  //typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
  //typedef hmpi::delegate_partitioned_graph
  //  <typename segment_manager_t::template allocator<void>::type> graph_type;
  
  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
  //havoqgt::havoqgt_init(&argc, &argv);
  havoqgt::init(&argc, &argv);
  {

  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  ///havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI Initialized With " << mpi_size << " Ranks." << std::endl;
    ///havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);  

  /////////////////////////////////////////////////////////////////////////////

  // parse command line
/*  std::string graph_input = argv[1];
  //std::string backup_filename; 

  // for pattern matching
  std::string vertex_data_input_filename = argv[2];
  //std::string pattern_input_filename = argv[3];
  std::string pattern_dir = argv[3];
  std::string vertex_rank_output_filename = argv[4];
  std::string backup_filename = argv[5];
  std::string result_dir = argv[6];
  bool use_degree_as_vertex_data = std::stoull(argv[7]); // 1 - yes, 0 - no 
  bool do_load_graph_from_backup_file = std::stoull(argv[8]); // 1 - yes, 0 - no // TODO: remove
*/
 
  std::string graph_input;
  std::string backup_graph_input;
  std::string vertex_metadata_input;
  std::string edge_metadata_input;
  std::string pattern_input;
  std::string result_output;  
  ///uint64_t tp_vertex_batch_size = havoqgt_env()->world_comm().size();   
  uint64_t tp_vertex_batch_size = comm_world().size();

  parse_cmd_line(argc, argv, graph_input, backup_graph_input, 
    vertex_metadata_input, edge_metadata_input, pattern_input, result_output, 
    tp_vertex_batch_size); 

  std::string pattern_dir = pattern_input; 
  std::string result_dir = result_output;

  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////

  // load graph

  if (mpi_rank == 0) {
    std::cout << "Loading Graph ... " << std::endl;
  }

  //if(do_load_graph_from_backup_file) { // TODO: remove do_load_graph_from_backup_file
  //  distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
  //}
  if (backup_graph_input.size() > 0) {
    distributed_db::transfer(backup_graph_input.c_str(), graph_input.c_str());
  }

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

  //segment_manager_t* segment_manager = ddb.get_segment_manager();
  //  bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  //graph_type *graph = segment_manager->
  //  find<graph_type>("graph_obj").first;
  //assert(graph != nullptr);  
  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  // edge data
  //if (mpi_rank == 0) {
    //std::cout << "Loading / Initializing Edge Data ... " << std::endl;
  //}

  typedef uint8_t edge_data_type; 
  // TODO: figure out a way to get it from graph_type
  // see edge_data_value_type in parallel_edge_list_reader.hpp

  typedef graph_type::edge_data<edge_data_type, 
    bip::allocator<edge_data_type, segment_manager_t>> edge_data_t;

  edge_data_t* edge_data_ptr = ddb.get_segment_manager()->
    find<edge_data_t>("graph_edge_data_obj").first;
//  assert(edge_data_ptr != nullptr); 

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Done Loading Graph." << std::endl;
  }

  //graph->print_graph_statistics(); // causes MPI error
  //MPI_Barrier(MPI_COMM_WORLD);  

  /////////////////////////////////////////////////////////////////////////////

  // pattern matching
  {

  // types used by the delegate partitioned graph
  typedef typename graph_type::vertex_iterator vitr_type;
  typedef typename graph_type::vertex_locator vloc_type;
  //typedef typename graph_type::edge_iterator eitr_type;
 
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData; // for string hash
  //typedef uint8_t VertexData; // for log binning 
  //typedef edge_data_type EdgeData;
  typedef uint8_t EdgeData; 
  
  typedef uint64_t VertexRankType; // TODO: delete

  static constexpr size_t max_bit_vector_size = 16; // TODO:
  static constexpr size_t max_template_vertex_count = 16;  
  typedef std::bitset<max_bit_vector_size> BitSet; // TODO: rename to TemplateVertexSet
  typedef BitSet TemplateVertexBitSet; 
  typedef uint16_t TemplateVertexType; // TODO: rename to TemplateVertexBitSetToUint
  
  typedef uint8_t Boolean; // TODO: replace all bool with Boolean?

  // TODO: mmap
  //typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetadata;
  //typedef graph_type::vertex_data<VertexRankType, SegmentAllocator<VertexRankType> > VertexRank;
  //typedef graph_type::vertex_data<bool, SegmentAllocator<bool> > VertexActive;
  //typedef graph_type::vertex_data<uint64_t, SegmentAllocator<uint64_t> > VertexIteration;

  typedef graph_type::vertex_data<VertexData, std::allocator<VertexData> > VertexMetadata;
  typedef graph_type::vertex_data<Boolean, std::allocator<Boolean> > VertexActive; // TODO: solution_graph // TODO: you are mixing bool and uint!
  typedef graph_type::vertex_data<TemplateVertexType, std::allocator<TemplateVertexType> > TemplateVertex; // TODO: solution_graph, rename to VertexTemplateVertexBitSetToUint

  typedef graph_type::vertex_data<uint64_t, std::allocator<uint64_t> > VertexIteration; // TODO: delete
  typedef graph_type::vertex_data<VertexRankType, std::allocator<VertexRankType> > VertexRank; // TODO: delete

  //typedef vertex_state<uint8_t> VertexState;
///  typedef prunejuice::vertex_state_generic<Vertex, VertexData, uint8_t, BitSet> VertexState;
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

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching ... " << std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////

  double time_start = MPI_Wtime();
  double time_end = MPI_Wtime();

  // per rank containers 
  VertexStateMap vertex_state_map; 

  // vertex containers 
   
  // TODO: need a new alloc_inst to use bip/mmap
//  VertexMetadata vertex_metadata(*graph, alloc_inst);
//  VertexRank vertex_rank(*graph, alloc_inst);
//  VertexActive vertex_active(*graph, alloc_inst);
//  VertexIteration vertex_iteration(*graph, alloc_inst);
 
  VertexMetadata vertex_metadata(*graph); 
  VertexActive vertex_active(*graph);
  TemplateVertex template_vertices(*graph);
  VertexUint8MapCollection vertex_active_edges_map(*graph);
  VertexSetCollection vertex_token_source_set(*graph); // per vertex set
//--  VertexSetCollection token_source_edge_set(*graph); // per vertex set // edge aware

//  VertexRank vertex_rank(*graph);
  uint8_t vertex_rank; // TODO: dummy
//  VertexIteration vertex_iteration(*graph);
  uint8_t vertex_iteration; // TODO: dummy  

  // edge containers
  //EdgeMetadata edge_metadata(*graph);
  ////EdgeActive edge_active(*graph);

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching | Allocated Vertex and Edge Containers" 
    << std::endl;
  }

  ///////////////////////////////////////////////////////////////////////////// 

  // application parameters // TODO: command line input

  // write vertex data to file
  bool do_output_vertex_data = false; // TODO: ?

  size_t token_passing_algo = 0; // TODO: ? 

  MPI_Barrier(MPI_COMM_WORLD);
 
  ///////////////////////////////////////////////////////////////////////////// 

  // build the distributed vertex data db
  time_start = MPI_Wtime();

  //if (use_degree_as_vertex_data) {
  if (vertex_metadata_input.size() > 0) {
    vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>
      //(graph, vertex_metadata, vertex_data_input_filename, 10000);      
      (graph, vertex_metadata, vertex_metadata_input, 10000);
      // TODO: each rank reads 10K lines from file at a time
  } else {
    vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>
      (graph, vertex_metadata);
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Vertex Data DB : " 
      << time_end - time_start << std::endl;
  }

  if (do_output_vertex_data) {
    std::string vertex_data_filename = result_dir + //"/" + std::to_string(0) + 
      "/all_ranks_vertex_data/vertex_data_" + std::to_string(mpi_rank);
    std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

    for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
      ++vitr) {
      vloc_type vertex = *vitr;
      vertex_data_file << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
        << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
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
  
  //MPI_Barrier(MPI_COMM_WORLD); // Test  
  //return 0; // Test

  ///////////////////////////////////////////////////////////////////////////// 
 
  // result
  std::string pattern_set_result_filename = result_dir + "/result_pattern_set";
  std::ofstream pattern_set_result_file;  
  if (mpi_rank == 0) {
    pattern_set_result_file = std::ofstream(pattern_set_result_filename, std::ofstream::out);
  }
  /////////////////////////////////////////////////////////////////////////////
  // TODO: setup pattern set 
  // a pattern set is a collection of directories containing pattern files 
  
  // TODO: code indentation
   
  // loop over pattern set
  for (size_t ps = 0; ps < 1; ps++) { // TODO: for now, only reading from pattern_dir/0 
  // beginning of the pattern set
   
  // setup pattern to search
  if(mpi_rank == 0) { 
    std::cout << "Setting up Pattern [" << ps << "] ... " << std::endl;
  }
   
  // setup pattern - for local constraint checking 
  //std::string pattern_input_filename = pattern_dir + "/" + std::to_string(ps) + "/pattern";
  std::string pattern_input_filename = pattern_dir + "/pattern";
  
  typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData, 
    EdgeData> PatternGraph;
  PatternGraph pattern_graph(
    pattern_input_filename + "_edge",
    pattern_input_filename + "_vertex",
    pattern_input_filename + "_vertex_data",
    pattern_input_filename + "_edge_data",
    pattern_input_filename + "_stat",
    false, false); // TODO: improve

  // TODO: can you create a graphical representation of the pattern 
  // using a 2D graphics library?
  // test print
  if(mpi_rank == 0) {
  std::cout << "Pattern Matching | Searching Pattern [" << ps
    << "] : " << std::endl;
  for (auto v = 0; v < pattern_graph.vertex_count; v++) {
    std::cout << v << " : off-set " << pattern_graph.vertices[v] 
    << " vertex_data " << pattern_graph.vertex_data[v] 
    << " vertex_degree " << pattern_graph.vertex_degree[v] << std::endl;
    std::cout << " neighbours : "; 
    for (auto e = pattern_graph.vertices[v]; e < pattern_graph.vertices[v + 1]; e++) {
      auto v_nbr = pattern_graph.edges[e];
      std::cout << v_nbr << ", " ;
    }
    std::cout << std::endl;
    std::cout << " neighbour vertex data count : ";
    for (auto& nd : pattern_graph.vertex_neighbor_data_count_map[v]) {
      std::cout << "(" << nd.first << ", " << nd.second << ") ";
    }
    std::cout << std::endl; 
  }
  std::cout << "diameter : " << pattern_graph.diameter << std::endl; 
  }
  // test print

  // TODO: remove from here    
  // setup pattern - for token passing
   
  //pattern_util<VertexData, Vertex> ptrn_util_two(pattern_input_filename, true); 
  //pattern_util<VertexData, Vertex> ptrn_util_two(pattern_input_filename + "_nem", true);
  
  //typedef pattern_util<VertexData, Vertex> PatternUtilities;
  
  //PatternUtilities ptrn_util_two(pattern_input_filename + "_nlc", 
  //  pattern_input_filename + "_non_local_constraint", true, true);

  //PatternUtilities ptrn_util_two(pattern_input_filename + "_nem", pattern_input_filename + "_tds", true, true);
  //pattern_util<VertexData, Vertex> ptrn_util_two(pattern_input_filename + "_pc", true);

  //PatternUtilities ptrn_util_two(pattern_input_filename + "_nlc",
  // pattern_input_filename + "_non_local_constraint", 
  // pattern_input_filename + "_aggregation", true, true);
  
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: ?
  
  typedef pattern_nonlocal_constraint<Vertex, Edge, VertexData, PatternGraph>
    PatternNonlocalConstraint;

  //PatternNonlocalConstraint ptrn_util_two(pattern_graph,
  //  pattern_input_filename + "_non_local_constraints",
  //  pattern_input_filename + "vertex_non_local_constraints");     

  PatternNonlocalConstraint ptrn_util_two(pattern_graph,
    //pattern_input_filename + "_nonlocal_constraint");
    pattern_dir + "/pattern_nonlocal_constraint"); 

  //auto pattern = std::get<0>(ptrn_util_two.input_patterns[0]); // TODO: remove
  //auto pattern_indices = std::get<1>(ptrn_util_two.input_patterns[0]); // TODO: remove

  // Test
  //if(mpi_rank == 0) {
  //  std::cout << "Token Passing | Agreegation Steps : " << std::endl;
  //  for (auto i : ptrn_util_two.aggregation_steps) {  
  //    pattern_util<uint8_t>::output_pattern(i);
  //  }
  //} 
  // Test

  //MPI_Barrier(MPI_COMM_WORLD); // Test
  //return 0; // Test  

  // initialization - per-pattern

  // initialize containers
  vertex_state_map.clear(); // important
  //vertex_rank.reset(0);
  vertex_active.reset(true); // initially all vertices are active
  //vertex_iteration.reset(0); // TODO: -1 ?
  vertex_active_edges_map.clear(); // important
  vertex_token_source_set.clear(); // clear all the sets on all the vertices
  //edge_metadata.reset(55); // Test
  //edge_active.reset(0); //  initially all edges are active / inactive

  // initialize application parameters  
  bool global_init_step = true; // TODO: Boolean 
  bool global_not_finished = false; // TODO: Boolean

  bool do_nonlocal_constraint_checking = true; // TODO: Boolean

  uint64_t global_itr_count = 0;
  uint64_t active_vertices_count = 0;
  uint64_t active_edges_count = 0;
  uint64_t message_count = 0; 

  // result
  std::string itr_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/result_iteration"; // TODO: improve
  std::ofstream itr_result_file(itr_result_filename, std::ofstream::out);

  std::string step_result_filename = result_dir + //"/" + std::to_string(ps) 
    "/result_step";
  std::ofstream step_result_file(step_result_filename, std::ofstream::out); 

  std::string superstep_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/result_superstep";
  std::ofstream superstep_result_file(superstep_result_filename, std::ofstream::out);

  std::string active_vertices_count_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/all_ranks_active_vertices_count/active_vertices_" + std::to_string(mpi_rank); 
  std::ofstream active_vertices_count_result_file(active_vertices_count_result_filename, std::ofstream::out);

  std::string active_vertices_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/all_ranks_active_vertices/active_vertices_" + std::to_string(mpi_rank);
  std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);

  std::string active_edges_count_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/all_ranks_active_edges_count/active_edges_" + std::to_string(mpi_rank);
  std::ofstream active_edges_count_result_file(active_edges_count_result_filename, std::ofstream::out);

  std::string active_edges_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/all_ranks_active_edges/active_edges_" + std::to_string(mpi_rank);
  std::ofstream active_edges_result_file(active_edges_result_filename, std::ofstream::out);

  //std::string paths_result_filename = result_dir + "/" +
  //  std::to_string(ps) + "/all_ranks_paths/paths_" + std::to_string(mpi_rank);
  //std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);

  std::string message_count_result_filename = result_dir + //"/" + std::to_string(ps) + 
    "/all_ranks_messages/messages_" + std::to_string(mpi_rank); // TODO:message_count
  std::ofstream message_count_result_file(message_count_result_filename, std::ofstream::out);

  MPI_Barrier(MPI_COMM_WORLD);

  double pattern_time_start = MPI_Wtime();

  /////////////////////////////////////////////////////////////////////////////

  // run application
///  do {
  
  if (mpi_rank == 0) {
    std::cout << "Running Constraint Checking ..." << std::endl;
  }

  global_not_finished = false;

  double itr_time_start = MPI_Wtime();

  /////////////////////////////////////////////////////////////////////////////

  // mark inactive edges
  //update_edge_state();

  /////////////////////////////////////////////////////////////////////////////
//#ifdef ENABLE_BLOCK
  // label propagation   
  double label_propagation_time_start = MPI_Wtime();

  // clone pattern matchng
  //fuzzy_pattern_matching(graph, vertex_metadata, pattern, pattern_indices, vertex_rank);

  // label propagation pattern matching 
//  label_propagation_pattern_matching<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    VertexRank, VertexActive, VertexIteration, VertexStateMap, PatternGraph>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph);

  // label propagation pattern matching bsp, iterative 
//  label_propagation_pattern_matching_bsp<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    /*VertexRank*/uint8_t, VertexActive, /*VertexIteration*/uint8_t, VertexStateMap, PatternGraph, EdgeActive, VertexSetCollection,
//    EdgeMetadata>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph, global_init_step, global_not_finished, 
//    global_itr_count, superstep_result_file, active_vertices_count_result_file, edge_active/**edge_data_ptr*/, edge_metadata);

 prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData, 
   graph_type, VertexMetadata, VertexStateMap, VertexActive, 
   VertexUint8MapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>
   (graph, vertex_metadata, vertex_state_map, vertex_active, 
   vertex_active_edges_map, template_vertices, pattern_graph, global_init_step, 
   global_not_finished, global_itr_count, superstep_result_file, 
   active_vertices_count_result_file, active_edges_count_result_file,
   message_count_result_file);

  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
  double label_propagation_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Local Constraint Checking : " 
      << label_propagation_time_end - label_propagation_time_start << std::endl;
  }

  // result
  if(mpi_rank == 0) {
    step_result_file << global_itr_count << ", LP, "
      << (label_propagation_time_end - label_propagation_time_start) << "\n";
  }

//#endif

  /////////////////////////////////////////////////////////////////////////////

  if (global_init_step) { // Important
    global_init_step = false;
  } 

  // global termination detection 
  //std::cout << "Fuzzy Pattern Matching | Global Not Finished status (local) : " << global_not_finished << std::endl; // Test
  // global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work
  ///global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD); 
  global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching | Global Finished Status : "; 
    if (global_not_finished) { 
      std::cout << "Continue" << std::endl;
    } else {
      std::cout << "Stop" << std::endl;
    } 
  }

  // global verification - are all vertex_state_maps empty
  // false - no active vertex left, true - active vertices left 
  bool global_active_vertex = true; //vertex_state_map.size() < 1 ? false : true;
  if (vertex_state_map.size() < 1) {
//    global_active_vertex = false;
  }

  // TODO: What is the issue here? Preventing stdout from this file beyond this point.  	
  // global_active_vertex = havoqgt::mpi::mpi_all_reduce(global_active_vertex, std::greater<uint8_t>(), MPI_COMM_WORLD); 
  // TODO: not working properly - why? // bool does not work
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

//  global_not_finished = global_active_vertex; // TODO: verify and fix

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching | Global Active Vertex Status : ";
    if (global_active_vertex) {
      std::cout << "Active vertices left." << std::endl; // TODO: ignore for now
    } else {
      std::cout << "No active vertex left." << std::endl;
    }
  }

  // Test
/*   uint64_t vertex_active_count = 0;
   uint64_t vertex_inactive_count = 0;
   for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
      ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_active[vertex]) {
        //std::cout << mpi_rank << ", l, " << graph->locator_to_label(vertex) 
        //  << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
        //  vertex_active_count++;
      } else { 
        vertex_inactive_count++;
      }  
    } 	
    	
    for(vitr_type vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_active[vertex]) {  
        if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank)) {
          //std::cout << mpi_rank << ", c, " << graph->locator_to_label(vertex) 
          //  << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
        } else {	
          //std::cout << mpi_rank << ", d, " << graph->locator_to_label(vertex) 
          //  << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n"; 
        }
        vertex_active_count++; 
      } else {
        vertex_inactive_count++;
      } 
    }

    //std::cout << mpi_rank << " vertex_active_count " << vertex_active_count << std::endl; 
    //std::cout << mpi_rank << " vertex_inactive_count " << vertex_inactive_count << std::endl;
    //std::cout << mpi_rank << " vertex_state_map size " << vertex_state_map.size() << std::endl;      
*/
  // Test

  /////////////////////////////////////////////////////////////////////////////
  
  // Test
  // forced token passing
  if (global_itr_count == 0) {
    global_not_finished = true; // TODO: for load balancing experiments, ?  
  }
  // Test  

//#ifdef ENABLE_BLOCK 
  // toekn passing
  double token_passing_time_start = MPI_Wtime(); 

  if (ptrn_util_two.input_patterns.size() < 1) {
    do_nonlocal_constraint_checking = false; 
  }   

  //if ((token_passing_algo == 0) && global_not_finished) { // do token passing ?
  if (do_nonlocal_constraint_checking && global_not_finished) { // TODO: do we need this? 

  global_not_finished = false;  

  //typedef std::unordered_map<Vertex, bool> TokenSourceMap; 
  //TokenSourceMap token_source_map;
  VertexUint8Map token_source_map; // per vertex state

  //std::vector<bool> pattern_found(ptrn_util_two.input_patterns.size(), false); 
  // TODO: bool does not work with mpi_all_reduce_inplace  
  //typedef std::vector<Boolean> VectorBoolean;  
  //std::vector<uint8_t> pattern_found(ptrn_util_two.input_patterns.size(), 0); 
  VectorBoolean pattern_found(ptrn_util_two.input_patterns.size(), 0); // per rank state
  VectorBoolean pattern_token_source_found(ptrn_util_two.input_patterns.size(), 0); // per rank state

  // TODO: code indentation
  // loop over the constraints and run token passing
  for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {

  // TODO: only output subgraphs when doing enumeration  
  // result
  std::string paths_result_filename = result_dir + //"/" +
    //std::to_string(ps) + "/all_ranks_paths/paths_" + 
    //std::to_string(ps) + 
    "/all_ranks_subgraphs/subgraphs_" +
    std::to_string(pl) + "_" + std::to_string(mpi_rank);
    std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);

  bool token_source_deleted = false;

  // TODO: This is actually a bad design. It should be one object per entry. 
  // ptrn_util_twois the object 
 
  // setup pattern 
  auto pattern_tp = std::get<0>(ptrn_util_two.input_patterns[pl]);
  auto pattern_indices_tp = std::get<1>(ptrn_util_two.input_patterns[pl]);
  auto pattern_cycle_length_tp = std::get<2>(ptrn_util_two.input_patterns[pl]); // uint
  auto pattern_valid_cycle_tp = std::get<3>(ptrn_util_two.input_patterns[pl]); // boolean
//  auto pattern_interleave_label_propagation_tp = std::get<4>(ptrn_util_two.input_patterns[pl]); // boolean
//--  auto pattern_seleted_edges_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean 
//  auto pattern_selected_vertices_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean
  
  auto pattern_selected_vertices_tp = 0; // TODO: remove
  
  auto pattern_is_tds_tp = std::get<4>(ptrn_util_two.input_patterns[pl]); // boolean
  auto pattern_interleave_label_propagation_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean

  auto pattern_enumeration_tp = ptrn_util_two.enumeration_patterns[pl]; 
  auto pattern_aggregation_steps_tp = ptrn_util_two.aggregation_steps[pl]; 

  // TODO: read from file / remove
  auto pattern_selected_edges_tp = false; // boolean
  auto pattern_mark_join_vertex_tp = false; // boolean
  auto pattern_ignore_join_vertex_tp = false; // boolean  
  size_t pattern_join_vertex_tp = 0; // TODO: 
  //bool do_tds_tp = false;

  message_count = 0;

  // Test
  // RDT_25
  //if (pl >= 4) {
  // /p/lscratchf/havoqgtu/reza2_tmp/WDC_patterns_12_tree_C_2 / C_4
  //if (pl >= 8) {
  // /p/lscratchf/havoqgtu/reza2_tmp/WDC_patterns_12_tree_C_2 - pruned graph
  //if (pl >= 0) {   
  // /p/lscratchf/havoqgtu/reza2_tmp/WDC_patterns_12_D
  //if (pl >= 4) {
  // /p/lscratchf/havoqgtu/reza2_tmp/WDC_patterns 
  //if (pl >= 8) {  
  // /p/lscratchf/havoqgtu/reza2_tmp/WDC_patterns_16/0
  //if (pl >= 4) {
  //if (pl >= 4) {
  // /p/lscratchf/havoqgtu/reza2_tmp/IMDB/graph/patterns_B_2/
  //if (pl >= 10 && pl <=27) {
  //if (pl >= 6) {
  //if (pl >= 0) {
  // RMAT_tree
  //if (pl >= 4) {
  //if (pl >= 99) { // Test
  //if (pattern_is_tds_tp) {
    //do_tds_tp = true;
    //if(mpi_rank == 0) {
    //  std::cout << "Token Passing [" << pl << "] | Template Driven Search " << std::endl;
    //}
  //}  
  // Test
  
  if(mpi_rank == 0) {
    std::cout << "Token Passing [" << pl << "] | Searching Subpattern : ";
    //pattern_util<VertexData>::output_pattern(pattern_tp);
    PatternNonlocalConstraint::output_pattern(pattern_tp);
    std::cout << "Token Passing [" << pl << "] | Vertices : ";
    //pattern_util<VertexData>::output_pattern(pattern_indices_tp);
    PatternNonlocalConstraint::output_pattern(pattern_indices_tp);
    std::cout << "Token Passing [" << pl << "] | Arguments : " 
      << pattern_cycle_length_tp << " " 
      << pattern_valid_cycle_tp << " " 
      << pattern_interleave_label_propagation_tp << " "
//--      << pattern_seleted_edges_tp << std::endl; // Test 		
      << pattern_selected_vertices_tp << std::endl; // Test

    std::cout << "Token Passing [" << pl << "] | Enumeration Indices : ";
    //pattern_util<Vertex>::output_pattern(pattern_enumeration_tp); 
    PatternNonlocalConstraint::output_pattern(pattern_enumeration_tp);
    std::cout << "Token Passing [" << pl << "] | Agreegation Steps : TODO" << std::endl;
    //PatternNonlocalConstraint::output_pattern(pattern_aggregation_steps_tp); // TODO:     
  }

  //return 0; // Test

  // initialize containers
  
  if (!pattern_selected_vertices_tp) {
    token_source_map.clear(); // Important
    vertex_token_source_set.clear(); // Important
  } else {
    // TODO: delete

    // delete invalid vertices from the token_source_map
    // reset the valid vertices in the token_source_map 

    /*for (auto itr = token_source_map.begin(); itr != token_source_map.end(); ) {
     if (!itr->second) {
       itr = token_source_map.erase(itr); // C++11
     } else {
       itr->second = 0;
       ++itr;
     }
    }*/
    token_source_map.clear(); // Important       
 
    // Test
    //if (token_source_map.size() > 0) {   
    //  std::cout << mpi_rank << " " << token_source_map.size() << std::endl;
    //}

    //for (auto& s : token_source_map) {
    //  if (pl == 3 && s.second) {
    //  std::cout << ">> >> Valid  " << pl << " " << mpi_rank << " " << s.first
    //    << " " << s.second << " " << vertex_metadata[graph->label_to_locator(s.first)]
    //    << " " << global_not_finished << std::endl;
    //  }
    //} 
    // Test

    // filter vertices, do not clear vertex_token_source_set for the vertices with label pattern_tp[0]      
    //if (mpi_rank == 0) { // Test
    //  std::cout << pattern_tp[0] << " " << pattern_indices_tp[0] << std::endl; // Test 
    //} // Test
    
    for (vitr_type vitr = graph->vertices_begin(); 
      vitr != graph->vertices_end(); ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_active[vertex] && 
        (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1])) {
        //std::cout << mpi_rank << " " << graph->locator_to_label(vertex) << " [l] " << vertex_token_source_set[vertex].size() << std::endl; // Test  
        continue; 
      } else {
        vertex_token_source_set[vertex].clear();
      }
    } 	
    	
    for(vitr_type vitr = graph->delegate_vertices_begin();
      vitr != graph->delegate_vertices_end(); ++vitr) {
      vloc_type vertex = *vitr;
      if (vertex_active[vertex] &&
        (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1])) {
        //std::cout << mpi_rank << " " << graph->locator_to_label(vertex) << " [d] " << vertex_token_source_set[vertex].size() << std::endl; // Test
        continue;
      } else {
        vertex_token_source_set[vertex].clear();
      } 
    }   
           
  } // if pattern_selected_vertices_tp
 
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
 
//--  //token_source_edge_set.clear(); // Important // edge aware
    
  time_start = MPI_Wtime();

  // token_passing_pattern_matching.hpp
  //token_passing_pattern_matching(graph, vertex_metadata, pattern_tp, 
  //  pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map, 
  //  token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp, pattern_found[pl], *edge_data_ptr, vertex_token_source_set);
  
//--  if (!pattern_seleted_edges_tp) { 
//--    token_source_edge_set.clear(); // edge aware 
//--  }

#ifdef TP_ASYNC
  // token_passing_pattern_matching_new.hpp
//  token_passing_pattern_matching(graph, vertex_metadata, pattern_tp,
//    pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map,
//    token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp, pattern_seleted_edges_tp, 
//    pattern_found[pl], edge_active, /**edge_data_ptr,*/ vertex_token_source_set, token_source_edge_set, vertex_active);
  //if (do_tds_tp) {
  if (pattern_is_tds_tp) {
    // token_passing_pattern_matching_nonunique_tds_batch_1.hpp
    // token_passing_pattern_matching_nonunique_iterative_tds_1.hpp
    // token_passing_pattern_matching_nonunique_iterative_tds_batch_1.hpp
    // token_passing_pattern_matching_nonunique_tds_batch_1.hpp 
    prunejuice::token_passing_pattern_matching<graph_type, Vertex, Edge, VertexData, 
      EdgeData, VertexMetadata, EdgeMetadata, VertexActive, 
      VertexUint8MapCollection, TemplateVertex, VertexStateMap, PatternGraph, 
      /*PatternUtilities*/PatternNonlocalConstraint, VertexUint8Map, VertexSetCollection, 
      DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>
      (graph, vertex_metadata, vertex_active, vertex_active_edges_map, 
       template_vertices, vertex_state_map, pattern_graph, ptrn_util_two, pl,
       token_source_map, vertex_token_source_set, 
       pattern_found[pl], tp_vertex_batch_size, paths_result_file, message_count);
    // token_passing_pattern_matching_nonunique_tds_1.hpp
    /*token_passing_pattern_matching<graph_type, Vertex, VertexMetaData, 
      decltype(pattern_tp), decltype(pattern_indices_tp), decltype(pattern_enumeration_tp),
      uint8_t, PatternGraph, VertexStateMap, VertexUint8Map, edge_data_t,
      VertexSetCollection, VertexActive, TemplateVertex, VertexUint8MapCollection, BitSet>
      (graph, vertex_metadata, pattern_tp, pattern_indices_tp, pattern_enumeration_tp, 
      vertex_rank, pattern_graph, vertex_state_map,
      token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
      pattern_found[pl], *edge_data_ptr, vertex_token_source_set, vertex_active,
      template_vertices, vertex_active_edges_map, pattern_selected_vertices_tp, 
      paths_result_file, message_count); // pass a boolean flag to indicate to use batching*/ 
  } else {     
    prunejuice::token_passing_pattern_matching<graph_type, VertexMetadata, decltype(pattern_tp), decltype(pattern_indices_tp), uint8_t, PatternGraph,
    VertexStateMap, VertexUint8Map, edge_data_t,
    VertexSetCollection, VertexActive, TemplateVertex, VertexUint8MapCollection, BitSet>(graph, vertex_metadata, pattern_tp,
    pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map,
    token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
    pattern_found[pl], *edge_data_ptr, vertex_token_source_set, vertex_active, 
    template_vertices, vertex_active_edges_map, pattern_selected_vertices_tp, //);
    pattern_selected_edges_tp, pattern_mark_join_vertex_tp,
    pattern_ignore_join_vertex_tp, pattern_join_vertex_tp, message_count);
  } 
#endif

// TODO: delete
//#ifdef TP_BATCH
  // token_passing_pattern_matching_batch.hpp
//  token_passing_pattern_matching(graph, vertex_metadata, pattern_tp,
//    pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map,
//    token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
//    pattern_found[pl], *edge_data_ptr, vertex_token_source_set, vertex_active, 
//    global_not_finished, token_source_deleted); 
//#endif 
 
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?    
  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Token Passing (Traversal) [" 
      << pl << "] : " << time_end - time_start << std::endl;
  }

  // result
//--  if(mpi_rank == 0) {
//--    superstep_result_file << global_itr_count << ", TP, "
//--      << pl << ", "
//--      << time_end - time_start << "\n";      
//--  }

//  } // for - loop over the patterns

  // pattren found ? // TODO: update for the new token passing loop
//  havoqgt::mpi::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
//  if(mpi_rank == 0) {   
//    for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
//      std::string s = pattern_found[pl] == 1 ? "true" : "false";
//      std::cout << "Token Passing [" << pl << "] | Found pattern : " << s << std::endl;
//    }
//  }

#ifdef TP_ASYNC

  // remove the invalid (token source) vertices from the vertex_state_map
  // for delegates, set vertex_active to false

  // TODO: In the case, a vertex is on multiple cycles/chains (not as the token source)
  // only invalidate it as a token source, but do not remove it from the vertex_state_map

  uint64_t remove_count = 0;      
  size_t token_source_pattern_indices_tp_index = 0; // TODO: this is confusing, update 
  
  bool is_token_source_map_not_empty = havoqgt::mpi_all_reduce
    (!token_source_map.empty(), std::greater<uint8_t>(), MPI_COMM_WORLD); // TODO: less did not work?
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

  pattern_token_source_found[pl] = is_token_source_map_not_empty;

//  if (pattern_selected_vertices_tp) {
    // TODO: delete
//    token_source_pattern_indices_tp_index = pattern_indices_tp.size() - 1; // Important     
    //std::cout << "token_source_map.size() " << token_source_map.size() << std::endl; // Test
//  }  

  for (auto& s : token_source_map) {
    if (!s.second) {
      //if (pl == 2) { // Test
      //  std::cout << "Invalid  " << pl << " " << mpi_rank << " " << s.first 
      //    << " " << s.second << " " << global_not_finished << std::endl;
      //} // Test 
      //if (pl <= 1) { // Test
      //  std::cout << "Invalid  " << pl << " " << mpi_rank << " " << s.first
      //    << " " << s.second << std::endl;  
      //} //Test

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
      //      //  std::cout << "MPI Rank: " << mpi_rank << " template vertices " << v_template_vertices << " global_not_finished " << global_not_finished << std::endl;
      //}  
    }  
    
    // Test
    //else { 
    //  if (pl == 2) {
    //    std::cout << "Valid  " << pl << " " << mpi_rank << " " << s.first
    //      << " " << s.second << " " << vertex_metadata[graph->label_to_locator(s.first)] 
    //      << " " << global_not_finished << std::endl;
    //  }   
    //}
    // Test
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

  vertex_active.all_min_reduce(); 
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

  // TODO: this is a temporary patch, forcing all the delegates to have no identity
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
  }
  MPI_Barrier(MPI_COMM_WORLD);

  template_vertices.all_max_reduce(); // ensure all the delegates have the same value as the controller
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

  // remove from vertex_state_map
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
 
      }
    }

    // Test 
    //if (pl == 2 && s.second) {
    //  std::cout << ">> Valid  " << pl << " " << mpi_rank << " " << s.first
    //    << " " << s.second << " " << vertex_metadata[graph->label_to_locator(s.first)] 
    //    << " " << global_not_finished << std::endl;
    //}
    // Test 
  }

//  if(mpi_rank == 0) {
//    std::cout << "Token Passing [" << pl << "] | MPI Rank [" << mpi_rank 
//      << "] | Removed " << remove_count << " vertices."<< std::endl; // TODO: not useful informtion
//  }

#endif // ifdef TP_ASYNC

  time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Token Passing [" << pl 
      << "] : " << time_end - time_start << std::endl;
  }

  // result
  if(mpi_rank == 0) {
    superstep_result_file << global_itr_count << ", TP, "
      << pl << ", "
      << time_end - time_start << "\n";
  }

  // Important : This may slow down things -only for presenting results
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

  // result
  
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

  if (is_token_source_map_not_empty) {

  //if (pl <= 1 && pattern_found[pl] == 1) { // Test
  //  std::cout << mpi_rank << " | Token Passing [" << pl << "] | Found Pattern : " << pattern_found[pl] << std::endl;             
  //} // Test    

  // pattren found ? // TODO: write output to file 
  ///havoqgt::mpi::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
  havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?   
  if(mpi_rank == 0) {   
//    for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
      std::string s = pattern_found[pl] == 1 ? "True" : "False";
      std::cout << "Token Passing [" << pl << "] | Found Subpattern : " << s << std::endl;
//    }
  }

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

  // verify global token source deleted status
  //token_source_deleted = havoqgt::mpi::mpi_all_reduce(token_source_deleted, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work  
  token_source_deleted = havoqgt::mpi_all_reduce(token_source_deleted, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 
  if(mpi_rank == 0) {
    std::cout << "Token Passing [" << pl << "] | Token Source Deleted Status : ";
    if (token_source_deleted) {
      std::cout << "Deleted" << std::endl;
    } else {
      std::cout << "Not Deleted" << std::endl;
    }
  } 

  } else {
    if(mpi_rank == 0) {
      std::cout << "Token Passing [" << pl << "] | No Token Source Found" << std::endl;   
    }
  } // is_token_source_map_not_empty

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  
  // interleave token passing with label propagation  
  if (token_source_deleted && pattern_interleave_label_propagation_tp) {

  bool global_not_finished_dummy = true; // TODO: do we need this?

  // lable propagation   
  label_propagation_time_start = MPI_Wtime();

  // label propagation pattern matching bsp, iterative 
//  label_propagation_pattern_matching_bsp<graph_type, VertexMetaData, VertexData, decltype(pattern), decltype(pattern_indices), 
//    /*VertexRank*/uint8_t, VertexActive, /*VertexIteration*/uint8_t, VertexStateMap, PatternGraph, EdgeActive, VertexSetCollection,
//    EdgeMetadata>
//    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active, 
//    vertex_iteration, vertex_state_map, pattern_graph, global_initstep, global_not_finished_dummy, 
//    global_itr_count, superstep_result_file, active_vertices_count_result_file, edge_active/**edge_data_ptr*/, edge_metadata);
#ifdef LCC_OLD 
  prunejuice::label_propagation_pattern_matching_bsp<graph_type, VertexMetadata, VertexData, decltype(pattern), decltype(pattern_indices),
    /*VertexRank*/uint8_t, VertexActive, /*VertexIteration*/uint8_t, VertexStateMap, PatternGraph, BitSet, TemplateVertex, 
    VertexUint8MapCollection>
    (graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active,
    vertex_iteration, vertex_state_map, pattern_graph, global_init_step, global_not_finished,
    global_itr_count, superstep_result_file, active_vertices_count_result_file, active_edges_count_result_file,
    template_vertices, vertex_active_edges_map, message_count_result_file);
#endif

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
    std::cout << "Fuzzy Pattern Matching Time | Local Constraint Checking (Interleaved) : " 
      << label_propagation_time_end - label_propagation_time_start << std::endl;
  }

  // result
  if(mpi_rank == 0) {
    step_result_file << global_itr_count << ", LP, "
      << (label_propagation_time_end - label_propagation_time_start) << "\n";
  }

  /////////////////////////////////////////////////////////////////////////////

  //if (global_init_step) {
  //  global_init_step = false;
  //} 

  // verify global termination condition
  //std::cout << "Fuzzy Pattern Matching | Global Not Finished status (local) : " << global_not_finished << std::endl; // Test
  // global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work
//--  global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD); 
//--  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

//--  if(mpi_rank == 0) {
//--    std::cout << "Fuzzy Pattern Matching | Global Finished Status : "; 
//--    if (global_not_finished) { 
//--      std::cout << "Continue" << std::endl;
//--    } else {
//--      std::cout << "Stop" << std::endl;
//--    } 
//--  }

  // global verification - are all vertex_state_maps empty
  // false - no active vertex left, true - active vertices left 
//--  global_active_vertex = true; //vertex_state_map.size() < 1 ? false : true;
//--  if (vertex_state_map.size() < 1) {
//    global_active_vertex = false;
//--  }

  // TODO: What is the issue here? Preventing stdout from this file beyond this point.  	
  // global_active_vertex = havoqgt::mpi::mpi_all_reduce(global_active_vertex, std::greater<uint8_t>(), MPI_COMM_WORLD); 
  // TODO: not working properly - why? // bool does not work
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

//  global_not_finished = global_active_vertex; // TODO: verify and fix

//--  if(mpi_rank == 0) {
//--    std::cout << "Fuzzy Pattern Matching | Global Active Vertex Status : ";
//--    if (global_active_vertex) {
//--      std::cout << "Active vertices left." << std::endl;
//--    } else {
//--      std::cout << "No active vertex left." << std::endl;
//--    }
//--  }
   	  
  } else { // if (token_source_deleted) // if (global_not_finished) - old code
    if(mpi_rank == 0) {
      std::cout << "Pattern Matching | Skipping Local Constraint Checking (Interleaved)." << std::endl;
    }        
  } // interleave token passing with label propagation

  // result
  paths_result_file.close();

  MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New 

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

  } // for - loop over the constraints and run token passing 

  // pattren found ? // TODO: write output to file 
  ///havoqgt::mpi::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
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

  /////////////////////////////////////////////////////////////////////////////

  // result
  // Important : This may slow down things -only for presenting results
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
      std::cout << "Pattern Matching | Skipping Token Passing." << std::endl;
    }
    global_not_finished = false;
  } // do token passing ?

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

//#endif 

  ///////////////////////////////////////////////////////////////////////////// 

  // global termination detection  
  //std::cout << "Fuzzy Pattern Matching | Global Not Finished status (local) : " << global_not_finished << std::endl; // Test
  //global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work  
  ///global_not_finished = havoqgt::mpi::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
  global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here 

  if(mpi_rank == 0) {
    std::cout << "Pattern Matching | Global Finished Status : ";
    if (global_not_finished) {
      std::cout << "Continue" << std::endl;
    } else {
      std::cout << "Stop" << std::endl;
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  double itr_time_end = MPI_Wtime(); 
//  if(mpi_rank == 0) { //TODO: sum of LCC and NLCC iterations
//    std::cout << "Pattern Matching Time | Pattern [" << ps 
//      << "] | Iteration [" << global_itr_count << "] : " 
//      << itr_time_end - itr_time_start << std::endl;  
//  }
 
  // result 
  if(mpi_rank == 0) {
    // iteration number, time
    itr_result_file << global_itr_count << ", "
      << (itr_time_end - itr_time_start) << "\n";
  }

  global_itr_count++; //TODO: sum of LCC and NLCC iterations

  // global termination
  // Test
  //if (global_itr_count > 0) {
  //  global_not_finished = false;
  //}
  //MPI_Barrier(MPI_COMM_WORLD);
  // Test

///  } while (global_not_finished); // application loop
 
  /////////////////////////////////////////////////////////////////////////////

  MPI_Barrier(MPI_COMM_WORLD);  
  double pattern_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Pattern Matching Time | Pattern [" << ps << "] : " 
      << pattern_time_end - pattern_time_start << std::endl;
  }
   
//  if(mpi_rank == 0) { //TODO: sum of LCC and NLCC iterations
//    std::cout << "Pattern Matching | Pattern [" << ps 
//      << "] | # Iterations : " << global_itr_count << std::endl;
//  }

  // TODO: state that if the pattern was found or not  
 
  /////////////////////////////////////////////////////////////////////////////

  // result
  if(mpi_rank == 0) {  
    // pattern set element ID, number of ranks, 
    // total number of iterations (lp+tp), time, 
    // #edges in the pattern, #vertices in the pattern,  
    // #token passing paths
    pattern_set_result_file << ps << ", " 
      << mpi_size << ", "
      << global_itr_count << ", " 
      << (pattern_time_end - pattern_time_start) << ", "
      << pattern_graph.edge_count << ", " 
      << pattern_graph.vertex_count << ", "
      << ptrn_util_two.input_patterns.size() << "\n"; 
  }

  // Important : This may slow things down -only for presenting results

  for (auto& v : vertex_state_map) {
    auto v_locator = graph->label_to_locator(v.first);
    if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
      BitSet v_template_vertices(template_vertices[v_locator]); 
      active_vertices_result_file << mpi_rank << ", "
        << v.first << ", " 
        << v.second.vertex_pattern_index << ", "
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
  
    } else if (!v_locator.is_delegate()) {
      BitSet v_template_vertices(template_vertices[v_locator]);
      active_vertices_result_file << mpi_rank << ", "
        << v.first << ", "
        << v.second.vertex_pattern_index << ", "
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

    }
  }

  // cleanup memeory
  //vertex_rank.clear(); // TODO: add clear() method to vertex_data.cpp   
  //vertex_active.clear(); // TODO: add clear() method to vertex_data.cpp
  //vertex_state_map.clear();

  /////////////////////////////////////////////////////////////////////////////

  // close files
  itr_result_file.close();
  step_result_file.close();
  superstep_result_file.close();
  active_vertices_count_result_file.close(); 
  active_vertices_result_file.close();
  active_edges_count_result_file.close();
  active_edges_result_file.close();
  //paths_result_file.close();
  message_count_result_file.close();

  // end of an elemnet in the pattern set
  } // for - loop over pattern set

  if (mpi_rank == 0) {
    pattern_set_result_file.close(); // close file
  }

  } // pattern matching  

  } // havoqgt_init
  ;
  // END Main MPI
  ///havoqgt::havoqgt_finalize();

  return 0;  
} // main
