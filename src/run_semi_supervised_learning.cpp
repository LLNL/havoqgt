#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/temporal_semi_supervised_learning.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/impl/vertex_data.hpp>
#include <havoqgt/wiki_link_metadata.hpp>
#include <havoqgt/segment_tree.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>

#define TIMER_START() std::chrono::high_resolution_clock::now()
#define TIMER_END()   std::chrono::high_resolution_clock::now()
#define TIME_DURATION_SECS() (end - start)


using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef graph_type::edge_data<wiki_link_metadata, segment_manager_t> edge_data_type;

typedef segment_tree<uint64_t, graph_type::edge_iterator> segment_tree_t;
typedef boost::interprocess::allocator<segment_tree_t, segment_manager_t> allocator_t;
//typedef graph_type::vertex_data<segment_tree_t, allocator_t> vertex_data_type;
typedef graph_type::vertex_data<segment_tree_t, std::allocator<segment_tree_t>> vertex_data_type;


/*struct vertex_hasher
      {
	std::size_t operator()(const typename graph_type::vertex_locator& other) const {
	  return other.hash();
	}
      };
*/
std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

typedef struct _random_walker_vertex_data random_walk_vertex_data;

void parse_cmd_line(int argc, char** argv, std::string& input_graph_filename, std::string& input_metadata_filename
		    , uint64_t& time_start, uint64_t& time_end, uint64_t& time_offset , uint64_t& num_of_walkers
		    , uint64_t& max_steps 
		    , uint64_t& restart_max, std::string& vertex_data_out_filename
		    , std::string& source_label_filename
		    , std::string& destination_label_filename
		    , double& gbyte_per_rank ) {

  if( havoqgt_env()->world_comm().rank() == 0 ) {
    std::cout <<  "CMD line:";
    for( int i = 0; i < argc; ++i ) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_graph_filename = false;
  bool found_input_metadata_filename = false;

  char c;
  while( (c = getopt(argc, argv, "i:s:e:o:m:n:c:r:x:f:h " )) != -1 ) {
    if( havoqgt_env()->world_comm().rank() == 0) {
      std::cout << "Processing " << c << " , Value : " << optarg << std::endl;
    }
    //defaults
    
    
    switch (c) {
    case 'i':
      found_input_graph_filename = true;
      input_graph_filename = optarg;
      break;
    case 'm':
      found_input_metadata_filename = true;
      input_metadata_filename = optarg;
      break;
    case 's':
      time_start = atoll(optarg);
      break;
    case 'e':
      time_end = atoll(optarg);
      break;
    case 'o':
      time_offset = atoll(optarg);
      break;
    case 'n':
      num_of_walkers = atoll(optarg);
      break;
    case 'c':
      max_steps = atoll(optarg);
      break;
    case 'r':
      vertex_data_out_filename = optarg;
      break;
    case 'x':
      restart_max = atoll(optarg);
      break;
    case 'f':
      gbyte_per_rank = atof(optarg);
      break;
    default:
      std::cout << "Unrecognized option: " <<  c <<", ignore. " << std::endl;
      break;
    }
  }

  if(!found_input_graph_filename || !found_input_metadata_filename ) {
    std::cout << "Necessary input file not provided" << std::endl;
    exit(-1);
  }
  
  source_label_filename = std::string(argv[optind]);
  destination_label_filename = std::string(argv[optind + 1]);
}

void read_label_from_file( std::string filename, std::vector<uint64_t>& labels) {
  std::ifstream ifs(filename);
  std::string line;
  while( std::getline(ifs, line)) {
    labels.push_back( std::atol( line.c_str() ) );
  }
  ifs.close();
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
    std::string input_graph_filename;
    std::string input_metadata_filename;
    std::string vertex_data_out_filename;

    std::string source_label_filename;
    std::string destination_label_filename;

    uint64_t time_start = 0;
    uint64_t time_end = 0;
    uint64_t time_offset = 0;
    uint64_t num_of_walkers = 1;
    uint64_t max_steps = 0;
    uint64_t restart_max = 0;
    std::vector<uint64_t> target_labels;
    std::vector<uint64_t> source_labels;
    
    double gbyte_per_rank;

    havoqgt::get_environment();
    
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();

      parse_cmd_line(argc, argv, input_graph_filename, input_metadata_filename, time_start, time_end, time_offset 
		     , num_of_walkers, max_steps, restart_max, vertex_data_out_filename, source_label_filename, destination_label_filename
		     , gbyte_per_rank);

      read_label_from_file( source_label_filename, source_labels);
      read_label_from_file( destination_label_filename, target_labels);

      std::vector<uint64_t> start_time_steps;
      for(uint64_t time = time_start; time <= time_end; time += time_offset) {
	start_time_steps.push_back(time);
      }

      havoqgt::distributed_db graph_ddb(havoqgt::db_open(), input_graph_filename.c_str());
      havoqgt::distributed_db metadata_ddb(havoqgt::db_open(), input_metadata_filename.c_str());

      MASTER_MSG("Initializing");
      
      graph_type *graph = graph_ddb.get_segment_manager()->find<graph_type>
	("graph_obj").first;
      assert( graph != NULL) ;
      
      havoqgt_env()->world_comm().barrier();
      
      edge_data_type *edge_data = metadata_ddb.get_segment_manager()->find<edge_data_type>
	("edge_data").first;
      assert( edge_data != NULL );
      
      //--------CREATE SEGMENT TREE------------------
      /*      MASTER_MSG("Creating Allocations for SEG_T");
      havoqgt::distributed_db segment_tree_ddb( havoqgt::db_create(), "/l/ssd/seg_tree", 4 * gbyte_per_rank );
      segment_manager_t *segment_manager= segment_tree_ddb.get_segment_manager();
      allocator_t alloc_inst( segment_manager );
      vertex_data_type *segment_tree_data = segment_manager->construct<vertex_data_type>("segment_tree")(*graph, alloc_inst);
      */
      vertex_data_type segment_tree_data(*graph);
      havoqgt_env()->world_comm().barrier();
      MASTER_MSG("Starting SEG_T construction");
      uint64_t missed_edges = 0;
      start = TIMER_START();
      for(auto itr = graph->vertices_begin(); itr != graph->vertices_end(); ++itr) {
	graph_type::vertex_locator vertex = *itr;
	if( vertex.owner() != mpi_rank) {
	  std::cout << "Shouldn't be here" << std::endl;
	  continue;
	}
	segment_tree_t tree(graph->edges_end(vertex) - graph->edges_begin(vertex));// = segment_tree_data[vertex];//( graph->edges_end(vertex) - graph->edges_begin(vertex));
	std::vector<uint64_t> points;
	uint64_t count = 0;
	for( auto edge_itr = graph->edges_begin(vertex); edge_itr != graph->edges_end(vertex); ++edge_itr ) {
	  wiki_link_metadata& metadata = (*edge_data)[edge_itr];
	  uint64_t start_point = metadata.start_time();
	  if( start_point == 0) { missed_edges++; continue; }
	  else count++;
	  uint64_t end_point = metadata.end_time();
	  if( end_point == 0 ) end_point = std::numeric_limits<uint64_t>::max();
	  points.push_back( start_point);
	  points.push_back( end_point);
	}
	if(count == 0) continue;
	
	tree.balanced_tree( std::less<uint64_t>(), points);

	for( auto edge_itr = graph->edges_begin(vertex); edge_itr != graph->edges_end(vertex); ++edge_itr ) {
	  auto metadata = (*edge_data)[edge_itr];
	  uint64_t start_point = metadata.start_time();
	  uint64_t end_point = metadata.end_time();
	  if( end_point == 0 ) end_point = std::numeric_limits<uint64_t>::max();
	  tree.add( std::make_pair( start_point, end_point), edge_itr );
	}

	(segment_tree_data)[vertex] = tree;
	
      }
      havoqgt_env()->world_comm().barrier();
      end = TIMER_END();
      uint64_t global_missed_edges = mpi_all_reduce( missed_edges, std::plus<int>(), havoqgt_env()->world_comm().comm());
      havoqgt_env()->world_comm().barrier();
      MASTER_DO(
		//REPORT
		std::chrono::duration<double> processing_time = TIME_DURATION_SECS();
		std::cout << "Preprocessing time  " 
		<< processing_time.count() << " secs";
		std::cout <<"  Total Missed Edges: " << global_missed_edges << std::endl;
		);
      havoqgt_env()->world_comm().barrier();

      //----------------------------------------------
      
      std::vector<graph_type::vertex_locator> targets;
      std::vector<graph_type::vertex_locator> sources;

      for( auto label: target_labels) {
	targets.push_back( graph->label_to_locator(label - 1) );
      }
      
      /*
      for( auto label: source_labels ) {
	graph_type::vertex_locator vertex = graph->label_to_locator(label - 1);
	if( vertex.owner() == mpi_rank ) { 
	  sources.push_back( vertex );
	}
	}*/

      std::string output_filename(vertex_data_out_filename);

      MASTER_DO(
		for(int i = 0; i < mpi_size; i++) {
		  std::stringstream ss;
		  ss << output_filename << "_" <<  i << "of" << mpi_size;
		  std::ofstream out_file( ss.str().c_str(), std::ofstream::trunc);
		  out_file << "Results" << std::endl;
		  out_file.flush();
		  out_file.close();
		}
		);
      std::stringstream ss;
      ss << output_filename << "_" << mpi_rank << "of" << mpi_size;
      output_filename = ss.str();
      
      /* ----------- Testing Segment Tree <starts>--------------
      havoqgt_env()->world_comm().barrier();
      std::ofstream out_file( output_filename.c_str(), std::ofstream::trunc);
      //execute test_cases
      for( auto &vertex : sources) {
	start = TIMER_START();
	uint64_t size_test = 0;
	for(int steps= 0; steps<num_of_walkers; steps++) {
	  for( auto &time : start_time_steps) {	    
	    size_test = 0;
	    (segment_tree_data[vertex]).query( time, size_test);
	    graph_type::edge_iterator itr;
	    if(size_test > 0)
	      (segment_tree_data[vertex]).query_ith( time, size_test - 1, itr);
	  }
	}
	end=TIMER_START();
	{
	  std::chrono::duration<double> processing_time = TIME_DURATION_SECS();
	
	  out_file << "Processing " << (graph->edges_end(vertex) - graph->edges_begin(vertex)) << " edges in time : "
		   << processing_time.count() << std::endl;
	  }
      }
      out_file.close();
      ---------------------Testing Segment Tree <ends>----------------*/

      havoqgt_env()->world_comm().barrier();
      start=TIMER_START();
      temporal_random_walk_simulation(graph, edge_data, &segment_tree_data, sources, targets, num_of_walkers, &start_time_steps, max_steps, restart_max, output_filename);
      havoqgt_env()->world_comm().barrier();
      end=TIMER_END();
      
       MASTER_DO(
		//REPORT
		std::chrono::duration<double> processing_time = TIME_DURATION_SECS();
		std::cout << "Total edges traversed in " 
		<< processing_time.count() << " secs"; 
		);
      havoqgt_env()->world_comm().barrier();
   }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
