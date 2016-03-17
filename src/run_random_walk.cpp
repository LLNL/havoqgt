#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/temporal_random_walk_simulation.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/impl/vertex_data.hpp>
#include <havoqgt/wiki_link_metadata.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <unordered_map>
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

std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

template<typename A, typename B>
std::ostream& operator<<(std::ostream& o, const std::pair< A, B>& out) {
  return o << out.first << " : " << out.second;
}

struct _random_walker_vertex_data {

  struct _value{
    uint64_t min;
    uint64_t max; 
    uint64_t sum;
    uint64_t count;

    _value() : min(std::numeric_limits<uint64_t>::max()), max(0), sum(0), count(0) {}

    friend std::ostream& operator<<(std::ostream& o, const _value& val ) {
      return o << " [ " << val.min <<", "<< val.max << " ] "
	       << val.sum << " " << val.count;
    }
  };

  typedef struct _value value;

  std::unordered_map<uint64_t, _value> map;
  
  value& at( uint64_t time ) {
    if(map.find(time) == map.end()) map[time] = value();
    return map[time];
  }

  uint64_t traversal_passed_through() {
    uint64_t count = 0;
    for( auto map_itr = map.begin(); map_itr != map.end(); ++map_itr) {
      count += map_itr->second.count;
    }
    return count;
  }

  friend std::ostream& operator<<(std::ostream& o, const _random_walker_vertex_data& rw_vertex_data) {
      std::copy( rw_vertex_data.map.cbegin()
		 , rw_vertex_data.map.cend()
		 , std::ostream_iterator<std::pair<uint64_t, value>>( o, " | "));
      return o;
    }
};

typedef struct _random_walker_vertex_data random_walk_vertex_data;
//using clock_t = std::chrono::high_resolution_clock;
//using time_point_t = std::chrono::time_point(clock_t);

void parse_cmd_line(int argc, char** argv, std::string& input_graph_filename, std::string& input_metadata_filename
		    , uint64_t& time_start, uint64_t& time_end, uint64_t& time_offset , uint64_t& num_of_walkers
		    , uint64_t& max_steps, std::string& vertex_data_out_filename ) {

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
  while( (c = getopt(argc, argv, "i:s:e:o:m:n:c:r:h " )) != -1 ) {
    if( havoqgt_env()->world_comm().rank() == 0) {
      std::cout << "Processing " << c << " , Value : " << optarg << std::endl;
    }
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
    case 'n':
      num_of_walkers = atoll(optarg);
      break;
    case 'c':
      max_steps = atoll(optarg);
      break;
    case 'r':
      vertex_data_out_filename = optarg;
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
}


int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
    std::string input_graph_filename;
    std::string input_metadata_filename;
    std::string vertex_data_out_filename;
    uint64_t time_start = 0;
    uint64_t time_end = 0;
    uint64_t time_offset = 0;
    uint64_t num_of_walkers = 0;
    uint64_t max_steps = 0;
    havoqgt::get_environment();
    
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();

      parse_cmd_line(argc, argv, input_graph_filename, input_metadata_filename, time_start, time_end, time_offset 
		     , num_of_walkers, max_steps, vertex_data_out_filename );

      std::vector<hmpi::time_point_t> start_time_steps;
      for(uint64_t time = time_start; time <= time_end; time += time_offset) {
	start_time_steps.push_back( hmpi::time_point_t(std::chrono::seconds(time)));
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

      graph_type::vertex_data<random_walk_vertex_data, std::allocator<random_walk_vertex_data>> rw_vertex_data(*graph);
      rw_vertex_data.reset(random_walk_vertex_data());

      graph_type::vertex_locator source, target;

      havoqgt_env()->world_comm().barrier();
      start=TIMER_START();
      temporal_random_walk_simulation(graph, edge_data, &rw_vertex_data, source, target, num_of_walkers, &start_time_steps, max_steps);
      havoqgt_env()->world_comm().barrier();
      end=TIMER_END();
      
      std::string output_filename(vertex_data_out_filename);
      MASTER_DO(
		for(int i = 0; i < mpi_size; i++) {
		  std::stringstream ss;
		  ss << output_filename << i << "of" << mpi_size;
		  std::ofstream out_file( ss.str().c_str(), std::ofstream::trunc);
		  out_file << "Results" << std::endl;
		  out_file.flush();
		  out_file.close();
		}
		);
      havoqgt_env()->world_comm().barrier();
      std::stringstream ss;
      ss << output_filename << mpi_rank << "of" << mpi_size;
      std::ofstream out_file( ss.str().c_str(), std::ofstream::app );

      MASTER_DO(
		std::cout << "Started writing now. :) " << std::endl;
		);
      uint64_t total_steps_rank = 0;
      uint64_t ver_count = 0;
      for( auto itr = graph->vertices_begin(); itr != graph->vertices_end(); ++itr) {
	if( itr->owner() == mpi_rank && rw_vertex_data[*itr].map.size() > 0 ) {
	  out_file << graph->locator_to_label(*itr) << " : " << rw_vertex_data[*itr] << "\n";
	  total_steps_rank += rw_vertex_data[*itr].traversal_passed_through();
	}
	ver_count++;
      }
      out_file.flush();
      out_file.close();
      havoqgt_env()->world_comm().barrier();
      uint64_t total_count = mpi_all_reduce( total_steps_rank, std::plus<uint64_t>(), havoqgt_env()->world_comm().comm());
      uint64_t total_ver_count = mpi_all_reduce( ver_count, std::plus<uint64_t>(), havoqgt_env()->world_comm().comm());
      havoqgt_env()->world_comm().barrier();
      MASTER_DO(
		//REPORT
		std::chrono::duration<double> processing_time = TIME_DURATION_SECS();
		std::cout << "Total edges required: " << total_ver_count * 10000 << std::endl;
		std::cout << "Total edges traversals: " << total_count << " in " 
		<< processing_time.count() <<"secs. Ratio = " << ((double)total_count)/((double)total_ver_count * 10000) 
		<< std::endl;
		std::cout << "Edges per seconds: " << ( ((double) total_count )/ ((double) processing_time.count()) ) << std::endl;
		);
      havoqgt_env()->world_comm().barrier();
   }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
