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

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef graph_type::edge_data<wiki_link_metadata, segment_manager_t> edge_data_type;

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

  friend std::ostream& operator<<(std::ostream& o, const _random_walker_vertex_data& rw_vertex_data) {
      std::copy( rw_vertex_data.map.cbegin()
		 , rw_vertex_data.map.cend()
		 , std::ostream_iterator<std::pair<uint64_t, value>>( o, "\n"));
      return o;
    }
};

typedef struct _random_walker_vertex_data random_walk_vertex_data;
//using clock_t = std::chrono::high_resolution_clock;
//using time_point_t = std::chrono::time_point(clock_t);

void parse_cmd_line(int argc, char** argv, std::string& input_graph_filename, std::string& input_metadata_filename
		    , uint64_t& time_start, uint64_t& time_end, uint64_t& time_offset , uint64_t& num_of_walkers) {

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
  while( (c = getopt(argc, argv, "i:s:e:o:m:n:h " )) != -1 ) {
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
    uint64_t time_start = 0;
    uint64_t time_end = 0;
    uint64_t time_offset = 0;
    uint64_t num_of_walkers = 0;
    havoqgt::get_environment();
    
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();

      parse_cmd_line(argc, argv, input_graph_filename, input_metadata_filename, time_start, time_end, time_offset 
		     , num_of_walkers );

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
      havoqgt_env()->world_comm().barrier();

      graph_type::vertex_locator source, target;

      temporal_random_walk_simulation(graph, edge_data, &rw_vertex_data, source, target, num_of_walkers, &start_time_steps);
      havoqgt_env()->world_comm().barrier();
      

      /*
      std::stringstream output_filename("/l/ssd/next_output_");
      output_filename << mpi_rank << "of" << mpi_size;
      std::ofstream out_file( output_filename.str(), std::ofstream::out );
      out_file << "Results goes here" << std::flush;
      MASTER_DO(
		std::cout << "Started writing now " << std::endl;
		);
      std::size_t i = 1;
      for( auto itr = graph->vertices_begin(); itr != graph->vertices_end(); ++itr) {
	out_file << i << " : " << rw_vertex_data[*itr] << "\n";
	i++;
      }
      out_file.flush();
      out_file.close();
      havoqgt_env()->world_comm().barrier();
      */

   }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
