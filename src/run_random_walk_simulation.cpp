#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#define POINT_RW_SIM

#include <havoqgt/environment.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/wiki_link_metadata.hpp>

#if defined INTERVAL_RW_SIM
#include <havoqgt/random_walkers/simulations/simple_rw_simulation.hpp>
#define SIM "Interval RW Simulation"
#elif defined POINT_RW_SIM
#include <havoqgt/random_walkers/simulations/uniform_time_rw_simulation.hpp>
#define SIM "Point RW Simulation"
#elif defined POINT_RW_SIM_V2
#include <havoqgt/random_walkers/simulations/uniform_time_rw_simulation_with_path.hpp>
#define SIM "Point RW Simulation with path traceback"
#endif


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

void parse_cmd_line(int argc, char** argv,
		    std::string& input_graph_filename, std::string& input_metadata_filename) {

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
  while( (c = getopt(argc, argv, "s:e:o:n:c:r:a:b:x:i:m:h " )) != -1 ) {
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
    default:
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

    havoqgt::get_environment();
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();
      if(mpi_rank == 0) {
	std::cout << "Executing " << SIM << std::endl;
      }
      //parse_cmd_line(argc, argv, input_graph_filename, input_metadata_filename);
      config_reader reader(argc, argv, "s:e:o:n:c:r:a:b:x:i:m:h ");
      auto identity =[] (std::string str) -> std::string {
	return str;
      };
      input_graph_filename = reader.get_value<std::string, decltype(identity)>('i', identity).second;
      input_metadata_filename = reader.get_value<std::string, decltype(identity)>('m', identity).second;

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

      simple_rw_simulation<graph_type, edge_data_type> simulation( graph, edge_data);
      simulation.parse_input_parameters( reader );
      start = TIMER_START();
      simulation.run();
      end   = TIMER_END();      
      havoqgt_env()->world_comm().barrier();
      MASTER_DO(
		std::chrono::duration<double> processing_time=TIME_DURATION_SECS();
		std::cout<<"Processing time: " <<processing_time.count() << " secs.\n";
		);
      havoqgt_env()->world_comm().barrier();
    }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}

