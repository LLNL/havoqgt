#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

#include <havoqgt/random_walkers/simulations/back_button_simulations.hpp>

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

std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);
  
  havoqgt::havoqgt_init(&argc, &argv);
  {
    std::string input_graph_filename;

    havoqgt::get_environment();
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();

      config_reader reader(argc, argv, "f:b:n:c:i:h ");
      auto identity =[] (std::string str) -> std::string {
	return str;
      };
      input_graph_filename = reader.get_value<std::string, decltype(identity)>('i', identity).second;

      havoqgt::distributed_db graph_ddb(havoqgt::db_open(), input_graph_filename.c_str());

      MASTER_MSG("Initializing");
      
      graph_type *graph = graph_ddb.get_segment_manager()->find<graph_type>
	("graph_obj").first;
      assert( graph != NULL) ;
      
      havoqgt_env()->world_comm().barrier();
      
      simple_rw_simulation<graph_type> simulation( graph);
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


