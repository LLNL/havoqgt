#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/sssp_in_time.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/impl/vertex_data.hpp>
#include <havoqgt/ingest_flow_edge_list.hpp>
#include <unistd.h>

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef graph_type::edge_data<flow, segment_manager_t> edge_data_type;

void parse_cmd_line(int argc, char** argv, std::string& input_graph_filename, std::string& input_metadata_filename
		    , uint64_t& source_vertex_label ) {

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
  while( (c = getopt(argc, argv, "i:s:m:h " )) != -1 ) {
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
      source_vertex_label = atoll(optarg);
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
    uint64_t source_vertex_label = 0;

    havoqgt::get_environment();
    
    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();
      
      havoqgt_env()->world_comm().barrier();

      parse_cmd_line(argc, argv, input_graph_filename, input_metadata_filename, source_vertex_label);

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
      
      havoqgt_env()->world_comm().barrier();

      graph_type::vertex_locator source = graph->label_to_locator( source_vertex_label );
      nano_time start;
      graph_type::edge_iterator edge_itr(graph);
      if( source.owner() == mpi_rank) {
	edge_itr = (graph->edges_begin(source)); //Run SSSP
      }
     
      do
      {
	if( source.owner() == mpi_rank ) {
	  if( (edge_itr) == graph->edges_end(source)) {
	    //send termination command as nano_time seconds negative
	    start = nano_time(-1, 0);
	  }else {
	    start = (*edge_data)[edge_itr].start_time();
	  }
	}
	mpi_bcast(start, source.owner(), havoqgt_env()->world_comm().comm());
	
	havoqgt_env()->world_comm().barrier();
	if(start.seconds == -1) break;
	if( source.owner() == mpi_rank ) edge_itr++;
	typedef typename graph_type::vertex_locator vertex_locator;
	typedef typename graph_type::vertex_iterator vertex_iterator;
	
	graph_type::vertex_data<nano_time, std::allocator<nano_time> >      sssp_arrival_time_data(*graph);
	sssp_arrival_time_data.reset(nano_time(std::numeric_limits<uint32_t>::max()
					       ,std::numeric_limits<uint32_t>::max()));

	graph_type::vertex_data<vertex_locator, std::allocator<vertex_locator>> sssp_parent_data(*graph);
	
	havoqgt_env()->world_comm().barrier();
	
	vertex_locator source = graph->label_to_locator( source_vertex_label );
	//nano_time start = (*edge_data)[edge_itr].start_time();
	MASTER_DO(
		  std::cout << start.to_string() << ";";
		  );
	sssp_in_time(graph, edge_data, &sssp_arrival_time_data, &sssp_parent_data, source, start);
	havoqgt_env()->world_comm().barrier();

	int connected_components = 0;
	for(vertex_iterator curr_vertex_itr = graph->vertices_begin();
	    curr_vertex_itr != graph->vertices_end();
	    curr_vertex_itr++ ) {
	  vertex_locator curr_vertex = *curr_vertex_itr;
	  if(sssp_arrival_time_data[curr_vertex] == nano_time(std::numeric_limits<uint32_t>::max(),
							      std::numeric_limits<uint32_t>::max() ) )
	    continue;

	  if(curr_vertex == source) continue;
#if 0	  
	  std::cout << graph->locator_to_label(curr_vertex) << " : " << "arrival time : "
		    << sssp_arrival_time_data[curr_vertex].to_string()
		    << "  parent label : "
		    << graph->locator_to_label(sssp_parent_data[curr_vertex])
		    << std::endl;
#endif
	  connected_components++;
	}
	havoqgt_env()->world_comm().barrier();
	
	int global_connected_components = mpi_all_reduce( connected_components, std::plus<int>(),
							  havoqgt_env()->world_comm().comm() );
	havoqgt_env()->world_comm().barrier();
	
	MASTER_DO(
		  std::cout << global_connected_components << std::endl;
        );
	
	//MASTER_MSG( "SSSP in time successfully completed" );
      }while(true);
      
    }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
