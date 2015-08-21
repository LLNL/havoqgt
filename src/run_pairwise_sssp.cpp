#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/impl/vertex_data.hpp>
#include <havoqgt/pairwise_sssp.hpp>
#include <havoqgt/child_count_visitor.hpp>
#include <havoqgt/betweenness_centrality.hpp>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef graph_type::edge_data<flow, segment_manager_t> edge_data_type;

void parse_cmd_line(int argc, char** argv, std::string& input_graph_filename
		    , std::string& input_metadata_filename
		    , std::vector<uint64_t>& vertices_list
		    , nano_time& start_time
		    , nano_time& waiting_time
		    , uint32_t& source_count) {
#if 0 //DEBUG CODE
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for( size_t i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }
#endif

  bool found_input_graph_filename = false;
  bool found_input_metadata_filename = false;

  char c;
  while( (c = getopt(argc, argv, "i:m:s:w:c:h ")) != -1 ) {
#if 0 //DEBUG CODE
    if( havoqgt_env()->world_comm().rank() == 0) {
      std::cout << "Processing " << c << " , Value : " << optarg << std::endl;
    }
#endif

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
      start_time = nano_time(optarg);
      break;
    case 'w':
      waiting_time = nano_time(optarg);
      break;
    case 'c':
      source_count = std::atoi(optarg);
      break;
    }
  }
    if( !found_input_graph_filename || !found_input_metadata_filename ) {
      std::cout << "Necessary input file not provided" << std::endl;
      exit(-1);
    }
  
    //all the other to the vector list
    for(int index = optind; index < argc; ++index) {
      vertices_list.push_back( atoll(argv[index]) );
    }
  
}

std::vector<uint64_t> random_vertices_list(graph_type* graph, int source_count) {
  std::vector<uint64_t> random_vertices;
  srand(time(NULL));
  int count = 0;
  do { //simple implementation - non-deterministic
  for(graph_type::vertex_iterator itr =  graph->vertices_begin();
      itr != graph->vertices_end() && count < source_count;
      itr++) {
    if( rand() % 2  == 0) {
      random_vertices.push_back( graph->locator_to_label(*itr));
      count++;
    }  
  }
  }while( count < source_count);
  
  return random_vertices;
}

clock_t start;
clock_t end;

void start_timer() { start = clock(); }
void end_timer() { end = clock(); }
double execution_time() { return ((double)(end - start))/CLOCKS_PER_SEC; }

int main( int argc, char **argv) {

  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
    std::string input_graph_filename;
    std::string input_metadata_filename;
    std::vector<uint64_t> vertices_list;
    nano_time waiting_time;
    nano_time start_time;
    uint32_t source_count;
    
    havoqgt::get_environment();

    {
      mpi_rank = havoqgt_env()->world_comm().rank();
      mpi_size = havoqgt_env()->world_comm().size();

      havoqgt_env()->world_comm().barrier();

      parse_cmd_line( argc, argv, input_graph_filename, input_metadata_filename, vertices_list, start_time, waiting_time, source_count);

      havoqgt::distributed_db graph_db(havoqgt::db_open(), input_graph_filename.c_str());
      havoqgt::distributed_db metadata_db(havoqgt::db_open(), input_metadata_filename.c_str());

      //MASTER_MSG("Initializing graph database");
      graph_type *graph = graph_db.get_segment_manager()->find<graph_type>("graph_obj").first;
      assert( graph != NULL );

      //MASTER_MSG("Initializing metadata database");
      edge_data_type *edge_data = metadata_db.get_segment_manager()->find<edge_data_type>("edge_data").first;
      assert( graph != NULL );

      havoqgt_env()->world_comm().barrier();

      if(vertices_list.size() == 0){
	/*for(uint64_t vertex_id = 0; vertex_id <= graph->max_global_vertex_id(); vertex_id++) {
	  vertices_list.push_back(vertex_id);
	  }*/
	if(source_count != 0)
	  vertices_list = random_vertices_list( graph, source_count);
	else{
	  for(graph_type::vertex_iterator itr = graph->vertices_begin();
	      itr != graph->vertices_end();
	      itr++ ) {
	    vertices_list.push_back( graph->locator_to_label(*itr));
	  }
	}
      }// else std::sort( vertices_list.begin(), vertices_list.end());
      MASTER_DO( std::cout << source_count <<","; );

      typedef vertex_state_map<graph_type, nano_time> vertex_state_type;
      vertex_state_type::traversal_count = 0;
      typedef graph_type::vertex_data< vertex_state_type, std::allocator<vertex_state_type> > vertex_data_type;
	
      typedef pairwise_sssp_visitor<graph_type, edge_data_type, vertex_data_type> pairwise_sssp_visitor_type;

      vertex_state_type::set_graph( graph );

      vertex_data_type pairwise_vertex_data(*graph);
      											      
      std::vector<graph_type::vertex_locator> vertices_locator_list;
      for( auto vertices_itr = vertices_list.begin(); vertices_itr != vertices_list.end(); ++vertices_itr ) {
	vertices_locator_list.push_back(graph->label_to_locator(*vertices_itr));
      }

      //MASTER_MSG("Starting pairwise SSSP Computation");
      havoqgt_env()->world_comm().barrier();
      MASTER_DO(
		start_timer();
		);

      //havoqgt_env()->world_comm().barrier();

      run_pairwise_sssp( graph, edge_data, &pairwise_vertex_data, vertices_locator_list, start_time, waiting_time);
      havoqgt_env()->world_comm().barrier();
      MASTER_DO(
		end_timer();
		);

      
      uint64_t total_count = mpi_all_reduce( vertex_state_type::traversal_count,
					     std::plus<uint64_t>() ,
					     havoqgt_env()->world_comm().comm() );

      //MASTER_DO( std::cout << "Total traversed vertices: " << total_count << std::endl; );
      MASTER_DO(std::cout << total_count << ",";);
      {
#if 0
	for(graph_type::vertex_iterator vert_itr = graph->vertices_begin();
	    vert_itr != graph->vertices_end();
	    vert_itr++) {
	  for(auto itr = pairwise_vertex_data[*vert_itr].source_begin();
	      itr != pairwise_vertex_data[*vert_itr].source_end();
	      itr++) {
	    std::cout << graph->locator_to_label(*vert_itr) << " : "
		      << pairwise_vertex_data[*vert_itr].get_state(graph->label_to_locator(itr->first)).to_string()
		      <<" : " << itr->first << " --> "
		      << graph->locator_to_label(itr->second) << std::endl;
	  }
	}
#endif
	//MASTER_MSG("Preliminary pairwise SSSP computation completed. " );
	/*MASTER_DO(
		  std::cout << "Time of Execution: Pairwise SSSP : " << execution_time() << std::endl;
		  );*/
	MASTER_DO(std::cout << execution_time() <<",";);
      }
      havoqgt_env()->world_comm().barrier();

      typedef vertex_state_map<graph_type, uint64_t> count_vertex_state_type;
      count_vertex_state_type::traversal_count=0;
      typedef graph_type::vertex_data< count_vertex_state_type, std::allocator<count_vertex_state_type>>
	count_vertex_data_type;

      count_vertex_state_type::set_graph(graph);
      count_vertex_data_type child_count_vertex_data(*graph);
      
      havoqgt_env()->world_comm().barrier();
      MASTER_DO( start_timer(); )
      count_child(graph, &child_count_vertex_data, &pairwise_vertex_data, vertices_locator_list);
      havoqgt_env()->world_comm().barrier();
      MASTER_DO( end_timer(); )
      total_count = mpi_all_reduce(count_vertex_state_type::traversal_count
				   , std::plus<uint64_t>()
				   , havoqgt_env()->world_comm().comm() );
      //MASTER_DO( std::cout << "Total traversed vertices: " << total_count << std::endl; );
      MASTER_DO(std::cout << total_count << ",";);
      //MASTER_MSG(" Counting child complete " );
      //MASTER_DO( std::cout << "Time of Execution: Child Count :" << execution_time() << std::endl; );
      MASTER_DO(std::cout << execution_time() << ","; );
            {
#if 0
	for( graph_type::vertex_iterator itr = graph->vertices_begin();
	     itr!= graph->vertices_end();
	     itr++ ) {
	  bool print = false;
	  std::stringstream msg;
	  
	  msg << graph->locator_to_label(*itr) << "  : ";
	  for(auto source_itr = pairwise_vertex_data[*itr].source_begin();
	      source_itr != pairwise_vertex_data[*itr].source_end();
	     source_itr++) {
	    
	    graph_type::vertex_locator vl = graph->label_to_locator(source_itr->first);
	    if( child_count_vertex_data[*itr].get_state(vl) == 0) continue;
	    print = true;
	    msg << source_itr->first << " ( " << child_count_vertex_data[*itr].get_state(vl) << ") | ";
	  
	  }
	  if(print) std::cout <<msg.str() <<  std::endl;
	}
#endif
	    }
      
      count_vertex_data_type shortest_path_vertex_data(*graph);
      count_vertex_state_type::traversal_count=0;
      havoqgt_env()->world_comm().barrier();
      MASTER_DO( start_timer(); )
      betweenness_centrality(graph, &child_count_vertex_data, &pairwise_vertex_data, &shortest_path_vertex_data, vertices_locator_list);

      havoqgt_env()->world_comm().barrier();
      MASTER_DO( end_timer(); )
      total_count = mpi_all_reduce( count_vertex_state_type::traversal_count
				    , std::plus<uint64_t>()
				    , havoqgt_env()->world_comm().comm() );
      //MASTER_DO( std::cout << "Total traversed vertices: " << total_count << std::endl; );
      MASTER_DO(std::cout << total_count << ",";);
      //MASTER_MSG(" Betweenness centrality computation completed. " );
      //MASTER_DO( std::cout << "Time of Execution : Betweenness Centrality : " << execution_time() << std::endl; );
      MASTER_DO(std::cout << execution_time() << std::endl;);
      typedef graph_type::vertex_data< uint64_t, std::allocator<uint64_t> > aggregate_vertex_data_type;
      aggregate_vertex_data_type shortest_path_aggregate( *graph );
      shortest_path_aggregate.reset( 0 );

      std::stringstream str;
      str << "betweenness_" << start_time.seconds << "_" << waiting_time.seconds << "_" <<  mpi_rank << "_of_" << mpi_size;
      
      std::ofstream out(str.str());
      
      for( auto itr = graph->vertices_begin(); itr != graph->vertices_end(); itr++) {
	for( auto source_itr = pairwise_vertex_data[*itr].source_begin();
	     source_itr != pairwise_vertex_data[*itr].source_end();
	     source_itr++ ) {
	  
	  shortest_path_aggregate[*itr] += (shortest_path_vertex_data[*itr].get_state(
						      graph->label_to_locator(source_itr->first)
									       ));
	}
	if(shortest_path_aggregate[*itr] != 0)
	  out << graph->locator_to_label(*itr) << ";" << shortest_path_aggregate[*itr] << ";" << std::endl;
   
      }
      out.flush();
      out.close();
      
      havoqgt_env()->world_comm().barrier();
      //MASTER_MSG(" Writing to file completed");
    }
    
  }

  havoqgt::havoqgt_finalize();
}
