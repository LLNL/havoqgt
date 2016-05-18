#ifndef BACK_BUTTON_SIMULATIONS
#define BACK_BUTTON_SIMULATIONS

#include <unistd.h>
#include <havoqgt/environment.hpp>

#include <havoqgt/non_temporal_rw_simulation.hpp>
#include <havoqgt/random_walkers/non_temporal_random_walker.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/fifo_queue.hpp>
#include <havoqgt/mpi.hpp>

#include <iostream>
#include <unordered_map>
#include <random>
#include <fstream>
#include <sstream>

namespace havoqgt { namespace mpi {
    /*********************************************************************************
     * Config Reader
     ********************************************************************************/
    class config_reader{
    public:
      std::string options;
      std::unordered_map<char, std::string> value_map;

      config_reader( int argc, char**argv , std::string _options) : options(_options) {
	char c;
	while( (c = getopt(argc, argv, options.c_str())) != -1 ) {
	  value_map[c]=optarg;
	  //	  std::cout << "inside config reader : " << c << " -> " << optarg << "\n";
	}
      }

      template<typename T, typename Op>
      std::pair<bool, T> get_value( char c, Op& op, T default_val = T() ) {
	if( value_map.find(c) == value_map.end()) {
	  return std::make_pair(false, default_val );
	}
	return std::make_pair(true, op(value_map[c]));
      }

      static uint64_t parse_to_long(std::string string) {
	return std::atoll(string.c_str());
      }
    };

    /*********************************************************************************
     * OutputIterator
     ********************************************************************************/
    template <typename Graph, typename random_walker_t>
    class output_iterator{
    public:
      uint64_t completed;
      output_iterator(Graph *_g ) : g(_g), completed(0) {
      }

      void operator()( const typename Graph::vertex_locator& vertex, const random_walker_t& random_walker) {
	uint64_t label = g->locator_to_label(vertex);
	if( hits.find(label) == hits.end() )
	  hits.insert( std::make_pair( 0, 0 ) );
	hits[ label ]++;
      }

      void close(std::string filename, int mpi_rank, int mpi_size) {
	std::stringstream ss;
	ss << filename << "_" << mpi_rank << "_of_" << mpi_size;
	std::ofstream ofs(ss.str(), std::ofstream::trunc );
	for( auto itr = hits.begin(); itr != hits.end(); ++itr) {
	  ofs << itr->first << " " << itr->second << "\n";
	}
	ofs.close();
      }

    private:
      std::unordered_map< uint64_t, uint64_t> hits;
      Graph* g;
    };

    /*********************************************************************************
     * Simulation & Input Iterator
     ********************************************************************************/
    template<typename Graph>
    class simple_rw_simulation{
    public:
  

      using random_walker_t = simple_random_walker<typename Graph::vertex_locator>;
      using output_iterator_t = output_iterator<Graph, random_walker_t>;
      typedef random_walk_simulation_visitor<Graph, output_iterator_t, random_walker_t > visitor_t;
      typedef std::unordered_set<typename Graph::vertex_locator, typename Graph::vertex_locator::hasher> vertex_locator_set;
      simple_rw_simulation(Graph *_graph) :
	graph(_graph) , itr(_graph){ 
	initialize();
      }

      void initialize() {
      }

      void read_local_vertices_from_file( std::string filename, vertex_locator_set& out){
	std::ifstream ifs( filename );
	std::string line;
	while( std::getline(ifs, line) ) {
	  // trim the line ??
	  uint64_t label = std::atoll( line.c_str() );
	  typename Graph::vertex_locator locator = graph->label_to_locator( label );
	  if( locator.owner() == havoqgt_env()->world_comm().rank())
	    out.insert(locator);
	}
      }

      bool parse_input_parameters(config_reader& reader) {
	//config_reader reader(argc, argv, "x:s:e:o:n:c:r:a:b:i:m:h ");

	//How to have them as variable of the config_reader class
	auto parse_to_long = [](std::string str) -> uint64_t {
	  return std::atoll(str.c_str());
	};
	auto identity = [] (std::string str) -> std::string {
	  return str;
	};
	auto parse_to_double = [](std::string str) -> double {
	  return std::atof( str.c_str());
	};

	using op_type = decltype(parse_to_long);
	using d_op_t = decltype( parse_to_double );

	num_of_walkers = reader.template get_value<uint64_t, op_type>( 'n', parse_to_long).second;
	max_steps = reader.template get_value<uint64_t, op_type>( 'c', parse_to_long).second;
	forward_jump = reader.template get_value<double, d_op_t>( 'f', parse_to_double).second;
	backup_jump = reader.template get_value<double, d_op_t>( 'b', parse_to_double).second;

	return true; // false means missing parameters
      }

      void run() {
	//	output_iterator_t itr(graph);

	visitor_t::set_output_iterator( &itr);
    
	random_walker_t::max_steps = max_steps;
	random_walker_t::forward_jump = forward_jump;
	random_walker_t::backup_jump = backup_jump;

	typedef visitor_queue<visitor_t, havoqgt::detail::fifo_queue, Graph> visitor_queue_type;
	BATCH_SIZE=100000;
	uint64_t sub_batch_processed = 0;
	uint64_t total_count = 0;
	uint64_t value = 0;
	for(uint64_t batch = 0; batch < num_of_walkers; batch++){
	  source_itr = graph->vertices_begin();
	  sub_batch_processed = 0;
	  do{
	    is_complete = false;
	    count = 0;
	    visitor_queue_type vq(graph);
	    //prepare states before execution
	    vq.template init_visitor_traversal_experiment<simple_rw_simulation<Graph>>(*this);
	    havoqgt_env()->world_comm().barrier();
	    uint64_t current_count = mpi_all_reduce( count, std::plus<uint64_t>(), havoqgt_env()->world_comm().comm() );
	    total_count += current_count;
	    if(havoqgt_env()->world_comm().rank() == 0) {
	      std::cout << "Batch " << ( batch + 1 ) << " [ Vertex Processed : " << ( total_count ) <<" ]\n";
	    }
	    sub_batch_processed++;
	    value = 0;
	    if( source_itr != graph->vertices_end()) value = 1;
	    value = mpi_all_reduce( value, std::plus<uint64_t>(), havoqgt_env()->world_comm().comm());
	  } while(value != 0 );
	}
	itr.close(std::string("/p/lscratchf/poudel1/logs/batch_back_button")
		  , havoqgt_env()->world_comm().rank(), havoqgt_env()->world_comm().size());
      }

      template< typename VisitorQueueHandle >
      void process( VisitorQueueHandle* queue ) {
	if(count == BATCH_SIZE || source_itr == graph->vertices_end()) {
	  is_complete = true;
	  return;
	}  
	random_walker_t rwalker(0);
	visitor_t visitor( *source_itr, rwalker);
	queue->queue_visitor( visitor );
	source_itr++;
	count++;
      }

      bool empty() {
	return is_complete;
      }
    private:
      /**Graph Stuffs**/
      Graph* graph;
      output_iterator_t itr;

      /**Experiment Related Constants**/
      uint64_t num_of_walkers;
      uint64_t max_steps;
      double forward_jump;
      double backup_jump;

      /**Experiment Progress Variables**/
      uint64_t BATCH_SIZE;
      uint64_t count;
      bool is_complete;
      typename Graph::vertex_iterator source_itr;
      
    };

  } }
#endif
