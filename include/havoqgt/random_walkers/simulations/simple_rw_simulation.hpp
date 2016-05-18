#ifndef _SIMPLE_RW_SIMULATION_HPP
#define _SIMPLE_RW_SIMULATION_HPP

#include <unistd.h>
#include <havoqgt/environment.hpp>

#include <havoqgt/temporal_random_walk_simulation.hpp>
#include <havoqgt/random_walkers/simple_random_walker.hpp>
#include <havoqgt/random_edge_container/interval_tree_random_edge_container.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/fifo_queue.hpp>
#include <iostream>
#include <unordered_set>
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
    class file_output_iterator{
    public:
      file_output_iterator(std::string filename,Graph *_g ) : g(_g) {
	ofs = new std::ofstream( filename, std::ofstream::trunc );
      }
      void operator()( const typename Graph::vertex_locator& vertex, const random_walker_t& random_walker) {
	(*ofs) << g->locator_to_label(vertex) << " "  
	       << g->locator_to_label(random_walker.start_from) << " " 
	       << random_walker << "\n";
      }

      void close() {
	ofs->close();
	delete(ofs);
      }

    private:
      std::ofstream* ofs;
      Graph* g;
    };

    /*********************************************************************************
     * Simulation & Input Iterator
     ********************************************************************************/
    template<typename Graph, typename EdgeMetaData>
    class simple_rw_simulation{
    public:
  
      using random_edge_container_t = random_edge_container<Graph, EdgeMetaData>;
      using random_walker_t = simple_random_walker<typename EdgeMetaData::value_type,
						   typename Graph::vertex_locator,
						   random_edge_container_t >;
      using output_iterator_t = file_output_iterator<Graph, random_walker_t>;
      typedef temporal_random_walk_simulation_visitor<Graph, EdgeMetaData, output_iterator_t, random_walker_t > visitor_t;
      typedef std::unordered_set<typename Graph::vertex_locator, typename Graph::vertex_locator::hasher> vertex_locator_set;
      simple_rw_simulation(Graph *_graph, EdgeMetaData* _edge_metadata) :
	graph(_graph), edge_metadata(_edge_metadata) { 
	initialize();
      }

      void initialize() {
	rwalker_dispatched = 0;
	is_complete = false;
	std::random_device r;
	mt.seed(r());
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

	using op_type = decltype(parse_to_long);
	time_start = reader.template get_value<uint64_t, op_type>( 's', parse_to_long).second;
	time_end   = reader.template get_value<uint64_t, op_type>( 'e', parse_to_long).second;
	time_offset = reader.template get_value<uint64_t, op_type>( 'o', parse_to_long).second;
	num_of_walkers = reader.template get_value<uint64_t, op_type>( 'n', parse_to_long).second;
	max_steps = reader.template get_value<uint64_t, op_type>( 'c', parse_to_long).second;

	output_file_prefix = reader.template get_value<std::string, decltype(identity)>('r', identity).second;
	source_label_file = reader.template get_value<std::string, decltype(identity)>('a', identity).second;
	target_label_file = reader.template get_value<std::string, decltype(identity)>('b', identity).second;

	//	std::cout << time_start<<" "<<time_end<<" "<<time_offset<<" "<<num_of_walkers<<" "<<max_steps
	//	  << " " << output_file_prefix << " " << source_label_file << " " << target_label_file << "\n";

	read_local_vertices_from_file( source_label_file, sources);
	read_local_vertices_from_file( target_label_file, targets);
	std::stringstream ss;
	ss << output_file_prefix << "_" 
	   << havoqgt_env()->world_comm().rank() << "_of_" 
	   << havoqgt_env()->world_comm().size();
	output_file_name = ss.str();
	return true; // false means missing parameters
      }

      void run() {
	output_iterator_t itr(output_file_name, graph);
	random_edge_container_t edge_container( graph, edge_metadata);

	visitor_t::set_output_iterator( &itr);
	visitor_t::set_edge_metadata( edge_metadata);
    
	random_walker_t::max_steps = max_steps;
	random_walker_t::random_edge_container = &edge_container;
	random_walker_t::local_targets = targets;

	typedef visitor_queue<visitor_t, havoqgt::detail::fifo_queue, Graph> visitor_queue_type;
	visitor_queue_type vq(graph);
    
	//prepare states before execution
	source_itr = sources.begin();

	vq.template init_visitor_traversal_experiment<simple_rw_simulation<Graph, EdgeMetaData>>(*this);

	itr.close();
      }

      std::pair<uint64_t, uint64_t> get_random_interval() {
	std::uniform_int_distribution<uint64_t> uniform_dist_start( time_start, time_end);
	uint64_t start = uniform_dist_start(mt);
	std::uniform_int_distribution<uint64_t> uniform_dist_end( start, time_end);
	uint64_t end = uniform_dist_end(mt);
	return std::make_pair( start, end);
      }

      std::pair<uint64_t, uint64_t> get_interval() {
	return std::make_pair( time_start, time_end);
      }

      template< typename VisitorQueueHandle >
      void process( VisitorQueueHandle* queue ) {
	while( source_itr == sources.end() &&
	       rwalker_dispatched < num_of_walkers ) {
	  source_itr = sources.begin();
	  rwalker_dispatched++;
	}
	if( rwalker_dispatched >= num_of_walkers ) {
	  is_complete = true;
	  return;
	}
	auto interval = get_interval();
	random_walker_t rw( rwalker_dispatched, *source_itr, interval, interval);    
	queue->queue_visitor( visitor_t(*source_itr, rw) );
	source_itr++;
      }

      bool empty() {
	return is_complete;
      }
    private:
      /**Graph Stuffs**/
      Graph* graph;
      EdgeMetaData* edge_metadata;

      /**Tools**/
      std::mt19937 mt;
      /**Experiment Related Constants**/
      uint64_t time_start;
      uint64_t time_end;
      uint64_t time_offset;
      uint64_t num_of_walkers;
      uint64_t max_steps;
      std::string output_file_prefix;
      std::string output_file_name;
      std::string source_label_file;
      std::string target_label_file;

      /**Experiment Progress Variables**/
      uint64_t rwalker_dispatched;
      bool is_complete;
      typename vertex_locator_set::iterator source_itr;
      vertex_locator_set sources;
      vertex_locator_set targets;
    };

  } }
#endif
