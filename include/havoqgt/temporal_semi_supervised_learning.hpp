#ifndef _TEMPORAL_RANDOM_WALK_SIMULATION_HPP
#define _TEMPORAL_RANDOM_WALK_SIMULATION_HPP

#define SEG_TREE 0


#include <chrono>
#include <random>
#include <tuple>
#include <sstream>
#include <iostream>
//#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/detail/fifo_queue.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/segment_tree.hpp>
#include <unordered_set>

namespace havoqgt { namespace mpi {


    using clock_t = std::chrono::high_resolution_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    /**********************************************************************
    Random Number Generator
    ***********************************************************************/
    struct random_number_generator {
    public:
      std::mt19937 en;
      
      static random_number_generator& get() {
	if(val == 0) {
	  _rng = random_number_generator();
	  val = 1;
	}
	return _rng;
      }

      template<typename distribution>
      typename distribution::result_type operator()(distribution& dist){
	return dist(en);
      }

    private:
      static random_number_generator _rng;
      static int val;
      random_number_generator(){ 
	std::random_device r;
	en.seed(r());
      }

    };

    typedef random_number_generator rng;
    rng rng::_rng;
    int rng::val = 0;

    /**********************************************************************
    Random Walker
    ***********************************************************************/
    template <typename Graph, typename EdgeMetaData, typename SegmentTreeData>
    class random_walker {
    public:
      
      struct vertex_hasher
      {
	std::size_t operator()(const typename Graph::vertex_locator& other) const {
	  return other.hash();
	}
      };

      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      typedef typename std::unordered_set<vertex_locator, vertex_hasher> set;

      template<std::size_t sz>
      class memory{
      public: 
	memory(): cur_size(0) {} 

	//TODO
	void push_back(vertex_locator vertex) {
	}

	std::size_t size() { return cur_size; }
	std::size_t cur_size;
	vertex_locator remembered_vertics[sz + 1];
      };
      
      random_walker() : id(0) {}

      random_walker( uint64_t _id, vertex_locator _start_from, uint64_t _cost
		     ,uint64_t _cur_time
		     ,uint64_t _started_at
		     , uint64_t _steps
		     , uint64_t _restart_count
		     )
	: id(_id), start_from(_start_from), cost(_cost), cur_time(_cur_time), started_at(_started_at),steps(_steps), target_hit(vertex_locator()), was_completed(false), restart_count(_restart_count) { }

      random_walker( uint64_t _id, vertex_locator _start_from, uint64_t _cost
		     ,uint64_t _cur_time
		     ,uint64_t _started_at
		     , uint64_t _steps
		     , uint64_t _restart_count
		     , bool _completed
		     , vertex_locator _target_hit
		     )
	: id(_id), start_from(_start_from), cost(_cost), cur_time(_cur_time), started_at(_started_at),steps(_steps), target_hit(_target_hit), was_completed(_completed), restart_count(_restart_count) { }

      bool is_complete(vertex_locator vertex) const {
	bool target_reached = ( targets.find(vertex) != targets.end());
	return target_reached;
      }

      random_walker register_complete(vertex_locator vertex) const{
	return random_walker(id, start_from, cost, cur_time, started_at, steps, restart_count,  true, vertex);
      }

      bool in_infinite_path(vertex_locator vertex) const {
	bool is_infinite = (steps >= random_walker::max_steps); 
	return is_infinite;
      }

      random_walker restart() const {
	return random_walker(id, start_from, 0, cur_time, started_at, 0, restart_count + 1); 
      }

      static void set_targets(std::unordered_set<vertex_locator>& _targets) {
	for( auto itr = _targets.begin(); itr != _targets.end();itr++) {
	  targets.insert(*itr);
	}
      }

      static void add_target(vertex_locator target) {
	targets.insert(target);
      }
      

      //Function returning the next Visitor
      std::tuple<bool, vertex_locator, random_walker> next(Graph& g, EdgeMetaData*& edges_metadata
							   , SegmentTreeData*& segment_tree_data, vertex_locator cur_vertex ) const{
#if SEG_TREE == 1
	std::cout << "Not here totally" << std::endl;
	auto &tree = (*segment_tree_data)[cur_vertex];
	uint64_t size = 0;
	tree.query( cur_time, size);
#else
	uint64_t size = g.edges_end(cur_vertex) - g.edges_begin(cur_vertex);
#endif
	if( size == 0 ) {
	  if( should_restart() ){
	    std::cout << "RW restarting" << std::endl;
	    random_walker restart_rw(id, start_from, 0, cur_time, started_at, 0 , restart_count + 1);
	    return std::make_tuple(true, start_from, restart_rw);
	  }
	  else return std::make_tuple( false, vertex_locator(),  random_walker() );

	}

	std::uniform_int_distribution<std::size_t> uniform_dist( 0, size - 1);

	std::size_t index;
	{
	  index = (rng::get())(uniform_dist);
	}
#if SEG_TREE == 1
	eitr_type itr = g.edges_end(cur_vertex);
	tree.query_ith( cur_time, index, itr);
#else
	eitr_type itr = g.edges_begin(cur_vertex);
	itr+=index;
#endif
	metadata_type& metadata = (*edges_metadata)[itr];
	random_walker next_rw(id, start_from, cost + (metadata.redirect ? 0 : 1 ), cur_time + 0, started_at, steps + 1, restart_count);
	return std::make_tuple( true , itr.target(), next_rw );
      }

      friend std::ostream& operator<<(std::ostream& o, const random_walker& rw) {
	return o << rw.started_at
		 << ";" << rw.steps
		 << ";" << rw.cost;
      }

      bool returned_to_mother(vertex_locator vertex) const {
	return was_completed && ( vertex == start_from );
      }

      bool should_restart() const {
	return restart_count < restart_max;
      }

      uint64_t id;
      vertex_locator start_from;
      static std::unordered_set<vertex_locator, vertex_hasher> targets;
      bool was_completed;
      vertex_locator target_hit;
      uint64_t cost;
      uint64_t steps;
      uint64_t cur_time;
      uint64_t started_at;
      memory<0> no_memory;
      uint64_t restart_count;
      static uint64_t max_steps;
      static uint64_t restart_max;
    };
    
    template<typename Graph, typename EdgeMetaData, typename SegmentTreeData>
    typename random_walker<Graph, EdgeMetaData, SegmentTreeData>::set random_walker<Graph, EdgeMetaData, SegmentTreeData>::targets = typename random_walker::set();
   
    template<typename Graph, typename EdgeMetaData, typename SegmentTreeData>
    uint64_t random_walker<Graph, EdgeMetaData, SegmentTreeData>::max_steps = 0;

    template<typename Graph, typename EdgeMetaData, typename SegmentTreeData>
    uint64_t random_walker<Graph, EdgeMetaData, SegmentTreeData>::restart_max = 0;

    /**********************************************************************
    Random Walker Visitor
    ***********************************************************************/
    template <typename Graph, typename EdgeMetaData, typename SegmentTreeData>
    class temporal_random_walk_simulation_visitor {

      enum TUPLE_INDEX {
	HAS_NEIGHBOR=0, NEIGHBOR, RANDOM_WALKER
      };

    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      using random_walker_t = random_walker<Graph, EdgeMetaData, SegmentTreeData>;

      temporal_random_walk_simulation_visitor() {}

      temporal_random_walk_simulation_visitor( vertex_locator _vertex)
	: vertex(_vertex) { }

      temporal_random_walk_simulation_visitor( vertex_locator _vertex, random_walker_t _rwalker)
	: vertex(_vertex), rwalker(_rwalker) { }
      
      template<typename VisitorQueueHandle>
      bool init_visit(Graph& g, VisitorQueueHandle vis_queue) const {
	uint64_t i = 1;
	for( auto x : (*start_time_steps()) ) {
	  random_walker_t _rwalker(i, vertex, 0, x, x, 0 , 0);
	  temporal_random_walk_simulation_visitor v( vertex, _rwalker);
	  vis_queue->queue_visitor(v);
	  i++;
	} 
      }

      bool pre_visit() const {
	return true;
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
	if( rwalker.returned_to_mother(vertex)) {
	  (*out_stream) << ( (rwalker.was_completed==true)?g.locator_to_label(rwalker.target_hit):0 ) << ";"
		        << g.locator_to_label(rwalker.start_from) << ";"
		        << rwalker << ";COMP_RET" << std::endl;
	}else {

	bool is_completed = rwalker.is_complete(vertex);
	if( is_completed ) {
	  /*	  temporal_random_walk_simulation_visitor neighbor( rwalker.start_from, rwalker.register_complete(vertex));
		  vis_queue->queue_visitor(neighbor); */
	  (*out_stream) << g.locator_to_label( vertex ) << ";"
			<< g.locator_to_label( rwalker.start_from) << ";"
			<< rwalker << ";COMP" << std::endl;

	} else if( rwalker.in_infinite_path(vertex) ){
	  if( rwalker.should_restart() ) {
	    temporal_random_walk_simulation_visitor neighbor( rwalker.start_from, rwalker.restart());
	    vis_queue->queue_visitor(neighbor);
	  } else {
	    (*out_stream) << "-1;"
		        << g.locator_to_label(rwalker.start_from) << ";"
		        << rwalker << ";DEL" << std::endl;
	  }
	}
	else {
	  std::tuple< bool, vertex_locator, random_walker_t> next = rwalker.next( g, edges_metadata(), segment_tree_data(), vertex);
	  if( std::get<HAS_NEIGHBOR>(next) ) {
	    temporal_random_walk_simulation_visitor neighbor( std::get<NEIGHBOR>(next), std::get<RANDOM_WALKER>(next));
	    vis_queue->queue_visitor(neighbor);
	  } 
	}
	}
	return true;
      }
      
      friend inline bool operator>(const temporal_random_walk_simulation_visitor& v1
				   , const temporal_random_walk_simulation_visitor& v2 ) {
	return v1.rwalker.steps > v2.rwalker.steps;
      }

      friend inline bool operator<(const temporal_random_walk_simulation_visitor& v1
				   , const temporal_random_walk_simulation_visitor& v2 ) {
	return v1.rwalker.steps < v2.rwalker.steps;
      }
      
      static void set_edge_metadata(EdgeMetaData* data) { edges_metadata() = data; }
      
      static EdgeMetaData*& edges_metadata() {
	static EdgeMetaData* data;
	return data;
      }

      static void set_segment_tree_data( SegmentTreeData* data) { segment_tree_data() = data; }
      
      static SegmentTreeData*& segment_tree_data() {
	static SegmentTreeData* data;
	return data;
      }

      static void set_start_time_steps(std::vector<uint64_t>* data) { start_time_steps() = data; }

      static std::vector<uint64_t>*& start_time_steps() {
	static std::vector<uint64_t>* data;
	return data;
     }

      static std::ofstream *out_stream;
      vertex_locator vertex;
      random_walker_t rwalker;
    };// __attribute__ ((packed));

    template<typename TGraph, typename EdgeMetaData, typename SegmentTreeData>
    std::ofstream *temporal_random_walk_simulation_visitor<TGraph, EdgeMetaData, SegmentTreeData>::out_stream;
    
    template< typename TGraph, typename EdgeMetaData, typename SegmentTreeData>
    void temporal_random_walk_simulation( TGraph* g,
					  EdgeMetaData* edge_metadata,
					  SegmentTreeData* segment_tree_data,
					  std::vector<typename TGraph::vertex_locator> sources,
					  std::vector<typename TGraph::vertex_locator> targets,
					  uint64_t num_of_walkers,
					  std::vector<uint64_t>* start_times,
					  uint64_t max_steps
					  , uint64_t restart_max,
					  std::string& out_filename) {
      typedef temporal_random_walk_simulation_visitor< TGraph, EdgeMetaData, SegmentTreeData> visitor_type;
      visitor_type::set_edge_metadata( edge_metadata);
      visitor_type::set_segment_tree_data( segment_tree_data);
      visitor_type::set_start_time_steps( start_times );
      visitor_type::random_walker_t::max_steps = max_steps;
      visitor_type::random_walker_t::restart_max = restart_max;
      visitor_type::random_walker_t::targets = typename visitor_type::random_walker_t::set();
      
      for( auto x:targets ) {
	visitor_type::random_walker_t::add_target(x);
      }
      
      std::ofstream out_file( out_filename, std::ofstream::trunc );
      visitor_type::out_stream = &out_file;
      
      typedef visitor_queue< visitor_type, detail::fifo_queue, TGraph> visitor_queue_type;
      visitor_queue_type vq(g);


	std::vector<visitor_type> visitor_list;

	uint64_t i = 1;
	//for( uint64_t c = 0; c < num_of_walkers; c++ ) {
	  //	  for( auto x : (*start_times) ) {
	  //for( auto vertex: sources ) {
	  //  typename visitor_type::random_walker_t _rwalker(i, vertex, 0, x, x, 0 , 0);
	  //  visitor_list.push_back( visitor_type( vertex, _rwalker) );	      
	  //  i++;
	  //}
	  // }
	  //}
	//      vq.init_visitor_traversal(visitor_list);
	vq.init_visitor_traversal_new(num_of_walkers);

	out_file.close();
    }

  }/*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif
