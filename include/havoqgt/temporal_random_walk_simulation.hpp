#ifndef _TEMPORAL_RANDOM_WALK_SIMULATION_HPP
#define _TEMPORAL_RANDOM_WALK_SIMULATION_HPP

#include <chrono>
#include <random>
#include <tuple>
#include <sstream>
#include <iostream>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>

namespace havoqgt { namespace mpi {

    using clock_t = std::chrono::high_resolution_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    /**********************************************************************
    Random Number Generator
    ***********************************************************************/
    struct random_number_generator {
    public:
      std::mt19937 en;
      random_number_generator(){ 
	std::random_device r;
	en.seed(r());
      }
      
      template<typename distribution>
      typename distribution::result_type operator()(distribution& dist){
	return dist(en);
      }
    };

    typedef random_number_generator rng;


    /**********************************************************************
    Random Walker
    ***********************************************************************/
    template <typename Graph, typename EdgeMetaData>
    class random_walker {
    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      
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
		     ,time_point_t _cur_time
		     ,time_point_t _started_at
		     , uint64_t _steps
		     ,vertex_locator _target = vertex_locator() )
	: id(_id), start_from(_start_from), cost(_cost), cur_time(_cur_time), started_at(_started_at),steps(_steps), target(_target) { }

      bool is_complete(vertex_locator vertex) const {
	return vertex == target || steps >= 100;
      }

      void set_target(vertex_locator vertex) {
	target = vertex;
      }
      

      //Function returning the next Visitor
      std::tuple<bool, vertex_locator, random_walker> next(Graph& g, EdgeMetaData*& edges_metadata,
					  vertex_locator cur_vertex ) const{

	std::vector<vertex_locator> adjacents;
	std::vector<uint32_t> cost_vec;

	for(eitr_type itr = g.edges_begin(cur_vertex); itr != g.edges_end(cur_vertex); ++itr ) {
	  metadata_type& metadata = (*edges_metadata)[itr];
	  if( metadata.start_time() <= cur_time && (metadata.end_time().time_since_epoch() == std::chrono::seconds(0) ||
						    metadata.end_time() >= cur_time )) {
	    adjacents.push_back(itr.target());
	    cost_vec.push_back(metadata.redirect?0:1);
	  }
	}
	if( adjacents.size() == 0 ) return std::make_tuple( false, vertex_locator(),  random_walker() );

	std::uniform_int_distribution<std::size_t> uniform_dist( 0, adjacents.size() - 1);
	// for memory -- uniform_dist(0, adjacents.size() + memory.size() - 1 );
	std::size_t index;
	{
	  rng _rng;
	  index = _rng(uniform_dist);
	}
	random_walker next_rw(id, start_from, cost + cost_vec[index], cur_time + std::chrono::seconds(3600), started_at, steps + 1, target);
	return std::make_tuple( true , adjacents[index], next_rw );
      }

      friend std::ostream& operator<<(std::ostream& o, const random_walker& rw) {
	return o << rw.id
		 << ";" << std::chrono::duration_cast<std::chrono::seconds>(rw.started_at.time_since_epoch()).count()
		 << ";" << rw.cost;
      }

      uint64_t id;
      vertex_locator start_from;
      vertex_locator target;
      uint64_t cost;
      uint64_t steps;
      time_point_t cur_time;
      time_point_t started_at;
      memory<0> no_memory;
    };

    /**********************************************************************
    Random Walker Visitor
    ***********************************************************************/
    template <typename Graph, typename EdgeMetaData, typename VertexData>
    class temporal_random_walk_simulation_visitor {

      enum TUPLE_INDEX {
	HAS_NEIGHBOR=0, NEIGHBOR, RANDOM_WALKER
      };

    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      using random_walker_t = random_walker<Graph, EdgeMetaData>;

      temporal_random_walk_simulation_visitor() {}

      temporal_random_walk_simulation_visitor( vertex_locator _vertex)
	: vertex(_vertex) { }

      temporal_random_walk_simulation_visitor( vertex_locator _vertex, random_walker_t _rwalker)
	: vertex(_vertex), rwalker(_rwalker) { }
      
      template<typename VisitorQueueHandle>
      bool init_visit(Graph& g, VisitorQueueHandle vis_queue) const {
	uint64_t i = 1;
	for( auto x : (*start_time_steps()) ) {
	  random_walker_t _rwalker(i, vertex, 0, x, x, 0);
	  temporal_random_walk_simulation_visitor v( vertex, _rwalker);
	  vis_queue->queue_visitor(v);
	  i++;
	}
      }

      bool pre_visit() const {
	bool is_complete = rwalker.is_complete(vertex);
	if( rwalker.steps != 0) { //not for the first time
	  (*vertex_data())[vertex].at(rwalker.id).sum += rwalker.cost;
	  (*vertex_data())[vertex].at(rwalker.id).min = std::min( (*vertex_data())[vertex].at(rwalker.id).min, rwalker.cost );
	  (*vertex_data())[vertex].at(rwalker.id).max = std::max( (*vertex_data())[vertex].at(rwalker.id).max, rwalker.cost );
	  (*vertex_data())[vertex].at(rwalker.id).count++;
	}
	return !is_complete;
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
	std::tuple< bool, vertex_locator, random_walker_t> next = rwalker.next( g, edges_metadata(), vertex);
	if( std::get<HAS_NEIGHBOR>(next) ) {
	  temporal_random_walk_simulation_visitor neighbor( std::get<NEIGHBOR>(next), std::get<RANDOM_WALKER>(next));
	  vis_queue->queue_visitor(neighbor);
	}
	return true;
      }
      
      friend inline bool operator>(const temporal_random_walk_simulation_visitor& v1
				   , const temporal_random_walk_simulation_visitor& v2 ) {
	return true;
      }

      friend inline bool operator<(const temporal_random_walk_simulation_visitor& v1
				   , const temporal_random_walk_simulation_visitor& v2 ) {
	return false;
      }
      
      static void set_edge_metadata(EdgeMetaData* data) { edges_metadata() = data; }
      
      static EdgeMetaData*& edges_metadata() {
	static EdgeMetaData* data;
	return data;
      }

      static void set_start_time_steps(std::vector<time_point_t>* data) { start_time_steps() = data; }

      static std::vector<time_point_t>*& start_time_steps() {
	static std::vector<time_point_t>* data;
	return data;
     }

      static void set_vertex_data(VertexData* data) { vertex_data() = data; }
      
      static VertexData*& vertex_data() {
	static VertexData* data;
	return data;
      }

      vertex_locator vertex;
      random_walker_t rwalker;
    } __attribute__ ((packed));
    
    template< typename TGraph, typename EdgeMetaData, typename VertexData>
    void temporal_random_walk_simulation( TGraph* g,
					  EdgeMetaData* edge_metadata,
					  VertexData* vertex_data,
					  typename TGraph::vertex_locator source,
					  typename TGraph::vertex_locator target,
					  uint64_t num_of_walkers,
					  std::vector<time_point_t>* start_times) {
      typedef temporal_random_walk_simulation_visitor< TGraph, EdgeMetaData, VertexData> visitor_type;
      visitor_type::set_edge_metadata( edge_metadata);
      visitor_type::set_vertex_data( vertex_data );
      visitor_type::set_start_time_steps( start_times );

      typedef visitor_queue< visitor_type, detail::visitor_priority_queue, TGraph> visitor_queue_type;
      visitor_queue_type vq(g);

      vq.init_visitor_traversal_new();
    }

  }/*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif
