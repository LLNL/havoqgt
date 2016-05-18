#ifndef _TEMPORAL_RANDOM_WALK_SIMULATION_HPP
#define _TEMPORAL_RANDOM_WALK_SIMULATION_HPP

#include <chrono>
#include <random>
#include <tuple>
#include <sstream>
#include <iostream>

namespace havoqgt { namespace mpi {
    /*
    using clock_t = std::chrono::high_resolution_clock;
    using time_point_t = std::chrono::time_point<clock_t>;
    */

    template <typename Graph, typename EdgeMetaData, typename OutputIterator, typename RandomWalker>
    class temporal_random_walk_simulation_visitor {

      enum TUPLE_INDEX {
	HAS_NEIGHBOR=0, NEIGHBOR, RANDOM_WALKER
      };

    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      using random_walker_t = RandomWalker;

      temporal_random_walk_simulation_visitor() {}

      temporal_random_walk_simulation_visitor( vertex_locator _vertex)
	: vertex(_vertex) { }

      temporal_random_walk_simulation_visitor( vertex_locator _vertex, random_walker_t _rwalker)
	: vertex(_vertex), rwalker(_rwalker) { }
      
      template<typename VisitorQueueHandle>
      bool init_visit(Graph& g, VisitorQueueHandle vis_queue) const {
	
      }

      bool pre_visit() const {
	bool is_complete = rwalker.is_complete(vertex);
	if( rwalker.steps != 0 )
	  (*output_iterator())(vertex, rwalker); // can do by state as well
	return !is_complete;
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
	std::tuple< bool, vertex_locator, random_walker_t> next = rwalker.next(vertex);
	if( std::get<HAS_NEIGHBOR>(next) ) {
	  temporal_random_walk_simulation_visitor neighbor( std::get<NEIGHBOR>(next), std::get<RANDOM_WALKER>(next));
	  vis_queue->queue_visitor(neighbor);
	}
	/*else {
	  (*output_iterator())( vertex, rwalker);
	  }*/
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
      
      static void set_output_iterator(OutputIterator* data) { output_iterator() = data; }

      static OutputIterator*& output_iterator() {
	static OutputIterator* data;
	return data;
      }

      vertex_locator vertex;
      random_walker_t rwalker;
    };// __attribute__ ((packed));
    /*    
    template< typename TGraph, typename EdgeMetaData, typename VertexData>
    void temporal_random_walk_simulation( TGraph* g,
					  EdgeMetaData* edge_metadata,
					  VertexData* vertex_data,
					  typename TGraph::vertex_locator source,
					  typename TGraph::vertex_locator target,
					  uint64_t num_of_walkers,
					  std::vector<time_point_t>* start_times,
					  uint64_t max_steps) {
      typedef temporal_random_walk_simulation_visitor< TGraph, EdgeMetaData, VertexData> visitor_type;
      visitor_type::set_edge_metadata( edge_metadata);
      visitor_type::set_vertex_data( vertex_data );
      visitor_type::set_start_time_steps( start_times );
      visitor_type::random_walker_t::max_steps = max_steps;

      typedef visitor_queue< visitor_type, detail::visitor_priority_queue, TGraph> visitor_queue_type;
      visitor_queue_type vq(g);

      vq.init_visitor_traversal_new();
    }
    */
  }/*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif
