#ifndef _TEMPORAL_RW_WITH_PATH_HISTORY_SIMULATION_HPP
#define _TEMPORAL_RW_WITH_PATH_HISTORY_SIMULATION_HPP

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

    template <typename Graph, typename EdgeMetaData, typename SimulationState, typename RandomWalker>
    class temporal_random_walk_simulation_visitor {

      enum TUPLE_INDEX {
	VALID_VALUES = 0, TARGET_VERTEX, CHOSEN_EDGE, RANDOM_WALKER
      };

    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
      typedef typename EdgeMetaData::value_type metadata_type;
      using random_walker_t = RandomWalker;

      temporal_random_walk_simulation_visitor() {}

      temporal_random_walk_simulation_visitor( vertex_locator _vertex)
	: vertex(_vertex) { }

      temporal_random_walk_simulation_visitor( vertex_locator _vertex, random_walker_t _rwalker
					       , vertex_locator _parent = vertex_locator(), bool _reverse = false)
	: vertex(_vertex), rwalker(_rwalker), parent(_parent), reverse(_reverse) { }
      
      template<typename VisitorQueueHandle>
      bool init_visit(Graph& g, VisitorQueueHandle vis_queue) const {
	
      }

      bool pre_visit() const {
	if(reverse) { return true; }
	else{
	  (*simulation_state()).register_parent( vertex, parent, rwalker);
	}
	return !rwalker.infinite_path(vertex);
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
	// if reached target  ( then send reverse as true )
	if( rwalker.success(vertex) && !reverse) {
	  if( parent == vertex_locator()){
	    std::cout << "Parent is default at" << g.locator_to_label( vertex) <<"\n";
	    std::cout << "Random walker started at" << g.locator_to_label( rwalker.start_from ) << "\n";
	    exit(1);
	  }
	  temporal_random_walk_simulation_visitor neighbor( parent, rwalker, vertex_locator(), true);
	  vis_queue->queue_visitor(neighbor); 
	  return true;
	}
	
	if(!reverse) {
	  std::tuple< uint32_t, vertex_locator, eitr_type, random_walker_t> next = rwalker.next(vertex);
	  if( std::get<VALID_VALUES>(next) != 0 /*HAS TARGET*/ ) {
	    //register the chosen edge
	    if( std::get<VALID_VALUES>(next) == 3 /* Edge Chosen */ ) {
	      //add vertex and edge
	      (*simulation_state()).register_edge( std::get< CHOSEN_EDGE >(next), rwalker );
	      temporal_random_walk_simulation_visitor neighbor( std::get<TARGET_VERTEX>(next), std::get<RANDOM_WALKER>(next), vertex);
	      vis_queue->queue_visitor(neighbor);
	    } else {
	      temporal_random_walk_simulation_visitor neighbor( std::get<TARGET_VERTEX>(next), std::get<RANDOM_WALKER>(next));
	      vis_queue->queue_visitor(neighbor);
	    }
	  }
	  /* else {
	    (*output_iterator())( vertex, rwalker);
	  }*/
	} else {
	  // Read the parent for this data for traversal ( current vertex --> rw_id --> steps --> parent )
	  // increment the edge chosen at this stage by 1 / path_length
	  (*simulation_state()).increment_edge(vertex, rwalker);
	  vertex_locator _parent = (*simulation_state()).get_parent(vertex, rwalker);
	  if(_parent != vertex_locator() && vertex != rwalker.start_from) {
	    temporal_random_walk_simulation_visitor neighbor(_parent, rwalker, vertex_locator(), true);
	    vis_queue->queue_visitor(neighbor);
	  } else {
	    //output the start time for the randomwalker id
	    (*simulation_state()).success_rw_id_time_map.insert( std::make_pair( rwalker.id, rwalker.started_at.first ) );
	  }
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
      
      static void set_simulation_state(SimulationState* data) { simulation_state() = data; }

      static SimulationState*& simulation_state() {
	static SimulationState* data;
	return data;
      }

      bool reverse;
      vertex_locator vertex;
      vertex_locator parent;
      random_walker_t rwalker;
    };// __attribute__ ((packed));
  }/*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif
