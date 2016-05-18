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

    template <typename Graph, typename OutputIterator, typename RandomWalker>
    class random_walk_simulation_visitor {

      enum TUPLE_INDEX {
	STATE=0, NEIGHBOR, RANDOM_WALKER
      };

    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;

      using random_walker_t = RandomWalker;

      random_walk_simulation_visitor() {}

      random_walk_simulation_visitor( vertex_locator _vertex)
	: vertex(_vertex) { }

      random_walk_simulation_visitor( vertex_locator _vertex, random_walker_t _rwalker, int _state = -1)
	: vertex(_vertex), rwalker(_rwalker), state(_state) { }
      
      template<typename VisitorQueueHandle>
      bool init_visit(Graph& g, VisitorQueueHandle vis_queue) const {
	
      }

      bool pre_visit() const {
	bool is_complete = rwalker.is_complete(vertex);
	if(state == 0 || state == 2)
	  (*output_iterator())(vertex, rwalker); // can do by state as well
	//(*output_iterator()).completed++;
	return !is_complete;
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
	std::tuple< int , vertex_locator, random_walker_t> next = rwalker.next(&g, vertex);
	random_walk_simulation_visitor neighbor( std::get<NEIGHBOR>(next), std::get<RANDOM_WALKER>(next), std::get<STATE>(next));
	  vis_queue->queue_visitor(neighbor);

	return true;
      }
      
      friend inline bool operator>(const random_walk_simulation_visitor& v1
				   , const random_walk_simulation_visitor& v2 ) {
	return true;
      }

      friend inline bool operator<(const random_walk_simulation_visitor& v1
				   , const random_walk_simulation_visitor& v2 ) {
	return false;
      }
      
      static void set_output_iterator(OutputIterator* data) { output_iterator() = data; }

      static OutputIterator*& output_iterator() {
	static OutputIterator* data;
	return data;
      }

      int state;
      vertex_locator vertex;
      random_walker_t rwalker;
    };// __attribute__ ((packed));

  }/*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif
