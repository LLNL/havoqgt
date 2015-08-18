#ifndef _CHILD_COUNT_VISITOR_HPP_INCLUDED
#define _CHILD_COUNT_VISITOR_HPP_INCLUDED

#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <vector>
#include <havoqgt/pairwise_sssp.hpp>

namespace havoqgt { namespace mpi {

    template<typename Graph, typename ChildCountData, typename ParentData>
    class child_count_visitor {
    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator edge_iterator;
      typedef ParentData value_type;
      child_count_visitor() { }

      child_count_visitor( vertex_locator _vertex)
	: vertex(_vertex), source( vertex_locator() ) { }
      
      child_count_visitor( vertex_locator _vertex, vertex_locator _source)
	: vertex(_vertex), source(_source) { }
       
      bool pre_visit() const {
	vertex_state_map<Graph, uint64_t>::traversal_count++;
	if( source != vertex_locator() ) {
	 
	  if( !(*child_count_data())[vertex].exists_with_state(source) ){
	    (*child_count_data())[vertex].set_state(source, 0);
	  }
	  (*child_count_data())[vertex].get_state(source)++;
	  return false;
	}
	return true;
      }

      template <typename VisitorQueueHandle>
      bool visit( Graph& g, VisitorQueueHandle queue) {
	//std::cout <<  "Visiting " << std::endl;
	if(source != vertex_locator() ) {
	  if((*parent_data())[vertex].exists(source)) {
	    queue->queue_visitor(child_count_visitor( (*parent_data())[vertex].get_parent(source), source));
	  }
	}
	else {
	  auto begin_itr = (*parent_data())[vertex].source_begin();
	  auto end_itr = (*parent_data())[vertex].source_end();
	  
	  for(auto itr = begin_itr; itr != end_itr; itr++) {
	    if(itr->second == vertex_locator()) continue; // for the vertices where source = current vertex
	    queue->queue_visitor(child_count_visitor( itr->second, g.label_to_locator(itr->first)));
	  }
	}
	return true;
      }

      friend bool operator>(const child_count_visitor& a,
			    const child_count_visitor& b) {
	if( b.vertex == a.vertex ) return false;
	return b.vertex < a.vertex;
      }

      static ChildCountData*& child_count_data() {
	static ChildCountData* data;
	return data;
      }

      static void set_child_count_data(ChildCountData* data) {
	child_count_data() = data;
      }

      static ParentData*& parent_data() {
	static ParentData* data;
	return data;
      }

      static void set_parent_data(ParentData* data) {
	parent_data() = data;
      }

      vertex_locator vertex;
      vertex_locator source;
    } __attribute__ ((packed));


    template<typename TGraph, typename ChildCountData, typename ParentData>
    void count_child(TGraph* g, ChildCountData *child_count_data, ParentData *parent_data, std::vector<typename TGraph::vertex_locator>& sources) {
      typedef child_count_visitor<TGraph, ChildCountData, ParentData> visitor_type;

      visitor_type::set_child_count_data(child_count_data);
      visitor_type::set_parent_data(parent_data);
      typedef visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue, TGraph> visitor_queue_type;
      visitor_queue_type vq(g);
      //std::cout << "Initiating the visitor traversal child " << std::endl;
      vq.init_visitor_traversal();
    }
    
  }
}


#endif  //  _CHILD_COUNT_VISITOR_HPP_INCLUDED
