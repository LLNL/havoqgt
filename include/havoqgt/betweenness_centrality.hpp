#ifndef _BETWEENNESS_CENTRALITY_HPP_INCLUDED
#define _BETWEENNESS_CENTRALITY_HPP_INCLUDED

#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/pairwise_sssp.hpp>

namespace havoqgt { namespace mpi {

    template<typename Graph, typename ChildCountData, typename ParentData, typename ShortestPathData>
    class betweenness_centrality_visitor {
    public:
      typedef typename Graph::vertex_locator vertex_locator;

      betweenness_centrality_visitor() {}

      betweenness_centrality_visitor( vertex_locator _vertex)
	: vertex(_vertex), source(vertex_locator()), count(0) {}
      
      betweenness_centrality_visitor( vertex_locator _vertex, vertex_locator _source, uint64_t _count )
	: vertex(_vertex), source(_source), count(_count) { }

      bool pre_visit() const {
	vertex_state_map<Graph, uint64_t>::traversal_count++;
	if( source != vertex_locator() ) {
	  if ( (*child_count_data())[vertex].get_state(source) == 0 ) {
	    return false;
	  }
	  (*shortest_path_data())[vertex].get_state(source) += count;
	  (*child_count_data())[vertex].get_state(source)--;
	  return ((*child_count_data())[vertex].get_state(source) == 0);
        } else {
	
	  bool should_visit = false;
	  auto begin_itr = (*parent_data())[vertex].source_begin();
	  auto end_itr = (*parent_data())[vertex].source_end();
	  
	  for(auto source_itr = begin_itr; source_itr != end_itr ; source_itr++) {
	    vertex_locator _source = vertex_state_map<Graph, uint64_t>::graph()->label_to_locator(source_itr->first);

	    if((*child_count_data())[vertex].get_state(_source) == 0 //I dont have child
	       && (*shortest_path_data())[vertex].get_state(_source) == 0 && //I dont have child with this source from very beginning
	       (*parent_data())[vertex].get_parent(_source) != vertex_locator() ) { //I do have a path from this source to myself
	      (*shortest_path_data())[vertex].get_state(_source) = 1;
	      should_visit = true;
	    } 	      
	  }
	  return should_visit;
	}
      }

      template<typename VisitorQueueHandle>
      bool visit(Graph& g, VisitorQueueHandle queue) const {
	if(source == vertex_locator() ) {
	  auto begin_itr = (*parent_data())[vertex].source_begin();
	  auto end_itr = (*parent_data())[vertex].source_end();
	  
	  for( auto source_itr = begin_itr; source_itr!=end_itr; source_itr++){
	    vertex_locator _source = g.label_to_locator(source_itr->first);
	    if ( (*child_count_data())[vertex].get_state(_source) == 0
		 && (*parent_data())[vertex].get_parent(_source) != vertex_locator() ) {
	      betweenness_centrality_visitor v((*parent_data())[vertex].get_parent(_source)
					    ,_source
					   ,(*shortest_path_data())[vertex].get_state(_source) );
	      queue->queue_visitor(v);
	  }
	  }
	}else {
	  if((*parent_data())[vertex].get_parent(source) != vertex_locator()) {
	    betweenness_centrality_visitor v((*parent_data())[vertex].get_parent(source),
					   source,
					   (*shortest_path_data())[vertex].get_state(source) );
	    queue->queue_visitor(v);
	  }
	}
      }

      friend bool operator>( const betweenness_centrality_visitor& a,
			     const betweenness_centrality_visitor& b ) {
	
	return (*child_count_data())[a.vertex].get_state(a.source) >
	  (*child_count_data())[b.vertex].get_state(b.source);
      }

      static ChildCountData*& child_count_data() {
	static ChildCountData* data;
	return data;
      }

      static void set_child_count_data(ChildCountData* data) { child_count_data() = data; }

      static ParentData*& parent_data() {
	static ParentData* data;
	return data;
      }
	  
      static void set_parent_data(ParentData* data) { parent_data() = data; }

      static ShortestPathData*& shortest_path_data() {
	static ShortestPathData* data;
	return data;
      }

      static void set_shortest_path_data(ShortestPathData* data) {
	shortest_path_data() = data;
      }

      static void set_sources(std::vector<vertex_locator>* _sources) {
	sources = _sources;
      }

      static std::vector<vertex_locator>* sources;
      vertex_locator vertex;
      vertex_locator source;
      int64_t count;
    } __attribute__ ((packed));

    template<typename TGraph, typename ChildCountData, typename ParentData, typename ShortestPathData>
    std::vector<typename TGraph::vertex_locator>* betweenness_centrality_visitor<TGraph, ChildCountData, ParentData, ShortestPathData>::sources;
    
    template<typename TGraph, typename ChildCountData, typename ParentData, typename ShortestPathData>
    void betweenness_centrality(TGraph* g,
				ChildCountData* child_count_data,
				ParentData* parent_data,
				ShortestPathData* shortest_path_data,
				std::vector<typename TGraph::vertex_locator>& sources) {
      typedef betweenness_centrality_visitor<TGraph, ChildCountData, ParentData, ShortestPathData> visitor_type;

      visitor_type::set_child_count_data(child_count_data);
      visitor_type::set_parent_data(parent_data);
      visitor_type::set_shortest_path_data(shortest_path_data);
      visitor_type::set_sources(&sources);

      typedef visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue, TGraph> visitor_queue_type;

      visitor_queue_type vq(g);
      vq.init_visitor_traversal();
    }
  }
}

#endif
