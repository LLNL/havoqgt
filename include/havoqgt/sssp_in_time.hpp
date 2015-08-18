#ifndef _SSSP_IN_TIME_HPP_INCLUDED
#define _SSSP_IN_TIME_HPP_INCLUDED

#include <havoqgt/nano_time.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>

namespace havoqgt { namespace mpi {

template <typename Graph, typename EdgeMetaData, typename ArrivalData, typename ParentData>
class sssp_in_time_visitor {

public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  typedef typename EdgeMetaData::value_type metadata_type;
  sssp_in_time_visitor() { }

  sssp_in_time_visitor(vertex_locator _vertex)
    : vertex(_vertex), arrival_time(nano_time()), parent(vertex_locator()) { }

  sssp_in_time_visitor(vertex_locator _vertex, nano_time _arrival_time, vertex_locator _parent)
    : vertex(_vertex), arrival_time(_arrival_time), parent(_parent) { }

  bool pre_visit() const {
    return true;
    //return ( arrival_time < (*arrival_data())[vertex] );
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    
    if(arrival_time < (*arrival_data())[vertex]) {
      (*arrival_data())[vertex] = arrival_time;
      (*parent_data())[vertex] = parent;
    }
    
    for(eitr_type itr = g.edges_begin(vertex); itr != g.edges_end(vertex); ++itr ) {
      metadata_type meta_data = (*edges_metadata())[itr];
      
      if( meta_data.end_time() < arrival_time // this edge has past
	  || meta_data.start_time() > ( arrival_time + nano_time( 900, 0) ) ) // this edge is in future
	continue;
      //for present edges
      sssp_in_time_visitor neighbour(itr.target(), meta_data.end_time(), vertex);
      vis_queue->queue_visitor(neighbour);
    }
    return true;
  }

  friend inline bool operator>(const sssp_in_time_visitor& v1
			       , const sssp_in_time_visitor& v2 ) {
    return v1.arrival_time > v2.arrival_time;
  }

  friend inline bool operator<(const sssp_in_time_visitor& v1
			       , const sssp_in_time_visitor& v2 ){
    return v1.arrival_time < v2.arrival_time;
  }

  static void set_edge_metadata(EdgeMetaData* data) { edges_metadata() = data; }

  static EdgeMetaData*& edges_metadata() {
    static EdgeMetaData* data;
    return data;
  }
    
  static void set_arrival_data(ArrivalData* data) { arrival_data() = data; }
    
  static ArrivalData*& arrival_data() {
    static ArrivalData* data;
    return data;
  }

  static void set_parent_data(ParentData* data) { parent_data() = data; }

  static ParentData*& parent_data() {
    static ParentData* data;
    return data;
  }

  vertex_locator vertex;
  nano_time arrival_time;
  vertex_locator parent;
} __attribute__ ((packed));

    template < typename TGraph, typename EdgeMetaData, typename ArrivalData, typename ParentData>
void sssp_in_time( TGraph* g,
		   EdgeMetaData* edge_metadata,
		   ArrivalData* arrival_data,
		   ParentData* parent_data,
		   typename TGraph::vertex_locator source,
		   nano_time start_time){
  typedef sssp_in_time_visitor<TGraph, EdgeMetaData, ArrivalData, ParentData> visitor_type;
  visitor_type::set_edge_metadata(edge_metadata);
  visitor_type::set_arrival_data(arrival_data);
  visitor_type::set_parent_data(parent_data);
  
  typedef visitor_queue<visitor_type, detail::visitor_priority_queue, TGraph> visitor_queue_type;

  visitor_queue_type vq(g);
  visitor_type source_visitor(source, start_time, typename TGraph::vertex_locator() ); 
  vq.init_visitor_traversal(source_visitor);
}
    
    
} /*namespace mpi ends*/ } /*namespace havoqgt ends*/

#endif  //  _SSSP_IN_TIME_HPP_INCLUDED
