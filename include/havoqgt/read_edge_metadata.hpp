#ifndef _EDGE_METADATA_VISITOR_HPP_INCLUDED
#define _EDGE_METADATA_VISITOR_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {
    
template<typename Graph, typename EdgeData, typename MetaData>
class edge_metadata_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  
  edge_metadata_visitor() {}
  edge_metadata_visitor(vertex_locator _vertex, vertex_locator _target, MetaData _meta_data)
    : vertex(_vertex)
    , target(_target)
    , meta_data(_meta_data) {
    
  }

  bool pre_visit() const {
    return true;
  }

  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    eitr_type *eitr_ptr = edge_iterator(vertex, target);
    if(eitr_ptr == NULL) (*edge_map())[std::make_pair(vertex, target)] = new eitr_type(g.edges_begin(vertex));
    for(eitr_type* eitr = edge_iterator(vertex, target); *eitr != g.edges_end(vertex); ++(*eitr)) {
      vertex_locator _target = eitr->target();
      if(_target == target) {
	(*edge_data())[(*eitr)] = meta_data;
	++(*eitr);
	count++;
	return true;
      }else {
	eitr_type *__eitr_ptr = edge_iterator(vertex, _target);
	if(__eitr_ptr == NULL) (*edge_map())[std::make_pair(vertex, _target)] = new eitr_type(*eitr);
      }
    }
    return false;
  }
  
  friend inline bool operator>(const edge_metadata_visitor& v1, const edge_metadata_visitor& v2) {
    return true;
  }
			   

  static void set_edge_data(EdgeData* _data) { edge_data() = _data; }

  static EdgeData*& edge_data() {
    static EdgeData* data;
    return data;
  }

  typedef boost::unordered_map<std::pair<vertex_locator, vertex_locator>, eitr_type*, boost::hash<std::pair<vertex_locator, vertex_locator>>> edge_to_eitr_map_type;
 
  
  static eitr_type* edge_iterator(vertex_locator _source, vertex_locator _target) {
    try{
      return edge_map()->at(std::make_pair(_source, _target));
    }catch(const std::out_of_range& ex){
      return NULL;
    }
  }

  static edge_to_eitr_map_type*& edge_map() {
    static edge_to_eitr_map_type* map = new edge_to_eitr_map_type();
    return map;
  }

  vertex_locator vertex;
  vertex_locator target;
  MetaData meta_data;
  static int count;
  //static edge_to_eitr_map_type* edge_to_eitr_map;
} __attribute__ ((packed));

template <typename TGraph, typename EdgeData>
void read_edge_metadata(TGraph *g,
			EdgeData* edge_data, ingest_flow_edge_list::flow_input_iterator flow_itr_begin, ingest_flow_edge_list::flow_input_iterator flow_itr_end, int& count) {
  typedef edge_metadata_visitor<TGraph, EdgeData, flow> visitor_type;
  visitor_type::set_edge_data( edge_data );
  
  typedef visitor_queue< visitor_type, detail::visitor_priority_queue, TGraph> visitor_queue_type;
  
  visitor_queue_type vq(g);
  vq.init_visitor_traversal_flow(flow_itr_begin, flow_itr_end);

  count = edge_metadata_visitor<TGraph, EdgeData, flow>::count;
}

}} //end namespace havoqgt::mpi

#endif   //  _EDGE_METADATA_VISITOR_HPP_INCLUDED
