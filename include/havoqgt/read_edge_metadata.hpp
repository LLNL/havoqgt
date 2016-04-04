#ifndef _EDGE_METADATA_VISITOR_HPP_INCLUDED
#define _EDGE_METADATA_VISITOR_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/fifo_queue.hpp>
#include <boost/range/algorithm/lower_bound.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace havoqgt { namespace mpi {

    template <typename Graph>
    class edge_list_iterator
      : public boost::iterator_facade< edge_list_iterator<Graph>,
			       typename Graph::edge_iterator,
			       boost::random_access_traversal_tag,
			       typename Graph::edge_iterator const& >
    {
    public:
      typedef std::ptrdiff_t difference;

      edge_list_iterator() { }

      edge_list_iterator(typename Graph::edge_iterator _iterator)
	: edge_itr(_iterator) { }

      edge_list_iterator & operator=(edge_list_iterator const& _iterator) {
	edge_itr = _iterator.edge_itr;
	return *this;
       }

      bool equal(const edge_list_iterator& x) const {
	return (edge_itr == x.edge_itr);
      }

      typename Graph::edge_iterator const& dereference() const {
	return edge_itr;
      }

      difference distance_to(const edge_list_iterator& z) const {
	return z.edge_itr - edge_itr;
      }

      void increment() { ++edge_itr; }

      void decrement() { --edge_itr; }

      void advance(int n) { edge_itr += n; }

      typename Graph::edge_iterator get_iterator() const {
	return edge_itr;
      }

    private:
      typename Graph::edge_iterator edge_itr;
    };

    template<typename Graph>
    class edge_list_per_vertex {
    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;

      typedef edge_list_iterator<Graph> const_iterator;
      typedef edge_list_iterator<Graph> iterator;

      edge_list_per_vertex(): g(NULL) { }
      
      edge_list_per_vertex(const Graph* const _g, vertex_locator _source)
	:g(_g), source(_source) {
      }

      edge_list_iterator<Graph> begin() const {
	return edge_list_iterator<Graph>( g->edges_begin(source) );
      }

      edge_list_iterator<Graph> end() const {
	return edge_list_iterator<Graph>( g->edges_end(source) );
      }

      void sanity_check() const {
	eitr_type first = g->edges_begin(source);
	if(first  == g->edges_end(source)){
	  //std::cout << "The list is empty " << std::endl;
	  return;
	}
	
	eitr_type second = ++first;
	while(second != g->edges_end(source)) {
	  if( second.target() < first.target() ) {
	    std::cerr << "The list here isn't sorted" << std::endl;
	    return;
	  }
	  ++first;
	  ++second;
	}
      }
            
    private:
      const Graph* const g;
      vertex_locator source;
      
    };

    template<typename Graph, typename EdgeData, typename MetaData>
    class edge_metadata_visitor;

    template<typename Graph, typename EdgeData, typename MetaData>
    class sort_predicate {
    public:
      typedef typename Graph::vertex_locator vertex_locator;
      typedef typename Graph::edge_iterator eitr_type;
            
      friend class edge_metadata_visitor<Graph, EdgeData, MetaData>;

      typedef edge_metadata_visitor<Graph, EdgeData, MetaData> edge_metadata_visitor_t;

      bool operator()(const eitr_type& e1, const vertex_locator &target) const {
	  if(target == e1.target()){
	    return (*(edge_metadata_visitor_t::edge_data()))[e1].is_recorded();
	  }else
	    return  e1.target() < target;
      }

    private:
      const Graph* const g;
    };
    
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
      
#if 0 // MAP STORAGE
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
	return false; //terminate the visit here
      }else {
	eitr_type *__eitr_ptr = edge_iterator(vertex, _target);
	if(__eitr_ptr == NULL) (*edge_map())[std::make_pair(vertex, _target)] = new eitr_type(*eitr);
      }
    }
    return (vertex.get_bcast() == 0 ) ? true : false ; // dont re broadcast
  }
#elseif 0 // FOR ALL VERTICES
  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator _target = eitr.target();
      if( _target == target && !((*edge_data())[eitr].is_recorded()) ) {
	(*edge_data())[eitr] = meta_data;
	(*edge_data())[eitr].register_recorded();
	++eitr;
	count++;
	return false;
      }
    }
    if( vertex.get_bcast() != 0 ) {
      return false;
    } else {
      return true;
    }
  }

#else // BINARY SEARCH
  template<typename VisitorQueueHandle>
  bool visit(Graph& g, VisitorQueueHandle vis_queue) const {
    edge_list_per_vertex<Graph> el(&g, vertex);
 #if 0 // Test if the edge list of the vertex is sorted   
    //  el.sanity_check();
 #endif

    if(g.edges_begin(vertex) == g.edges_end(vertex)){
      if(vertex.is_delegate()){
	return true;
      }
      else{
	std::cerr << "Logic Problem!!" << std::endl;
	return false;
      }
    }
    
    typedef sort_predicate<Graph, EdgeData, MetaData> sort_predicate;
    edge_list_iterator<Graph> found_itr = boost::lower_bound(el, target, sort_predicate());

    eitr_type found = found_itr.dereference();
    
    if(found_itr == el.end() || !( target == found.target() /* this should have been covered by first statement*/) ) {
      assert(vertex.is_delegate());
      return true;
    }else {
#if 0 // Test result of the Binary Search
      assert(sort_predicate()(found, target) == false);
      assert(target == found.target());
      assert((*edge_data())[found].is_recorded() == false );
#endif
      count++;
      (*edge_data())[found] = meta_data;
      (*edge_data())[found].register_recorded();
      return false;
    }
  }
#endif
 
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
  static uint64_t count;
  //static edge_to_eitr_map_type* edge_to_eitr_map;
} __attribute__ ((packed));

    /*
template <typename TGraph, typename EdgeData>
void read_edge_metadata(TGraph *g,
			EdgeData* edge_data, ingest_flow_edge_list::flow_input_iterator flow_itr_begin, ingest_flow_edge_list::flow_input_iterator flow_itr_end, int& count) {
  typedef edge_metadata_visitor<TGraph, EdgeData, flow> visitor_type;
  visitor_type::set_edge_data( edge_data );
  
  typedef visitor_queue< visitor_type, detail::visitor_priority_queue, TGraph> visitor_queue_type;
  
  visitor_queue_type vq(g);
  vq.init_visitor_traversal_flow(flow_itr_begin, flow_itr_end);

  count = edge_metadata_visitor<TGraph, EdgeData, flow>::count;
  }*/

  template<typename TGraph, typename EdgeData, typename metadata_itr_t, typename metadata_t>
  void generic_read_edge_metadata( TGraph *g, EdgeData* edge_data, metadata_itr_t metadata_itr_begin
				   , metadata_itr_t metadata_itr_end, uint64_t& count ) {
    typedef edge_metadata_visitor<TGraph, EdgeData, metadata_t> visitor_type;
    visitor_type::set_edge_data( edge_data);
    
    typedef visitor_queue< visitor_type, detail::fifo_queue, TGraph> visitor_queue_type;
    
    visitor_queue_type vq(g);
    vq.template init_visitor_traversal_metadata<metadata_itr_t, metadata_t>( metadata_itr_begin, metadata_itr_end );

    count = edge_metadata_visitor<TGraph, EdgeData, metadata_t>::count;
  }

}} //end namespace havoqgt::mpi

#endif   //  _EDGE_METADATA_VISITOR_HPP_INCLUDED
