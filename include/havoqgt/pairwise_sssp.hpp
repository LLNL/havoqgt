#ifndef _PAIRWISE_SSSP_IN_TIME_HPP_INCLUDED
#define _PAIRWISE_SSSP_IN_TIME_HPP_INCLUDED

#include <havoqgt/nano_time.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <unordered_map>

namespace hdetail = havoqgt::detail;

//namespace havoqgt { namespace mpi {
template<typename Graph, typename T>
class vertex_state_map {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  
  std::unordered_map<uint64_t, vertex_locator> parent_map;
  std::unordered_map<uint64_t, T> state_map;

  static uint64_t traversal_count;
  
  static Graph*& graph() {
    static Graph* graph;
    return graph;
  }

  static void set_graph(Graph* g) {
    graph() = g;
  }
  
  bool exists(vertex_locator source) const {
    auto source_idx = parent_map.find( graph()->locator_to_label(source) );
    return (source_idx != parent_map.end());
  }

  bool exists_with_state(vertex_locator source) const {
    auto source_idx = state_map.find( graph()->locator_to_label(source) );
    return (source_idx != state_map.end());
  }
 
  T& get_state(const vertex_locator& source) {
    return state_map[ graph()->locator_to_label(source) ];
  }

  void set_state(const vertex_locator& source, T state) {
    state_map.insert(std::make_pair(graph()->locator_to_label(source), state));
  }

  vertex_locator& get_parent(const vertex_locator& source) {
    return parent_map[graph()->locator_to_label(source)];
  }

  void set_parent(const vertex_locator& source, const vertex_locator& parent) {
    parent_map.insert(std::make_pair(graph()->locator_to_label(source), parent));
  }

  typename std::unordered_map<uint64_t, vertex_locator>::const_iterator source_begin() const{
    return parent_map.cbegin();
  }

  typename std::unordered_map<uint64_t, vertex_locator>::const_iterator source_end() const {
    return parent_map.cend();
  }

};

template<typename Graph, typename T>
uint64_t vertex_state_map<Graph, T>::traversal_count = 0;
/*
template <typename Graph, typename T>
class vertex_state {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  
  static Graph*& graph(){
    static Graph* graph;
    return graph;
  }

  static void set_graph(Graph* g) {
    graph() = g;
  }

  static size_t& number_of_sorted_source() {
    static size_t sz = 0;
    return sz;
  }

  static void set_number_of_sorted_source(size_t sz) { number_of_sorted_source() = sz; }

  static uint64_t*& sorted_source_labels() {
    static uint64_t* labels;
    return labels;
  }

  static void set_sorted_source_labels(uint64_t* labels) { sorted_source_labels() = labels; }
  
    static size_t get_source_label_index(uint64_t label_index) {
      //return binary search index
      size_t left = 0, right = number_of_sorted_source() - 1;
      while( left < right ) {
	size_t mid = left + ( (right - left) >> 1 );	
	if( sorted_source_labels()[mid] > label_index ) right = mid - 1;
	else if( sorted_source_labels()[mid] < label_index ) left = mid + 1;
	else return mid;
      }
      return left;
    }

    vertex_state() {}
    
    vertex_state(int _n, T val)
      : n(_n) {
      labels = new T[_n]; // make an array for each of the sources
      parents = new vertex_locator[_n];
      reset(val);
    }

    void reset(T val) {
      for( size_t i = 0; i < n; i++) {
	labels[i] = val;
      }
    }

  bool exists(vertex_locator locator) const {
    uint64_t label = graph()->locator_to_label(locator);
    size_t index = vertex_state::get_source_label_index(label);
    if(index < 0 || index >= number_of_sorted_source() || sorted_source_labels()[index] != label) return false;
    return true;
  }

    T& get_state(vertex_locator locator) const {
      uint64_t label = graph()->locator_to_label(locator);
      return labels[ vertex_state::get_source_label_index(label) ];
    }

  vertex_locator& get_parent(const vertex_locator locator) const {
    return parents[vertex_state::get_source_label_index( graph()->locator_to_label(locator) ) ];
  }
    vertex_locator& set_parent(const vertex_locator& locator) const{
      uint64_t label = graph()->locator_to_label(locator);
      return parents[ vertex_state::get_source_label_index(label) ];
    }
    
  private:
    uint64_t n;
    T *labels;
    vertex_locator *parents;
    
  };
*/
namespace havoqgt { namespace mpi {

template<typename Graph, typename EdgeMetaData, typename VertexData>
class pairwise_sssp_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator edge_iterator;

  pairwise_sssp_visitor() { }

  pairwise_sssp_visitor( vertex_locator _vertex, vertex_locator _source, nano_time _arrival_time, vertex_locator _parent)
    : vertex( _vertex ), source( _source ), arrival_time( _arrival_time), parent(_parent) { }

  bool pre_visit() const {
    vertex_state_map<Graph, nano_time>::traversal_count++;
    if( !((*vertex_data())[vertex].exists_with_state(source)) ) {
      (*vertex_data())[vertex].set_state(source, nano_time( std::numeric_limits<uint32_t>::max(),
							    std::numeric_limits<uint32_t>::max() ) );
      (*vertex_data())[vertex].set_parent(source, parent);
    }
    
    if( arrival_time < (*vertex_data())[vertex].get_state(source)) {
      (*vertex_data())[vertex].get_state(source) = arrival_time;
      (*vertex_data())[vertex].set_parent(source, parent);
      return true;
    }
    return false;
  }

  template<typename VisitorQueueHandle>
  bool visit( Graph& g, VisitorQueueHandle vis_queue) const {
    //Boiler plate stuffs  ( we always loop through the out edges and filter out the edges not available in time )
    typedef typename EdgeMetaData::value_type metadata_type;
    
    for(edge_iterator edge_itr = g.edges_begin(vertex);
	edge_itr != g.edges_end(vertex);
	++edge_itr ) {
      metadata_type metadata = (*(edge_metadata()))[ edge_itr ];
      if( arrival_time > metadata.end_time() || ( arrival_time + waiting_time ) < metadata.start_time() ) continue;

      //else queue new visitor
      vis_queue->queue_visitor(pairwise_sssp_visitor( edge_itr.target(), source, metadata.end_time(), vertex ) );
    }
    return true;
  }

  friend bool operator<( const pairwise_sssp_visitor& a, const pairwise_sssp_visitor& b) {
    return  a.arrival_time < b.arrival_time;
  }

  friend bool operator>(const pairwise_sssp_visitor& a, const pairwise_sssp_visitor& b) {
    return a.arrival_time > b.arrival_time;
  }

  static void set_edge_metadata( EdgeMetaData *data) {
    edge_metadata() = data;
  }

  static EdgeMetaData*& edge_metadata(){
    static EdgeMetaData* data;
    return data;
  }

  static void set_vertex_data( VertexData *data) {
    vertex_data() = data;
  }

  static VertexData*& vertex_data() {
    static VertexData* data;
    return data;
  }

  static nano_time waiting_time;

  
  nano_time arrival_time;
  vertex_locator vertex;
  vertex_locator source;
  vertex_locator parent;
} __attribute__ ((packed));

template<typename TGraph, typename EdgeMetadata, typename VertexData>
nano_time pairwise_sssp_visitor<TGraph, EdgeMetadata, VertexData>::waiting_time = nano_time(0, 0);

template <typename TGraph, typename EdgeMetadata, typename VertexData>
void run_pairwise_sssp( TGraph* g,
			EdgeMetadata* edge_metadata,
			VertexData* vertex_data,
			std::vector<typename TGraph::vertex_locator>& vertex_list,
			nano_time& start_time,
			nano_time& waiting_time) {

  typedef pairwise_sssp_visitor<TGraph, EdgeMetadata, VertexData> visitor_type;

  visitor_type::set_edge_metadata( edge_metadata );
  visitor_type::set_vertex_data( vertex_data );
  visitor_type::waiting_time = waiting_time;
  
  typedef visitor_queue<visitor_type, hdetail::visitor_priority_queue, TGraph> visitor_queue_type;

  visitor_queue_type vq(g);
  
  std::vector<visitor_type> visitor_list;

  for(auto vertex_list_itr = vertex_list.begin(); vertex_list_itr != vertex_list.end(); vertex_list_itr++) {
    visitor_list.push_back( visitor_type ( *vertex_list_itr, *vertex_list_itr, start_time, typename TGraph::vertex_locator())); 
  }

  vq.init_visitor_traversal( visitor_list );  
}
  
} }

#endif  //  _PAIRWISE_SSSP_IN_TIME_HPP_INCLUDED
 
