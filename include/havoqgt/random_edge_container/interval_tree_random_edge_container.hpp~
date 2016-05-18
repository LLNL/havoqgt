#ifndef _INTERVAL_TREE_RANDOM_EDGE_CONTAINER_HPP
#define _INTERVAL_TREE_RANDOM_EDGE_CONTAINER_HPP

#include <random>
#include <havoqgt/interval_tree.hpp>
#include <havoqgt/impl/vertex_data.hpp>
/**********************************************************************
    Random Number Generator
***********************************************************************/
namespace havoqgt{
  namespace mpi{
struct random_number_generator {
public:
  std::mt19937 en;
      
  static random_number_generator& get_rng() {
    if(val == 0) {
      _rng = random_number_generator();
      val = 1;
    }
    return _rng;
  }

  template<typename distribution>
  typename distribution::result_type operator()(distribution& dist){
    return dist(en);
  }

private:
  static random_number_generator _rng;
  static int val;
  random_number_generator(){ 
    std::random_device r;
    en.seed(r());
  }

};
typedef random_number_generator rng;
rng rng::_rng;
int rng::val = 0;

template<typename Graph, typename EdgeMetaData>
class random_edge_container{
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  typedef typename EdgeMetaData::value_type metadata_type;

  typedef interval_tree<uint64_t, eitr_type, EdgeMetaData> interval_tree_t;
  typedef typename Graph::template vertex_data<interval_tree_t, std::allocator<interval_tree_t>> vertex_data_type;

  vertex_data_type interval_tree_data;

  random_edge_container(Graph* _g, EdgeMetaData *_edge_metadata) : g(_g), edge_metadata(_edge_metadata)
								 , interval_tree_data(*_g) {
    //    std::cout << "Initializing tree" << std::endl;
    for( auto v_itr = g->vertices_begin(); v_itr != g->vertices_end(); ++v_itr) {
      vertex_locator vertex = *v_itr;
      interval_tree_t& tree = interval_tree_data[ vertex ];
      tree.data_ref( _edge_metadata);

      for( auto e_itr = g->edges_begin(vertex); e_itr != g->edges_end(vertex); ++e_itr) {
	auto& metadata = (*_edge_metadata)[e_itr];
	tree.add_keys( metadata.start_time() );
	tree.add_keys( metadata.end_time()   );
      }

      tree.create_tree( g->edges_begin(vertex), g->edges_end(vertex) );
    }
    //std::cout << "Finalizing tree"<<std::endl;
  }

  template<typename interval>
  std::pair<bool, interval> find_intersection( const interval& cur_time, const metadata_type& metadata) const {
    interval intersect;
    intersect.first = std::max( cur_time.first, metadata.start_time());
    intersect.second = cur_time.second;
    if( metadata.end_time() != 0) {
      intersect.second = std::min( cur_time.second, metadata.end_time());
    }
    return std::make_pair( intersect.first <= intersect.second, intersect);
  }

  template<typename operation, typename interval>
  std::pair<bool, interval> get_random_weighted_edge(operation& op, interval cur_time, vertex_locator cur_vertex) {
    auto& tree = interval_tree_data[cur_vertex];
    std::size_t size = 0;
    tree.query_size( cur_time.first, cur_time.second, size);
    if(size == 0) return std::make_pair(false, cur_time);
    std::uniform_int_distribution<std::size_t> uniform_dist( 0, size - 1);
    std::size_t index;
    {
      index = (rng::get_rng())(uniform_dist);
    }
    eitr_type itr = g->edges_begin(cur_vertex);
    tree.query_ith(cur_time.first, cur_time.second, index, itr);
    op( (*edge_metadata)[itr], itr.target());    
    return find_intersection(cur_time, (*edge_metadata)[itr]);
  }

private:
  Graph *g;
  EdgeMetaData *edge_metadata;
};
  }}
#endif
