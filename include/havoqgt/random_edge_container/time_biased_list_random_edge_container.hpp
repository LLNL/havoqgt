#ifndef LIST_RANDOM_EDGE_CONTAINER_HPP
#define LIST_RANDOM_EDGE_CONTAINER_HPP

#include <random>

template<typename Graph, typename EdgeMetaData>
class random_edge_container{
public:

  /**********************************************************************
    Random Number Generator
  ***********************************************************************/
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


  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  typedef typename EdgeMetaData::value_type metadata_type;
 
  random_edge_container(Graph* _g, EdgeMetaData *_edge_metadata) : g(_g), edge_metadata(_edge_metadata) { }

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
  std::pair<bool, interval> get_random_weighted_edge(operation& op, interval cur_time, vertex_locator cur_vertex) const {
    std::vector<eitr_type> adjacents;
    std::vector<double> weights;
    for(eitr_type itr = g->edges_begin(cur_vertex); itr != g->edges_end(cur_vertex); ++itr ) {
      metadata_type& metadata = (*edge_metadata)[itr];
      std::pair<bool, interval> intersect_data = find_intersection( cur_time, metadata);
      if(intersect_data.first) {
	adjacents.push_back(itr);
	weights.push_back( (double)1.0/((double)( cur_time.first - metadata.start_time()) + 0.1) ); 
      }
    }
    if( adjacents.size() == 0) return std::make_pair(false, cur_time);

    std::discrete_distribution<uint32_t> discrete_dist(weights.begin(), weights.end());
    uint32_t index;
    {
      index = (rng::get_rng())(discrete_dist);
    }
    auto& itr = adjacents[index];
    op( (*edge_metadata)[itr], itr.target());    
    return find_intersection(cur_time, (*edge_metadata)[itr]);
  }

private:
  Graph *g;
  EdgeMetaData *edge_metadata;
};
template<typename G, typename M>
typename random_edge_container<G, M>::rng random_edge_container<G, M>::rng::_rng;

template<typename G, typename M>
int random_edge_container<G, M>::rng::val = 0;

#endif
