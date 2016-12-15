#ifndef HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED
#define HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class fpm_queue {

public:
  fpm_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    data.pop_back();
  }
 
  Visitor const& top() {
    return data.back();
  } 
  
  size_t size() const {
    return data.size();;
  }

  bool empty() const {
    return data.empty();
  }

  void clear() {
    data.clear();
  }

protected:
  std::vector<Visitor> data;

};

template<typename Graph>
class fpm_walker {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  fpm_walker() {}

  fpm_walker(vertex_locator _vertex) :
    vertex(_vertex),
    current_pattern_index(0) {
  }

  fpm_walker(vertex_locator _vertex, size_t _current_pattern_index) :
    vertex(_vertex), 
    current_pattern_index(_current_pattern_index) {
  }

  ~fpm_walker() {
    //std::cout << "Terminating walker." << std::endl;
  } 

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);
    
    if (vertex_data != pattern[current_pattern_index]) {
      //std::cout << "Pre-vist: Fail" << std::endl; 
      return false;
    } else { 
      //std::cout << "Pre-vist: Pass" << std::endl; 
      return true;
    }
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);

    if (vertex_data != pattern[current_pattern_index]) {
      //std::cout << "Fail" << std::endl;
      return false;
    } else {
      //std::cout << "Pass" << std::endl;
      if (current_pattern_index == pattern.size() - 1) {
        std::cout << "Found" << std::endl;
        return false;
      }
    }

    size_t next_pattern_index = current_pattern_index + 1;

    //std::cout << g.locator_to_label(vertex) << " " << vertex_data << " " << pattern[3] << std::endl; 
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      auto neighbor_data = std::get<0>(alg_data)[neighbor]; 

      // create clones
      fpm_walker new_visitor(neighbor, next_pattern_index);
      vis_queue->queue_visitor(new_visitor);  
    }
    return true;
  }

  friend inline bool operator>(const fpm_walker& v1, const fpm_walker& v2) {
    return false;
  }

  friend inline bool operator<(const fpm_walker& v1, const fpm_walker& v2) {
    return false;
  }

  vertex_locator vertex;
  //vertex_locator parent;
  size_t current_pattern_index;
} __attribute__ ((packed));

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices>
void  fuzzy_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices) {
  std::cout << "fuzzy_pattern_matching.hpp" << std::endl;

  typedef fpm_walker<TGraph> visitor_type;
  size_t current_pattern_index = 0;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);
  vq.init_visitor_traversal_new();  
} 

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED
