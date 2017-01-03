#ifndef HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_HPP_INCLUDED
#define HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class lppm_queue {

public:
  lppm_queue() {}

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

// label propagation pattern matching visitor class
template<typename Graph, typename VertexData>
class lppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  lppm_visitor() : 
    itr_count(0) {}

  lppm_visitor(vertex_locator _vertex, uint64_t _itr_count = 0) : 
    vertex(_vertex), 
    itr_count(_itr_count) {}

  lppm_visitor(vertex_locator _vertex, vertex_locator _parent, VertexData _parent_vertex_data, uint64_t _itr_count) :
    vertex(_vertex), 
    parent(_parent),
    parent_vertex_data(_parent_vertex_data),
    parent_pattern_index(0), 
    itr_count(_itr_count) {}

  ~lppm_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);
    return true;
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);


    // does vertex_data match any entry in the query pattern
    bool match_found = false;
    size_t vertex_pattern_index = 0;
    for (; vertex_pattern_index < pattern.size(); vertex_pattern_index++) { 
      if (pattern[vertex_pattern_index] == vertex_data) {
        match_found = true;
        break;  
      }       
    }

    if (!match_found) {
      std::get<4>(alg_data)[vertex] = false;
      return false;
    } 

    std::cout << "visiting " << itr_count << std::endl;

    if (itr_count == 0 && match_found) {
      // send metadata and which vertex it is (in pattern_indices) to all the neighbours
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        // do this for all the pattern indices with the same pattern label 
        lppm_visitor new_visitor(neighbor, vertex, vertex_data, itr_count + 1); 
        vis_queue->queue_visitor(new_visitor); 
      }
    }

    if ((itr_count >= 1 && itr_count < pattern.size()) && match_found) {
      // heard from a neighbour (parent) whose vertex_data match an entry in the query pattern
      
      // verify if parent_vertex_data meet constrains for current vertex
      // what pattern label the vartex_data of this vertex corresponds to  
      // what are the valid parent labels for this vertex`   
      if (vertex_pattern_index == 0 && 
        !(parent_vertex_data == pattern[vertex_pattern_index + 1])) {
        return false;
      }
   
      if (!(parent_vertex_data == pattern[vertex_pattern_index - 1] || 
        parent_vertex_data == pattern[vertex_pattern_index + 1])) {
        return false; 
      }

      if (vertex_pattern_index == pattern.size() - 1 && 
        !(parent_vertex_data == pattern[vertex_pattern_index - 1])) {
        return false;
      }  

      lppm_visitor new_visitor(parent, vertex, vertex_data, itr_count + 1);
      vis_queue->queue_visitor(new_visitor); 
    }
 
    return true;
  }

  friend inline bool operator>(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator parent;
  VertexData parent_vertex_data; 
  size_t parent_pattern_index;
  uint64_t itr_count;
}; 

template <typename TGraph, typename VertexMetaData, typename VertexData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename VertexActive, 
  typename VertexIteration>
void label_propagation_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, VertexRank& vertex_rank,
  VertexActive& vertex_active, VertexIteration& vertex_iteration) {
  std::cout << "label_propagation_pattern_matching.hpp" << std::endl;

  typedef lppm_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank,
    vertex_active, vertex_iteration);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);
  vq.init_visitor_traversal_new(); 
  MPI_Barrier(MPI_COMM_WORLD);
  vertex_rank.all_reduce();
  MPI_Barrier(MPI_COMM_WORLD);
}  

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_HPP_INCLUDED 
