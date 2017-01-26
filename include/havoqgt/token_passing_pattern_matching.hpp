#ifndef HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_HPP_INCLUDED
#define HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename Visitor>
class tppm_queue {

public:
  tppm_queue() {}

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

// token passing pattern matching visitor class
template<typename Graph>
class tppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  tppm_visitor() : 
  itr_count(0), 
  do_pass_token(false), 
  is_init_step(true), 
  source_index_pattern_indices(0), 
  parent_pattern_index(0) {}

  tppm_visitor(vertex_locator _vertex) :  
    vertex(_vertex), 
    itr_count(0),
    do_pass_token(false), 
    is_init_step(true), 
    source_index_pattern_indices(0), 
    parent_pattern_index(0) {}
   
  tppm_visitor(vertex_locator _vertex, vertex_locator _target_vertex, 
    size_t _itr_count, size_t _max_itr_count, size_t _source_index_pattern_indices, 
    size_t _parent_pattern_index, 
    bool _expect_target_vertex = true, bool _do_pass_token = true, 
    bool _is_init_step = false) : 
    vertex(_vertex),
    target_vertex(_target_vertex), 
    itr_count(_itr_count), 
    max_itr_count(_max_itr_count), 
    expect_target_vertex(_expect_target_vertex), 
    do_pass_token(_do_pass_token), 
    is_init_step(_is_init_step), 
    source_index_pattern_indices(_source_index_pattern_indices), 
    parent_pattern_index(_parent_pattern_index) {}  

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    // TODO: pre-visit on local vertices 
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    // TODO: verify if this vertex is alive
    // TODO: add parent and verify if the originating parent is valid in terms of label and pattern index
    
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);   
    auto& pattern_graph = std::get<4>(alg_data); 

    if (!do_pass_token && is_init_step) {
      // create visitors only for the source vertices
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();

        if (g.locator_to_label(neighbor) != 46) { // TODO: put all the valid vertices form label propagation in a queue and pop here and create visitors like this 
          continue;    
        }  
        // 65 66 67 66 67 69 68 66 67 66 67 : 0 1 2 3 4 5 6 3 4 1 2 
        tppm_visitor new_visitor(neighbor, g.label_to_locator(46), 0, 7, 1, pattern_indices[1], true, true, true); 
        vis_queue->queue_visitor(new_visitor);
      }
      return true;
 
    } else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) && itr_count == 0 && is_init_step) {
      // initiate token passing from the source vertex
      std::cout << "found source vertex " << g.locator_to_label(vertex) << " vertex_data " << vertex_data << std::endl; // Test
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        tppm_visitor new_visitor(neighbor, target_vertex, itr_count, max_itr_count, 
          source_index_pattern_indices, parent_pattern_index, true, true, false);
        vis_queue->queue_visitor(new_visitor);
      }  
      return true;
 
    } else if (!is_init_step) { // else if      

      bool do_forward_token = false;
      auto new_itr_count = itr_count + 1;
      auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index  
 
      auto vertex_pattern_index = pattern_indices[source_index_pattern_indices + new_itr_count]; // TODO: read from the map !!!!!        

      if (max_itr_count > itr_count) {

        // are vertex_data and vertex_pattern_index valid
        if (vertex_data == pattern[pattern_indices[next_pattern_index]] && 
          vertex_pattern_index == pattern_indices[next_pattern_index]) {
          // verify if received from a valid parent
          for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
            e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
            if (pattern_graph.edges[e] == parent_pattern_index) {
              do_forward_token = true; 
              break; 
            }  
          } // for      
        } // if    
 
      } else if (max_itr_count <= itr_count) {
        // is this the target vertex
        // TODO: also verify parent_pattern_index 
        if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex)) {
          // found loop
          std::cout << "found loop - vertex " <<  " | parent_pattern_index " <<  parent_pattern_index <<  " | " << g.locator_to_label(vertex) << " itr " << itr_count << std::endl; // Test
          return false; // TODO: true ?		
        } else {
          // reached max iteration but did not find the target vertex or a loop
          //std::cout << "At " << g.locator_to_label(vertex) <<  ", did not find target " 
          //<< g.locator_to_label(target_vertex) <<  " after " << itr_count 
          //<< " iterations" <<std::endl; // Test
          return false;  
        }   
      } else {
        std::cerr << "Error: wrong code branch." << std::endl;    
        return false;
      }   

      if (!do_forward_token) {
        return false;
      }

      // forward along
      std::cout << g.locator_to_label(vertex) << " " << new_itr_count << " forwarding ... " << g.locator_to_label(target_vertex) << std::endl; // Test

      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        tppm_visitor new_visitor(neighbor, target_vertex, new_itr_count, max_itr_count, 
          source_index_pattern_indices, vertex_pattern_index); // vertex_pattern_index = parent_pattern_index for the neighbours 
        vis_queue->queue_visitor(new_visitor);
      }
      return true;

      // else if
    } else {
      return false;
    }
  
  }

  friend inline bool operator>(const tppm_visitor& v1, const tppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const tppm_visitor& v1, const tppm_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator target_vertex;
  size_t itr_count; // TODO: change type
  size_t max_itr_count; // equal to diameter - 1 of the pattern as itr_count is initialized to 0 // TODO: change type
  bool expect_target_vertex;
  bool do_pass_token;
  bool is_init_step;
  size_t source_index_pattern_indices; // index of the token source in the pattern_indices container 
  size_t parent_pattern_index; // TODO: change to the same type as in the pattern_graph
};

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename PatternGraph>
void token_passing_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, 
  VertexRank& vertex_rank, PatternGraph& pattern_graph) {
  //std::cout << "token_passing_pattern_matching.hpp" << std::endl;
  typedef tppm_visitor<TGraph> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank, pattern_graph);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);
  vq.init_visitor_traversal_new();
  MPI_Barrier(MPI_COMM_WORLD);
}

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_HPP_INCLUDED

