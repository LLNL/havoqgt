#pragma once

#include <unordered_map>

namespace prunejuice {

enum class ApproximateQuery {None, EditDistance, OptionalEdge};

template<typename Vertex, typename VertexData, typename BitSet>
class vertex_state {
  public :
    vertex_state() :
    //vertex_pattern_index(0),
    is_active(false), new_global_vertex_ID(0) {}
  
    BitSet template_vertices;
    BitSet template_neighbors;

    //std::unordered_map<VertexData, IntegralType>
    //  template_neighbor_metadata_count_map;
 
    //size_t vertex_pattern_index; // TODO: dummy, to be removed
    
    bool is_active; 		 

    // remapped pruned graph 
    Vertex new_global_vertex_ID;    
    std::unordered_set<Vertex> new_global_neighbor_IDs; // TODO: can it be a vector?
};
  
} // end namespace prunejuice 
