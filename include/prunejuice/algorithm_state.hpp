// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>

namespace prunejuice {

template<typename Vertex, typename VertexData, typename BitSet>
class vertex_state {
  public :
    vertex_state() :
    vertex_pattern_index(0),
    is_active(false) {}
  
    BitSet template_vertices;
    BitSet template_neighbors;
    //std::unordered_map<VertexData, IntegralType>
    //  template_neighbor_metadata_count_map;
 
    size_t vertex_pattern_index; // TODO: dummy, to be removed
    bool is_active; 		  
};
  
} // end namespace prunejuice 
