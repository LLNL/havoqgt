// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
//#include <unordered_map>

#include <boost/container/scoped_allocator.hpp>
#include <boost/container/vector.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

namespace prunejuice { namespace container {

template<typename Vertex, typename VertexData, typename Uint = uint64_t, 
  typename allocator_t = std::allocator<char>>
class motif_container {

  public :
    
    //motif_container() {} 

    //motif_container() : 
    //max_motif_count(0) {}

    explicit motif_container(allocator_t allocator = allocator_t())
      : motif_map(allocator),
        motif_unique_map(allocator),  
        vertex_motif_map(allocator) {} 

    //~motif_container() {}

    using unordered_set_allocator_t = 
      typename std::allocator_traits<allocator_t>::template rebind_alloc<Uint>;

    using unordered_set_t = 
      boost::unordered_set<Uint, std::hash<Uint>, std::equal_to<Uint>, 
      unordered_set_allocator_t>;

    using unordered_map_allocator_t = 
      boost::container::scoped_allocator_adaptor
      <typename std::allocator_traits<allocator_t>::
      template rebind_alloc<std::pair<const Uint, Uint>>>; // TODO: Uint, string
  
    using unordered_map_t = boost::unordered_map<Uint, Uint, 
      std::hash<Uint>, std::equal_to<Uint>, unordered_map_allocator_t>;

    using unordered_motif_set_map_allocator_t =
      boost::container::scoped_allocator_adaptor
      <typename std::allocator_traits<allocator_t>::
      template rebind_alloc<std::pair<const Vertex, unordered_set_t>>>;

    using unordered_motif_set_map_t = 
      boost::unordered_map<Vertex, unordered_set_t, std::hash<Vertex>, 
      std::equal_to<Vertex>, unordered_motif_set_map_allocator_t>;

    void insert(Vertex key, Vertex value) {
      //m_map[key].emplace_back(value);
      motif_map.emplace(key, value);
    }

    size_t size() {
      return motif_map.size(); 
    }  

    unordered_map_t motif_map;
    
    unordered_map_t motif_unique_map; 

    unordered_motif_set_map_t vertex_motif_map;

    /*size_t max_motif_count; 
  
    std::unordered_map<std::string, Uint> motif_map;
    
    typedef std::unordered_set<Uint> UintSet;
    std::unordered_map<Vertex, UintSet>vertex_motif_set_map;

    Uint getMotifID() {
      return 0;
    }*/
    
    //std::string collection_name;
    //std::string motif_name;   
	
};

}} // end namespace prunejuice::container 
 
