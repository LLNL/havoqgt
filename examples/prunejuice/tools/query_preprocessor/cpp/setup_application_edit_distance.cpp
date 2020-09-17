//#include <atomic>
//#include <cmath>
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "omp.h"

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/binomial.hpp>
//#include <boost/random.hpp>
//#include <boost/random/random_device.hpp>
//#include <boost/random/discrete_distribution.hpp>

//#include "graph.hpp"
//#include "random_walk.hpp"
#include "template_constraint.hpp"
#include "template_graph.hpp"
#include "combination.hpp"
#include "graph_algorithm.hpp"
#include "edge_dfs.hpp"
#include "file_utilities.hpp"
#include "utilities.hpp"

//enum DISTRIBUTION {UNIFORM, RANDOM, BIASED, ALL};

////////////////////////////////////////////////////////////////////////////////

// TODO: move them to hpp files

////////////////////////////////////////////////////////////////////////////////

/**
 * For hash-based containers
 */ 
template <typename T> 
bool has_common_element(T& a, T& b) {
  bool found_common_element = false;
  for (auto& i : b) {
    auto find_key = a.find(i.first);
    if (find_key != a.end()) {
      found_common_element = true;
      break; 
    } 
  }  
  return found_common_element;
}

////////////////////////////////////////////////////////////////////////////////
 
template <typename Vertex, typename Edge, typename EdgeListMap, 
  typename EdgeSet>
bool is_subset_edgeset(EdgeListMap& a, EdgeSet& b) {
  bool is_subset = true;
  for(auto it = b.begin(); it != b.end(); ++it) {
    auto find_it = a.find(it->first);  
    if (find_it == a.end()) {
      is_subset = false;
      break;   
    }
  } // for
  return is_subset;  
}  

////////////////////////////////////////////////////////////////////////////////

constexpr size_t VERTEX_COUNT = 64;
typedef std::bitset<VERTEX_COUNT> VertexBitSet;

template <typename Vertex, typename VertexBitSet>
void update_vertex_bitset(Vertex s, Vertex t, VertexBitSet& vertex_bitset) {
  vertex_bitset.set(static_cast<size_t>(s));
  vertex_bitset.set(static_cast<size_t>(t));       
}

template <typename Vertex, typename Edge>
Edge edge_hash_vertex_bitset(Vertex s, Vertex t) {
  VertexBitSet vertex_bitset; 
  vertex_bitset.set(static_cast<size_t>(s));
  vertex_bitset.set(static_cast<size_t>(t));
  return static_cast<Edge>(vertex_bitset.to_ullong());   
}

////////////////////////////////////////////////////////////////////////////////

/**
 * All template constraints
 */
template <typename Vertex, typename Edge, typename EdgeSet, 
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths, typename EdgeListTuple,
  typename EdgeListTupleVector, typename EdgeListTupleVectors,
  typename TemplateConstraintVector>
void populate_template_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique, 
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  GraphNonLocalPropertiesUnique& graph_tds_cycles_unique,
  EdgeListTupleVectors& graph_tds_cycles_all_paths,
  EdgeListTuple& graph_enumeration_unique,   
  TemplateConstraintVector& template_constraints) {
  
  for (auto it = graph_cycles_unique.begin(); 
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first, 
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));    
  } // for
  
  for (auto it = graph_tds_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      graph_tds_cycles_all_paths[it->first], // EdgeListTupleVector 
      template_constraints.size(), TemplateConstraint::TDS));
  } // for 

  EdgeListTupleVector graph_enumeration_path(0); // TODO: improve 
    graph_enumeration_path.push_back(graph_enumeration_unique);

  template_constraints.push_back(TemplateConstraint(0, 
    graph_enumeration_path, // EdgeListTupleVector //  
    template_constraints.size(), TemplateConstraint::ENUMERATION));  
} 

template <typename Vertex, typename Edge, typename EdgeSet, 
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths, typename EdgeListTupleVectors,
  typename TemplateConstraintVector>
void populate_template_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique, 
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  GraphNonLocalPropertiesUnique& graph_tds_cycles_unique,
  EdgeListTupleVectors& graph_tds_cycles_all_paths,   
  TemplateConstraintVector& template_constraints) {
  
  for (auto it = graph_cycles_unique.begin(); 
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first, 
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));    
  } // for
  
  for (auto it = graph_tds_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      graph_tds_cycles_all_paths[it->first], // EdgeListTupleVector 
      template_constraints.size(), TemplateConstraint::TDS));
  } // for  
} 

template <typename Vertex, typename Edge, typename EdgeSet,
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique,
  typename GraphNonLocalPropertiesAllPaths,
  typename TemplateConstraintVector>
void populate_cycle_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique,
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  TemplateConstraintVector& template_constraints) {

  for (auto it = graph_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));
  } // for
}

////////////////////////////////////////////////////////////////////////////////

template <typename Vertex, typename VertexSet, typename EdgeSet>
VertexSet edgeset_to_vertexset(EdgeSet& edgeset) {
  VertexSet vertexset(0);
  for (auto e : edgeset) {
    auto find_s = vertexset.find(std::get<0>(e.second)); 
    if (find_s == vertexset.end()) {
      vertexset.insert(std::get<0>(e.second));
    }
    
    auto find_t = vertexset.find(std::get<1>(e.second));  
    if (find_t == vertexset.end()) {
      vertexset.insert(std::get<1>(e.second));
    }
  } 
  return vertexset; 
}

/**
 * Generate TDS cycle constraints
 */
template <typename Vertex, typename Edge, typename EdgeSet,
  typename GraphNonLocalPropertiesUnique, typename  EdgeListTuple, 
  typename EdgeListTupleVector, typename EdgeListTupleVectors>
void generate_tds_cycle_constraints
  (GraphNonLocalPropertiesUnique& graph_tds_cycles_unique, 
  EdgeListTupleVectors& graph_tds_cycles_all_paths, bool do_print = false) {
 
  typedef std::unordered_set<Vertex> VertexSet;

  GraphNonLocalPropertiesUnique new_graph_tds_cycles_unique(0);

  size_t edgeset_hash = 0;

  for (auto it = graph_tds_cycles_unique.begin(); 
    it != graph_tds_cycles_unique.end(); ++it) {

    //std::cout << it->first << "(x) : ";
    //for (auto& j : it->second) { // it->second - EdgeSet
    //  std::cout << "(" << j.first << ", "
    //  << std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", "
    //  << std::get<2>(j.second) << "(x)), ";
    //}
    //std::cout << std::endl;  

    // it->second is of type EdgeSet
    VertexSet vertexset = edgeset_to_vertexset<Vertex, VertexSet>(it->second);   

    EdgeListTupleVector edgeset_all_paths(0); 

    for (auto v : vertexset) {
      //std::cout << v << ", ";
      edgeset_all_paths.push_back(prunejuice::pattern::graphalgorithm::edge_dfs
        <Vertex, Edge, EdgeListTuple, EdgeSet>(it->second, v));

      // Test 
      assert(edgeset_all_paths.size() > 0);   
      if (do_print) {
        for (auto i : edgeset_all_paths[edgeset_all_paths.size() - 1]) {
          std::cout << "(" << std::get<0>(i) << ", " << std::get<1>(i) << "), ";
        }
        std::cout << std::endl;
      }       
      // Test
       
    } // for 
    
    // Test
    if (do_print) { 
      std::cout << std::endl;
    } 
    // Test
 
    graph_tds_cycles_all_paths.push_back(edgeset_all_paths);

    // update graph_tds_cycles_unique to have new edgeset_hash
    // graph_tds_cycles_all_paths[0] contains the paths for edgeset_hash == 0 
    new_graph_tds_cycles_unique.insert({edgeset_hash, it->second});
    edgeset_hash++;
  
    //std::cout << "graph_tds_cycles_all_paths.size(): " 
    //  << graph_tds_cycles_all_paths.size() << std::endl; // Test

  } // for

  graph_tds_cycles_unique = new_graph_tds_cycles_unique; // Important:

} 

/**
 * Generate TDS cycles from all the unique cycles
 */ 
template <typename TDSEdgeSet, typename EdgeSetIterator>
void generate_tds_cycles_recursive(TDSEdgeSet& tds_edge_set,
  EdgeSetIterator it, EdgeSetIterator it_end) {
  if (it == it_end) {
    return;
  } else {
    // add edge
    auto find_edge_hash = tds_edge_set.find(it->first);
    if (find_edge_hash == tds_edge_set.end()) {
      //tds_edge_set.insert({it->first, 
      //  std::forward_as_tuple(std::get<0>(it->second),std::get<1>(it->second),
      //  std::get<2>(it->second))});
      tds_edge_set.insert({it->first, it->second});  
    }
    generate_tds_cycles_recursive(tds_edge_set, ++it, it_end);
  }  
}

template <typename Veretx, typename Edge, 
  typename GraphNonLocalPropertiesUnique>
void generate_tds_cycles(GraphNonLocalPropertiesUnique graph_cycles_unique, // copy
  GraphNonLocalPropertiesUnique& tds_cycles) {
  
  for (auto it = graph_cycles_unique.begin(); it != graph_cycles_unique.end();) {
   
    // add to tds_cycles  
    auto find_it = tds_cycles.find(it->first); 
    if (find_it == tds_cycles.end()) {
      auto insert_status = tds_cycles.insert({it->first, it->second});

      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map. " << std::endl;
        return;
      } 

      auto& tds_cycles_it = insert_status.first; // iterator 
      it = graph_cycles_unique.erase(it); // C++11 // remove from graph_cycles_unique      
    
      //std::cout << "tds_cycles size: " << tds_cycles.size() << " " << std::endl; // Test
      //  std::cout << "graph_cycles_unique size: " << graph_cycles_unique.size() << std::endl; // Test

      bool is_tds_cycle = false;

      // find intersection of cycle constraints
      for (auto it_nest = it; it_nest != graph_cycles_unique.end();) {
 
        // if intersection is found, it_nest is merged with tds_cycles_it 
        // and removed from graph_cycles_unique  
        bool found_intersection = has_common_element
          (insert_status.first->second, it_nest->second);  

        if (found_intersection) {

          // merege cycles
          generate_tds_cycles_recursive(insert_status.first->second, 
            it_nest->second.begin(), it_nest->second.end());

          if (it == it_nest) { // Important: must update it
            it_nest = graph_cycles_unique.erase(it_nest); // C++11 // remove it_nest from graph_cycles_unique
            it = it_nest;
          } else {
            it_nest = graph_cycles_unique.erase(it_nest); // C++11 // remove it_nest from graph_cycles_unique
          }

          is_tds_cycle = true;

        } else { // no intersection found
          it_nest++;
        } 
      } // for  

      // nothing was merged with this cycle, so not a TDS cycle; 
      // therefore, remove from the list of TDS cycles 
      if (!is_tds_cycle) {
        //std::cout << "Removing " << tds_cycles_it->first << " from tds_cycles." << std::endl; // Test
        tds_cycles.erase(tds_cycles_it);          
      }          

    } else {
      std::cerr << "Error: unexpected item in the map. " << std::endl;  
    }

  } // for 

  //std::cout << std::endl; // Test 
  //std::cout << "> tds_cycles size: " << tds_cycles.size() 
  //  << " " << std::endl; // Test 
  //std::cout << "> graph_cycles_unique size: " << graph_cycles_unique.size() 
  //  << std::endl; // Test
}

//------------------------------------------------------------------------------

/**
 * Not used  
 */
template <typename TDSEdgeSet, typename EdgeSetIterator>
void generate_tds_cycles_recursive_2(TDSEdgeSet& tds_edge_set, 
  EdgeSetIterator it, EdgeSetIterator it_end) {
  if (it == it_end) {
    return;
  } else {
    // add edge
    auto find_edge_hash = tds_edge_set.find(it->first);
    if (find_edge_hash == tds_edge_set.end()) { 
      tds_edge_set.insert({it->first, std::forward_as_tuple
        (std::get<0>(it->second),std::get<1>(it->second))});
    } 
    generate_tds_cycles_recursive_2(tds_edge_set, ++it, it_end);
  }
}

template <typename Veretx, typename Edge, typename GraphNonLocalPropertiesUnique, typename TDSEdgeSet>
void generate_tds_cycles(GraphNonLocalPropertiesUnique& graph_cycles_unique, TDSEdgeSet& tds_edge_set) {
  /*for (auto it_a = graph_cycles_unique.begin(), it_b = graph_cycles_unique.begin(); it_a != graph_cycles_unique.end(); ++it_a) {
    if (it_b != graph_cycles_unique.end()) { 
      ++it_b;
         
      // process
      std::cout << "a" << it_a->first;
      std::cout << " -- b" << it_b->first;
 
    } else {
      break;  
    }   
    std::cout << std::endl;
  }*/ 

  //auto it = graph_cycles_unique.begin();

  for (auto i :  graph_cycles_unique) {
    //auto it = i.begin();
    generate_tds_cycles_recursive_2(tds_edge_set, i.second.begin(), i.second.end());
  }  

  // add the first item (edge by edge) to tds_edge_set and then ++it
  
  //graph_tds_cycles_unique.insert({it->first, it->second});

  //generate_tds_cycles_recursive(tds_edge_set, ++it, graph_cycles_unique.end());  
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Find all the unique cycles in the graph
 */
template <typename Vertex, typename Edge, 
  typename VertexNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesUnique>
void find_unique_cycles(
  VertexNonLocalPropertiesUnique& vertex_cycles_unique, 
  GraphNonLocalPropertiesUnique& graph_cycles_unique) {

  // Note: preserving edge order in a path is not required here   

  for (size_t i = 0; i < vertex_cycles_unique.size(); i++) { // vector
    for (auto& j : vertex_cycles_unique[i]) { // j is a map
      // j.first // path hash // Note: path hash is unique for a cycle
      // j.second // path map         
      auto find_path = graph_cycles_unique.find(j.first); 
      if (find_path == graph_cycles_unique.end()) {
        graph_cycles_unique.insert({j.first, j.second});
      } 
      //for (auto& k : j.second) { // k is a map
        // k.first // edge hash
        // std::get<0>(k.second), std::get<1>(k.second) // edge - s, t 
        // std::get<2>(k.second) // walk step
      //} // for      
    } // for  
  } // for 
}

template <typename Vertex, typename Edge,  
  typename VertexNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths>
void find_unique_cycles(
  VertexNonLocalPropertiesUnique& vertex_cycles_unique, 
  GraphNonLocalPropertiesUnique& graph_cycles_unique,
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths) {

  // Note: preserving edge order in a path is not required here   

  for (size_t i = 0; i < vertex_cycles_unique.size(); i++) { // vector
    for (auto& j : vertex_cycles_unique[i]) { // j is a map
      // j.first // path hash // Note: path hash is unique for a cycle
      // j.second // path map         
      auto find_path = graph_cycles_unique.find(j.first); 
      if (find_path == graph_cycles_unique.end()) {
        graph_cycles_unique.insert({j.first, j.second});
      }

      auto find_paths = graph_cycles_all_paths.find(j.first);
      graph_cycles_all_paths[j.first].push_back(j.second);        
         
      //for (auto& k : j.second) { // k is a map
        // k.first // edge hash
        // std::get<0>(k.second), std::get<1>(k.second) // edge - s, t 
        // std::get<2>(k.second) // walk step
      //} // for      
    } // for  
  } // for 
}

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique,
  typename OrderedPath, typename OrderedPathEdges, typename UnorderedPathEdges>
void find_cycles_recursive(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges,
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,  
  Vertex source_vertex, Vertex v, OrderedPath walk_history, 
  OrderedPathEdges walk_history_edges, 
  UnorderedPathEdges walk_history_edges_unordered, size_t r) {

    // v is the current vertex    
    // r must be greater than 2 for it to be considered a cycle    
    // do not forward a reference to walk_history, make a new copy

    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];

      /*std::cout << " -> " << v_nbr << " <" << v << "> : " ;
      for (auto& i : walk_history) {
        std::cout << "(" << i.first << ", " << i.second << ") ";
      }
      std::cout << std::endl;*/

      if (vertex_degree[v_nbr] < 2) {
        //continue;
      }
      
      if (v_nbr == source_vertex && r >= 3) { 
        // cycle found
         
        //-std::cout << source_vertex << " : ";   
        //-for (auto& i : walk_history) {
        //-  std::cout << "(" << i.first << ", " << i.second << ") ";
        //-}    

//        OrderedPathEdges new_walk_history_edges(walk_history_edges);
//        new_walk_history_edges.push_back(std::forward_as_tuple(v, v_nbr)); 
   
//        for (auto& i : new_walk_history_edges) {
//          std::cout << "(" << std::get<0>(i) << " -- " << std::get<1>(i) << "), ";  
//        }      

        // tuple - vector, size_t 
        OrderedPathEdges new_walk_history_edges(walk_history_edges);
        std::get<0>(new_walk_history_edges).push_back(std::forward_as_tuple(v, v_nbr));
        //std::get<1>(new_walk_history_edges) = 
        //  edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);
        update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges));  
 
        //-for (auto& i : std::get<0>(new_walk_history_edges)) {   
        //-  std::cout << "(" << std::get<0>(i) << " -- " << std::get<1>(i) << "), ";  
        //-}
        //std::cout << ": " << std::get<1>(new_walk_history_edges);   
        //-std::cout << ": " << std::get<1>(new_walk_history_edges).to_ullong(); 

        //-std::cout << " : Found Cycle." << std::endl;

        // tuple - map, size_t  
        UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);

        Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

        auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
        if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
          std::get<0>(new_walk_history_edges_unordered).
            insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered));
        }

        // add path to vertex_cycles
 
        //----------------------------------------------------------------------       
  
        // vertex_cycles_unique - vector based walk history edges     
 
        auto find_path = vertex_cycles[source_vertex].find(std::get<1>(new_walk_history_edges).to_ullong());  
        if (find_path ==  vertex_cycles[source_vertex].end()) {    
          vertex_cycles[source_vertex].
            insert({std::get<1>(new_walk_history_edges).to_ullong(), std::get<0>(new_walk_history_edges)});
        } //else {
          //std::cerr << "Error: unexpected item in the map." << std::endl; 
        //}

        //----------------------------------------------------------------------

        // vertex_cycles_unique - map based walk history edges    
        //walk_history_edges_unordered

        auto find_path_3 = vertex_cycles_unique[source_vertex].find(std::get<1>(new_walk_history_edges_unordered).to_ullong());
        if (find_path_3 ==  vertex_cycles_unique[source_vertex].end()) {  
          vertex_cycles_unique[source_vertex].
            insert({std::get<1>(new_walk_history_edges_unordered).to_ullong(), std::get<0>(new_walk_history_edges_unordered)});    
        }

        //return; 
      } else if (v_nbr == source_vertex && r < 3) {
        // path length is invalid, ignore
        //return; 
        //continue;        
      } else if (v_nbr != source_vertex && r > 0) { 
        // forward  
        auto find_v_nbr = walk_history.find(v_nbr); 
        if (find_v_nbr == walk_history.end()) {
          //walk_history.insert({v_nbr, r});           

          // vertex -     
          OrderedPath new_walk_history(walk_history);              
          new_walk_history.insert({v_nbr, r + 1}); 
 
          //--------------------------------------------------------------------

//          OrderedPathEdges new_walk_history_edges(walk_history_edges);  
//          new_walk_history_edges.push_back(std::forward_as_tuple(v, v_nbr));

          // tuple - vector, size_t
          OrderedPathEdges new_walk_history_edges(walk_history_edges);
          std::get<0>(new_walk_history_edges).push_back(std::forward_as_tuple(v, v_nbr));
          //std::get<1>(new_walk_history_edges) = 
          //  edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges)); 

          //-------------------------------------------------------------------- 

          // tuple - map, size_t 
          UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
             
          Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

          auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
          if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
            std::get<0>(new_walk_history_edges_unordered).
              insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});  
            update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered)); 
          }      

          //--------------------------------------------------------------------  

          //std::cout << " : " << v_nbr << " -> " << std::endl;
  
          // forward
          find_cycles_recursive<Vertex, Edge, VertexList, EdgeList,
            VertexNonLocalProperties, VertexNonLocalPropertiesUnique,
            OrderedPath, OrderedPathEdges, UnorderedPathEdges>
            (vertices, vertex_count, vertex_degree, edges, vertex_cycles,
            vertex_cycles_unique,  
            source_vertex, v_nbr, new_walk_history, new_walk_history_edges, 
            new_walk_history_edges_unordered,
            r + 1);

        } else { 
          //return; // repeated vertex, ignore   
          //continue;
        }
      } // if
 
    } // for   
} 

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique>
void find_cycles_parallel(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges, 
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique) {

  { 
  #pragma omp parallel for
  for (Vertex v = 0; v < vertex_count; v++) {
    //Vertex v = 0; // Test   

    if (vertex_degree[v] < 2) {
      //continue;
    }

    /*std::cout << v << " : ";
    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e]; 
      std::cout << v_nbr << ", "; 
    }
    std::cout << std::endl;*/
    //continue;

    size_t r = 1; // walk step, initialized to 1

    // map based walk history of vertices with step order
    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0); 
  
    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map. " << std::endl;
    }

    // vector based walk history of edges
    //typedef std::vector<std::tuple<Vertex, Vertex>> OrderedPathEdges; 
    typedef std::tuple<std::vector<std::tuple<Vertex, Vertex>>, VertexBitSet> OrderedPathEdges;
    //OrderedPathEdges walk_history_edges(0);     
 
    //walk_history_edges.push_back(std::forward_as_tuple(v, v));
    
    std::vector<std::tuple<Vertex, Vertex>> wh_vector(0);
    //wh_vector.push_back(std::forward_as_tuple(v, v)); // not required 
    VertexBitSet vertex_bitset;
    //update_vertex_bitset(v, v, vertex_bitset); // not required
    OrderedPathEdges walk_history_edges(wh_vector, vertex_bitset);

    // map based walk history of edges with step order
    typedef std::tuple< std::unordered_map< Edge, std::tuple<Vertex, Vertex, size_t> >, VertexBitSet> 
      UnorderedPathEdges;    
    std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> wh_unordered_map(0);
    UnorderedPathEdges walk_history_edges_unordered(wh_unordered_map, vertex_bitset); 

    find_cycles_recursive<Vertex, Edge, VertexList, EdgeList, 
      VertexNonLocalProperties, VertexNonLocalPropertiesUnique, 
      OrderedPath, OrderedPathEdges, UnorderedPathEdges>
      (vertices, vertex_count, vertex_degree, edges, vertex_cycles, 
      vertex_cycles_unique, 
      v, v, walk_history, walk_history_edges, walk_history_edges_unordered, r); 
     
  } // for
  } 
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Template prototype generation - up to k-edit distance
 */

template <typename Edge, typename TemplateGraph, typename O, typename S> 
void generate_edge_combinations(TemplateGraph& input_template, 
  size_t k_edit_distance, O& edge_combinations, S& edge_set_combinations) {
  //typedef typename TemplateGraph::Edge_t Edge_t; 

  // setup input
  typedef std::vector<Edge> NCollection;
  NCollection n_items(0); 

  for (auto& e : input_template.edgelist_unique) {
    auto edge_uint = e.first;
    n_items.push_back(static_cast<Edge>(edge_uint));
  }

  std::cout << "Sorting unique edges ..." << std::endl;
  std::stable_sort(n_items.begin(), n_items.end(),
    [](const Edge& a, const Edge& b) -> bool {
         return a < b;
       });
  
  // Test
  //for (auto& i : n_items) {
  //  std::cout << i << std::endl;
  //}
  // Test

  // ~setup input  

  std::cout << "Calculating edge combinaions ..." << std::endl;

  //size_t n_collection = 7; // Test
  size_t n_collection = n_items.size();
  size_t k_selection = k_edit_distance;

  assert(n_collection >= k_selection);

  //size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
  //  (static_cast<size_t>(n_collection.size()), k_selection));
  size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
    (n_collection, k_selection));

  std::cout << "n : " << n_collection << ", k : " << k_selection
    << ", nCk : " << nCk << std::endl;

  std::cout << "Total edge combinations: " << nCk << std::endl;

  //typedef std::vector<std::vector<size_t>> O; // indices // TODO
  //typedef std::vector<std::unordered_set<Edge_t>> S; // edge_uint // TODO

  //O edge_combinations(0);
  //S edge_set_combinations(0);

  std::cout << "Generating edge combinaions ..." << std::endl;

  // k is the edit distance in this context

  // generate all combinations of uinque unsigned integers
  prunejuice::pattern::combination::generate_combinations_recursive
    <NCollection, Edge, O, S>
    (n_collection, k_selection, n_items, edge_combinations,
    edge_set_combinations); // sequential routine

  assert(nCk == edge_combinations.size());
  assert(nCk == edge_set_combinations.size());

  // Test
  /*for (auto& e : edge_combinations) {
    for (auto& j : e) {
      std::cout << n_items[j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << "--- --- --- ---" << std::endl;
  for (auto& e : edge_set_combinations) {
    for (auto& j : e) {
      std::cout << j << ", ";
    }
    std::cout << std::endl;
  }*/
  // Test 
  
}

template <typename TemplateGraph, typename TemplateGraphVector, typename S,
  typename TemplateConstraint, typename TemplateConstraintVector>
void generate_k_edit_distance_prototypes_parallel(TemplateGraph& input_template, 
  TemplateGraphVector& prototypes, S& edge_set_combinations, 
  TemplateConstraintVector& template_constraints) {

  {
  #pragma omp parallel for 
  for (size_t i = 0; i < edge_set_combinations.size(); i++) {
    //size_t thread_ID = omp_get_thread_num();
    //std::cout << "Thread# " << thread_ID << " "
    //  << " Combination# " << i << std::endl;

    typedef typename TemplateGraph::Vertex_t Vertex_t; 
    typedef typename TemplateGraph::Edge_t Edge_t;
    
    typedef typename TemplateGraph::EdgeListTuple_t EdgeListTuple_t;
    typedef typename TemplateGraph::EdgeListMap_t EdgeListMap_t;
    typedef typename TemplateGraph::EdgeList_t EdgeList_t;    
    typedef typename TemplateGraph::VertexList_t VertexList_t;

    typedef typename TemplateGraph::EdgeListTupleVector_t EdgeListTupleVector_t;
    typedef std::vector<EdgeListTupleVector_t> EdgeListTupleVectors_t; 

    typedef std::unordered_map<Edge_t, std::tuple<Vertex_t, Vertex_t, size_t>> EdgeSet_t; // TODO: ?
    typedef std::vector<EdgeSet_t> EdgeSetVector_t; // TODO: ?

    EdgeListTuple_t new_edge_list = TemplateGraph::filter_edge_list
      //<Edge_t, EdgeListTuple_t, std::unordered_set<Edge_t>> // TODO: define type
      (input_template.edgelist, edge_set_combinations[i]);

    //std::cout << new_edge_list.size() << std::endl;

    // create a new template graph - a prototype
    TemplateGraph new_template(new_edge_list, input_template.vertex_data, 
      input_template.max_vertex);  
    // Important: must use the max_vertex from input_template.
    // Here the number of vertices is fixed. 

    if (prunejuice::pattern::graph_algorithm::is_connected_component
      <TemplateGraph>(new_template)) {
      /*std::cout << "Connected component." << std::endl;
      std::cout << "vertex_count " << new_template.vertex_count << std::endl;  
      std::cout << "max_vertex " << new_template.max_vertex << std::endl;
      std::cout << "max_degree " << new_template.max_degree << std::endl;
      std::cout << "edge_count " << new_template.edge_count << std::endl;
      std::cout << "vertices.size() " << new_template.vertices.size() << std::endl;
      std::cout << "vertex_degree.size() " << new_template.vertex_degree.size() << std::endl;
      std::cout << "vertex_data.size() " << new_template.vertex_data.size() << std::endl;
      std::cout << "edges.size() " << new_template.edges.size() << std::endl;
      std::cout << "edgelist.size() " << new_template.edgelist.size() << std::endl;       
      std::cout << "edgelist_unique.size() " << new_template.edgelist_unique.size() << std::endl 
        << std::endl;*/   
      
      // identify template constraints for the prototype 
      /*for (size_t j = 0; j < template_constraints.size(); j++) {
        if(TemplateGraph::is_subset_edgeset
          (new_template.edgelist_unique, template_constraints[j].edgeset)) {
          assert(j == template_constraints[j].constraint_ID);
          new_template.template_constraints.
            push_back(template_constraints[j].constraint_ID);
        } 
      } // for*/

      // Note: only cycle constraints are generated form the input template,
      // path and TDS constraints are generated for each prototype 
      // (in parallel). template_constraints only contains cycle constraints.

      // identify cycle constraints for the prototype  
      typedef std::unordered_map<size_t, EdgeSet_t> GraphNonLocalPropertiesUnique;
      GraphNonLocalPropertiesUnique prototype_cycles_unique(0);
   
      typedef std::unordered_map<size_t, EdgeSetVector_t> GraphNonLocalPropertiesAllPaths;
      GraphNonLocalPropertiesAllPaths prototype_cycles_all_paths(0);

      GraphNonLocalPropertiesUnique prototype_tds_cycles_unique(0);

      for (size_t j = 0; j < template_constraints.size(); j++) {
        if (template_constraints[j].constraint_type == TemplateConstraint::CYCLE) { 
          if(TemplateGraph::is_subset_edgeset
            (new_template.edgelist_unique, template_constraints[j].edgeset)) {
            assert(j == template_constraints[j].constraint_ID);
            new_template.template_constraints.
              push_back(template_constraints[j].constraint_ID);

            {
            auto find_edgeset_hash = prototype_cycles_unique.
              find(template_constraints[j].edgeset_hash); 
            if (find_edgeset_hash ==  prototype_cycles_unique.end()) {
              prototype_cycles_unique.insert
                ({template_constraints[j].edgeset_hash, 
                template_constraints[j].edgeset});                   
            } else {
              std::cerr << "Error: unexpected item in the map." << std::endl;
              //return; // TODO: graceful exit; return causes OpenMP error   
            }
            }    

            { 
            auto find_edgeset_hash = prototype_cycles_all_paths.
              find(template_constraints[j].edgeset_hash); 
            if (find_edgeset_hash ==  prototype_cycles_all_paths.end()) {
              prototype_cycles_all_paths.insert
                ({template_constraints[j].edgeset_hash, 
                template_constraints[j].edgeset_vector});                   
            } else {
              std::cerr << "Error: unexpected item in the map." << std::endl;
              //return; // TODO: graceful exit; return causes OpenMP error   
            }
            }           
 
          }
        }
      } 

      // TODO: remove template_constraints from TemplateGraph ?     

      // generate prototype constraints 
       
      generate_tds_cycles<Vertex_t, Edge_t>(prototype_cycles_unique, 
        prototype_tds_cycles_unique);
  
      EdgeListTupleVectors_t prototype_tds_cycles_all_paths(0);

      generate_tds_cycle_constraints<Vertex_t, Edge_t, EdgeSet_t,
        GraphNonLocalPropertiesUnique, EdgeListTuple_t, EdgeListTupleVector_t,
        EdgeListTupleVectors_t>(prototype_tds_cycles_unique, 
        prototype_tds_cycles_all_paths);

      EdgeListTuple_t prototype_enumeration_unique = 
        prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
          (new_template, static_cast<Vertex_t>(0)); 
      
      TemplateConstraintVector protoype_constraints(0);
     
      //populate_template_constraints<Vertex_t, Edge_t, EdgeSet_t, 
      //  TemplateConstraint>(prototype_cycles_unique, prototype_cycles_all_paths, 
      //  prototype_tds_cycles_unique, // TODO: prototype_tds_cycles_all_paths
      //  protoype_constraints);               

      populate_template_constraints<Vertex_t, Edge_t, EdgeSet_t,
        TemplateConstraint, GraphNonLocalPropertiesUnique, 
        GraphNonLocalPropertiesAllPaths, EdgeListTuple_t, 
        EdgeListTupleVector_t, EdgeListTupleVectors_t,
        TemplateConstraintVector>(prototype_cycles_unique, 
        prototype_cycles_all_paths, prototype_tds_cycles_unique, 
        prototype_tds_cycles_all_paths, prototype_enumeration_unique,
        protoype_constraints);

      new_template.template_nonlocal_constraints = std::move(protoype_constraints); // Important:

      // add the template prototype to prototypes 
      {
      #pragma omp critical(prototypes)
      {
        //size_t thread_ID = omp_get_thread_num();
        //std::cout << "Thread# " << thread_ID  
        //  << " entered the critical section" << std::endl;

        prototypes.push_back(new_template); // TODO: add prototype_template_constraints

      } // #pragma omp citical 
      } 

    } else {
      //std::cout << "Not a connected component." << std::endl;
    }  

  } // for
  } // #pragma omp parallel for

  if (prototypes.size() < 1) {
    std::cout << "No prototype was generated" << std::endl;
    //std::cerr << "Error: No prototype was generated. "
    //  << "Aborting ... " << std::endl;
    //return -1;
  } else {
    std::cout << "Number of prototypes (each one is a single connected component): " 
      << prototypes.size() << std::endl;  
  } 

  // Test

  /*size_t nlc_count = 0; // Test

  for (size_t i = 0; i < prototypes.size(); i++) {
    //std::cout << "P " << i << " : ";
    //std::cout << prototypes[i].edgelist_unique.size() << std::endl;
    //std::cout << prototypes[i].template_constraints.size() << std::endl;
    if (prototypes[i].template_constraints.size() > 0) {          
      //~std::cout << "P " << i << " : ";
      //std::cout << prototypes[i].template_constraints.size() << std::endl; 
      //for (auto& e : prototypes[i].edgelist) {
      //  std::cout << "(" << std::get<0>(e) << " - " << std::get<1>(e) << "), ";  
      //}
      for (auto& c : prototypes[i].template_constraints) {
        std::cout << c << ", ";
        std::cout << template_constraints[c].constraint_ID << ", "; 
        std::cout << template_constraints[c].edgeset_vector.size() << std::endl;   
      } 
      //~std::cout << std::endl;
      nlc_count++;
    }
    //std::cout << std::endl; 
  }

  std::cout << "Number of prototypes with NLCs: " << nlc_count << std::endl;*/
  // Test
  
}

template <typename TemplateGraph, typename TemplateGraphVector, 
  typename TemplateConstraint, typename TemplateConstraintVector>
void generate_k_edit_distance_prototypes
  (TemplateGraph& input_template, size_t k_edit_distance, 
  TemplateGraphVector& prototypes, 
  TemplateConstraintVector& template_constraints) {  

  std::cout << "Generating edit distance " << k_edit_distance 
    << " template prototypes ..." << std::endl;

  typedef typename TemplateGraph::Edge_t Edge_t;
  typedef std::vector<std::vector<size_t>> O; // indices // TODO
  typedef std::vector<std::unordered_set<Edge_t>> S; // edge_uint // TODO: define this type in template_graph class?

  O edge_combinations(0);
  S edge_set_combinations(0);

  generate_edge_combinations<Edge_t, TemplateGraph, O, S>(input_template, 
    k_edit_distance, edge_combinations, edge_set_combinations);

  generate_k_edit_distance_prototypes_parallel<TemplateGraph, 
    TemplateGraphVector, S, TemplateConstraint, TemplateConstraintVector>
    (input_template, prototypes, edge_set_combinations, template_constraints);

  std::cout << std::endl;
}

template <typename TemplateGraph, typename TemplateGraphVector, 
  typename TemplateGraphVectors, typename TemplateConstraint,
  typename TemplateConstraintVector>
void generate_up_to_k_edit_distance_prototypes 
  (TemplateGraph& input_template, size_t max_k_edit_distance, 
  TemplateGraphVectors& k_edit_distance_prototypes,
  TemplateConstraintVector& template_constraints) {

   assert(k_edit_distance_prototypes.size() > 0); 
   // the first element in k_edit_distance_prototypes is the input template
  
   // TODO: make this loop parallel? // Note: this step can be sequential
   for (size_t k = 1; k <= max_k_edit_distance; k++) {
     generate_k_edit_distance_prototypes<TemplateGraph, TemplateGraphVector, 
       TemplateConstraint, TemplateConstraintVector>
       (input_template, k, k_edit_distance_prototypes[k], template_constraints);      
   } 
}

////////////////////////////////////////////////////////////////////////////////

/**
 * File operations
 */

template <typename TemplateGraphVectors, typename FileUtilities>
void create_template_prototype_directories
  (TemplateGraphVectors& k_edit_distance_prototypes, 
  std::string prototype_dir_name) {
 
  for (size_t k = 0;  k < k_edit_distance_prototypes.size(); k++) {
    if (k_edit_distance_prototypes[k].size() < 1) {      
      continue;
    }
  
    {
    //#pragma omp parallel for // TODO: ?
    for (size_t p = 0; p < k_edit_distance_prototypes[k].size(); p++) {
      //std::cout << "Creating directory " << k << "_" << p << std::endl; 
      FileUtilities::create_template_directories(prototype_dir_name,  
        std::to_string(k) + "_" + std::to_string(p));        
    } // for
    }

    FileUtilities::create_graph_directory(prototype_dir_name, 
      std::to_string(k) + "_" + 
      std::string(FileUtilities::PRUNED_GRAPH_UNION_DIR_NAME));
    
  } // for 
}

//------------------------------------------------------------------------------

template <typename TemplateGraph, typename FileUtilities>
void write_template_edge_file(size_t k, size_t p, 
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name, 
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR), 
    std::string(FileUtilities::PATTERN_EDGE_FILE));       

  for (size_t e = 0; e < template_graph.edgelist.size(); e++) {
    file_utilities.output_file 
      << std::get<0>(template_graph.edgelist[e]) << " "
      << std::get<1>(template_graph.edgelist[e]) << "\n";
  }  
}

template <typename TemplateGraph, typename FileUtilities>
void write_template_edge_data_file(size_t k, size_t p, 
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name, 
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR), 
    std::string(FileUtilities::PATTERN_EDGE_DATA_FILE));       

  for (size_t e = 0; e < template_graph.edgelist.size(); e++) {
    file_utilities.output_file 
      << std::get<0>(template_graph.edgelist[e]) << " "
      << std::get<1>(template_graph.edgelist[e]) << " "
      << 1 << " " << 1 << "\n"; // TODO: dummy
  }  
}

// TODO: read vertex data from file
template <typename TemplateGraph, typename FileUtilities>
void write_template_vertex_data_file(size_t k, size_t p,
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name,
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR),
    std::string(FileUtilities::PATTERN_VERTEX_DATA_FILE));

  typedef typename TemplateGraph::Vertex_t Vertex_t;

  for (Vertex_t v = 0; v < template_graph.vertex_data.size(); v++) {
    file_utilities.output_file
      << v << " "
      << template_graph.vertex_data[v] << "\n";
  }
}

// TODO: read from file or compute
template <typename TemplateGraph, typename FileUtilities>
void write_template_stat_file(size_t k, size_t p,
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name,
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR),
    std::string(FileUtilities::PATTERN_STAT_FILE));

    file_utilities.output_file << "diameter : 10" << "\n";  // TODO: dummy 
}

// TODO:
template <typename TemplateGraph, typename FileUtilities>
void write_template_lc_file(size_t k, size_t p,
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name,
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR),
    std::string(FileUtilities::PATTERN_LC_FILE));

    file_utilities.output_file << "1111" << "\n";  // TODO: dummy 
}

//------------------------------------------------------------------------------

// write_template_nlc_file helper functions

template <typename Vertex, typename EdgeSet, typename OrderedPathVertices>
void edgeset_to_ordered_path_vertices(EdgeSet& edgeset, 
  OrderedPathVertices& ordered_path_vertices) {
  
  for (auto it = edgeset.begin(); it != edgeset.end(); ++it) {
    //Edge edge_hash = it->first; //ignore
    Vertex s = std::get<0>(it->second);
    Vertex t = std::get<1>(it->second);
    size_t hop = std::get<2>(it->second);
  
    if (hop == 1) {
      ordered_path_vertices[0] = s;
      ordered_path_vertices[hop] = t;
    } else {
      ordered_path_vertices[hop] = t;
    }   
  } 
}

// TODO: for now, it assumes unique labels
template <typename TemplateGraph, typename TemplateConstraint, 
  typename FileUtilities>
void write_cycle_constraint(TemplateConstraint& template_constraint, 
  FileUtilities& file_utilities) {
  assert(template_constraint.constraint_type == TemplateConstraint::CYCLE);

  typedef typename TemplateGraph::Vertex_t Vertex_t;  
  typedef std::vector<Vertex_t> OrderedPathVertices;

  for (size_t i = 0; i < template_constraint.edgeset_vector.size(); i++) {

    std::string path_string = "";

    size_t path_length = template_constraint.edgeset_vector[i].size() + 1; // Important:   
    OrderedPathVertices ordered_path_vertices(path_length);

    edgeset_to_ordered_path_vertices<Vertex_t>(template_constraint.edgeset_vector[i], 
      ordered_path_vertices);     
    assert(ordered_path_vertices[0] == 
      ordered_path_vertices[ordered_path_vertices.size() - 1]); 

    // path
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {  
      path_string += std::to_string(ordered_path_vertices[j]) + " ";   
    }

    path_string += ": ";

    // path indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the path indices accordingly.  
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      if (j == ordered_path_vertices.size() - 1) {
        path_string += std::to_string(0) + " ";  
      } else {
        path_string += std::to_string(j) + " ";
      }
    }
 
    path_string += ": ";

    // aggregation indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the aggregation vertices accordingly. 
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {    
      path_string += std::to_string(0) + " ";
    }

    // is cyclic : is TDS : do invoke LCC : constraint_ID 
    // TODO: use ConstraintType, set the parameters dynamically 
    // - is cyclic is 1 when source == destination, no need to ack the source
    // - is TDS is 1 if a cycle or a path has duplicate intermediate vertices 
    // (other than the source)
    // TODO: constraint_ID, remove ? 
    std::string path_parameter_string = ": 1 : 0 : 1 : " + 
      std::to_string(template_constraint.constraint_ID);  
    //path_string += ": 1 : 0 : 1";
    path_string += path_parameter_string; 

    // write to the file
    
    file_utilities.output_file << path_string << "\n"; 
    
  } // for        
}

template <typename Vertex, typename EdgeListTuple, typename OrderedPathVertices>
void edgelist_to_ordered_path_vertices(EdgeListTuple& edgelist, 
  OrderedPathVertices& ordered_path_vertices, 
  OrderedPathVertices& ordered_path_indices) {

  std::unordered_map<Vertex, size_t> vertex_path_index_map(0);  
  Vertex current_vertex = 0;
 
  for (size_t i = 0; i < edgelist.size(); i++) {
    Vertex s = std::get<0>(edgelist[i]);
    Vertex t = std::get<1>(edgelist[i]);

    if (i == 0) {
      ordered_path_vertices.push_back(s);

      auto find_s = vertex_path_index_map.find(s);
      if (find_s == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({s, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1); 
      }

      ordered_path_vertices.push_back(t);
      current_vertex = t; 

      auto find_t = vertex_path_index_map.find(t);
      if (find_t == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({t, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
      }

      continue;
    } else {
      if (s != current_vertex) {
        ordered_path_vertices.push_back(s);         
       
        auto find_s = vertex_path_index_map.find(s);
        if (find_s == vertex_path_index_map.end()) {
          vertex_path_index_map.insert({s, ordered_path_vertices.size() - 1});

          ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
        } else {
          ordered_path_indices.push_back(find_s->second);  
        }
           
      }
      ordered_path_vertices.push_back(t);
      current_vertex = t;

      auto find_t = vertex_path_index_map.find(t);
      if (find_t == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({t, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
      } else {
        ordered_path_indices.push_back(find_t->second); 
      } 

      continue;
    } 
  } // for

  assert(ordered_path_vertices.size() == ordered_path_indices.size()); 
}  

template <typename TemplateGraph, typename TemplateConstraint, 
  typename FileUtilities>
void write_tds_constraint(TemplateConstraint& template_constraint, 
  FileUtilities& file_utilities) {
  //assert(template_constraint.constraint_type == TemplateConstraint::TDS);
  if (!(template_constraint.constraint_type == TemplateConstraint::TDS || 
    template_constraint.constraint_type == TemplateConstraint::ENUMERATION)) {
    assert(false);
  }

  typedef typename TemplateGraph::Vertex_t Vertex_t;  
  typedef std::vector<Vertex_t> OrderedPathVertices;

  for (size_t i = 0; i < template_constraint.edgelist_vector.size(); i++) {

    std::string path_string = "";

    OrderedPathVertices ordered_path_vertices(0);
    OrderedPathVertices ordered_path_indices(0);

    edgelist_to_ordered_path_vertices<Vertex_t>
      (template_constraint.edgelist_vector[i], ordered_path_vertices, 
      ordered_path_indices);

    // path
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      path_string += std::to_string(ordered_path_vertices[j]) + " ";
    }

    path_string += ": ";

    // path indices
    for (size_t j = 0; j < ordered_path_indices.size(); j++) {
      path_string += std::to_string(ordered_path_indices[j]) + " ";
    }

    path_string += ": ";

    // aggregation indices // TODO: set the aggregation vertices 
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      path_string += std::to_string(1) + " ";
    }

    path_string += ": ";

    bool is_cyclic = ordered_path_vertices[0] == 
      ordered_path_vertices[ordered_path_vertices.size() - 1] ? true : false;   

    path_string += std::to_string(is_cyclic) + " "; 

    path_string += ": ";

    bool is_TDS = true; // TODO: constraint type ?

    path_string += std::to_string(is_TDS) + " ";

    path_string += ": "; 

    bool do_lcc = true; // TODO: ?

    path_string += std::to_string(do_lcc) + " "; 

    path_string += ": ";

    path_string +=  std::to_string(template_constraint.constraint_ID); // TODO: constraint_ID, remove ?
     
    /*size_t path_length = template_constraint.edgelist_vector[i].size() + 1; // Important:   
    OrderedPathVertices ordered_path_vertices(path_length);

    edgeset_to_ordered_path_vertices<Vertex_t>(template_constraint.edgeset_vector[i], 
      ordered_path_vertices);     
    assert(ordered_path_vertices[0] == 
      ordered_path_vertices[ordered_path_vertices.size() - 1]); 

    // path
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {  
      path_string += std::to_string(ordered_path_vertices[j]) + " ";   
    }

    path_string += ": ";

    // path indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the path indices accordingly.  
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      if (j == ordered_path_vertices.size() - 1) {
        path_string += std::to_string(0) + " ";  
      } else {
        path_string += std::to_string(j) + " ";
      }
    }
 
    path_string += ": ";

    // aggregation indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the aggregation vertices accordingly. 
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {    
      path_string += std::to_string(0) + " ";
    }

    // is cyclic : is TDS : do invoke LCC : constraint_ID 
    // TODO: use ConstraintType, set the parameters dynamically 
    // - is cyclic is 1 when source == destination, no need to ack the source
    // - is TDS is 1 if a cycle or a path has duplicate intermediate vertices 
    // (other than the source)
    // TODO: constraint_ID, remove ? 
    std::string path_parameter_string = ": 1 : 0 : 1 : " + 
      std::to_string(template_constraint.constraint_ID);  
    //path_string += ": 1 : 0 : 1";
    path_string += path_parameter_string;*/ 

    // write to the file
    
    file_utilities.output_file << path_string << "\n"; 
    
  } // for        
}

//------------------------------------------------------------------------------
 
template <typename TemplateGraph, typename TemplateConstraint, 
  typename TemplateConstraintVector, typename FileUtilities>
void write_template_nlc_file(size_t k, size_t p,
  TemplateGraph& template_graph, std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name,
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR),
    std::string(FileUtilities::PATTERN_NLC_FILE));

    for (size_t i = 0; i < template_graph.template_nonlocal_constraints.size(); 
      i++) {
      //std::cout <<  
      //  template_graph.template_nonlocal_constraints[i].constraint_type 
      //  << ", ";  
      if (template_graph.template_nonlocal_constraints[i].
        constraint_type == TemplateConstraint::CYCLE) {
        write_cycle_constraint<TemplateGraph, TemplateConstraint, FileUtilities>
          (template_graph.template_nonlocal_constraints[i], 
          file_utilities); 
      } else if (template_graph.template_nonlocal_constraints[i].
        constraint_type == TemplateConstraint::TDS) {
        write_tds_constraint<TemplateGraph, TemplateConstraint, FileUtilities>
          (template_graph.template_nonlocal_constraints[i],
          file_utilities); 
      } else if (template_graph.template_nonlocal_constraints[i].
        constraint_type == TemplateConstraint::ENUMERATION) {
        write_tds_constraint<TemplateGraph, TemplateConstraint, FileUtilities>
          (template_graph.template_nonlocal_constraints[i],
          file_utilities);
      }
           
    } // for
    //std::cout << std::endl; 
    //
}

// Note: not used
template <typename TemplateGraph, typename TemplateConstraint, 
  typename TemplateConstraintVector, typename FileUtilities>
void write_template_nlc_file(size_t k, size_t p,
  TemplateGraph& template_graph, TemplateConstraintVector& template_constraints,
  std::string prototype_dir_name) {

  FileUtilities file_utilities(prototype_dir_name,
    std::to_string(k) + "_" + std::to_string(p),
    std::string(FileUtilities::PATTERN_DIR),
    std::string(FileUtilities::PATTERN_NLC_FILE));

    for (size_t i = 0; i < template_graph.template_constraints.size(); 
      i++) {
      //std::cout <<  
      //  template_constraints[template_graph.template_constraints[i]].constraint_type 
      //  << ", ";  
      if (template_constraints[template_graph.template_constraints[i]].
        constraint_type == TemplateConstraint::CYCLE) {
        write_cycle_constraint<TemplateGraph, TemplateConstraint, FileUtilities>
          (template_constraints[template_graph.template_constraints[i]], 
          file_utilities); 
      } 
    } // for
    //std::cout << std::endl; 
    //
}

//------------------------------------------------------------------------------

template <typename TemplateGraph, typename TemplateGraphVectors, 
  typename TemplateConstraint, typename TemplateConstraintVector, 
  typename FileUtilities>
void write_template_prototype_input_files
  (TemplateGraphVectors& k_edit_distance_prototypes,
  //TemplateConstraintVector& template_constraints, 
  std::string prototype_dir_name) {

  for (size_t k = 0;  k < k_edit_distance_prototypes.size(); k++) {
    if (k_edit_distance_prototypes[k].size() < 1) {
      continue;
    }
   
    {
    #pragma omp parallel for
    for (size_t p = 0; p < k_edit_distance_prototypes[k].size(); p++) {
      //std::cout << "Writing to file " << k << "_" << p << std::endl;
  
      // write edge 
      write_template_edge_file<TemplateGraph, FileUtilities>(k,  p, 
        k_edit_distance_prototypes[k][p], prototype_dir_name);

      // write edge_data
      write_template_edge_data_file<TemplateGraph, FileUtilities>(k,  p,
        k_edit_distance_prototypes[k][p], prototype_dir_name);
     
      // write vertex_data 
      write_template_vertex_data_file<TemplateGraph, FileUtilities>(k,  p,
        k_edit_distance_prototypes[k][p], prototype_dir_name);  

      // write stat
      write_template_stat_file<TemplateGraph, FileUtilities>(k,  p,
        k_edit_distance_prototypes[k][p], prototype_dir_name);
        
      // write local_constraint
      write_template_lc_file<TemplateGraph, FileUtilities>(k,  p,
        k_edit_distance_prototypes[k][p], prototype_dir_name);      
 
      // write nonlocal_constraint
      //write_template_nlc_file<TemplateGraph, TemplateConstraint, 
      //  TemplateConstraintVector, FileUtilities>(k,  p,
      //  k_edit_distance_prototypes[k][p], template_constraints, 
      //  prototype_dir_name);

      write_template_nlc_file<TemplateGraph, TemplateConstraint, 
        TemplateConstraintVector, FileUtilities>(k,  p,
        k_edit_distance_prototypes[k][p], prototype_dir_name);

    } // for 
    } // #pragma omp parallel for

  } // for
}  

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

/**
 * DFS walk on the TDS constraint
 */
template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename OrderedPath, typename Visited>
void dfs_recursive(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges,
  Vertex source_vertex, Vertex v, OrderedPath& walk_history, 
  Visited& visited, size_t r) {

  visited[static_cast<size_t>(v)] = 1;  
  //std::cout << v << ", "; // Test 

  for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
    Vertex v_nbr = edges[e];

    if (visited[static_cast<size_t>(v_nbr)] == 0) { 
      std::cout << "(" << v << " - " << v_nbr << "), "; // Test
      dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
        (vertices, vertex_count, vertex_degree, edges, source_vertex, v_nbr, 
        walk_history, visited, r + 1);

    } else {
      std::cout << "[" << v << " - " << v_nbr << "], "; // Test
    }
    
  } // for 
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList>
void dfs_parallel(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges) {
  
  {
  //#pragma omp parallel for schedule(static, 1)
  for (Vertex v = 0; v < vertex_count; v++) 
  {

    std::vector<uint8_t> visited(vertex_count);
    //for (auto& i : visited) {
    //  i = 0;
    //}

    //Vertex v = 0;

    size_t r = 1; // walk step, initialized to 1

    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0);      

    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map." << std::endl;
    }

    dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
      (vertices, vertex_count, vertex_degree, edges, v, v, walk_history, 
      visited, r);
     
    std::cout << std::endl; // Test
  } // for
  } 
} 
 
template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList>
void dfs_single_source(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges, Vertex v) {
  
  {
  //#pragma omp parallel for schedule(static, 1)
  for (Vertex v = 0; v < vertex_count; v++) 
  {

    std::vector<uint8_t> visited(vertex_count);
    //for (auto& i : visited) {
    //  i = 0;
    //}

    //Vertex v = 0;    
    
    if (vertex_degree[v] < 1) {
      continue;
    }

    size_t r = 1; // walk step, initialized to 1

    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0);      

    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map." << std::endl;
    }

    dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
      (vertices, vertex_count, vertex_degree, edges, v, v, walk_history, 
      visited, r);
     
    std::cout << std::endl; // Test

    break; // only require one successful walk from one source
  } // for
  } 
} 

////////////////////////////////////////////////////////////////////////////////

/**
 * This routine is called by multiple graphs is parallel, so the routine itself
 * is siquential 
 */
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList> 
bool is_connected_component(VertexList& vertices, Vertex vertex_count, 
  EdgeList& edges) {

  std::vector<uint8_t> visited(vertex_count);
  for (auto& v : visited) {
    v = 0;
  }

  Vertex source = 0;
  visited[static_cast<size_t>(source)] = 1;  
  bool finished = false;
  bool is_cc = true;

  // TODO: just count the number of visited vertices?

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
      //std::cout << v << " : ";
      if (visited[static_cast<size_t>(v)] == 1) { 
        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          //std::cout << v_nbr << ", ";
          if (visited[static_cast<size_t>(v_nbr)] == 0) {  
            visited[static_cast<size_t>(v_nbr)] = 1;
            finished = false; 
          } // if
        } // for 
      } // if  
      //std::cout << std::endl;
    } // for
  } while(!finished);

  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cout << v << " : " 
    //  << static_cast<uint64_t>(visited[static_cast<size_t>(v)]) << std::endl;       
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      break;
    }  
  }

  return is_cc;
}  

////////////////////////////////////////////////////////////////////////////////

/**
 * Generate combinations
 */
//template <typename N, typename K, typename NCK>
template <typename NCollection, typename ItemType, typename O, 
  typename CombinationBuffer, typename S>
void generate_combinations_recursive(size_t s, size_t n, size_t k, size_t r, 
  CombinationBuffer& buffer, O& combinations, NCollection& n_items, 
  S& combinations_items) {
  if (r == k) {

    //for (auto& t : buffer) {  
    //  std::cout << t; 
    //} 
    //std::cout << std::endl;

    combinations.push_back(buffer);
    
    std::unordered_set<ItemType> buffer_set(0); 
    for (auto i : buffer) {
      buffer_set.insert(n_items[i]);    
    }
    assert(buffer.size() == buffer_set.size());
    combinations_items.push_back(buffer_set);

    return; 
  }

  for (size_t i = s; i < n ; i++) {
    buffer[r] = i;
    generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
      (i + 1, n, k, r + 1, buffer, combinations, n_items, combinations_items);
  } // for  
}

//template <typename N, typename K, typename NCK>
template <typename NCollection, typename ItemType, typename O, typename S>
void generate_combinations_recursive(size_t n, size_t k, 
  NCollection& n_items, O& combinations, S& combinations_items) {

  typedef std::vector<size_t> CombinationBuffer;
  CombinationBuffer buffer(k); // size_t - array index type

  generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
    (0, n, k, 0, buffer, combinations, n_items, combinations_items);
    // size_t s, size_t n, size_t k, size_t r, ...	
}

//------------------------------------------------------------------------------

/**
 * Not used
 */
template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel_3(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) {

  std::cout << "n : " << _n << ", k : " << _k 
    << ", nCk : " << _nCk << std::endl;   

  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  assert(nCk == combinations.size()); 

  {
  //#pragma omp parallel for schedule(static, 1)
  for (size_t i = 0; i < n; ++i) {

    //size_t thread_ID = omp_get_thread_num();
    //std::cout << "Thread# " << thread_ID << " " 
    //  << " Combination# " << i << std::endl;     
       
    //if (i > n - k) {
    //  continue;
    //}
    std::cout << i << ", ";        
    for (size_t j = i + 1; j < n; j++) {      
      std::cout << j << ", ";  
      //for (size_t f = j + 1; f < (j + 1) + k ; f++) {
      //  assert(f < n); 
      //  std::cout << f << ", "; 
      //} // for   
      std::cout<<std::endl;  
    } // for
    std::cout<<std::endl;  
    
  } // for // omp parallel
 
  } 
   
} 

template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel_2(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) {

  std::cout << "n : " << _n << ", k : " << _k 
    << ", nCk : " << _nCk << std::endl;   

  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  size_t s = n - k; // number of items to remove 

  assert(nCk == combinations.size()); 

  {
  //#pragma omp parallel for schedule(static, 1)
  for (size_t i = 0; i < nCk; ++i) {

    size_t thread_ID = omp_get_thread_num();
    //std::cout << "Thread# " << thread_ID << " " 
    //  << " Combination# " << i << std::endl;  
    
    //size_t r = i % n; // r is always smaller than n
    //assert(r < n);

    std::cout << i << ", ";
    for (size_t j = i*s; j < i*s + s; j++) {
      size_t r = j % n;
      
      assert(r < n);    
      assert(combinations[i].test(r) == false);
      combinations[i].set(r);
      std::cout << j << ", " << r << ", ";    
    }
    std::cout << std::endl;   
    
  } // for 
 
  } 
   
} 

template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) {

  std::cout << "n : " << _n << ", k : " << _k 
    << ", nCk : " << _nCk << std::endl;   

  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  assert(nCk == combinations.size()); 

  size_t s = n - k; // number of items to remove 

  {
  #pragma omp parallel for schedule(static, 1)
  for (size_t i = 0; i < nCk; ++i) {

    size_t thread_ID = omp_get_thread_num();
    //std::cout << "Thread# " << thread_ID << " " 
    //  << " Combination# " << i << std::endl;  
    
    size_t r = i % n; // r is always smaller than n
    assert(r < n);
    //if (!(r < n)) {
    //  std::cerr << "Error: value error, r." << std::endl;
    //}
    size_t offset = i / n; // integer division 
    //if (offset >= n) {
    //  std::cerr << "i: " << i << ", r:" << r
    //    << ", offset: " << offset 
    //    << ", n: " << n
    //    << std::endl; // Test
    //}
    assert(offset < n);
    
    //if (offset >= n) {
    //  offset = offset - n;
    //  assert(offset < n);   
    //}  

    // r is the first item
    assert(combinations[i].test(r) == false); 
    combinations[i].set(r);
    
    for (size_t j = 0; j < s - 1; j++) {

      size_t offset_j = offset + j + 1;

      /*if (offset_j >= n) {
        std::cerr << "i: " << i << ", r:" << r << ", j: " << j
          << ", offset: " << offset << ", offset_j: " << offset_j
          << ", n: " << n
          << std::endl; // Test
      }*/ 
      //assert (offset_j < n);

      if (offset_j > (n - 1)) {
        offset_j = offset_j - n;
        assert(offset_j < n);
      }  

      size_t next_item = r + offset_j; 
      // next_item can be > n
 
      if (next_item > (n - 1)) {
        next_item = next_item - n;

        /*if (next_item >= n) {  
          std::cerr << "i: " << i << ", r:" << r << ", j: " << j 
            << ", offset: " << offset << ", offset_j: " << offset_j
            << ", next_item: " << next_item << ", n: " << n  
            << std::endl; // Test 
        }*/ 
        assert (next_item < n);  
      }
      combinations[i].set(next_item);    
    } // for

 
  } // for 
 
  } 
   
} 

////////////////////////////////////////////////////////////////////////////////

/**
 * Graph construction
 */
template <typename EdgeListTuple>
void write_edge_list_file(EdgeListTuple& edge_list, size_t unique_ID, 
  std::string output_filepath) {
  std::string output_filename = output_filepath + "/edgelist_" + 
    std::to_string(unique_ID);
  std::ofstream output_file(output_filename, std::ofstream::out);
  for (size_t e = 0; e < edge_list.size(); e++) {    
    output_file << std::get<0>(edge_list[e]) << " " 
      << std::get<1>(edge_list[e]) << "\n";  
  }  
  output_file.close();
} 

template <typename Vertex, typename VertexData, typename VertexDataList>
void read_vertex_data_file(const std::string vertex_data_input_filename, 
  VertexDataList& vertex_data) {
  std::ifstream vertex_data_input_file(vertex_data_input_filename,
    std::ifstream::in);
  std::string line;
  while (std::getline(vertex_data_input_file, line)) {
    std::istringstream iss(line);
    Vertex v_source(0);
    VertexData v_data(0);
    iss >> v_source >> v_data;
    vertex_data.push_back(v_data);
  }
  vertex_data_input_file.close();
}

template <typename Vertex, typename Edge, typename EdgeListTuple>
Edge read_edge_list_file(const std::string input_filename, EdgeListTuple& edge_list,
  Vertex& max_vertex, Edge skip_lines_count = 0) {
  std::ifstream input_file(input_filename, std::ifstream::in);
  Edge edge_count(0);
  Vertex s(0), t(0);
  std::string line;

  while(std::getline(input_file, line)) {
    std::istringstream iss(line);
    edge_count++;
    //if (edge_count >= skip_lines_count || skip_lines_count == 0) {
      iss >> s >> t;
      edge_list.push_back(std::forward_as_tuple(s, t));

      auto tmp_max_vertex = s >= t ? s : t;

      if (max_vertex < tmp_max_vertex) {
        max_vertex = tmp_max_vertex;
      }

    //}
  }
  input_file.close();
  edge_count-= skip_lines_count;
  assert(edge_count > 0);
  return edge_count;
}

template <typename EdgeListTuple, typename EdgeList>
void generate_edge_list(EdgeListTuple& edge_list, EdgeList& edges) {
  for (size_t e = 0; e < edge_list.size(); e++) {
    edges.push_back(std::get<1>(edge_list[e]));
  }
}

template <typename EdgeListTuple>
EdgeListTuple directed_to_undirected_edge_list
  (EdgeListTuple& directed_edge_list) {  
  EdgeListTuple undirected_edge_list(directed_edge_list);
  assert(directed_edge_list.size() == undirected_edge_list.size());
  for (size_t e = 0; e < directed_edge_list.size(); e++) {
    undirected_edge_list.
      push_back(std::forward_as_tuple(std::get<1>(directed_edge_list[e]), 
      std::get<0>(directed_edge_list[e]))); 
  } // for
  return undirected_edge_list; 
} 

template <typename Vertex, typename Edge, typename EdgeListTuple, 
  typename EdgeListMap>
void generate_unique_edge_list(EdgeListTuple& edge_list, 
  EdgeListMap& edge_list_unique) {

  for (auto edge : edge_list) {
    auto s = std::get<0>(edge);
    auto t = std::get<1>(edge);

    //if (s > t) {
    //  std::swap(s , t);  
    //}

    // unique edge identifier
    std::bitset<VERTEX_COUNT> edge_bitset; // TODO: this is actually vertex bitset / edge hash? 
    edge_bitset.set(static_cast<size_t>(s));
    edge_bitset.set(static_cast<size_t>(t));      
 
    Edge edge_uint = static_cast<Edge>(edge_bitset.to_ullong());

    edge_list_unique.insert({edge_uint, edge});
     
//    edge_list_unique.insert({edge_hash_vertex_bitset<Vertex, Edge>
//      (std::get<0>(edge), std::get<1>(edge)), edge});

    //std::cout << s << " - " << t << " | " << edge_bitset << " | " 
    //  << edge_uint << std::endl;
  }      
}  

// retun a new edge list without the edges in edge_list_filter 
template <typename Edge, typename EdgeListTuple, typename EdgeListFilter>
EdgeListTuple filter_edge_list(EdgeListTuple& edge_list, 
  EdgeListFilter& edge_list_filter) {

  EdgeListTuple new_edge_list(0);

  for (auto edge : edge_list) {
    auto s = std::get<0>(edge);
    auto t = std::get<1>(edge);
   
    // unique edge identifier
    std::bitset<VERTEX_COUNT> edge_bitset;
    edge_bitset.set(static_cast<size_t>(s));
    edge_bitset.set(static_cast<size_t>(t));

    Edge edge_uint = static_cast<Edge>(edge_bitset.to_ullong());

    auto find_edge = edge_list_filter.find(edge_uint); // ignore if found
    if (find_edge == edge_list_filter.end()) { 
      //std::cerr << "Error: did not find the expected item in the set." 
      //  << std::endl;
      //return NULL;     
      // edge not found, add it to the new_edge_list
      new_edge_list.push_back(edge);   
    }
      
  } // for 
  
  return new_edge_list; 
}   

template <typename EdgeListTuple, typename VertexList, typename Vertex>          
Vertex generate_vertex_list(EdgeListTuple& edge_list, VertexList& vertices,
  VertexList& vertex_degree, const Vertex max_vertex) {

  //std::ofstream vertex_file("vertex_file_tmp", std::ofstream::out);

  Vertex vertex_count = 0;
  //Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
  Vertex l = 0; // edge list index
  Vertex degree = 0;
  Vertex current_vertex = vertex_count;
  Vertex source;
  //Vertex target;

  do {
    auto edge = edge_list[l];
    source = std::get<0>(edge);
    //target = std::get<1>(edge);
    if (source == current_vertex) {
      degree++;
      l++;
    } else {
      //std::cout << current_vertex << std::endl;
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
      vertex_degree.push_back(degree);

      //VertexData v_data = get_random_uint(rnd_a, rnd_b, rnd_eng);
      //vertex_data.push_back(v_data);

      // write vertex info to file
      //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " " << "\n";

      //vertex_data_file << current_vertex << " " << v_data << "\n";

      // update vertices array
      degree = 0;
      vertex_count++;
      current_vertex = vertex_count;
    }
  } while(current_vertex <= max_vertex);

  // add the last dummy vertex to the vertex list
  vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));
  //vertex_file << current_vertex << " " << degree << " " << vertices[vertices.size() - 1] << " "
  //  << "\n";
  //vertex_file.close();
  //vertex_data_file.close();

  return vertex_count;
}

/**
 * Generate CSR graph from a given edgelist
 */
// TODO: this function assumes the input edgelist is undirected? 
// Note: it also works for a directed edgelist
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeListTuple, typename EdgeListMap>
void generate_graph(Vertex& vertex_count, VertexList& vertices, 
  VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
  EdgeListMap& edge_list_unique, Vertex& max_vertex) {
  
  //Vertex max_vertex;
  //Edge edge_count = edge_list.size();
  
  // sort edges by source
//  std::cout << "Sorting edges by source vertex ..." << std::endl;
  std::stable_sort(edge_list.begin(), edge_list.end(),
    [](const std::tuple<Vertex, Vertex>& a,
       const std::tuple<Vertex, Vertex>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });
  
  // identify the max vertex
  
  // wrong  
  //Vertex tmp_max_vertex = std::get<0>(edge_list[edge_list.size() - 1]); 
  //if (max_vertex < tmp_max_vertex ) {
  //  std::swap(max_vertex, tmp_max_vertex);
  //}
  
  //max_vertex = std::get<0>(edge_list[edge_list.size() - 1]); 
  // assuming edge_list is sorted by the source  

  // TODO: make sure this is not the issue again 
  // Important: max_vertex > 0 means, max_vertex is defined already
  if (max_vertex < 1) {  
    for (auto& e : edge_list) {
      auto s = std::get<0>(e);
      auto t = std::get<1>(e);
  
      auto tmp_max_vertex = s >= t ? s : t;

      if (max_vertex < tmp_max_vertex) {
          max_vertex = tmp_max_vertex;
      }
    } 
  } // if

  // generate vetex list
  //std::cout << "Generating vertex list ..." << std::endl;
//  std::cout << "Creating CSR row pointers ..." << std::endl;
  vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree, 
    max_vertex);
//  std::cout << "Size of vertex list: " << vertices.size() << std::endl;
//  std::cout << "Size of vertex degree list: " << vertex_degree.size() << std::endl;

  // sort targets in increasing order // TODO: move this inside generate_edge_list  
//  std::cout << "Sorting neighbors ..." << std::endl;
  {
  #pragma omp parallel for
  for (size_t v = 0; v < vertices.size() - 1; v++) {
    size_t start = vertices[v];
    size_t end = vertices[v + 1];
    std::stable_sort(edge_list.begin() + start,
                     edge_list.begin() + end,
       [](const std::tuple<Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
  }
  }

  // generate edge list
//  std::cout << "Generating edge list ..." << std::endl;
  generate_edge_list(edge_list, edges);
//  std::cout << "Size of edge list: " << edges.size() << std::endl;

  if (edge_list_unique.size() < 1) {
    //std::cout << "Generating unique edges ..." << std::endl;
    generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
      (edge_list, edge_list_unique);
    //std::cout << "Number of unique edges: " << edge_list_unique.size()
    //  << std::endl;
  }   
     
//  std::cout << "CSR Graph generation completed." << std::endl;
//  std::cout << "Number of vertices: " << vertex_count << std::endl;
//  std::cout << "Number of edges: " << edge_count << std::endl;
//  std::cout << "Max vertex: " << max_vertex << std::endl;         
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeListTuple>
void generate_graph(Vertex& vertex_count, VertexList& vertices, 
  VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
  Vertex& max_vertex) {

  // TODO: imporve
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> T;  
  T dummy_edge_list_unique;
  dummy_edge_list_unique.insert({0, std::forward_as_tuple(0, 0)});
 
  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple, T>
    (vertex_count, vertices, vertex_degree, edges, edge_list, 
    dummy_edge_list_unique, max_vertex); 
} 

////////////////////////////////////////////////////////////////////////////////

/**
 * UI
 */
void usage(uint64_t walker_count, uint64_t walker_hop_count, 
  size_t thread_count) {  
  std::cerr << "Usage: -k <uint> -t <uint> [edgelist_input_filename vertex_data_input_filename output_filepath]\n"
    //<< "Usage: -d <string> -w <uint> -s <uint> -t <uint> [edgelist_input_filename walks_output_filepath]\n"
    //[vertex_input_file edge_input_file vertex_data_input_file \
    //pattern_input_file vertex_rank_output_file]\n"
    //<< " -d <string>   - distribution type (uniform, random, or biased, Default is " << dist_opt << ")\n"
    << " -h            - print help and exit\n"
    //<< " -w <uint>     - number of walkers (nonzero unsigned integer, Default is " << walker_count << ")\n"
    //<< " -s <uint>     - walker max step count (nonzero unsigned integer, Default is " << walker_hop_count << ")\n"
    << " -t <uint>     - number of CPU threads (nonzero unsigned integer, Default is " << thread_count << ")\n"
    << " -k <uint>     - value of k for generating up to k-edit distance template prototypes (nonzero unsigned integer, Default is  1)\n"  
    << "[file ...] - list of input and output files (required)\n\n";
}

void parse_cmd_line(int argc, char** argv, 
  std::string& vertex_input_filename, std::string& edge_input_filename, 
  std::string& vertex_data_input_filename, std::string& pattern_input_filename, 
  std::string& vertex_rank_output_filename, std::string& walks_output_filename, 
  uint64_t& walker_count, uint64_t& walker_hop_count, size_t& thread_count, 
  size_t default_thread_count, size_t& k_input) {

  std::cout << "CMD line: ";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;
  std::cout << std::endl;

  //std::vector<std::string> dist_opts = {"uniform", "random", "biased", "all"};
  //std::string dist_opt = dist_opts[wlkr_dist];
  //uint64_t walker_count_d = walker_count;
  //uint64_t walker_hop_count_d = walker_hop_count;

  if (argc < 3) { 
  //if (argc < 4) {
    std::cerr << "Too few arguments. " <<std::endl;
    usage(walker_count, walker_hop_count, default_thread_count);
    exit(-1);
  }

  bool prn_help = false;
  char c;

  while ((c = getopt(argc, argv, "d:w:s:t:k:h ")) != -1) {
     switch (c) {
       case 'h':
         prn_help = true;
         break;
       case 'd':
         prn_help = true;
         /*for (size_t i = 0; i < dist_opts.size(); i++) {
           if(boost::iequals(dist_opts[i], optarg)) {
             wlkr_dist = static_cast<DISTRIBUTION>(i);
             prn_help = false;
             break;
           }
         }*/
         break;
      case 'w':
         walker_count = std::stoull(optarg);
         break;
      case 's':
         walker_hop_count = std::stoull(optarg);
         break;
      case 't':
         thread_count = std::stoull(optarg);
         break;
      case 'k':
         k_input = std::stoull(optarg);
         if (k_input < 1) {
           prn_help = true;              
         }
         break;
      default:
         std::cerr << "Unrecognized option: " << c << "." <<std::endl;
         prn_help = true;
         break;
     }
   }

   if (prn_help) {
     usage(walker_count, walker_hop_count, default_thread_count);
     exit(-1);
   }

   // optind is initialized to 1
   /*vertex_input_filename = argv[optind];
   edge_input_filename = argv[optind + 1];
   vertex_data_input_filename = argv[optind + 2];
   //pattern_input_filename = argv[optind + 3];
   //vertex_rank_output_filename = argv[optind + 4];*/

   edge_input_filename = argv[optind];
   vertex_data_input_filename = argv[optind + 1];
   walks_output_filename = argv[optind + 2]; 

   /*std::cout << "\n" //<< vertex_input_filename << " "
     << edge_input_filename << " "
     //<< vertex_data_input_filename << " "
     //<< pattern_input_filename << " "
     //<< vertex_rank_output_filename << " "
     << walks_output_filename << " "
     //<< walker_count << " "
     //<< wlkr_dist 
     << "\n" << std::endl;*/
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Main
 */
int main(int argc, char** argv) {
 
  std::chrono::time_point<std::chrono::steady_clock> start_time;
  std::chrono::time_point<std::chrono::steady_clock> end_time;
  double elapsed_time;
 
  // parse commandline input
   
  std::string vertex_input_filename = "dummy";
  std::string edge_input_filename;
  std::string vertex_data_input_filename;
  std::string edge_data_input_filename = "dummy"; 
  std::string pattern_input_filename = "dummy";
  std::string vertex_rank_output_filename = "dummy";

  size_t thread_count = 0;

  // randomwalk  
  //std::string walks_output_filename; // filepath
  uint64_t walker_count = 10; // dummy
  uint64_t walker_hop_count = 5; // dummy
  //DISTRIBUTION wlkr_dist = RANDOM;
 
  // k-edit distance  
  size_t k_input = 1;
   
  std::string prototype_dir_name;

  // TODO: mandatory or optional edges as input?

  size_t host_thread_count = 0;

  {
  #pragma omp parallel
  { 
    host_thread_count = omp_get_num_threads();
  }
  //std::cout << "Number of available threads: " << host_thread_count 
  //  << std::endl;
  }

  parse_cmd_line(argc, argv, vertex_input_filename, edge_input_filename,
    vertex_data_input_filename, pattern_input_filename,
    vertex_rank_output_filename, /*walks_output_filename,*/ prototype_dir_name,
    walker_count, walker_hop_count, thread_count, host_thread_count, k_input); 
 
  //std::cout << "Application ... " << std::endl;

  // host info
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  std::cout << "Hostname: " << hostname << std::endl;
   
  //size_t host_thread_count = 0;
   
  //{
  //#pragma omp parallel
  //{ 
  //  host_thread_count = omp_get_num_threads();
  //}
  std::cout << "Number of available threads on " << hostname << ": " 
    << host_thread_count << std::endl;
  //}
   
  if (thread_count > 0) {
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);
    {
    #pragma omp parallel
    {
      host_thread_count = omp_get_num_threads();
    }
    std::cout << "Number of threads to be used: " << host_thread_count << std::endl;
    } 
  }

  std::cout << std::endl;
 
  //////////////////////////////////////////////////////////////////////////////

  // type definitions 

  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData;
  typedef uint64_t EdgeData;

  typedef std::vector<std::tuple<Vertex, Vertex>> EdgeListTuple;
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> EdgeListMap;
  typedef std::vector<Vertex> EdgeList;
  typedef std::vector<Edge> VertexList;
  typedef std::vector<VertexData> VertexDataList;

  typedef std::vector<EdgeListTuple> EdgeListTupleVector;
  typedef std::vector<EdgeListTupleVector> EdgeListTupleVectors;

  typedef prunejuice::pattern::template_graph<Vertex, Edge, VertexData, 
    VertexList, EdgeList, EdgeListTuple, EdgeListMap, VertexDataList> TemplateGraph;
  typedef std::vector<TemplateGraph> TemplateGraphVector;
  typedef std::vector<TemplateGraphVector> TemplateGraphVectors; 

  // edge hash, s, t, hop ID - edge hash is used by cycles, hop ID is used by cycles
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> EdgeSet;
  typedef std::vector<EdgeSet> EdgeSetVector;

  typedef prunejuice::pattern::template_constraint<Vertex, Edge, EdgeSet, 
    EdgeSetVector, EdgeListTupleVector> TemplateConstraint;
  typedef std::vector<TemplateConstraint> TemplateConstraintVector;

  typedef prunejuice::pattern::file_utilities<TemplateGraph> FileUtilities;

  typedef std::unordered_map<VertexData, std::unordered_map<Edge, 
    std::tuple<Vertex, Vertex>>> VetexDataVertexPairs;

  //////////////////////////////////////////////////////////////////////////////

  // graph construction

  EdgeListTuple edge_list(0);
  EdgeListMap edge_list_unique(0);
  EdgeList edges(0);
  VertexList vertices(0);
  VertexList vertex_degree(0);
  VertexDataList vertex_data(0);

  Edge edge_count = 0;
  Vertex vertex_count = 0;
  Vertex max_vertex = 0;
  Edge max_degree = 0;
 
  // build in-memory CSR graph (undirected) 

  Edge skip_lines_count = 0; // TODO: read from the commandline 

  std::cout << "Building in-memory CSR graph (undirected) ..." << std::endl;

  std::chrono::time_point<std::chrono::steady_clock> global_start_time =
    std::chrono::steady_clock::now();

  // read an edgelist file
  std::cout << "Reading edgelist file ..." << std::endl;
  edge_count = read_edge_list_file<Vertex, Edge, EdgeListTuple>
    (edge_input_filename, edge_list, max_vertex, skip_lines_count);
  std::cout << "Size of edge list: " << edge_list.size() << std::endl;
  std::cout << "Max vertex: " << max_vertex << std::endl;

  // sort edges by source
  std::cout << "Sorting edges by source vertex ..." << std::endl;
  std::stable_sort(edge_list.begin(), edge_list.end(),
    [](const std::tuple<Vertex, Vertex>& a,
       const std::tuple<Vertex, Vertex>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });

  // generate vetex list
  //std::cout << "Generating vertex list ..." << std::endl;
  std::cout << "Creating CSR row pointers ..." << std::endl;
  vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree, max_vertex);
  std::cout << "Size of vertex list: " << vertices.size() << std::endl;
  std::cout << "Size of vertex degree list: " << vertex_degree.size() << std::endl;

  // sort targets in increasing order
  std::cout << "Sorting neighbors ..." << std::endl;
  {
  #pragma omp parallel for
  for (size_t v = 0; v < vertices.size() - 1; v++) {
    size_t start = vertices[v];
    size_t end = vertices[v + 1];
    std::stable_sort(edge_list.begin() + start,
                     edge_list.begin() + end,
       [](const std::tuple<Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
  }
  }

  // generate edge list
  std::cout << "Generating edge list ..." << std::endl;    
  generate_edge_list(edge_list, edges);
  std::cout << "Size of edge list: " << edges.size() << std::endl;

  std::cout << "Generating unique edges ..." << std::endl;
  generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
   (edge_list, edge_list_unique);
  std::cout << "Number of unique edges: " << edge_list_unique.size() 
    << std::endl;

  std::cout << "CSR Graph generation completed." << std::endl;
  std::cout << "Number of vertices: " << vertex_count << std::endl;
  std::cout << "Number of edges: " << edge_count << std::endl;  
  std::cout << "Max vertex: " << max_vertex << std::endl;

  // read vertex data
  std::cout << "Reading vertex data file ..." << std::endl;
  read_vertex_data_file<Vertex, VertexData, VertexDataList>
    (vertex_data_input_filename, vertex_data);
  std::cout << "Size of vertex data list: " << vertex_data.size() << std::endl;

  // create a graph object for the input template
  //TemplateGraph input_template(edge_list, edge_list_unique, edges, edge_count, 
  //  vertices, vertex_count, max_vertex, vertex_degree);
  TemplateGraph input_template(edge_list, vertex_data, 0);  

  // Test
  /*std::cout << "Number of vertices: " << input_template.vertex_count << std::endl;
  std::cout << "Number of edges: " << input_template.edge_count << std::endl;
  std::cout << "Max vertex: " << input_template.max_vertex << std::endl;*/
  // Test
  
  // Test
  /*std::cout << "vertex_count " << input_template.vertex_count << std::endl;  
  std::cout << "max_vertex " << input_template.max_vertex << std::endl;
  std::cout << "max_degree " << input_template.max_degree << std::endl;
  std::cout << "edge_count " << input_template.edge_count << std::endl;
  std::cout << "vertices.size() " << input_template.vertices.size() << std::endl;
  std::cout << "vertex_degree.size() " << input_template.vertex_degree.size() << std::endl;
  std::cout << "vertex_data.size() " << input_template.vertex_data.size() << std::endl; 
  std::cout << "edges.size() " << input_template.edges.size() << std::endl;
  std::cout << "edgelist.size() " << input_template.edgelist.size() << std::endl;       
  std::cout << "edgelist_unique.size() " << input_template.edgelist_unique.size() 
    << std::endl << std::endl;*/
  // Test

  // Test
  /*for (Vertex v = 0; v < vertex_count; v++) {
    std::cout << v << " " << vertex_degree[v] << " " 
      << vertices[v] << std::endl; 
    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      uint64_t v_nbr = edges[e];
      assert(v_nbr <= max_vertex);
    }  
  }
  std::cout << std::endl;*/ 
  // Test
 
  // Test
  /*for (auto& e : edge_list_unique) {
    auto edge_uint = e.first;
    auto edge = e.second;
    auto s =  std::get<0>(edge);
    auto t = std::get<1>(edge);
    std::cout << edge_uint << " | " << s << " - " << t << std::endl;  
  }*/
  // Test 
  
  // Test
  std::cout << "Identical Vertex Data Vertex Pairs: " << std::endl; 
  for (auto i : input_template.vertex_data_vertex_pairs) {
    std::cout << i.first << " : " << i.second.size() << " : "; //<< std::endl;
    for (auto j : i.second) {       
      auto u = std::get<0>(j.second);
      auto v = std::get<1>(j.second);  
      std::cout << "(" << std::get<0>(j.second) << ", " 
        << std::get<1>(j.second) << "), ";
      //prunejuice::pattern::graph_algorithm::get_shortest_path
      //  <Vertex, Edge, TemplateGraph, EdgeSet>(input_template, u, v);
      //prunejuice::pattern::graph_algorithm::get_shortest_path
      //  <Vertex, Edge, TemplateGraph, EdgeSet>(input_template, v, u); 
    }   
    std::cout << std::endl;
  } 
  // Test 

  std::chrono::time_point<std::chrono::steady_clock> global_end_time =
    std::chrono::steady_clock::now();
  double global_elapsed_time = getElapsedTimeSecond(global_start_time, 
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds." 
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  
  // build in-memory CSR graph (directed) // Note: not used
  
  std::cout << "Building in-memory CSR graph (directed) ..." << std::endl;

  global_start_time =
    std::chrono::steady_clock::now();

  EdgeListTuple edge_list_directed(0);
  EdgeList edges_directed(0);
  VertexList vertices_directed(0);
  VertexList vertex_degree_directed(0);

  Vertex vertex_count_directed = 0;
  Vertex max_vertex_directed = 0; 

  for (auto& e : edge_list_unique) {
    edge_list_directed.push_back(std::forward_as_tuple(std::get<0>(e.second), 
      std::get<1>(e.second)));
  }

  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
    (vertex_count_directed, vertices_directed, vertex_degree_directed, 
    edges_directed, edge_list_directed, max_vertex_directed);

  std::cout << "Number of vertices: " << vertex_count_directed << std::endl;
  std::cout << "Number of edges: " << edges_directed.size() << std::endl;
  std::cout << "Max vertex: " << max_vertex_directed << std::endl;
 
  // Test
  /*for (Vertex v = 0; v < vertex_count_directed; v++) {
    std::cout << v << " " << vertex_degree_directed[v] << " "
      << vertices_directed[v] << std::endl;
    for (Edge e = vertices_directed[v]; e < vertices_directed[v + 1]; e++) {
      uint64_t v_nbr = edges_directed[e];
      assert(v_nbr <= max_vertex_directed);
    }
  }
  std::cout << std::endl;*/
  // Test
 
  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // generate template constraints // TODO: move this to a class / function

  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Generating template constraints ..." << std::endl << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // generate cycle constrains

  std::cout << "Generating cycle constraints ..." << std::endl;

  // vertex_cycles // Note: not used
  //typedef std::vector< std::vector < std::unordered_map <Vertex, size_t> > > VertexNonLocalProperties;
  //typedef std::vector<std::tuple<std::vector<std::tuple<Vertex, Vertex>>, VertexBitSet>> VertexNonLocalProperties;
  //VertexNonLocalProperties vertex_cycles(vertex_count);
   
  typedef std::vector< std::unordered_map<size_t, std::vector<std::tuple<Vertex, Vertex>> > > 
    VertexNonLocalProperties;
  VertexNonLocalProperties vertex_cycles(vertex_count);

  // vertex_cycles_unique
  //typedef std::vector< std::unordered_map<size_t, std::unordered_map<Edge, 
  //  std::tuple<Vertex, Vertex, size_t> > > > VertexNonLocalPropertiesUnique; 
  typedef std::vector< std::unordered_map<size_t, EdgeSet> > VertexNonLocalPropertiesUnique;
  VertexNonLocalPropertiesUnique vertex_cycles_unique(vertex_count);
 
  // find cycles 
  //TODO: we are using vertex_cycles_unique, remove vertex_cycles? 
 
  find_cycles_parallel<Vertex, Edge, VertexList, EdgeList>(vertices, 
    vertex_count, vertex_degree, edges, vertex_cycles, vertex_cycles_unique);

  // Test
  // VertexNonLocalProperties 
  // vertex_cycles
  /*std::cout << "Per-vertex cycles (vector-based path)" << std::endl; 
  for (size_t i = 0; i < vertex_cycles.size(); i++) { // vector 
    //std::cout << vertex_cycles[i].size() << std::endl;
    for (auto& j : vertex_cycles[i]) { // j is a map
      std::cout << i << " | " << j.first << " | ";
      for (auto& k : j.second) { // k is a vector
        std::cout << "(" << std::get<0>(k) << " -- " << std::get<1>(k) << "), "; 
      } // for  
      std::cout << std::endl;  
    } // for
    std::cout << std::endl;
  } // for*/
  // Test   
  
  // Test
  // VertexNonLocalPropertiesUnique
  // vertex_cycles_unique
  std::cout << std::endl; 
  std::cout << "Per-vertex cycles (map-based path)" << std::endl;
  for (size_t i = 0; i < vertex_cycles_unique.size(); i++) { // vector
    //std::cout << vertex_cycles_unique[i].size() << std::endl; 
    for (auto& j : vertex_cycles_unique[i]) { // j is a map
      std::cout << i << " | " << j.first << " | ";
      for (auto& k : j.second) { // k is a map
        std::cout << "(" << k.first << ", " <<  
        std::get<0>(k.second) << " -- " << std::get<1>(k.second) << ", " << 
        std::get<2>(k.second) << "), ";  
      } // for 
      std::cout << std::endl; 
    } // for   
    std::cout << std::endl;
  } // for  
  // Test

  //----------------------------------------------------------------------------
 
  // identify unique cycles in the graph

  // graph_cycles
  
  // vector based path
  typedef std::unordered_map<size_t, std::vector<std::tuple<Vertex, Vertex>> > 
    NonLocalProperties; 
  NonLocalProperties graph_cycles(0);

  // identify unique cycles // TODO: graph_cycles is not used, remove? 
  for (size_t i = 0; i < vertex_cycles.size(); i++) {
    for (auto& j : vertex_cycles[i]) {
      //std::cout << j.first << std::endl; 
      auto find_path = graph_cycles.find(j.first);
      if (find_path == graph_cycles.end()) {
        graph_cycles.insert({j.first, j.second});   
      }    
    } // for       
  } // for 
 
  // Test
  /*std::cout << "Graph cycles (vector-based path)" << std::endl; 
  for (auto& i : graph_cycles) {
    std::cout << i.first << " : ";
    for (auto& j : i.second) {
      std::cout << "(" << std::get<0>(j) << " -- " << std::get<1>(j) << "), ";
    } 
    std::cout << std::endl;  
  }
  std::cout << std::endl;*/
  // Test
 
  // Note: graph_cycles is not used

  // graph_cycles_unique
    
  // map based path
  //typedef std::unordered_map<size_t, std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>>> 
  //  GraphNonLocalPropertiesUnique; 
  typedef std::unordered_map<size_t, EdgeSet> GraphNonLocalPropertiesUnique; 
  GraphNonLocalPropertiesUnique graph_cycles_unique(0);

  typedef std::unordered_map<size_t, EdgeSetVector> GraphNonLocalPropertiesAllPaths;
  GraphNonLocalPropertiesAllPaths graph_cycles_all_paths(0); // grouped by identical path / constraint 

  // find unique cycles in the graph
  //find_unique_cycles<Vertex, Edge>(vertex_cycles_unique, graph_cycles_unique);
  find_unique_cycles<Vertex, Edge>(vertex_cycles_unique, graph_cycles_unique, 
    graph_cycles_all_paths);
 
  // Test
  std::cout << "Graph cycles (map-based path / edge set)" << std::endl; 
  for (auto& i : graph_cycles_unique) {
    std::cout << i.first << " | "; // path hash
    for (auto& j : i.second) { // iterate over the EdgeSet - map<Edge, tuple<Vertex, Vertex, size_t>> 
      std::cout << "(" << j.first << ", " <<
        std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", " <<
        std::get<2>(j.second) << "), "; 
    } 
    std::cout << std::endl;  
  } 
  std::cout << std::endl;
  // Test
  
  // Test
  std::cout << "Graph cycles (map-based path / edge set), all paths - grouped by identical path" << std::endl; 
  for (auto& i : graph_cycles_all_paths) {
    //std::cout << i.first << " | "; // path hash
    for (auto& k : i.second) { // iterate over the EdgeSetVector
      std::cout << i.first << " | "; // path hash
      for (auto& j : k) { // iterate over the EdgeSet - map<Edge, tuple<Vertex, Vertex, size_t>> 
        std::cout << "(" << j.first << ", " <<
        std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", " <<
        std::get<2>(j.second) << "), ";
      } // for
      std::cout << std::endl;  
    } // for 
    std::cout << std::endl; 
  } // for 
  std::cout << std::endl;
  // Test 

  std::cout << "Number of unique cycles in the graph: " << graph_cycles_unique.size() << std::endl;
  
  //----------------------------------------------------------------------------

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  
  // generate duplicate label constraints
  
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Generating duplicate label constraints ..." << std::endl;

  // TODO: move it to a function 
  for (auto i : input_template.vertex_data_vertex_pairs) {
    std::cout << i.first << " : " << i.second.size() << " : "; //<< std::endl;  
    for (auto j : i.second) {
      Vertex s = std::get<0>(j.second);
      Vertex t = std::get<1>(j.second);   
      std::cout << "(" << std::get<0>(j.second) << ", "
        << std::get<1>(j.second) << "), ";

      EdgeSet forward_path = TemplateGraph::edgelisttuple_to_edgeset
        ( prunejuice::pattern::graph_algorithm::get_shortest_path
        <Vertex, Edge, EdgeListTuple, TemplateGraph>(input_template, s, t) );  

      EdgeSet reverse_path = TemplateGraph::edgelisttuple_to_edgeset
        ( prunejuice::pattern::graph_algorithm::get_shortest_path
        <Vertex, Edge, EdgeListTuple, TemplateGraph>(input_template, t, s) );

      // Test
      std::cout << "Graph vertex pair paths (map-based path / edge set)" << std::endl; 
      //for (auto& i : graph_cycles_unique) {
        //std::cout << i.first << " | "; // path hash
        for (auto& j : forward_path) { // iterate over the EdgeSet - map<Edge, tuple<Vertex, Vertex, size_t>> 
          std::cout << "(" << j.first << ", " <<
            std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", " <<
            std::get<2>(j.second) << "), "; 
        } 
        std::cout << std::endl;  
      //}   
      std::cout << std::endl;
      // Test 
 
    } // for 
  
    std::cout << std::endl; 
    
  } // for 
  // TODO: move it to a function

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  
  // generate TDS constraints
  
  global_start_time = std::chrono::steady_clock::now();

  // generate TDS cycle constraints

  std::cout << "Generating TDS cycle constraints ..." << std::endl;

  // assuming a single TDS cycle constraint

  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> TDSEdgeSet;  
  TDSEdgeSet tds_cycles_edge_set;  
  generate_tds_cycles<Vertex, Edge>(graph_cycles_unique, tds_cycles_edge_set); // single edgelist 
  // Note: not used

  // Test
  //std::cout << std::endl; 
  //std::cout << "Graph TDS cycles (map-based path)" << std::endl; 

  /*std::cout << "Full template" << std::endl;
  for (auto& i : tds_cycles_edge_set) {
    std::cout << "(" << i.first << ", " << std::get<0>(i.second) << " -- " << 
      std::get<1>(i.second) << "), " << std::endl;
  } 
  std::cout << std::endl;*/
  // Test

  // multiple TDS cycle constraints
  // recursively fuse cycles
 
  //typedef std::unordered_map<size_t, std::unordered_map<Edge, 
  //  std::tuple<Vertex, Vertex, size_t>>> GraphNonLocalPropertiesUnique;
  GraphNonLocalPropertiesUnique graph_tds_cycles_unique(0); // map of edge sets
 
  generate_tds_cycles<Vertex, Edge>(graph_cycles_unique, graph_tds_cycles_unique);

  // Test
  // Note: it->first and j.second are not useful for TDS edge sets 
  std::cout << "Graph TDS cycles (map-based edge set) " << std::endl; 
  for (auto it = graph_tds_cycles_unique.begin(); it != graph_tds_cycles_unique.end(); ++it) {
    std::cout << it->first << "(x) : ";   
    for (auto& j : it->second) {
      std::cout << "(" << j.first << ", " 
      << std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", " 
      << std::get<2>(j.second) << "(x)), ";  
    }
    std::cout << std::endl;         
  }
  // Test
 
  //----------------------------------------------------------------------------

  // generate TDS cycle constraints - all paths

  // Note: graph_tds_cycles_unique only contain the edges, not the ordered paths
  // required by the NLCC routine  
 
  std::cout << "Graph TDS cycles (vector-based path), all paths - grouped by identical path" << std::endl; // Test

  //GraphNonLocalPropertiesAllPaths graph_tds_cycles_all_paths(0);
  EdgeListTupleVectors graph_tds_cycles_all_paths(0);  

  generate_tds_cycle_constraints<Vertex, Edge, EdgeSet, 
    GraphNonLocalPropertiesUnique, EdgeListTuple, EdgeListTupleVector, 
    EdgeListTupleVectors>(graph_tds_cycles_unique, graph_tds_cycles_all_paths, 
    true);

  // Test
  // Note: it->first and j.second are not useful for TDS edge sets 
  std::cout << "Graph TDS cycles (map-based edge set - updated edgeset_hash) " << std::endl; 
  for (auto it = graph_tds_cycles_unique.begin(); it != graph_tds_cycles_unique.end(); ++it) {
    std::cout << it->first << " : ";   
    for (auto& j : it->second) {
      std::cout << "(" << j.first << ", " 
      << std::get<0>(j.second) << " -- " << std::get<1>(j.second) << ", " 
      << std::get<2>(j.second) << "(x)), ";  
    }
    std::cout << std::endl;         
  }
  // Test

  //---------------------------------------------------------------------------- 
 
  std::cout << "Number of TDS cycles in the graph: " << graph_tds_cycles_unique.size() << std::endl;   
  
  std::cout << std::endl;

  //----------------------------------------------------------------------------
#ifdef ENABLE_BLOCK 
 
  // output TDS cycle constraints // Note: not used 
 
  // map-based edgelist to vector-based edgelist  
   
  EdgeListTuple tds_cycles_new_directed_edge_list(0); 
  for (auto& i : tds_cycles_edge_set) {
    tds_cycles_new_directed_edge_list.push_back(std::forward_as_tuple(std::get<0>(i.second), std::get<1>(i.second))); 
  } // not used
  
  typedef std::vector<EdgeListTuple> VectorEdgeListTuple; 
  VectorEdgeListTuple graph_tds_cycles_unique_directed_edge_list(graph_tds_cycles_unique.size());

  size_t graph_tds_cycles_unique_index = 0;
  for (auto it = graph_tds_cycles_unique.begin(); it != graph_tds_cycles_unique.end(); ++it) {        
    for (auto& j : it->second) {
      graph_tds_cycles_unique_directed_edge_list[graph_tds_cycles_unique_index].push_back(std::forward_as_tuple
        (std::get<0>(j.second), std::get<1>(j.second)));
    }
    graph_tds_cycles_unique_index++;  
  }

  // Test
  /*for (size_t i = 0; i < graph_tds_cycles_unique_directed_edge_list.size(); i++) {
    auto& k = graph_tds_cycles_unique_directed_edge_list[i];
 
    for (size_t j = 0; j < k.size(); j++) {    
      std::cout << std::get<0>(k[j]) << " -- " << 
      std::get<1>(k[j]) << ", ";  
    }

    std::cout << std::endl; 
  }  
  std::cout << "graph_tds_cycles_unique_directed_edge_list size " 
    << graph_tds_cycles_unique_directed_edge_list.size() << std::endl;*/
  std::cout << "Graph TDS cycle constraints: " << std::endl; 
  // Test
 
  { // this part can be sequential, usually not many constraints 
  //#pragma omp parallel for 
  for (size_t i = 0; i < graph_tds_cycles_unique_directed_edge_list.size(); i++) 
  {

  // construct undirected edgelist
  //EdgeListTuple tds_cycles_new_edge_list = 
  //  directed_to_undirected_edge_list(tds_cycles_new_directed_edge_list); // not used
  //EdgeListTuple tds_cycles_new_edge_list = 
  //  directed_to_undirected_edge_list(graph_tds_cycles_unique_directed_edge_list[i]);  

  // directed edgelist    
  //EdgeListTuple& tds_cycles_new_edge_list = tds_cycles_new_directed_edge_list; // not used
  EdgeListTuple& tds_cycles_new_edge_list = graph_tds_cycles_unique_directed_edge_list[i]; 

  // Test
  //for (auto& i : tds_cycles_new_edge_list) {
  //  std::cout << "(" << std::get<0>(i) << " -- " << std::get<1>(i) << "), " << std::endl;
  //}  
  // Test 
  
  EdgeList tds_cycles_new_edges(0);
  VertexList tds_cycles_new_vertices(0);
  VertexList tds_cycles_new_vertex_degree(0);   

  Vertex tds_cycles_new_vertex_count = 0;
  Vertex tds_cycles_new_max_vertex = 0;
 
  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
    (tds_cycles_new_vertex_count, tds_cycles_new_vertices, 
    tds_cycles_new_vertex_degree, tds_cycles_new_edges,
    tds_cycles_new_edge_list, tds_cycles_new_max_vertex);

  // Note: the generated graph contains all the vertices upto max_vertex,   
  // so the is_connected_component routine won't work for the constraints.
  // It works for prototypes.

  // is connected component     
  /*if (is_connected_component<Vertex, Edge, VertexList, EdgeList> 
    (tds_cycles_new_vertices, tds_cycles_new_vertex_count, 
    tds_cycles_new_edges)) {
    std::cout << "Connected component." << std::endl; 
    //write_edge_list_file(edge_list, edge_set_combinations.size(), prototype_dir_name); // output
  } else {  
    //std::cout << "Not a connected component." << std::endl;
    //std::cerr << "Error: Input (sub)template is not a single connected component. "
    //  << "Aborting ... " << std::endl;
    //  return -1;
  }*/

  // generate TDS cycle constraints
  dfs_parallel<Vertex, Edge>(tds_cycles_new_vertices, 
    tds_cycles_new_vertex_count, tds_cycles_new_vertex_degree, 
    tds_cycles_new_edges);
  // TODO: use edge_dfs

  } // for
  } //#pragma omp parallel for
#endif
  //----------------------------------------------------------------------------

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  
  // generate enumeration constraint
  
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Generating template enumeration constraint ..." << std::endl;

  std::cout << "Full template enumeration constraint: " << std::endl;
#ifdef ENABLE_BLOCK
  dfs_single_source<Vertex, Edge>(vertices, vertex_count, vertex_degree, edges, 
    static_cast<Vertex>(0));    
  //dfs_parallel<Vertex, Edge>(vertices, vertex_count, vertex_degree, edges);
  // TODO: use edge_dfs
  
  //dfs_single_source
  //dfs_parallel<Vertex, Edge>(vertices_directed, vertex_count_directed,
  //  vertex_degree_directed, edges_directed); //, static_cast<Vertex>(0));
#endif

  EdgeListTuple graph_enumeration_unique = 
    prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
      (input_template, static_cast<Vertex>(0));
  
  // Test
  for (Vertex v = 0; v < input_template.vertex_count; v++) {
    EdgeListTuple graph_enumeration_all_paths = 
      prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
        (input_template, v);
    
    for (auto i : graph_enumeration_all_paths) {
      std::cout << "(" << std::get<0>(i) << ", " << std::get<1>(i) << "), ";
    }   
    std::cout << std::endl;
  }
  // Test
   
  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
 
  // all constraints
   
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Combining template constraints ..." << std::endl;
 
  //size_t graph_constraint_count = graph_cycles_unique.size() + 
  //  graph_tds_cycles_unique.size(); 

  TemplateConstraintVector template_constraints(0);  
  TemplateConstraintVector cycle_constraints(0);

  //populate_template_constraints<Vertex, Edge, EdgeSet, TemplateConstraint>
  //  (graph_cycles_unique, graph_cycles_all_paths, graph_tds_cycles_unique, 
  //  graph_tds_cycles_all_paths, template_constraints); 

  populate_template_constraints<Vertex, Edge, EdgeSet, TemplateConstraint, 
    GraphNonLocalPropertiesUnique, GraphNonLocalPropertiesAllPaths, 
    EdgeListTuple, EdgeListTupleVector, EdgeListTupleVectors,
    TemplateConstraintVector>
    (graph_cycles_unique, graph_cycles_all_paths, graph_tds_cycles_unique, 
    graph_tds_cycles_all_paths, graph_enumeration_unique, 
    template_constraints); 

  // Note: only cycle constraints are generated form the input template,
  // path and TDS constraints are generated for each prototype (in parallel).
  populate_cycle_constraints<Vertex, Edge, EdgeSet, TemplateConstraint>
    (graph_cycles_unique, graph_cycles_all_paths, cycle_constraints);

  //---------------------------------------------------------------------------- 

  // input template constraints  

  input_template.template_nonlocal_constraints = template_constraints; // Important: 

  // identify template constraints for the input_template 
  for (size_t j = 0; j < template_constraints.size(); j++) {
    if (TemplateGraph::is_subset_edgeset
      (input_template.edgelist_unique, template_constraints[j].edgeset)) {
      assert(j == template_constraints[j].constraint_ID);
      input_template.template_constraints.
        push_back(template_constraints[j].constraint_ID);
    }
  } // for  

  //----------------------------------------------------------------------------

  // Test
  /*for (auto& i : template_constraints) {
    std::cout << i.edgeset_hash << ", " << i.edgeset.size() << ", "
      << i.edgeset_vector.size() << std::endl;
  }  
  std::cout << std::endl;*/
  // Test

  std::cout << "Number of template constraints: " << template_constraints.size() 
    << std::endl;

  std::cout << "Number of cycle constraints in the input template: " 
    << cycle_constraints.size() << std::endl;

  // TODO: sort entries in template_constraints by the length of template_constraint.edgeset  

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;
   
  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  
  // verify if the input template is a connected component
   
  if (prunejuice::pattern::graph_algorithm::is_connected_component
    <TemplateGraph>(input_template)) {
    std::cout << "Input template is a connected component."
      << " Continuing ..." << std::endl;     
  } else {
    std::cerr << "Error: Input template is not a single connected component. "
      << "Aborting ..." << std::endl;
    return -1;  
  } 

  std::cout << std::endl; 

  //////////////////////////////////////////////////////////////////////////////

  // end of constraint generation for the input template

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // generate up to k-edit distance template prototypes
  
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Generating k-edit distance template prototypes ..." 
    << std::endl << std::endl;  

  TemplateGraphVectors k_edit_distance_prototypes(k_input + 1); 

  // the first element in k_edit_distance_prototypes is the input template
  k_edit_distance_prototypes[0].push_back(input_template);
//#ifdef ENABLE_BLOCK
  generate_up_to_k_edit_distance_prototypes<TemplateGraph, TemplateGraphVector, 
    TemplateGraphVectors, TemplateConstraint, TemplateConstraintVector>
    (input_template, k_input, k_edit_distance_prototypes, cycle_constraints); //template_constraints);
//#endif
  // Test
  for (size_t k = 0;  k < k_edit_distance_prototypes.size(); k++) {
    std::cout << "Edit distance k = " << k << " has " 
      << k_edit_distance_prototypes[k].size() << " prototypes" << std::endl;
  }
  // Test

  // TODO: produce a input file to help the python launch script

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////
 
  // create input and output directories / write template prototypes to files
   
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Writing template prototypes to files ..." 
    << std::endl << std::endl;  
//#ifdef ENABLE_BLOCK  
  create_template_prototype_directories<TemplateGraphVectors, FileUtilities>
    (k_edit_distance_prototypes, prototype_dir_name);
//#endif
//#ifdef ENABLE_BLOCK
  write_template_prototype_input_files<TemplateGraph, TemplateGraphVectors,  
    TemplateConstraint, TemplateConstraintVector, 
    FileUtilities>
    //(k_edit_distance_prototypes, template_constraints, prototype_dir_name);
    (k_edit_distance_prototypes, prototype_dir_name);
//#endif
  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  return 0;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   
  //////////////////////////////////////////////////////////////////////////////

  // generate template prototypes
 
  global_start_time = std::chrono::steady_clock::now();

  std::cout << "Generating template prototypes ..." << std::endl;

  // generate edge combinations 

  // setup input
 
  typedef std::vector<Edge> NCollection;
  //NCollection n_items {0, 1, 2, 3, 4, 5, 6}; // Test
  NCollection n_items(0);  
  
  for (auto& e : edge_list_unique) {
    auto edge_uint = e.first;
    n_items.push_back(static_cast<Edge>(edge_uint)); 
  } 

  std::cout << "Sorting unique edges ..." << std::endl;
  std::stable_sort(n_items.begin(), n_items.end(),
    [](const Edge& a, const Edge& b) -> bool {
         return a < b;
       }); 

  // Test
  //for (auto& i : n_items) {
  //  std::cout << i << std::endl;
  //} 
  // Test

  // ~setup input   

  std::cout << "Calculating edge combinaions ..." << std::endl;
 
  //size_t n_collection = 7; // Test
  size_t n_collection = n_items.size();  
  size_t k_selection = k_input;  

  assert(n_collection >= k_selection); 
 
  //size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
  //  (static_cast<size_t>(n_collection.size()), k_selection));
  size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
    (n_collection, k_selection)); 

  std::cout << "n : " << n_collection << ", k : " << k_selection
    << ", nCk : " << nCk << std::endl;

  std::cout << "Total edge combinations: " << nCk << std::endl;

  // setup output

  //constexpr size_t bitset_size = 16; 
  // bitset_size must be greater or equal to n
  //typedef std::bitset<bitset_size> CombinationBitSet;
  //typedef std::vector<CombinationBitSet> CombinationCollection;     

  //CombinationCollection edge_combinations(nCk);   
  //for (auto& i : edge_combinations) {
  //  i.reset();
  //} 

  //assert(bitset_size >= n_collection);
  //assert(nCk == edge_combinations.size());
  
  //generate_combintions_parallel<size_t, size_t, size_t, CombinationBitSet, 
  //  CombinationCollection>(n_collection, k_selection, nCk, combinations);
  
  //generate_combintions_parallel_2<size_t, size_t, size_t, CombinationBitSet,
  //  CombinationCollection>(n_collection, k_selection, nCk, combinations);  
  
  //generate_combintions_parallel_3<size_t, size_t, size_t, CombinationBitSet,
  //  CombinationCollection>(n_collection, k_selection, nCk, combinations);  
 
  // Test 
/*  for (size_t i = 0; i < combinations.size(); ++i) {
    std::cout << "[" << i << "] : " << combinations[i] << " : " ;
    for (size_t j = 0; j < items.size(); j++) {
      if (!combinations[i].test(j)) {
        std::cout << items[j] << ", ";   
      } 
    } // for
    std::cout << std::endl;
  } // for */
  // Test 
 
  typedef std::vector<std::vector<size_t>> O; // indices // TODO   
  typedef std::vector<std::unordered_set<Edge>> S; // edge_uint // TODO  

  O edge_combinations(0); 
  S edge_set_combinations(0);

  std::cout << "Generating edge combinaions ..." << std::endl;

  // k is the edit distance in this context

  // generate all combinations of uinque unsigned integers
  generate_combinations_recursive<NCollection, Edge, O, S>
    (n_collection, k_selection, n_items, edge_combinations, 
    edge_set_combinations); // sequential routine
  
  assert(nCk == edge_combinations.size());
  assert(nCk == edge_set_combinations.size()); 
 
  // Test 
  /*for (auto& e : edge_combinations) {
    for (auto& j : e) {
      std::cout << n_items[j] << ", ";  
    }    
    std::cout << std::endl;
  }
  std::cout << "--- --- --- ---" << std::endl;
  for (auto& e : edge_set_combinations) {
    for (auto& j : e) {
      std::cout << j << ", "; 
    }      
    std::cout << std::endl;
  }*/
  // Test
   
  //----------------------------------------------------------------------------
  
  // is the maximal-edge template is a connected component
    
  if (is_connected_component<Vertex, Edge, VertexList, EdgeList> 
    (vertices, vertex_count, edges)) {
    //std::cout << "Connected component." << std::endl; 
//~~    write_edge_list_file(edge_list, edge_set_combinations.size(), 
//~~      prototype_dir_name); // output
  } else {  
    //std::cout << "Not a connected component." << std::endl;
    std::cerr << "Error: Input (maximal-edge) template is not a single connected component. "
      << "Aborting ... " << std::endl;
    return -1;
  } 

  //----------------------------------------------------------------------------  
    
  // identify the edgelists (prototypes) that form a single connected component
  
  TemplateGraphVector prototypes(0);
   
  {
  #pragma omp parallel for
  for (size_t i = 0; i < edge_set_combinations.size(); i++) {
      //size_t thread_ID = omp_get_thread_num();
      //std::cout << "Thread# " << thread_ID << " "
      //  << " Combination# " << i << std::endl;
      
      EdgeListTuple new_edge_list = filter_edge_list
        <Edge, EdgeListTuple, std::unordered_set<size_t>>
        (edge_list, edge_set_combinations[i]); 

      //std::cout << new_edge_list.size() << std::endl;
      
      EdgeListMap new_edge_list_unique(0); 
      
      // generate new graph from new_edge_list
      EdgeList new_edges(0);
      VertexList new_vertices(0);
      VertexList new_vertex_degree(0);

      Vertex new_vertex_count = 0;

      generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
        (new_vertex_count, new_vertices, new_vertex_degree, new_edges, 
        new_edge_list, new_edge_list_unique, max_vertex); 
        // Important: must use the max_vertex in the maximal-edge template.
        // Here the number of vertices is fixed. 
     
      //assert(new_vertices.size() > 2 && new_edges.size() > 2);

      if (is_connected_component<Vertex, Edge, VertexList, EdgeList>
        (new_vertices, new_vertex_count, new_edges)) {

        TemplateGraph prototype(new_edge_list, new_edge_list_unique, 
          new_edges, new_vertices, new_vertex_degree);

        // identify template constraint matches 
        for (size_t j = 0; j < template_constraints.size(); j++) {
          if(is_subset_edgeset<Vertex, Edge, EdgeListMap, EdgeSet>
            (new_edge_list_unique, template_constraints[j].edgeset)) {
            assert(j == template_constraints[j].constraint_ID);
            prototype.template_constraints.push_back(template_constraints[j].constraint_ID);
          }            
        } 
 
        // add to prototypes
        {
        #pragma omp critical(prototypes)
        {
          //size_t thread_ID = omp_get_thread_num();
          //std::cout << "Thread# " << thread_ID  
          //  << " entered the critical section" << std::endl;
             
          prototypes.push_back(prototype);
          
        } // #pragma omp citical 
        }

        //std::cout << "Connected component." << std::endl;
//~~        write_edge_list_file(new_edge_list, i, prototype_dir_name); // output
      } else {
        //std::cout << "Not a connected component." << std::endl;
      }      
  } // for 
  }

  std::cout << "Number of prototypes (each one is a single connected component): " << prototypes.size() << std::endl; 

  if (prototypes.size() < 1) {
    std::cerr << "Error: No prototype was generated. "
      << "Aborting ... " << std::endl;
    return -1;
  } 

  // Test
   
  /*size_t nlc_count = 0; // Test

  for (size_t i = 0; i < prototypes.size(); i++) {
    //std::cout << "P " << i << " : ";
    //std::cout << prototypes[i].edgelist_unique.size() << std::endl;
    //std::cout << prototypes[i].template_constraints.size() << std::endl;
    if (prototypes[i].template_constraints.size() > 0) {          
      std::cout << "P " << i << " : ";
      //std::cout << prototypes[i].template_constraints.size() << std::endl; 
      //for (auto& e : prototypes[i].edgelist) {
      //  std::cout << "(" << std::get<0>(e) << " - " << std::get<1>(e) << "), ";  
      //}
      for (auto& c : prototypes[i].template_constraints) {
        std::cout << c << ", ";  
      } 
      std::cout << std::endl;
      nlc_count++;
    }
    //std::cout << std::endl; 
  }

  std::cout << "Number of prototypes with NLCs: " << nlc_count << std::endl;*/ 
  // Test
 
  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;
 
  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // prototype-constraint graph (pcg) 
 
  global_start_time = std::chrono::steady_clock::now();
 
  std::cout << "Generating prototype-constraint graph ..." << std::endl;

  // prototype vertex range: 0 - prototype_max_vertex (=prototype.size() - 1)
  // constraint vertex range: prototype_max_vertex + 1 - constraint_max_vertex (=template_constraints.size() - 1)
  // the offset in the template_constraints vector is -(prototype_max_vertex + 1)
  
  Vertex prototype_min_vertex = static_cast<Vertex>(0);
  Vertex prototype_max_vertex = static_cast<Vertex>(prototypes.size() - 1);  
  Vertex constraint_min_vertex = prototype_max_vertex + static_cast<Vertex>(1);
  Vertex constraint_max_vertex = constraint_min_vertex + (static_cast<Vertex>(template_constraints.size()) - 1); 

  Vertex constraint_vertex_offset = constraint_min_vertex;  

  // build edgelist

  EdgeListTuple pcg_edge_list(0);

  for (size_t i = static_cast<size_t>(prototype_min_vertex); 
    i <= static_cast<size_t>(prototype_max_vertex); i++) {
    Vertex s = static_cast<Vertex>(i);
 
    if (prototypes[i].template_constraints.size() < 1) {
      continue;
    }  

    for (size_t j = 0; j < prototypes[i].template_constraints.size(); j++) {
      Vertex t = static_cast<Vertex>(prototypes[i].template_constraints[j]) + constraint_vertex_offset;
      assert(t >= constraint_min_vertex);
      assert(t <= constraint_max_vertex);

      pcg_edge_list.push_back(std::forward_as_tuple(s, t));
      pcg_edge_list.push_back(std::forward_as_tuple(t, s));
    
    } // for
  } // for

  if (pcg_edge_list.size() < 1) { 
    std::cerr << "Error: Prototype-constraint edgelist is empty. "
      << "Aborting ... " << std::endl;
    return -1; 
  } 

  // generate graph

  //EdgeListMap pcg_edge_list_unique(0); 
  EdgeList pcg_edges(0);
  VertexList pcg_vertices(0);
  VertexList pcg_vertex_degree(0);
   
  Vertex pcg_vertex_count = 0;
  Vertex pcg_max_vertex = constraint_max_vertex; // Important: do not update it later 

  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
    (pcg_vertex_count, pcg_vertices, pcg_vertex_degree, pcg_edges,
    pcg_edge_list, pcg_max_vertex);

  std::cout << "Number of vertices: " << pcg_vertex_count << std::endl;
  std::cout << "Number of edges: " << pcg_edges.size() << std::endl;
  std::cout << "Max vertex: " << pcg_max_vertex << std::endl;
  std::cout << "Prototype vertices: " << prototype_min_vertex << " - " 
    << prototype_max_vertex << std::endl;
  std::cout << "Constraint vertices: " << constraint_min_vertex << " - " 
    << constraint_max_vertex << std::endl; 

  assert(pcg_vertex_count > 0);

  // populate vertex metadata  

  const Vertex PROTOTYPE_VERTEX = static_cast<Vertex>('P');   
  const Vertex CONSTRAINT_VERTEX = static_cast<Vertex>('C');

  VertexList pcg_vertex_metadata(pcg_vertex_count);
 
  { 
  #pragma omp parallel for 
  for (Vertex v = 0; v < pcg_vertex_count; v++) {
    if (v >= prototype_min_vertex && v <= prototype_max_vertex) {
      pcg_vertex_metadata[static_cast<size_t>(v)] = PROTOTYPE_VERTEX;
    } else if (v >= constraint_min_vertex && v <= constraint_max_vertex) {
      pcg_vertex_metadata[static_cast<size_t>(v)] = CONSTRAINT_VERTEX; 
    }    
  }
  }   

  // Test
  /*std::cout << std::endl;
  std::cout << "Prototype-constraint edgelist:" << std::endl;
  for (auto& e : pcg_edge_list ) {
    std::cerr << std::get<0>(e) << " " << std::get<1>(e) << std::endl;
  }*/ 
  // Test
  
  // Test
  /*std::cout << std::endl;
  std::cout << "Prototype-constraint graph vertex metadata" << std::endl;
  for (Vertex v = 0; v < pcg_vertex_count; v++) {
    std::cerr << v << " " << pcg_vertex_metadata[v] << std::endl;
  }*/ 
  // Test  

  // Test  
  /*std::cout << std::endl;  
  std::cout << "Prototype-constraint graph:" << std::endl;
  for (Vertex v = 0; v < pcg_vertex_count; v++) {
    std::cout << v << " | " << pcg_vertex_degree[v] << " | "
    << pcg_vertex_metadata[v] << " | ";
      //<< pcg_vertices[v] << std::endl;
    for (Edge e = pcg_vertices[v]; e < pcg_vertices[v + 1]; e++) {
      Vertex v_nbr = pcg_edges[e];
      assert(v_nbr <= pcg_max_vertex);
      std::cout << v_nbr << ", ";
    }
    std::cout << std::endl; 
  }
  std::cout << std::endl;*/ 
  // Test 
  
  //---------------------------------------------------------------------------- 

  std::cout << std::endl;
  std::cout << "Generating prototype-constraint graph (with prototype-prototype edges) ..." << std::endl;

  // create an edge between the prototype vertices that have at least one common 
  // constraint   
   
  typedef std::vector<std::unordered_set<Vertex>> AdjacencySetVector;
  AdjacencySetVector vertices_adjacency_set(pcg_vertex_count); 
  
  {
  #pragma omp parallel for 
  for (Vertex p = prototype_min_vertex; p <= prototype_max_vertex; p++) {
    for (Edge ep = pcg_vertices[p]; ep < pcg_vertices[p + 1]; ep++) {
      Vertex c = pcg_edges[ep];
      // c is a neighbor of v and a constraint type vertex
      assert(pcg_vertex_metadata[static_cast<size_t>(c)] == CONSTRAINT_VERTEX); 
      for (Edge ec = pcg_vertices[c]; ec < pcg_vertices[c + 1]; ec++) {
        Vertex q = pcg_edges[ec];
        // q is a neighbor of c and a prototype type vertex
        assert(pcg_vertex_metadata[static_cast<size_t>(q)] == PROTOTYPE_VERTEX);
        if (p != q) {
          auto find_q = vertices_adjacency_set[p].find(q);
          if (find_q == vertices_adjacency_set[p].end()) {
            vertices_adjacency_set[p].insert(q);    
          }
        }    
      } // for 
      
    } // for
  } // #pragma omp parallel for
  }

  // Test 
  /*for (Vertex v = 0; v < vertices_adjacency_set.size(); v++) { 
    std::cout << v << " | ";
    for (auto& v_nbr : vertices_adjacency_set[v]) {
      std::cout << v_nbr << ", ";   
      assert(v != v_nbr); 
    } 
    std::cout << std::endl;   
  }*/
  // Test
 
  //---------------------------------------------------------------------------

  // create edge between the constraint vertices: 
  // u --> v, where v's edgeset is a subset of u  

  AdjacencySetVector constrint_vertices_adjacency_set(pcg_vertex_count);

  //---------------------------------------------------------------------------
 
  // build edgelist (updated with prototype-prototype edges)
  
  for (Vertex v = 0; v < vertices_adjacency_set.size(); v++) {
    if (pcg_vertex_metadata[static_cast<size_t>(v)] != PROTOTYPE_VERTEX) {
      continue;
    }   
    for (auto& v_nbr_new : vertices_adjacency_set[static_cast<size_t>(v)]) {
       assert(pcg_vertex_metadata[static_cast<size_t>(v_nbr_new)] == 
         PROTOTYPE_VERTEX); 
       assert(v != v_nbr_new); 
       // Important: only add the out-edge 
       pcg_edge_list.push_back(std::forward_as_tuple(v, v_nbr_new));  
    }  
  } 
  
  if (pcg_edge_list.size() < 1) {
    std::cerr << "Error: Prototype-constraint edgelist is empty. "
      << "Aborting ... " << std::endl;
    return -1;
  } else {
    //std::cout << "Updated prototype-constraint edgelist size: " << pcg_edge_list.size() << std::endl; 
  } 

  // generate graph (updated with prototype-prototype edges)

  //EdgeListMap pcg_edge_list_unique(0); 
  //EdgeList pcg_edges(0);
  //VertexList pcg_vertices(0);
  //VertexList pcg_vertex_degree(0);
  
  pcg_edges.clear();
  pcg_vertices.clear();
  pcg_vertex_degree.clear();
 
  //Vertex pcg_vertex_count = 0;
  //Vertex pcg_max_vertex = constraint_max_vertex; // Important: do not update it later 

  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
    (pcg_vertex_count, pcg_vertices, pcg_vertex_degree, pcg_edges,
    pcg_edge_list, pcg_max_vertex);

  std::cout << "Number of vertices: " << pcg_vertex_count << std::endl;
  std::cout << "Number of edges: " << pcg_edges.size() << std::endl;
  std::cout << "Max vertex: " << pcg_max_vertex << std::endl;
  std::cout << "Prototype vertices: " << prototype_min_vertex << " - " 
    << prototype_max_vertex << std::endl;
  std::cout << "Constraint vertices: " << constraint_min_vertex << " - " 
    << constraint_max_vertex << std::endl; 

  assert(pcg_vertex_count > 0);

  // Test
  /*std::cout << std::endl;
  std::cout << "Prototype-constraint edgelist:" << std::endl;
  for (auto& e : pcg_edge_list ) {
    std::cerr << std::get<0>(e) << " " << std::get<1>(e) << std::endl;
  }*/ 
  // Test
  
  // Test
  /*std::cout << std::endl;
  std::cout << "Prototype-constraint graph vertex metadata" << std::endl;
  for (Vertex v = 0; v < pcg_vertex_count; v++) {
    std::cerr << v << " " << pcg_vertex_metadata[v] << std::endl;
  }*/ 
  // Test  

  // Test  
  std::cout << std::endl;  
  std::cout << "Prototype-constraint graph: " << std::endl;
  for (Vertex v = 0; v < pcg_vertex_count; v++) {
    std::cout << v << " | " << pcg_vertex_degree[v] << " | "
    << pcg_vertex_metadata[v] << " | ";
      //<< pcg_vertices[v] << std::endl;
    for (Edge e = pcg_vertices[v]; e < pcg_vertices[v + 1]; e++) {
      Vertex v_nbr = pcg_edges[e];
      assert(v_nbr <= pcg_max_vertex);
      std::cout << v_nbr << ", ";
    }
    std::cout << std::endl; 
  }
  std::cout << std::endl; 
  // Test 

  //--------------------------------------------------------------------------- 

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // TODO:
  // produce all constraints (walks starting from all vertices - cycle, TDS cycle for now)  
  // order constraints based on path lengths
  // reorder vertices in a path based on label frequency

  return 0;
 
} // end of main
