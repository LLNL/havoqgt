//#include <atomic>
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "omp.h"

#include <assert.h>
#include <unistd.h>

//#include <boost/algorithm/string.hpp>
//#include <boost/math/special_functions/binomial.hpp>
//#include <boost/random.hpp>
//#include <boost/random/random_device.hpp>
//#include <boost/random/discrete_distribution.hpp>

////////////////////////////////////////////////////////////////////////////////

constexpr size_t VERTEX_STATE_SIZE = 16;

const uint64_t vertex_log_degree_max = 3.0;
const uint64_t r_max = std::numeric_limits<uint64_t>::max();
static size_t r_max_depth = 0;

/**
 * Returns elapsed (monotonic) time.
 */
double getElapsedTimeSecond(
  std::chrono::time_point<std::chrono::steady_clock> startTime,
  std::chrono::time_point<std::chrono::steady_clock> endTime) {
  return ((endTime - startTime).count()) *
    std::chrono::steady_clock::period::num /
    static_cast<double>(std::chrono::steady_clock::period::den);
}

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
#ifdef SPOKE
      edgeset_all_paths.push_back(prunejuice::pattern::graphalgorithm::edge_dfs
        <Vertex, Edge, EdgeListTuple, EdgeSet>(it->second, v));
#endif
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
  typename OrderedPath, typename OrderedPathEdges, typename UnorderedPathEdges,
  typename VertexNonLocalPropertiesVector>
void find_cycles_recursive(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges,
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,  
  Vertex source_vertex, Vertex v, OrderedPath walk_history, 
  OrderedPathEdges walk_history_edges, 
  UnorderedPathEdges walk_history_edges_unordered, size_t r, 
  VertexNonLocalPropertiesVector& vertex_cycles_vector) {

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

      // Test 
      if (r > UINT64_MAX - 1) { //10) {
        if (v_nbr == source_vertex) {
          std::cerr << source_vertex << " : "; 
          for (auto& i : walk_history) {
            std::cerr << "(" << i.first << ", " << i.second << ") ";
          }
          std::cerr << "(" << v_nbr << ", " << r << ") " << std::endl;  
        }
        //continue; 
        return;
      }
      // Test

      if (vertex_degree[v_nbr] < 2) {
        //continue;
      }
      
      if (v_nbr == source_vertex && r >= 3) { 
        // cycle found
        
        // Test  
        std::cerr << source_vertex << " : ";   
        for (auto& i : walk_history) {
          std::cerr << "(" << i.first << ", " << i.second << ") ";
        }   
        std::cerr << "(" << v_nbr << ", " << r << ") " << std::endl; 
        // Test  

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
#ifdef SPOKE
        update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges));  
#endif 
        //-for (auto& i : std::get<0>(new_walk_history_edges)) {   
        //-  std::cout << "(" << std::get<0>(i) << " -- " << std::get<1>(i) << "), ";  
        //-}
        //std::cout << ": " << std::get<1>(new_walk_history_edges);   
        //-std::cout << ": " << std::get<1>(new_walk_history_edges).to_ullong(); 

        //-std::cout << " : Found Cycle." << std::endl;

        // tuple - map, size_t  
        UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
#ifdef SPOKE
        Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

        auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
        if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
          std::get<0>(new_walk_history_edges_unordered).
            insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered));
        }
#endif
        // add path to vertex_cycles
 
        //----------------------------------------------------------------------       
  
        // vertex_cycles_unique - vector based walk history edges     
#ifdef SPOKE 
        auto find_path = vertex_cycles[source_vertex].find(std::get<1>(new_walk_history_edges).to_ullong());  
        if (find_path ==  vertex_cycles[source_vertex].end()) {    
          vertex_cycles[source_vertex].
            insert({std::get<1>(new_walk_history_edges).to_ullong(), std::get<0>(new_walk_history_edges)});
        } //else {
          //std::cerr << "Error: unexpected item in the map." << std::endl; 
        //}
#endif

        // spoke
        //auto find_path = vertex_cycles[source_vertex].find(source_vertex);
        //if (find_path ==  vertex_cycles[source_vertex].end()) {
        //  vertex_cycles[source_vertex].
        //    insert({source_vertex, std::get<0>(new_walk_history_edges)});
        //} //else {
        //std::cerr << "Error: unexpected item in the map." << std::endl;
        //        //}
       
        vertex_cycles_vector[source_vertex].push_back(std::get<0>(new_walk_history_edges));     
        
        // spoke       

        //----------------------------------------------------------------------

        // vertex_cycles_unique - map based walk history edges    
        //walk_history_edges_unordered
#ifdef SPOKE
        auto find_path_3 = vertex_cycles_unique[source_vertex].find(std::get<1>(new_walk_history_edges_unordered).to_ullong());
        if (find_path_3 ==  vertex_cycles_unique[source_vertex].end()) {  
          vertex_cycles_unique[source_vertex].
            insert({std::get<1>(new_walk_history_edges_unordered).to_ullong(), std::get<0>(new_walk_history_edges_unordered)});    
        }
#endif
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
#ifdef SPOKE            
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges)); 
#endif
          //-------------------------------------------------------------------- 

          // tuple - map, size_t 
          UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
#ifdef SPOKE             
          Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

          auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
          if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
            std::get<0>(new_walk_history_edges_unordered).
              insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});  
            update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered)); 
          }      
#endif
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
            r + 1, vertex_cycles_vector);

        } else { 
          //return; // repeated vertex, ignore   
          //continue;
        }
      } // if
 
    } // for   
} 

template <typename Vertex, typename Edge, typename VertexData, typename VertexList,
  typename EdgeList, typename VertexDataList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique, typename VertexNonLocalPropertiesVector>
void find_cycles_parallel(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges, VertexDataList& vertex_data, 
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,
  VertexNonLocalPropertiesVector& vertex_cycles_vector) {

  { 
  #pragma omp parallel for
  for (Vertex v = 0; v < vertex_count; v++) {
    //Vertex v = 0; // Test   

    if (vertex_degree[v] < 2) {
      //continue;
    }

    //if (v < 100 && vertex_data[v] == 1546094233688003844) {	
    //if (v == 0 && vertex_data[v] == 1) {
      std::cout << v << " " ; // << vertex_data[v] << std::endl; 
    //} else {
    //  continue;
    //}    

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

    if (vertex_data[v] == 1546094233688003844) {  

    find_cycles_recursive<Vertex, Edge, VertexList, EdgeList, 
      VertexNonLocalProperties, VertexNonLocalPropertiesUnique, 
      OrderedPath, OrderedPathEdges, UnorderedPathEdges, 
      VertexNonLocalPropertiesVector>
      (vertices, vertex_count, vertex_degree, edges, vertex_cycles, 
      vertex_cycles_unique, 
      v, v, walk_history, walk_history_edges, walk_history_edges_unordered, r, 
      vertex_cycles_vector);

    } // if 

    //} // if (vertex_data == filter) 
     
  } // for
  } 
}

////////////////////////////////////////////////////////////////////////////////

// spoke

// find cycles in a directed graph

template <typename Vertex, typename Edge, typename Boolean, typename VertexList,
  typename EdgeList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique,
  typename OrderedPath, typename OrderedPathEdges, typename UnorderedPathEdges,
  typename VertexNonLocalPropertiesVector, typename VertexBooleanList>
void find_cycles_recursive(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges,
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,  
  Vertex source_vertex, Vertex v, OrderedPath walk_history, 
  OrderedPathEdges walk_history_edges, 
  UnorderedPathEdges walk_history_edges_unordered, size_t r, 
  VertexNonLocalPropertiesVector& vertex_cycles_vector, 
  VertexBooleanList& vertex_active) {

    // v is the current vertex    
    // r must be greater than 2 for it to be considered a cycle    
    // do not forward a reference to walk_history, make a new copy

    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];

      if (vertex_active[v_nbr] == static_cast<Boolean>(0)) { 
        continue;
      }    

      /*std::cout << " -> " << v_nbr << " <" << v << "> : " ;
      for (auto& i : walk_history) {
        std::cout << "(" << i.first << ", " << i.second << ") ";
      }
      std::cout << std::endl;*/

      // Test 
      if (r > r_max) { // 4) { 
        if (v_nbr == source_vertex) {
          //std::cerr << source_vertex << " : "; 
          std::string output_string = "" + std::to_string(source_vertex) + " : " + std::to_string(walk_history.size()) + " : ";
          for (auto& i : walk_history) {
            //std::cerr << "(" << i.first << ", " << i.second << ") ";
            output_string = output_string + "" + std::to_string(i.first) + " " + std::to_string(i.second) + ", "; 
          }
          //std::cerr << "(" << v_nbr << ", " << r << ") " << std::endl; 
          {
          #pragma omp critical(cerr) 
          { 
            std::cerr << output_string << "" << v_nbr << " " << r << ", " << std::endl; 
          }
          }
        }
        //continue; 
        return;
      }
      // Test

      if (vertex_degree[v_nbr] < 2) {
        //continue;
      }
      
      if (v_nbr == source_vertex && r >= 3) { 
        // cycle found
        
        // Test 
        //std::cerr << source_vertex << " : ";   
        std::string output_string = "" + std::to_string(source_vertex) + " : " + std::to_string(walk_history.size()) + " : " ;
        for (auto& i : walk_history) {
          //std::cerr << "(" << i.first << ", " << i.second << ") ";
          output_string = output_string + "" + std::to_string(i.first) + " " + std::to_string(i.second) + ", "; 
        }   
        //std::cerr << "(" << v_nbr << ", " << r << ") " << std::endl;
        {
        #pragma omp critical(cerr)
        { 
          std::cerr << output_string << "" << v_nbr << " " << r << ", " << std::endl;
        }
        }   
        // Test  

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
#ifdef SPOKE
        update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges));  
#endif 
        //-for (auto& i : std::get<0>(new_walk_history_edges)) {   
        //-  std::cout << "(" << std::get<0>(i) << " -- " << std::get<1>(i) << "), ";  
        //-}
        //std::cout << ": " << std::get<1>(new_walk_history_edges);   
        //-std::cout << ": " << std::get<1>(new_walk_history_edges).to_ullong(); 

        //-std::cout << " : Found Cycle." << std::endl;

        // tuple - map, size_t  
        UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
#ifdef SPOKE
        Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

        auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
        if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
          std::get<0>(new_walk_history_edges_unordered).
            insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered));
        }
#endif
        // add path to vertex_cycles
 
        //----------------------------------------------------------------------       
  
        // vertex_cycles_unique - vector based walk history edges     
#ifdef SPOKE 
        auto find_path = vertex_cycles[source_vertex].find(std::get<1>(new_walk_history_edges).to_ullong());  
        if (find_path ==  vertex_cycles[source_vertex].end()) {    
          vertex_cycles[source_vertex].
            insert({std::get<1>(new_walk_history_edges).to_ullong(), std::get<0>(new_walk_history_edges)});
        } //else {
          //std::cerr << "Error: unexpected item in the map." << std::endl; 
        //}
#endif

        // spoke
        //auto find_path = vertex_cycles[source_vertex].find(source_vertex);
        //if (find_path ==  vertex_cycles[source_vertex].end()) {
        //  vertex_cycles[source_vertex].
        //    insert({source_vertex, std::get<0>(new_walk_history_edges)});
        //} //else {
        //std::cerr << "Error: unexpected item in the map." << std::endl;
        //        //}
       
        vertex_cycles_vector[source_vertex].push_back(std::get<0>(new_walk_history_edges));     
        
        // spoke       

        //----------------------------------------------------------------------

        // vertex_cycles_unique - map based walk history edges    
        //walk_history_edges_unordered
#ifdef SPOKE
        auto find_path_3 = vertex_cycles_unique[source_vertex].find(std::get<1>(new_walk_history_edges_unordered).to_ullong());
        if (find_path_3 ==  vertex_cycles_unique[source_vertex].end()) {  
          vertex_cycles_unique[source_vertex].
            insert({std::get<1>(new_walk_history_edges_unordered).to_ullong(), std::get<0>(new_walk_history_edges_unordered)});    
        }
#endif
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
#ifdef SPOKE            
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges)); 
#endif
          //-------------------------------------------------------------------- 

          // tuple - map, size_t 
          UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
#ifdef SPOKE             
          Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

          auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
          if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
            std::get<0>(new_walk_history_edges_unordered).
              insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});  
            update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered)); 
          }      
#endif
          //--------------------------------------------------------------------  

          //std::cout << " : " << v_nbr << " -> " << std::endl;
//#ifdef SPOKE  
          // forward
          find_cycles_recursive<Vertex, Edge, Boolean, VertexList, EdgeList,
            VertexNonLocalProperties, VertexNonLocalPropertiesUnique,
            OrderedPath, OrderedPathEdges, UnorderedPathEdges>
            (vertices, vertex_count, vertex_degree, edges, vertex_cycles,
            vertex_cycles_unique,  
            source_vertex, v_nbr, new_walk_history, new_walk_history_edges, 
            new_walk_history_edges_unordered,
            r + 1, vertex_cycles_vector, vertex_active);
//#endif
        } else { 
          //return; // repeated vertex, ignore   
          //continue;
        }
      } // if
 
    } // for   
} 

template <typename Vertex, typename Edge, typename Boolean, typename VertexData, typename VertexList,
  typename EdgeList, typename VertexDataList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique, typename VertexNonLocalPropertiesVector,
  typename VertexBooleanList>
void find_cycles_parallel(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges, VertexDataList& vertex_data, 
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,
  VertexNonLocalPropertiesVector& vertex_cycles_vector, 
  VertexBooleanList& vertex_active) {

  { 
  #pragma omp parallel for
  for (Vertex v = 0; v < vertex_count; v++) {
    //Vertex v = 0; // Test   

    if (vertex_active[v] == static_cast<Boolean>(0)) {
      continue;
    }

    if (vertex_degree[v] < 2) {
      //continue;
    }

    //if (v < 100 && vertex_data[v] == 1546094233688003844) {	
    //if (v == 0 && vertex_data[v] == 1) {
    //  std::cout << v << " " ; // << vertex_data[v] << std::endl; 
    //} else {
    //  continue;
    //}    

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

    if (vertex_data[v] == 1546094233688003844) {  

    //std::cout << v << " " ;

//#ifdef SPOKE
    find_cycles_recursive<Vertex, Edge, Boolean, VertexList, EdgeList, 
      VertexNonLocalProperties, VertexNonLocalPropertiesUnique, 
      OrderedPath, OrderedPathEdges, UnorderedPathEdges, 
      VertexNonLocalPropertiesVector, VertexBooleanList>
      (vertices, vertex_count, vertex_degree, edges, vertex_cycles, 
      vertex_cycles_unique, 
      v, v, walk_history, walk_history_edges, walk_history_edges_unordered, r, 
      vertex_cycles_vector, vertex_active);
//#endif
    } // if 

    //} // if (vertex_data == filter) 
     
  } // for
  } 
}

////////////////////////////////////////////////////////////////////////////////

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

  if (r > r_max_depth) {
    r_max_depth = r;
  }

  visited[static_cast<size_t>(v)] = 1;  
  //std::cout << v << ", "; // Test 

  for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
    Vertex v_nbr = edges[e];

    if (visited[static_cast<size_t>(v_nbr)] == 0) { 
      //std::cout << "(" << v << " - " << v_nbr << "), "; // Test
      std::cerr << "(" << v << " - " << v_nbr << " - " << r << "), "; // Test
      dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
        (vertices, vertex_count, vertex_degree, edges, source_vertex, v_nbr, 
        walk_history, visited, r + 1);

    } else {
      //std::cout << "[" << v << " - " << v_nbr << "], "; // Test
      std::cerr << "[" << v << " - " << v_nbr << " - " << r << "], "; // Test
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
  
  std::cout << "DFS source: " << v << std::endl;

  {
  //#pragma omp parallel for schedule(static, 1)
  //for (Vertex v = 0; v < vertex_count; v++) 
  {

    std::vector<uint8_t> visited(vertex_count);
    //for (auto& i : visited) {
    //  i = 0;
    //}

    //Vertex v = 0;    
    
    if (vertex_degree[v] < 1) {
      //continue;      
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

    //break; // only require one successful walk from one source
    
    size_t visited_count = 0;

    for (auto i : visited) {
      if (i == static_cast<uint8_t>(1)) {
        visited_count++; 
      } 
    }

    std::cout << "DFS visited_count: " << visited_count << std::endl;       

    std::cout << "DFS r_max_depth: " << r_max_depth << std::endl;  

  } // for
  } 
} 


////////////////////////////////////////////////////////////////////////////////

/**
 * Level synchrnous BFS
 */
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeStateList> 
bool bfs_level_synchronous_single_source_0(VertexList& vertices, 
  Vertex vertex_count, EdgeList& edges, Vertex source, size_t max_depth,
   
  EdgeStateList& edge_state_list, size_t edge_state_index) {

  //std::cout << "BFS source: " << source << std::endl;
  
  size_t visited_count = 0; 

  std::vector<uint8_t> visited(vertex_count);
  for (auto& v : visited) {
    v = 0;
  }

  std::vector<size_t> level(vertex_count);
  for (auto& v : level) {
    v = std::numeric_limits<uint64_t>::max();
  }

  //Vertex source = 0;
  visited[static_cast<size_t>(source)] = 1;  
  level[static_cast<size_t>(source)] = 0;
  size_t current_level = level[static_cast<size_t>(source)];

  bool finished = false;
  bool is_cc = true;

  // TODO: count the number of visited

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
      //std::cout << v << " : ";
      if (visited[static_cast<size_t>(v)] == 1 && 
        level[static_cast<size_t>(v)] == current_level) {

        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          //std::cout << v_nbr << ", ";
          
          if (v_nbr == v) {
            continue; 
          }    
          
          if (level[static_cast<size_t>(v_nbr)] == 
            std::numeric_limits<uint64_t>::max()) {
            visited[static_cast<size_t>(v_nbr)] = 1;
            level[static_cast<size_t>(v_nbr)] = current_level + 1;

            // update edge state
            edge_state_list[e][edge_state_index] = static_cast<uint8_t>(1);

            if (finished == true) {
              finished = false; 
            }
          } // if

        } // for

        visited[static_cast<size_t>(v)] = 0;
 
      } // if  
      //std::cout << std::endl;
    } // for

    //std::cout << current_level << std::endl;  
 
    if (current_level == max_depth) {
      finished == true; 
    } else {
      current_level++;
    }

  } while(!finished);

  Vertex max_level_vertex = 0;
  size_t max_level = 0;
  size_t not_visited_count = 0; 

  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cerr << v << " " << level[static_cast<size_t>(v)] << std::endl;
    if (level[static_cast<size_t>(v)] == std::numeric_limits<uint64_t>::max()) { 
      not_visited_count++; 
    } else {
      if (level[static_cast<size_t>(v)] > max_level) {
        max_level = level[static_cast<size_t>(v)];
        max_level_vertex = v; 
      } 
      visited_count++;
    } 
  }

  //visited_count = visited.size() - not_visited_count;  

  std::cout << source << " " << max_level_vertex << " " << max_level << " " 
    << visited_count << std::endl;    
  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  /*not_visited_count = 0;
  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cout << v << " : " 
    //  << static_cast<uint64_t>(visited[static_cast<size_t>(v)]) << std::endl;       
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      //break;
      not_visited_count++;
    }  
  }*/

  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  return is_cc;
}

// filter according to vertex_active

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeStateList, typename Boolean,
  typename VertexBooleanList, typename BFSLevelList> 
bool bfs_level_synchronous_single_source_1(VertexList& vertices, 
  Vertex vertex_count, EdgeList& edges, Vertex source, size_t max_depth, 
  BFSLevelList& bfs_level, EdgeStateList& edge_state_list, size_t edge_state_index, 
  VertexBooleanList& vertex_active) {

  //std::cout << "BFS source: " << source << std::endl;
  
  //size_t visited_count = 0; 

  ///std::vector<uint8_t> visited(vertex_count); // TODO: remove this state? 
  ///for (auto& v : visited) {
  ///  v = 0;
  ///}

  ///std::vector<size_t> bfs_level(vertex_count, std::numeric_limits<uint64_t>::max()); // TODO: make it global, similar to edge_state_list
  //for (auto& v : level) {
  //  v = std::numeric_limits<uint64_t>::max();
  //}

  //Vertex source = 0;
  ///visited[static_cast<size_t>(source)] = 1;  
  bfs_level[static_cast<size_t>(source)] = 0;
  size_t current_level = bfs_level[static_cast<size_t>(source)];

  bool finished = false;
  bool is_cc = true; // TODO: change it to found source/sink/target 

  // TODO: count the number of visited

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
 
      if (vertex_active[v] == static_cast<Boolean>(0)) {
        continue;
      }  

      //std::cout << v << " : ";
      if ///(visited[static_cast<size_t>(v)] == 1 && 
        (bfs_level[static_cast<size_t>(v)] == current_level) {

        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          //std::cout << v_nbr << ", ";
          
          if (v_nbr == v || vertex_active[v_nbr] == static_cast<Boolean>(0)) {
            continue; 
          }

          // update edge state; an edge is visited only once 
          edge_state_list[e][edge_state_index] = static_cast<uint8_t>(1);    
          
          if (bfs_level[static_cast<size_t>(v_nbr)] == 
            std::numeric_limits<uint64_t>::max()) {
            ///visited[static_cast<size_t>(v_nbr)] = 1;
            bfs_level[static_cast<size_t>(v_nbr)] = current_level + 1;

            // update edge state
            //edge_state_list[e][edge_state_index] = static_cast<uint8_t>(1);

            if (finished == true) {
              finished = false; 
            }
          } // if

        } // for

        ///visited[static_cast<size_t>(v)] = 0;
 
      } // if  
      //std::cout << std::endl;
    } // for

    //std::cout << current_level << std::endl;  
 
    if (current_level == max_depth) {
      finished == true;         
    } else {
      current_level++;
    }

    // TODO: check if target(s) is found; if yes, stop early    

  } while(!finished);

  // Test

  Vertex max_level_vertex = 0;
  size_t max_level = 0;

  size_t visited_count = 0;
  //size_t not_visited_count = 0; 

  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cerr << v << " " << level[static_cast<size_t>(v)] << std::endl;
    if (bfs_level[static_cast<size_t>(v)] == std::numeric_limits<uint64_t>::max()) { 
      //not_visited_count++; 
    } else {
      if (bfs_level[static_cast<size_t>(v)] > max_level) {
        max_level = bfs_level[static_cast<size_t>(v)];
        max_level_vertex = v; 
      } 
      visited_count++;
    } 
  }

  //visited_count = visited.size() - not_visited_count;  

  std::cout << source << " " << max_level_vertex << " " << max_level << " " 
    << visited_count << std::endl;    
  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  /*not_visited_count = 0;
  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cout << v << " : " 
    //  << static_cast<uint64_t>(visited[static_cast<size_t>(v)]) << std::endl;       
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      //break;
      not_visited_count++;
    }  
  }*/

  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  return is_cc;
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeStateList, typename Boolean,
  typename VertexBooleanList, typename BFSLevelList> 
bool bfs_level_synchronous_single_source(VertexList& vertices, 
  Vertex vertex_count, EdgeList& edges, Vertex source, size_t max_depth, 
  BFSLevelList& bfs_level, EdgeStateList& edge_state_list, size_t state_index, 
  VertexBooleanList& vertex_active) {

  //std::cout << "BFS source: " << source << std::endl;
  
  //size_t visited_count = 0; 

  ///std::vector<uint8_t> visited(vertex_count); // TODO: remove this state? 
  ///for (auto& v : visited) {
  ///  v = 0;
  ///}

  ///std::vector<size_t> bfs_level(vertex_count, std::numeric_limits<uint64_t>::max()); // TODO: make it global, similar to edge_state_list
  //for (auto& v : level) {
  //  v = std::numeric_limits<uint64_t>::max();
  //}

  //Vertex source = 0;
  ///visited[static_cast<size_t>(source)] = 1;  
  bfs_level[static_cast<size_t>(source)][state_index] = 0;
  size_t current_level = bfs_level[static_cast<size_t>(source)][state_index];

  bool finished = false;
  bool is_cc = true; // TODO: change it to found source/sink/target 

  // TODO: count the number of visited

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
 
      if (vertex_active[v] == static_cast<Boolean>(0)) {
        continue;
      }  

      //std::cout << v << " : ";
      if ///(visited[static_cast<size_t>(v)] == 1 && 
        (bfs_level[static_cast<size_t>(v)][state_index] == current_level) {

        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          //std::cout << v_nbr << ", ";
          
          if (v_nbr == v || vertex_active[v_nbr] == static_cast<Boolean>(0)) {
            continue; 
          }

          // update edge state; an edge is visited only once 
          edge_state_list[e][state_index] = static_cast<uint8_t>(1);    
          
          if (bfs_level[static_cast<size_t>(v_nbr)][state_index] == 
            std::numeric_limits<uint64_t>::max()) {
            ///visited[static_cast<size_t>(v_nbr)] = 1;
            bfs_level[static_cast<size_t>(v_nbr)][state_index] = 
              current_level + 1;

            // update edge state
            //edge_state_list[e][state_index] = static_cast<uint8_t>(1);

            if (finished == true) {
              finished = false; 
            }
          } // if

        } // for

        ///visited[static_cast<size_t>(v)] = 0;
 
      } // if  
      //std::cout << std::endl;
    } // for

    //std::cout << current_level << std::endl;  
 
    if (current_level == max_depth) {
      finished == true;         
    } else {
      current_level++;
    }

    // TODO: check if target(s) is found; if yes, stop early    

  } while(!finished);

  // Test

  Vertex max_level_vertex = 0;
  size_t max_level = 0;

  size_t visited_count = 0;
  //size_t not_visited_count = 0; 

  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cerr << v << " " << level[static_cast<size_t>(v)] << std::endl;
    if (bfs_level[static_cast<size_t>(v)][state_index] ==  
      std::numeric_limits<uint64_t>::max()) { 
      //not_visited_count++; 
    } else {
      if (bfs_level[static_cast<size_t>(v)][state_index] > max_level) {
        max_level = bfs_level[static_cast<size_t>(v)][state_index];
        max_level_vertex = v; 
      } 
      visited_count++;
    } 
  }

  //visited_count = visited.size() - not_visited_count;  

  std::cout << source << " " << max_level_vertex << " " << max_level << " " 
    << visited_count << std::endl;    
  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  /*not_visited_count = 0;
  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cout << v << " : " 
    //  << static_cast<uint64_t>(visited[static_cast<size_t>(v)]) << std::endl;       
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      //break;
      not_visited_count++;
    }  
  }*/

  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  return is_cc;
}

////////////////////////////////////////////////////////////////////////////////

/**
 * This routine is called by multiple graphs is parallel, so the routine itself
 * is siquential 
 */
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeStateList> 
bool bfs_single_source(VertexList& vertices, Vertex vertex_count, 
  EdgeList& edges, Vertex source, EdgeStateList& edge_state_list, 
  size_t edge_state_index) {

  //std::cout << "BFS source: " << source << std::endl;
  
  size_t visited_count = 0; 

  std::vector<uint8_t> visited(vertex_count);
  for (auto& v : visited) {
    v = 0;
  }

  std::vector<size_t> level(vertex_count);
  for (auto& v : level) {
    v = std::numeric_limits<uint64_t>::max();
  }

  //Vertex source = 0;
  visited[static_cast<size_t>(source)] = 1;  
  level[static_cast<size_t>(source)] = 0;
  bool finished = false;
  bool is_cc = true;
  size_t current_max_level = level[static_cast<size_t>(source)];

  // TODO: just count the number of visited vertices?

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
      //std::cout << v << " : ";
      if (visited[static_cast<size_t>(v)] == 1) {

        auto v_level = level[static_cast<size_t>(v)];   

        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          //std::cout << v_nbr << ", ";
          if (visited[static_cast<size_t>(v_nbr)] == 0) {  
            visited[static_cast<size_t>(v_nbr)] = 1;

            // update edge state
            edge_state_list[e][edge_state_index] = static_cast<uint8_t>(1);

            if (finished == true) {
              finished = false; 
            }
          } // if

          if (level[static_cast<size_t>(v_nbr)] > (v_level + 1)) {
            level[static_cast<size_t>(v_nbr)] = v_level + 1;
            if (finished == true) {
              finished = false;
            }  

            // Test
            // not useful - this code is not level synchrnous 
            //if (current_max_level < level[static_cast<size_t>(v_nbr)]) { 
            //  current_max_level = level[static_cast<size_t>(v_nbr)];
            //}
            // Test

          } 

        } // for 
      } // if  
      //std::cout << std::endl;
    } // for
    //std::cout << current_max_level << std::endl; 
  } while(!finished);

  Vertex max_level_vertex = 0;
  size_t max_level = 0;
  size_t not_visited_count = 0; 

  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cerr << v << " " << level[static_cast<size_t>(v)] << std::endl;
    if (level[static_cast<size_t>(v)] == std::numeric_limits<uint64_t>::max()) { 
      not_visited_count++; 
    } else {
      if (level[static_cast<size_t>(v)] > max_level) {
        max_level = level[static_cast<size_t>(v)];
        max_level_vertex = v; 
      } 
    } 
  }

  visited_count = visited.size() - not_visited_count;

  std::cout << source << " " << max_level_vertex << " " << max_level << " " 
    << visited_count << std::endl;    
  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  /*not_visited_count = 0;
  for (Vertex v = 0; v < vertex_count; v++) { 
    //std::cout << v << " : " 
    //  << static_cast<uint64_t>(visited[static_cast<size_t>(v)]) << std::endl;       
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      //break;
      not_visited_count++;
    }  
  }*/

  //std::cout << "BFS not visited count: " << not_visited_count << std::endl;

  return is_cc;
}

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList>
void bfs_parallel(VertexList& vertices, Vertex vertex_count,
  EdgeList& edges) {

  //#pragma omp parallel for
  for (Vertex v = 0; v < vertex_count; v++) {
    //bfs_single_source<Vertex, Edge, VertexList, EdgeList>
    //  (vertices, vertex_count, edges, v); 
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
void read_vertex_data_file_0(const std::string vertex_data_input_filename, 
  VertexDataList& vertex_data) {
  std::ifstream vertex_data_input_file(vertex_data_input_filename,
    std::ifstream::in);
  std::string line;
  while (std::getline(vertex_data_input_file, line)) {
    //std::istringstream iss(line);
    //Vertex v_source(0);
    //VertexData v_data("");
    //iss >> v_source; //>> v_data;
    //vertex_data.push_back(v_data);
    vertex_data.push_back(line);	
  }
  vertex_data_input_file.close();
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
    //vertex_data.push_back(v_data);
    vertex_data[v_source] = v_data; 
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

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename VertexEdgeList>
void generate_transpose_graph(Vertex vertex_count, VertexList& vertices, 
  EdgeList& edges, VertexList& vertex_in_degree, 
  VertexEdgeList& vertex_in_edges) {

  vertex_in_degree.resize(vertex_count);
  vertex_in_edges.resize(vertex_count);
 
  for (Vertex v = 0; v < vertex_count; v++) {
    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e]; 
      vertex_in_degree[v_nbr]++;    		
      vertex_in_edges[v_nbr].push_back(v);
    }
  }
} 

template <typename Vertex, typename Edge, typename Boolean, typename VertexList,
  typename EdgeList, typename VertexEdgeList, typename VertexBooleanList>
void eliminate_non_cycle_vertices(Vertex vertex_count, VertexList& vertices, 
    VertexList& vertex_degree, VertexList& vertex_in_degree, EdgeList& edges, 
    VertexEdgeList& vertex_in_edges, VertexBooleanList& vertex_active) {

    vertex_active.resize(vertex_count, static_cast<Boolean>(1)); 

    //vertex_active[3] = static_cast<Boolean>(0); // Test

    // Test	

    // filter by vertex degree

    double max_vertex_log_degree = 0.0; 

    {
    #pragma omp parallel for 
    for (Vertex v = 0; v < vertex_count; v++) {
      double v_log_degree = ceil(log2(vertex_degree[v] + 1));
      if (v_log_degree > vertex_log_degree_max) {
        vertex_active[v] = static_cast<Boolean>(0);
      }   
      if (v_log_degree > max_vertex_log_degree) {
        max_vertex_log_degree = v_log_degree;
      }
    }
    }     

    std::cout << "Max out degree (log2): " << max_vertex_log_degree 
      << std::endl; 

    max_vertex_log_degree = 0.0; 

    {
    #pragma omp parallel for 
    for (Vertex v = 0; v < vertex_count; v++) {
      double v_log_degree = ceil(log2(vertex_in_degree[v] + 1));
      if (v_log_degree > vertex_log_degree_max) {
        vertex_active[v] = static_cast<Boolean>(0);
      }   
      if (v_log_degree > max_vertex_log_degree) {
        max_vertex_log_degree = v_log_degree;
      }
    }
    }     

    std::cout << "Max in degree (log2): " << max_vertex_log_degree 
      << std::endl; 

    // Test

    bool finished = true;

    size_t iteration_count = 0; 

    do {

      finished = true; 
      size_t deleted_vertex_count = 0; 

      for (Vertex v = 0; v < vertex_count; v++) {
 
        if (vertex_active[v] == static_cast<Boolean>(0)) {
         continue;
        } 
        if (vertex_degree[v] < 1 || vertex_in_degree[v] < 1) {
          vertex_active[v] = static_cast<Boolean>(0);
          finished = false;   
          deleted_vertex_count++;
          continue;
        } else {
        
          bool has_active_out_nbr = false;
          bool has_active_in_nbr = false;
   
          // out neighbors 
          for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
            Vertex v_nbr = edges[e];
            if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {
              if (!has_active_out_nbr) {
                has_active_out_nbr = true;
                break;
              }
            } 
          } // for 

          // in neighbors
          for (size_t i = 0; i < vertex_in_edges[v].size(); i++) {   
            Vertex v_nbr = vertex_in_edges[v][i];
            if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {
              if (!has_active_in_nbr) {
                has_active_in_nbr = true;
                break;
              }
            } 
          } // for
  
          if (!has_active_out_nbr || !has_active_in_nbr) {
            vertex_active[v] = static_cast<Boolean>(0);   
            finished = false;   
            deleted_vertex_count++;
          }    

        } // else 

      } // for

      std::cout << "Iteration: " << iteration_count 
        << ", deleted vertex count: " << deleted_vertex_count << std::endl;
      iteration_count++; 

    } while (!finished);
 
} 

template <typename Vertex, typename Edge, typename Boolean, typename VertexList,
  typename EdgeList, typename VertexBooleanList, typename VertexSet>
void eliminate_non_path_vertices_undirected(Vertex vertex_count, 
    VertexList& vertices, VertexList& vertex_degree, EdgeList& edges, 
    VertexBooleanList& vertex_active, VertexSet& uv_vertices) {

    vertex_active.resize(vertex_count, static_cast<Boolean>(1)); // TODO: initialize in main ? 

    //vertex_active[3] = static_cast<Boolean>(0); // Test

    // Test	

    // filter by vertex degree

    double max_vertex_log_degree = 0.0; 

    {
    //#pragma omp parallel for 
    for (Vertex v = 0; v < vertex_count; v++) {
      double v_log_degree = ceil(log2(vertex_degree[v] + 1));
      //if (v_log_degree > vertex_log_degree_max) {
      //  vertex_active[v] = static_cast<Boolean>(0);
      //}   
      if (v_log_degree > max_vertex_log_degree) {
        max_vertex_log_degree = v_log_degree;
      }
    }
    }     

    std::cout << "Max out degree (log2): " << max_vertex_log_degree 
      << std::endl; 

    /*max_vertex_log_degree = 0.0; 

    {
    //#pragma omp parallel for 
    for (Vertex v = 0; v < vertex_count; v++) {
      double v_log_degree = ceil(log2(vertex_in_degree[v] + 1));
      if (v_log_degree > vertex_log_degree_max) {
        vertex_active[v] = static_cast<Boolean>(0);
      }   
      if (v_log_degree > max_vertex_log_degree) {
        max_vertex_log_degree = v_log_degree;
      }
    }
    }     

    std::cout << "Max in degree (log2): " << max_vertex_log_degree 
      << std::endl;*/ 

    // Test

    bool finished = true;

    size_t iteration_count = 0; 

    do {

      finished = true; 
      size_t deleted_vertex_count = 0; 

      for (Vertex v = 0; v < vertex_count; v++) {
 
        if (vertex_active[v] == static_cast<Boolean>(0)) {
         continue;
        } 

        bool is_uv_vertex = false;
        auto find_uv_veretx = uv_vertices.find(v); 
        if (find_uv_veretx != uv_vertices.end()) {
          is_uv_vertex = true;  
        } 

        //if (vertex_degree[v] < 1 || vertex_in_degree[v] < 1) {
        if (vertex_degree[v] < 2 && !is_uv_vertex) { 
          vertex_active[v] = static_cast<Boolean>(0);
          finished = false;   
          deleted_vertex_count++;
          continue;
        } else if (vertex_degree[v] < 1 && is_uv_vertex) {
          vertex_active[v] = static_cast<Boolean>(0);
          finished = false;
          deleted_vertex_count++;
          std::cerr << v << " has no neighbor" << std::endl;
          continue;
        } else {
        
          //bool has_active_out_nbr = false;
          //bool has_active_in_nbr = false;
          
          size_t active_out_nbr_count = 0;  
   
          // out neighbors 
          for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
            Vertex v_nbr = edges[e];
            if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {
              //if (!has_active_out_nbr) {
                //has_active_out_nbr = true;
                //break;
              //}
              
              //bool is_uv_nbr = false;
              //auto find_uv_nbr = uv_vertices.find(v_nbr);
              //if (find_uv_nbr != uv_vertices.end()) {
              //  is_uv_nbr = true;
              //}

              active_out_nbr_count++; // TODO: break early 

            } 
          } // for 

          // in neighbors
          /*for (size_t i = 0; i < vertex_in_edges[v].size(); i++) {   
            Vertex v_nbr = vertex_in_edges[v][i];
            if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {
              if (!has_active_in_nbr) {
                has_active_in_nbr = true;
                break;
              }
            } 
          } // for*/
  
          //if (!has_active_out_nbr || !has_active_in_nbr) {
          //if (!has_active_out_nbr) {
          if (active_out_nbr_count < 2 && !is_uv_vertex)  {
            if (vertex_active[v] == static_cast<Boolean>(1)) {
              vertex_active[v] = static_cast<Boolean>(0);   
              finished = false;   
              deleted_vertex_count++;
            }  
          } else if (active_out_nbr_count < 1 && is_uv_vertex) {
            if (vertex_active[v] == static_cast<Boolean>(1)) {
              vertex_active[v] = static_cast<Boolean>(0);
              finished = false;
              deleted_vertex_count++;
            }
          }    


        } // else 

      } // for

      std::cout << "Iteration: " << iteration_count 
        << ", deleted vertex count: " << deleted_vertex_count << std::endl;
      iteration_count++; 

    } while (!finished); 

} 

// vertex_active - vertices active in the BFS intersection tree
template <typename Vertex, typename Edge, typename Boolean, typename VertexList,
  typename EdgeList, typename VertexBooleanList, typename VertexSet, 
  typename BFSLevelList>
void eliminate_non_path_vertices_from_bfs_intersection_tree_undirected
    (Vertex vertex_count, 
    VertexList& vertices, VertexList& vertex_degree, EdgeList& edges, 
    VertexBooleanList& vertex_active, VertexSet& uv_vertices, 
    BFSLevelList& bfs_level) {

    //vertex_active.resize(vertex_count, static_cast<Boolean>(1)); // TODO: initialize in main ? 

    //vertex_active[3] = static_cast<Boolean>(0); // Test

    // Test	

    // filter by vertex degree

    //double max_vertex_log_degree = 0.0; 

    //{
    //#pragma omp parallel for 
    //for (Vertex v = 0; v < vertex_count; v++) {
      //double v_log_degree = ceil(log2(vertex_degree[v] + 1));
      //if (v_log_degree > vertex_log_degree_max) {
      //  vertex_active[v] = static_cast<Boolean>(0);
      //}   
      //if (v_log_degree > max_vertex_log_degree) {
        //max_vertex_log_degree = v_log_degree;
      //}
    //}
    //}     

    //std::cout << "Max out degree (log2): " << max_vertex_log_degree 
    //  << std::endl; 

    /*max_vertex_log_degree = 0.0; 

    {
    //#pragma omp parallel for 
    for (Vertex v = 0; v < vertex_count; v++) {
      double v_log_degree = ceil(log2(vertex_in_degree[v] + 1));
      if (v_log_degree > vertex_log_degree_max) {
        vertex_active[v] = static_cast<Boolean>(0);
      }   
      if (v_log_degree > max_vertex_log_degree) {
        max_vertex_log_degree = v_log_degree;
      }
    }
    }     

    std::cout << "Max in degree (log2): " << max_vertex_log_degree 
      << std::endl;*/ 

    // Test    

    bool finished = true;

    size_t vertex_state_size = uv_vertices.size();

    //std::bitset<VERTEX_STATE_SIZE> v_successor_bitset; 
    //v_successor_bitset.reset();

    //std::bitset<VERTEX_STATE_SIZE> v_predecessor_bitset;
    //v_predecessor_bitset.reset();     
 
    size_t iteration_count = 0; 

    do {

      finished = true; 
      size_t deleted_vertex_count = 0; 

      for (Vertex v = 0; v < vertex_count; v++) {
 
        if (vertex_active[v] == static_cast<Boolean>(0)) {
         continue;
        } 

        //bool is_uv_vertex = false;
        auto find_uv_veretx = uv_vertices.find(v); 
        if (find_uv_veretx != uv_vertices.end()) {
          //is_uv_vertex = true;  
          continue; 
        } 

        // v (not in uv_vertices) must have at least two neighbors 
        // one with larger bfs level, one with smaller bfs level
        // in all of the bfs searches          

        size_t active_nbr_count = 0; 
 
        //bool has_successor = false;
        bool has_predecessor = false;       

        //bool is_v_valid = false;
        //v_successor_bitset.reset();
        //v_predecessor_bitset.reset();          
 
        //assert(v_successor_bitset.count() == static_cast<size_t>(0));  
        //assert(v_predecessor_bitset.count() == static_cast<size_t>(0));
 
        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {

          Vertex v_nbr = edges[e];

          if (vertex_active[v_nbr] == static_cast<Boolean>(0)) {
            continue;
          }

          //for (size_t i = 0; i < vertex_state_size; i++) {
          for (size_t i = 0; i < static_cast<size_t>(1); i++) { 
            //size_t i = 0;

            //if (bfs_level[v_nbr][i] > bfs_level[v][i]) {                 
              //v_successor_bitset.set(i);
              //is_v_valid = true;
              //break; 
              //has_successor = true; 
            //} 

            if (bfs_level[v_nbr][i] < bfs_level[v][i]) {
              //v_predecessor_bitset.set(i);   
              has_predecessor = true;
              break;  
            } 

            //if (v_successor_bitset.all() && v_predecessor_bitset.all()) {
            //  is_v_valid = true;     
            //  break;  
            //}      

            //if ((v_successor_bitset.count() == vertex_state_size) && 
            //  (v_predecessor_bitset.count() == vertex_state_size)) {                
            //  std::cout << v << " " << v_successor_bitset.count() << std::endl;  
            //  is_v_valid = true;
            //  break;   
            //}   

          } // for
 
          active_nbr_count++;
         
          //if (is_v_valid) {  
          //if (has_successor && has_predecessor) {            
            //break;
          //}   

          if (has_predecessor && (active_nbr_count > 1)) {
            break;
          }    

        } // for

        //if (!is_v_valid) {         
        //if (!has_successor || !has_predecessor) {
        if (!(has_predecessor && (active_nbr_count > 1))) {
          vertex_active[v] = static_cast<Boolean>(0); 
          finished = false;
          deleted_vertex_count++; 
        }   

      } // for

      std::cout << "Iteration: " << iteration_count 
        << ", deleted vertex count: " << deleted_vertex_count << std::endl;
      iteration_count++; 

    } while (!finished); 

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

template <typename VertexIDList>
void parse_cmd_line(int argc, char** argv, 
  std::string& vertex_input_filename, std::string& edge_input_filename, 
  std::string& vertex_data_input_filename, std::string& pattern_input_filename, 
  std::string& vertex_rank_output_filename, std::string& walks_output_filename, 
  uint64_t& walker_count, uint64_t& walker_hop_count, size_t& thread_count, 
  size_t default_thread_count, size_t& k_input, VertexIDList& u_vertices, 
  VertexIDList& v_vertices, size_t& bfs_max_depth) {

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

  while ((c = getopt(argc, argv, "d:w:s:t:k:u:v:h ")) != -1) {
     switch (c) {
       case 'h':
         prn_help = true;
         break;
       case 'd':
         //prn_help = true;
         /*for (size_t i = 0; i < dist_opts.size(); i++) {
           if(boost::iequals(dist_opts[i], optarg)) {
             wlkr_dist = static_cast<DISTRIBUTION>(i);
             prn_help = false;
             break;
           }
         }*/
         bfs_max_depth = std::stoull(optarg);
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
      case 'u' :
         //std::cout << optarg << std::endl;
         u_vertices.push_back(std::stoull(optarg));
         break;
      case 'v' :
         //std::cout << optarg << std::endl;
         v_vertices.push_back(std::stoull(optarg));
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

  typedef uint64_t Vertex;
 
  typedef std::vector<Vertex> VertexIDList;

  std::chrono::time_point<std::chrono::steady_clock> start_time;
  std::chrono::time_point<std::chrono::steady_clock> end_time;
  double elapsed_time;
 
  // parse command line input
   
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
   
  // bfs
  size_t bfs_max_depth = std::numeric_limits<size_t>::max();  

  std::string prototype_dir_name;

  VertexIDList u_vertices(0);
  VertexIDList v_vertices(0); 

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

  parse_cmd_line<VertexIDList>(argc, argv, vertex_input_filename, 
    edge_input_filename, vertex_data_input_filename, pattern_input_filename,
    vertex_rank_output_filename, /*walks_output_filename,*/ prototype_dir_name,
    walker_count, walker_hop_count, thread_count, host_thread_count, k_input,
    u_vertices, v_vertices, bfs_max_depth); 
 
  //std::cout << "Application ... " << std::endl;

  // host info
  size_t _HOST_NAME_MAX = 256;
  char hostname[_HOST_NAME_MAX];
  gethostname(hostname, _HOST_NAME_MAX);
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
 
  //return 0;

  //////////////////////////////////////////////////////////////////////////////

  // type definitions 

  //typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData;
  typedef uint64_t EdgeData;

  typedef uint8_t Boolean;

  typedef std::vector<std::tuple<Vertex, Vertex>> EdgeListTuple;
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> EdgeListMap;
  typedef std::vector<Vertex> EdgeList; // contiguous edges
  typedef std::vector<Edge> EdgeIDList;

  typedef std::vector<Edge> VertexList;
  typedef std::vector<VertexData> VertexDataList;  
  typedef std::vector<std::string> VertexDataStringList;

  typedef std::vector<EdgeListTuple> EdgeListTupleVector;
  typedef std::vector<EdgeListTupleVector> EdgeListTupleVectors;
  typedef std::vector<EdgeList> VertexEdgeList; // non-contiguous edges 
  // TODO: replace with an actual transpose graph
  
  typedef std::vector<Boolean> VertexBooleanList; 
#ifdef SPOKE
  typedef prunejuice::pattern::template_graph<Vertex, Edge, VertexData, 
    VertexList, EdgeList, EdgeListTuple, EdgeListMap, VertexDataList> TemplateGraph;
  typedef std::vector<TemplateGraph> TemplateGraphVector;
  typedef std::vector<TemplateGraphVector> TemplateGraphVectors; 
#endif
  // edge hash, s, t, hop ID - edge hash is used by cycles, hop ID is used by cycles
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> EdgeSet;
  typedef std::vector<EdgeSet> EdgeSetVector;
#ifdef SPOKE
  typedef prunejuice::pattern::template_constraint<Vertex, Edge, EdgeSet, 
    EdgeSetVector, EdgeListTupleVector> TemplateConstraint;
  typedef std::vector<TemplateConstraint> TemplateConstraintVector;

  typedef prunejuice::pattern::file_utilities<TemplateGraph> FileUtilities;
#endif
  typedef std::unordered_map<VertexData, std::unordered_map<Edge, 
    std::tuple<Vertex, Vertex>>> VetexDataVertexPairs;

  typedef std::vector<std::uint64_t> VectorUint64;  
  typedef std::vector<std::vector<uint8_t>> VectorVectorUint8;
  typedef std::vector<VectorUint64> VectorVectorUint64;
  
  typedef VectorVectorUint8 EdgeStateList;  

  typedef std::unordered_set<Vertex> VertexSet;

  //////////////////////////////////////////////////////////////////////////////

  // graph construction

  EdgeListTuple edge_list(0);
  EdgeListMap edge_list_unique(0);
  EdgeList edges(0); // out edges
  EdgeIDList edge_ID_list(0); // consistent with edges

  VertexList vertices(0);
  VertexList vertex_degree(0); // out degree for directed graphs

  VertexDataList vertex_data(0);
  VertexDataStringList vertex_data_string(0);

  VertexEdgeList vertex_in_edges(0); // in edges for directed graphs
  VertexList vertex_in_degree(0); // in degree for directed graphs

  VertexBooleanList vertex_active(0);

  EdgeStateList edge_state_list(0);

  Edge edge_count = 0;
  Vertex vertex_count = 0;
  Vertex max_vertex = 0;
  Edge max_degree = 0;
  Edge max_in_degree = 0;
 
  // build in-memory CSR graph (undirected) 

  Edge skip_lines_count = 0; // TODO: read from the commandline 

  std::cout << "Building in-memory CSR graph (undirected) ..." << std::endl;
  //std::cout << "Building in-memory CSR graph (directed) ..." << std::endl; // spoke

  std::chrono::time_point<std::chrono::steady_clock> global_start_time =
    std::chrono::steady_clock::now();

  // read an edge list file
  std::cout << "Reading edge list file ..." << std::endl;
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

#ifdef SPOKE
  std::cout << "Generating unique edges ..." << std::endl;
  generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
   (edge_list, edge_list_unique);
  std::cout << "Number of unique edges: " << edge_list_unique.size() 
    << std::endl;
#endif

  std::cout << "CSR Graph generation completed." << std::endl;
  std::cout << "Number of vertices: " << vertex_count << std::endl;
  std::cout << "Number of edges: " << edge_count << std::endl;  
  std::cout << "Max vertex: " << max_vertex << std::endl;

  // read vertex data
#ifdef SPOKE 
  std::cout << "Reading vertex data file ..." << std::endl;
  read_vertex_data_file_0<Vertex, std::string, VertexDataStringList>
    ("/p/lustre2/havoqgtu/Spoke/spoke-20201130/vertex_data_parts/vertex_data_2_rec_cmp", 
    vertex_data_string);
  std::cout << "Size of vertex data list: " << vertex_data_string.size() << std::endl;
#endif

  std::cout << "Reading vertex data file ..." << std::endl;
  vertex_data.resize(vertex_count); // TODO: improve
  read_vertex_data_file<Vertex, VertexData, VertexDataList>
    (vertex_data_input_filename, vertex_data);
  std::cout << "Size of vertex data list: " << vertex_data.size() << std::endl;

  // Test 
  //for (auto& v : vertex_data) {
  //  std::cout << v << std::endl;
  //}
  // Test
  
  // edge state list
  edge_state_list.resize(edge_count); // TODO: define size here 

  size_t edge_state_size = u_vertices.size() + v_vertices.size();

  for (auto& e : edge_state_list) { 
    e.resize(edge_state_size);
  }  

  std::cout << "Size of edge state list: " << edge_state_list.size() << ", " << 
    edge_state_list[0].size() << std::endl;

#ifdef SPOKE
  // Test 
  std::ofstream test_output_file("/p/lustre2/havoqgtu/Spoke/spoke-20201130/edge_list_2_3_undir.totem", std::ofstream::out);
  test_output_file << "#Nodes:" << vertex_count << "\n";
  test_output_file << "#Edges:" << edge_count << "\n";
  test_output_file << "#Undirected" << "\n"; 
  for (Vertex v = 0; v < vertex_count; v++) {
    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      uint64_t v_nbr = edges[e];
      test_output_file << v << " " << v_nbr << "\n";   
    }
  }
  test_output_file.close();  
  // Test
#endif
   
  //return 0;

  //////////////////////////////////////////////////////////////////////////////

  // spoke
#ifdef SPOKE
  // generate transpose graph 
  
  std::cout << "Generating transpose graph ..." << std::endl;
  
  generate_transpose_graph<Vertex, Edge>(vertex_count, vertices, edges, 
     vertex_in_degree, vertex_in_edges);

  // Test
  Edge tmp_max_in_degree = 0;  
  for (Vertex v = 0; v < vertex_count; v++) {
    assert(vertex_in_edges[v].size() == vertex_in_degree[v]);    
    if (vertex_in_degree[v] > tmp_max_in_degree) { 
      tmp_max_in_degree = vertex_in_degree[v];
    }
  }  
  std::cout << "Max in degree: " << tmp_max_in_degree << std::endl;
  // Test
#endif

#ifdef ENABLE_BLOCK
  // Test
  { 
  std::ofstream test_output_file("/p/lustre2/havoqgtu/Spoke/spoke-20201130/vertex_data_parts/vertex_data_6_rec_cmp", 
    std::ofstream::out);
  for (Vertex v = 0; v < vertex_count; v++) {
    test_output_file << v << " " << vertex_degree[v] << " " << vertex_in_degree[v] 
      << " " << vertex_data_string[v] << "\n";
  }
  test_output_file.close(); 
  }
  // Test
#endif
  //////////////////////////////////////////////////////////////////////////////
#ifdef SPOKE
  // create a graph object for the input template
  //TemplateGraph input_template(edge_list, edge_list_unique, edges, edge_count, 
  //  vertices, vertex_count, max_vertex, vertex_degree);
  TemplateGraph input_template(edge_list, vertex_data, 0);  
#endif
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
    std::cerr << v << " " << vertex_degree[v] << " "  
      << vertices[v] << " " << vertex_data[v] << std::endl; 
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
#ifdef SPOKE 
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
#endif
  std::chrono::time_point<std::chrono::steady_clock> global_end_time =
    std::chrono::steady_clock::now();
  double global_elapsed_time = getElapsedTimeSecond(global_start_time, 
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds." 
    << std::endl;

  std::cout << std::endl;
#ifdef SPOKE
  return 0;	
#endif
  //////////////////////////////////////////////////////////////////////////////
  
  // build in-memory CSR graph (directed) // Note: not used
#ifdef SPOKE
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
#endif

#ifdef SPOKE
  return 0;
#endif

  //////////////////////////////////////////////////////////////////////////////

  // eliminate non-path vertices 

  std::cout << "Eliminating non-path vertices ... " << std::endl;

  global_start_time =
    std::chrono::steady_clock::now();

  VertexSet uv_vertices(u_vertices.size() + v_vertices.size());

  std::copy(u_vertices.begin(), u_vertices.end(), std::inserter(uv_vertices, 
    uv_vertices.end()));
  std::copy(v_vertices.begin(), v_vertices.end(), std::inserter(uv_vertices, 
    uv_vertices.end()));
  std::cout << uv_vertices.size() << std::endl;
 
  for (auto& v: uv_vertices) {
    std::cout << v << " ";
  }
  std::cout << std::endl;  
  
  // vertex_active must not be initialized
  eliminate_non_path_vertices_undirected<Vertex, Edge, Boolean>(vertex_count,
    vertices, vertex_degree, edges, vertex_active, uv_vertices);  

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  //std::cout << std::endl;

  // Test
  {  
  std::cout << "Size of vertex active list: " << vertex_active.size() 
    << std::endl;

  size_t tmp_active_vertex_count = 0;
  size_t tmp_active_edge_count = 0;

  for (Vertex v = 0; v < vertex_count; v++) {
    if (vertex_active[v] == static_cast<Boolean>(1)) {
      tmp_active_vertex_count++;
    } else {
      continue;
    }

    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];
      if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {
        tmp_active_edge_count++;
        //std::cerr << v << " " << v_nbr << std::endl; // TODO: uncomment to output the edge list
        //std::cerr << v_nbr << " " << v << std::endl;	 
      }
    }

  } // for

  std::cout << "Active vertex count: " << tmp_active_vertex_count << std::endl;
  std::cout << "Active edge count: " << tmp_active_edge_count << std::endl;
  std::cout << std::endl; 
  }
  // Test    
 
  //return 0;  
  
  //////////////////////////////////////////////////////////////////////////////
 
  // s-t path search

  std::cout << "Searching paths ... " << std::endl;

  std::cout << "Source vertices: ";
  for (auto& v: u_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";   
  }
  std::cout << std::endl;

  std::cout << "Target vertices: ";
  for (auto& v: v_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";
  }
  std::cout << std::endl;
 
  // bfs 

  //typedef VectorUint64 BFSLevelList;
  //BFSLevelList bfs_level(vertex_count, std::numeric_limits<uint64_t>::max());

  typedef VectorVectorUint64 BFSLevelList;
  VectorUint64 vertex_bfs_levels(u_vertices.size() + v_vertices.size(), 
    std::numeric_limits<uint64_t>::max());
  BFSLevelList bfs_level(vertex_count, vertex_bfs_levels);

  //std::cout << bfs_level[0].size() << " " << 
  //  static_cast<size_t>(bfs_level[0][0]) << std::endl; // Test 

  bfs_max_depth = bfs_max_depth - 1; // the code computes current_level + 1 neighbors 
  bool is_cc = false;

  for (size_t i = 0; i < u_vertices.size(); i++) {

    global_start_time = std::chrono::steady_clock::now();   

    is_cc = bfs_level_synchronous_single_source<Vertex, Edge, VertexList, 
      EdgeList, EdgeStateList, Boolean, VertexBooleanList, BFSLevelList>
      (vertices, vertex_count, edges, static_cast<Vertex>(u_vertices[i]), 
      bfs_max_depth, bfs_level, edge_state_list, i, 
      vertex_active); 
 
    global_end_time =  std::chrono::steady_clock::now();
    global_elapsed_time = getElapsedTimeSecond(global_start_time,
      global_end_time);
    std::cout << "Time: " << global_elapsed_time << " seconds."
      << std::endl;
  }
 
  size_t v_edge_state_index_begin = u_vertices.size();
 
  for (size_t i = 0; i < v_vertices.size(); i++) {

    global_start_time = std::chrono::steady_clock::now();   

    is_cc = bfs_level_synchronous_single_source<Vertex, Edge, VertexList, 
      EdgeList, EdgeStateList, Boolean, VertexBooleanList, BFSLevelList>
      (vertices, vertex_count, edges, static_cast<Vertex>(v_vertices[i]), 
      bfs_max_depth, bfs_level, edge_state_list, i + v_edge_state_index_begin, 
      vertex_active);  

    global_end_time =  std::chrono::steady_clock::now();
    global_elapsed_time = getElapsedTimeSecond(global_start_time,
      global_end_time);
    std::cout << "Time: " << global_elapsed_time << " seconds."
      << std::endl;
  }

  // Test
  std::cout << std::endl; 
  for (size_t i = 0; i < edge_state_size; i++) {
    std::cout << u_vertices[0] << " " << 
      bfs_level[static_cast<size_t>(u_vertices[0])].size() << " " <<
      static_cast<size_t>(bfs_level[static_cast<size_t>(u_vertices[0])][i]) <<
      std::endl;
  }
  
  for (size_t i = 0; i < edge_state_size; i++) {
    std::cout << v_vertices[0] << " " << 
      bfs_level[static_cast<size_t>(v_vertices[0])].size() << " " <<
      static_cast<size_t>(bfs_level[static_cast<size_t>(v_vertices[0])][i]) << 
      std::endl;   
  }
  std::cout << std::endl;
  // Test

  /*global_start_time = std::chrono::steady_clock::now();

  std::cout << "BFS" << std::endl;

  is_cc = bfs_single_source<Vertex, Edge, VertexList, EdgeList, EdgeStateList>
    //(vertices, vertex_count, edges, static_cast<Vertex>(3)); // reaction - compound
    (vertices, vertex_count, edges, static_cast<Vertex>(u_vertices[0]), 
    edge_state_list, 0); // all 
    // 2294390 1771176 	

  //std::cout << "Connected component: " << is_cc << std::endl;
  
  //bfs_parallel<Vertex, Edge, VertexList, EdgeList>
    //(vertices, vertex_count, edges);

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;*/

  //return 0; 
  
  //////////////////////////////////////////////////////////////////////////////

  // reset states
  
  std::cout << "Reset states ... " << std::endl;

  std::fill(vertex_active.begin(), vertex_active.end(),  
    static_cast<Boolean>(0));

  std::cout << "Source vertices: ";
  for (auto& v: u_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";
  }
  std::cout << std::endl;

  std::cout << "Target vertices: ";
  for (auto& v: v_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";
  }
  std::cout << std::endl;
 
  std::cout << std::endl;

  //////////////////////////////////////////////////////////////////////////////

  // intersection of edges in the BFS trees (for undirected graph) 

  std::cout << "BFS tree intersection ... " << std::endl;
 
  global_start_time = std::chrono::steady_clock::now();

  size_t edges_processed = 0; 

  {
  //#pragma omp parallel for schedule (guided) 
  for (size_t e = 0; e < edge_count; e++) {
    auto s = std::get<0>(edge_list[e]);
    auto t = std::get<1>(edge_list[e]); 

    if (s > t) {
      continue;
    } //else if (s == t) {      
      //std::cerr << s << " " << t << " " << vertex_data[s] 
      //  << std::endl; // self edge check
      //continue;
    //}  

    edges_processed++; 

    // neighbors of t
    for (Edge f = vertices[t]; f < vertices[t + 1]; f++) {
      auto t_nbr = edges[f];   
      if (t_nbr == s) {

        bool is_valid_edge = true; 

        for (size_t i = 0; i < edge_state_size; i++ ) {
          auto edge_state_max = std::max(edge_state_list[e][i], 
            edge_state_list[f][i]); 

          // all u,v states must be set
          if (edge_state_max < static_cast<uint8_t>(1)) { 
            is_valid_edge = false;
            break;
          }   

          // TODO: remove ?
          if (edge_state_list[e][i] < edge_state_max) {
            edge_state_list[e][i] = edge_state_max;
          } 
          if (edge_state_list[f][i] < edge_state_max) {
            edge_state_list[f][i] = edge_state_max;
          }

        } // for
 
        //for (size_t i = 0; i < edge_state_size; i++ ) {
        //  if (edge_state_list[e][i] != static_cast<uint8_t>(1) || 
        //    edge_state_list[f][i] != static_cast<uint8_t>(1)) {
        //    is_valid_edge = false;
        //    break;
        //  }   
        //}  

        if (is_valid_edge) { // TODO: output to file
          //std::cerr << s << " " << t << std::endl; 
          //std::cerr << t << " " << s << std::endl;
          
          if (vertex_active[static_cast<size_t>(s)] == static_cast<Boolean>(0)) {
            {
            //#pragma omp atomic write 
            vertex_active[static_cast<size_t>(s)] = static_cast<Boolean>(1);
            }  
          } 
 
          if (vertex_active[static_cast<size_t>(t)] == static_cast<Boolean>(0)) {
            {
            //#pragma omp atomic write
            vertex_active[static_cast<size_t>(t)] = static_cast<Boolean>(1);       
            }
          }

        }   
 
        break;  
      } 
    } // for    
   
  } // for	
  }

  // Test
  std::cout << "Edges processed: " << edges_processed 
    << " | " << (edge_count / 2) << std::endl;

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl; 
 
  // Test
  {  
  std::cout << "Source vertices: ";
  for (auto& v: u_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";
  }
  std::cout << std::endl;

  std::cout << "Target vertices: ";
  for (auto& v: v_vertices) {
    std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", ";
  }
  std::cout << std::endl;

  size_t tmp_active_vertex_count = 0;
  for (Vertex v = 0; v < vertex_count; v++) { 
    if (static_cast<size_t>(vertex_active[v]) == static_cast<Boolean>(1)) {
      //std::cout << v << " " << static_cast<size_t>(vertex_active[v]) << ", "; 
      tmp_active_vertex_count++;
    }  
  }
  std::cout << "Active vertex count: " << tmp_active_vertex_count << std::endl;

  std::cout << std::endl;
  }
  // Test

  //return 0;

  //////////////////////////////////////////////////////////////////////////////

  // filter

  std::cout << "Filter ... " << std::endl; 

  {
  #pragma omp parallel for schedule (guided) 
  for (Vertex v = 0; v < vertex_count; v++) {
    auto find_uv_veretx = uv_vertices.find(v);
    if (find_uv_veretx != uv_vertices.end()) {          
      continue;
    } else {

      // TODO: user input, put them in a set   
      if ((vertex_data[v] == static_cast<VertexData>(1546094233688003844)) || 
        (vertex_data[v] == static_cast<VertexData>(13996733750220222916U))) {
        // do nothing 
      } else {
         vertex_active[static_cast<size_t>(v)] = static_cast<Boolean>(0);
         continue;
      }  

      if (vertex_active[static_cast<size_t>(v)] == static_cast<Boolean>(1)) {
        for (size_t i = 0; i < edge_state_size; i++) {
          if (bfs_level[v][i] >= (bfs_max_depth + 1)) {
            vertex_active[static_cast<size_t>(v)] = static_cast<Boolean>(0);
            break; 
          }   
        } // for 
      }    
    }
  } // for 
  }

  //////////////////////////////////////////////////////////////////////////////

  // eliminate invalid edges 

  std::cout << "Eliminate invalid edges from the BFS intersection tree ... " 
    << std::endl;

  global_start_time = std::chrono::steady_clock::now();
 
  eliminate_non_path_vertices_from_bfs_intersection_tree_undirected
    <Vertex, Edge, Boolean>(vertex_count,
    vertices, vertex_degree, edges, vertex_active, uv_vertices, bfs_level);
 
  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  // output  
 
  size_t active_vertex_count = 0;
  size_t active_edge_count = 0;  
  for (Vertex v = 0; v < vertex_count; v++) {
    if (vertex_active[v] == static_cast<Boolean>(0)) {
      continue;
    } else {
      active_vertex_count++; 
    } 
    
    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];
      if (vertex_active[v_nbr] == static_cast<Boolean>(0)) {
        continue;
      } else {
        //if (v < v_nbr) {
          active_edge_count++;           
          std::cerr << v << " " << v_nbr << std::endl;        
          //std::cerr << v << ", " << bfs_level[v][0] << ", " << bfs_level[v][1] 
          //  << " | " << v_nbr << ", " << bfs_level[v_nbr][0] << ", " 
          //  << bfs_level[v_nbr][1] << std::endl;
        //}
      }  
    } // for	
  } // for

  std::cout << "Active vertex count: " << active_vertex_count << std::endl;
  std::cout << "Active edge count: " << active_edge_count << std::endl;
  std::cout << std::endl; 
 
  return 0; 

  //////////////////////////////////////////////////////////////////////////////

  // spoke
#ifdef SPOKE
  // eliminate non-cycle edges

  std::cout << "Eliminating non-cycle vertices ... " << std::endl;

  global_start_time =
    std::chrono::steady_clock::now();
 
  eliminate_non_cycle_vertices<Vertex, Edge, Boolean>(vertex_count, vertices, 
    vertex_degree, vertex_in_degree,  edges, vertex_in_edges, vertex_active);

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl; 

  std::cout << std::endl;

  // Test
  {  
  std::cout << "Size of vertex active list: " << vertex_active.size() << std::endl;

  size_t tmp_active_vertex_count = 0;
  size_t tmp_active_edge_count = 0; 

  for (Vertex v = 0; v < vertex_count; v++) {
    if (vertex_active[v] == static_cast<Boolean>(1)) {
      tmp_active_vertex_count++;  
    } else {
      continue;
    }

    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];
      if (vertex_active[v_nbr] == static_cast<Boolean>(1)) {       
        tmp_active_edge_count++;
        //std::cerr << v << " " << v_nbr << std::endl; // TODO: uncomment to output the edge list
        //std::cerr << v_nbr << " " << v << std::endl;  
      }  
    }
 
  }

  std::cout << "Active vertex count: " << tmp_active_vertex_count << std::endl;
  std::cout << "Active edge count: " << tmp_active_edge_count << std::endl;
  }
  // Test
#endif
    
#ifdef SPOKE
  return 0;
#endif

  //////////////////////////////////////////////////////////////////////////////

  // dfs

  //std::cout << "DFS" << std::endl;

  //dfs_single_source<Vertex, Edge>(vertices, vertex_count, vertex_degree, edges,
  //  static_cast<Vertex>(3));  
    //static_cast<Vertex>(17638));
    //static_cast<Vertex>(1579)); 

  //return 0;

  //////////////////////////////////////////////////////////////////////////////

  global_start_time = std::chrono::steady_clock::now();

  std::cout << "BFS" << std::endl;

//  is_cc = bfs_single_source<Vertex, Edge, VertexList, EdgeList>
    //(vertices, vertex_count, edges, static_cast<Vertex>(3)); // reaction - compound
//    (vertices, vertex_count, edges, static_cast<Vertex>(2294390)); // all 
    // 2294390 1771176 	

  std::cout << "Connected component: " << is_cc << std::endl;
  
  //bfs_parallel<Vertex, Edge, VertexList, EdgeList>
    //(vertices, vertex_count, edges);

  global_end_time =  std::chrono::steady_clock::now();
  global_elapsed_time = getElapsedTimeSecond(global_start_time,
    global_end_time);
  std::cout << "Time: " << global_elapsed_time << " seconds."
    << std::endl;

  std::cout << std::endl;

  return 0;

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

  // spoke
  typedef std::vector< std::vector<std::vector<std::tuple<Vertex, Vertex>> > > 
    VertexNonLocalPropertiesVector;
  VertexNonLocalPropertiesVector  vertex_cycles_vector(vertex_count); 
  // spoke


  // vertex_cycles_unique
  //typedef std::vector< std::unordered_map<size_t, std::unordered_map<Edge, 
  //  std::tuple<Vertex, Vertex, size_t> > > > VertexNonLocalPropertiesUnique; 
  typedef std::vector< std::unordered_map<size_t, EdgeSet> > VertexNonLocalPropertiesUnique;
  VertexNonLocalPropertiesUnique vertex_cycles_unique(vertex_count);

  //return 0; // Test
 
  // find cycles 
  //TODO: we are using vertex_cycles_unique, remove vertex_cycles? 
#ifdef SPOKE 
  find_cycles_parallel<Vertex, Edge, VertexData, VertexList, EdgeList, VertexDataList>(vertices,  
    vertex_count, vertex_degree, edges, vertex_data, vertex_cycles, vertex_cycles_unique, vertex_cycles_vector);
#endif
  find_cycles_parallel<Vertex, Edge, Boolean, VertexData, VertexList, EdgeList, 
    VertexDataList>
    (vertices, vertex_count, vertex_degree, edges, vertex_data, vertex_cycles, 
    vertex_cycles_unique, vertex_cycles_vector, vertex_active);

  // Test
  // VertexNonLocalProperties 
  // vertex_cycles
  /*std::cout << "Per-vertex cycles (vector-based path)" << std::endl; 
  for (size_t i = 0; i < vertex_cycles.size(); i++) { // vector 
    //std::cout << vertex_cycles[i].size() << std::endl;
    for (auto& j : vertex_cycles[i]) { // j is a map
      std::cout << i << " | " << j.first << " | ";
      for (auto& k : j.second) { // k is a tuple
        std::cout << "(" << std::get<0>(k) << " -- " << std::get<1>(k) << "), "; 
      } // for  
      std::cout << std::endl;  
    } // for
    std::cout << std::endl;
  } // for*/

  /*std::cout << "Per-vertex cycles (vector-based path)" << std::endl;
  for (size_t i = 0; i < vertex_cycles_vector.size(); i++) { // vector
    for (auto& j : vertex_cycles_vector[i]) { // j is a vector
      if (j.size() > 0) {
        std::cerr << i << " : " << j.size() << " : ";
      }
      for (auto& k : j) { // k is a tuple
        std::cerr << "" << std::get<0>(k) << " " << std::get<1>(k) << ", ";
      }	 // for
      if (j.size() > 0) {   
        std::cerr << std::endl;	
      }
    } // for
    //if (vertex_cycles_vector[i].size() > 0) {  
    //  std::cerr << std::endl;
    //}
  } // for*/ 

  // Test   
  
  return 0;

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


  //////////////////////////////////////////////////////////////////////////////

  return 0;

 
} // end of main
