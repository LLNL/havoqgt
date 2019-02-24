#pragma once

#include <bitset>
#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>

namespace prunejuice { namespace pattern { namespace graphalgorithm {

  constexpr size_t VERTEX_COUNT = 64;
  typedef std::bitset<VERTEX_COUNT> VertexBitSet;
  
  template <typename Vertex, typename Edge>
  Edge edge_hash_vertex_bitset(Vertex s, Vertex t) {
    VertexBitSet vertex_bitset;
    vertex_bitset.set(static_cast<size_t>(s));
    vertex_bitset.set(static_cast<size_t>(t));
    return static_cast<Edge>(vertex_bitset.to_ullong());
  }

  template <typename Vertex, typename EdgeDFSGraph, 
    typename NeighborSetIterator>
  NeighborSetIterator neighbor_set_iterator(EdgeDFSGraph& edge_dfs_graph, 
    Vertex vertex) {
    auto find_vertex = edge_dfs_graph.find(vertex);
    if (find_vertex == edge_dfs_graph.end()) {
      std::cerr << "Error: did not find the expected item in the map." 
        << std::endl;
      return NeighborSetIterator(); 
    } else {    
      return find_vertex->second.begin();
    } 
  }

  // assuming edgelist is undirected
  template <typename EdgeListTuple, typename NeighborSet, typename EdgeDFSGraph>
  EdgeDFSGraph edgelist_to_edge_dfs_graph(EdgeListTuple& edgelist) {

    EdgeDFSGraph edge_dfs_graph;

    for (auto edge : edgelist) {
      auto s = std::get<0>(edge);
      auto t = std::get<1>(edge);

      auto find_s = edge_dfs_graph.find(s);
      if (find_s == edge_dfs_graph.end()) {
        auto insert_status = edge_dfs_graph.insert({s, NeighborSet()});
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the map." 
            << std::endl;
          return EdgeDFSGraph();
        } else {
          find_s = insert_status.first;
          find_s->second.insert(t);  
        } 
      } else {
        auto find_t = find_s->second.find(t);
        if (find_t == find_s->second.end()) {
          find_s->second.insert(t);
        }
      }   
    } // for     
 
    return edge_dfs_graph; 
  }

  template <typename EdgeListTuple, typename EdgeSet, typename NeighborSet, 
    typename EdgeDFSGraph>
  EdgeDFSGraph edgeset_to_edge_dfs_graph(EdgeSet& edgeset) {

    //EdgeDFSGraph edge_dfs_graph;
     
    EdgeListTuple edgelist(0);

    for (auto e : edgeset) {
      auto s = std::get<0>(e.second);
      auto t = std::get<1>(e.second);
 
      // directed EdgeDFSGraph
      /*auto find_s = edge_dfs_graph.find(s);
      if (find_s == edge_dfs_graph.end()) {
        auto insert_status = edge_dfs_graph.insert({s, NeighborSet()});
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the map." 
            << std::endl;
          return EdgeDFSGraph();
        } else {
          find_s = insert_status.first;
          find_s->second.insert(t);  
        } 
      } else {
        auto find_t = find_s->second.find(t);
        if (find_t == find_s->second.end()) {
          find_s->second.insert(t);
        }
      }*/      

      // undirected EdgeDFSGraph 
      // edgeset is undirected    
      edgelist.push_back(std::forward_as_tuple(s,t));
      edgelist.push_back(std::forward_as_tuple(t,s));

    } // for

    return edgelist_to_edge_dfs_graph<EdgeListTuple, NeighborSet, EdgeDFSGraph>
      (edgelist);     

    //return edge_dfs_graph;
  }

  template <typename EdgeDFSGraph> 
  void print_edge_dfs_graph(EdgeDFSGraph& edge_dfs_graph) {
    for (auto i : edge_dfs_graph) { 
      std::cout << i.first << " : ";
      auto edge_it = i.second.begin(); 
      while (edge_it != i.second.end()) {
        std::cout << *edge_it << ", "; 
        edge_it = std::next(edge_it);  
      }
      std::cout << std::endl; 
    }
  }

  template <typename T>
  void print_container(T& container) {
    for (auto i : container) {
      std::cout << i << ", ";
    }
    std::cout << std::endl; 
  } 

  template <typename Vertex, typename Edge, typename EdgeListTuple, 
    typename NeighborSet, typename EdgeDFSGraph>
  void edge_dfs_algorithm(EdgeDFSGraph& edge_dfs_graph, Vertex source_vertex, 
    EdgeListTuple& edge_dfs_output) {
    //std::cout << "Edge DFS Algorithm ... " << std::endl; 

    //typedef typename Graph::EdgeListTuple_t EdgeListTuple_t;

    //typedef std::set<Vertex> NeighborSet;        
    //typedef std::unordered_map<Vertex, NeighborSet> EdgeDFSGraph;
 
    typedef typename NeighborSet::iterator NeighborSetIterator;
    typedef std::tuple<Vertex, NeighborSetIterator> 
      VertexNeighborSetIteratorTuple;     
  
    //EdgeDFSGraph edge_dfs_graph = edgelist_to_edge_dfs_graph
    //  <EdgeListTuple_t, NeighborSet, EdgeDFSGraph>(graph.edgelist);   
   
    //print_edge_dfs_graph(edge_dfs_graph); 
 
    std::unordered_set<Vertex> visited_vertices;
    std::unordered_set<Edge> visited_edges;

    std::unordered_map<Vertex, Vertex> edge_map;
    
    //std::stack<Vertex> vertex_stack;
    //vertex_stack.push(source_vertex);

    Vertex current_vertex = source_vertex;     
  
    std::stack<VertexNeighborSetIteratorTuple> vertex_stack;
   
    vertex_stack.push(std::forward_as_tuple(source_vertex, 
      neighbor_set_iterator<Vertex, EdgeDFSGraph, NeighborSetIterator>
      (edge_dfs_graph, source_vertex)));
    
    while(!vertex_stack.empty()) {

      //std::cout << " > current_vertex " << current_vertex << std::endl;
      //print_container(visited_vertices);
      //print_container(visited_edges); 

      //current_vertex = vertex_stack.top();
      current_vertex = std::get<0>(vertex_stack.top());
     
      auto find_visited_vertex = visited_vertices.find(current_vertex);
      if (find_visited_vertex == visited_vertices.end()) {
        visited_vertices.insert(current_vertex);

      } else {
            
        auto find_current_vertex = edge_dfs_graph.find(current_vertex);
        if (find_current_vertex == edge_dfs_graph.end()) {
          std::cerr << "Error: did not find the expected item in the map."
            << std::endl;
          return;
        } 

        // visit the next neighbor of current_vertex         
        NeighborSetIterator edge_it = std::get<1>(vertex_stack.top());
        if (edge_it == find_current_vertex->second.end()) {
          // all neighbors of current_vertex has been visited
          vertex_stack.pop();  

        } else {

          // TODO: remove 64 vertex limit
          Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>
            (current_vertex, *edge_it);           

          auto find_edge_hash = visited_edges.find(edge_hash);
          if (find_edge_hash == visited_edges.end()) {  
            visited_edges.insert(edge_hash);

            std::get<1>(vertex_stack.top()) = std::next(edge_it); // update 
        
            vertex_stack.push(std::forward_as_tuple(*edge_it,
              neighbor_set_iterator<Vertex, EdgeDFSGraph, NeighborSetIterator>
              (edge_dfs_graph, *edge_it)));

            //std::cout << "(" << current_vertex << ", " << *edge_it << "), ";
            //std::cout << " > vertex_stack.size() " << vertex_stack.size() 
            //  << std::endl;
            edge_dfs_output.push_back(std::forward_as_tuple(current_vertex, 
              *edge_it)); 
  
          } else { 
            std::get<1>(vertex_stack.top()) = std::next(edge_it); // update
          }   

        } // if 

      } // if       
 
    } // while 

    //std::cout << std::endl; 
 
  } // edge_dfs

  template <typename Graph>
  typename Graph::EdgeListTuple_t edge_dfs(Graph& graph, 
    typename Graph::Vertex_t source_vertex) {

    //std::cout << "Edge DFS Algorithm ... " << std::endl; 

    typedef typename Graph::Vertex_t Vertex;
    typedef typename Graph::Edge_t Edge; 
    typedef typename Graph::EdgeListTuple_t EdgeListTuple;

    typedef std::set<Vertex> NeighborSet;        
    typedef std::unordered_map<Vertex, NeighborSet> EdgeDFSGraph;
 
    //typedef typename NeighborSet::iterator NeighborSetIterator;
    //typedef std::tuple<Vertex, NeighborSetIterator> 
    //  VertexNeighborSetIteratorTuple;     
  
    EdgeDFSGraph edge_dfs_graph = edgelist_to_edge_dfs_graph
      <EdgeListTuple, NeighborSet, EdgeDFSGraph>(graph.edgelist);   
   
    //if (source_vertex == 0) {
    //  print_edge_dfs_graph(edge_dfs_graph); 
    //}

    EdgeListTuple edge_dfs_output(0);

    edge_dfs_algorithm<Vertex, Edge, EdgeListTuple, NeighborSet, EdgeDFSGraph>
      (edge_dfs_graph, source_vertex, edge_dfs_output);   

    return edge_dfs_output;
 
  } // edge_dfs

  template <typename Vertex, typename Edge, typename EdgeListTuple, 
    typename EdgeSet>
  EdgeListTuple edge_dfs(EdgeSet& edgeset, Vertex source_vertex) {

    //std::cout << "Edge DFS Algorithm ... " << std::endl;
    
    typedef std::set<Vertex> NeighborSet;
    typedef std::unordered_map<Vertex, NeighborSet> EdgeDFSGraph;

    EdgeListTuple edge_dfs_output(0);

    EdgeDFSGraph edge_dfs_graph = edgeset_to_edge_dfs_graph
      <EdgeListTuple, EdgeSet, NeighborSet, EdgeDFSGraph>(edgeset);

    //if (source_vertex == 0) {
    //print_edge_dfs_graph(edge_dfs_graph);
    //}    

    edge_dfs_algorithm<Vertex, Edge, EdgeListTuple, NeighborSet, EdgeDFSGraph>
      (edge_dfs_graph, source_vertex, edge_dfs_output);

    return edge_dfs_output;
  
  } // edge_dfs  

}}} // end namespace prunejuice::pattern::graph_algorithm
