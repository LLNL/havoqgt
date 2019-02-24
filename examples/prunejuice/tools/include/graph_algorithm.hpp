#pragma once

namespace prunejuice { namespace pattern {

class graph_algorithm {

public:

  /**
   * Is the given graph a connected component
   * This routine is called by multiple graphs is parallel, so the routine 
   * itself is siquential 
   */
  template <typename Graph>
  static bool is_connected_component(Graph& graph) {
    return is_connected_component<typename Graph::Vertex_t, 
      typename Graph::Edge_t, typename Graph::VertexList_t, 
      typename Graph::EdgeList_t> (graph.vertices, graph.vertex_count, 
      graph.edges); 
  } 

  template <typename Vertex, typename Edge, typename VertexList,
    typename EdgeList>
  static bool is_connected_component(VertexList& vertices, Vertex vertex_count,
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
    } // for

    return is_cc;
  }

  template <typename Vertex, typename Edge, typename Graph, 
    typename OrderedPathEdges>
  static void find_path_recursive(Graph& graph, Vertex source, Vertex target, 
    uint64_t target_bfs_level, OrderedPathEdges& shortest_path_edges, Vertex v, 
    OrderedPathEdges walk_history_edges, uint64_t r, bool& finished) {

    if (r > target_bfs_level || finished) {
      return;
    }

    for (Edge e = graph.vertices[v]; e < graph.vertices[v + 1]; e++) {
      Vertex v_nbr = graph.edges[e];

      if (v_nbr == source) {
        //return; 
      } else if ((v_nbr != source || v_nbr != target) && r < target_bfs_level) {  
        // forward
        OrderedPathEdges new_walk_history_edges(walk_history_edges);
        new_walk_history_edges.push_back(std::forward_as_tuple(v, v_nbr));        

        find_path_recursive<Vertex, Edge, Graph, OrderedPathEdges>(graph, source,
          target, target_bfs_level, shortest_path_edges, v_nbr, 
          new_walk_history_edges, r + 1, finished);  
 
      } else if (v_nbr == target && r != target_bfs_level) {
        //return;  
      } else if (v_nbr == target && r == target_bfs_level) {
        // path found

        OrderedPathEdges new_walk_history_edges(walk_history_edges);
        new_walk_history_edges.push_back(std::forward_as_tuple(v, v_nbr));

        //std::cout << "Found Path." << std::endl;   
        //for (auto i : new_walk_history_edges) {
        //  std::cout << "(" << std::get<0>(i) << ", " << std::get<1>(i) << "), "; 
        //}
        //std::cout << std::endl; 

        shortest_path_edges = std::move(new_walk_history_edges);

        finished = true;

        //return;
      }   
 
    } // for 
  } 

  template <typename Vertex, typename Edge, typename OrderedPathEdges, 
    typename Graph>
  static OrderedPathEdges find_path(Graph& graph, Vertex source, Vertex target, 
    uint64_t target_bfs_level) {
  
    uint64_t r = 1; // walk step, initialized to 1

    // vector based walk history of edges
    //typedef std::vector<std::tuple<Vertex, Vertex>> OrderedPathEdges;
    OrderedPathEdges walk_history_edges(0);
    OrderedPathEdges shortest_path_edges(0); 
 
    bool finished = false;

    find_path_recursive<Vertex, Edge, Graph, OrderedPathEdges>(graph, source, 
      target, target_bfs_level, shortest_path_edges, source, walk_history_edges, 
      r, finished);

    assert (shortest_path_edges.size() == target_bfs_level); 

    return shortest_path_edges;  
  } 

  template <typename Vertex, typename Edge, typename EdgeListTuple, 
    typename Graph>
  static EdgeListTuple get_shortest_path(Graph& graph, Vertex source, 
    Vertex target) {

    assert(source < graph.vertex_count);
    assert(target < graph.vertex_count);

    std::vector<uint64_t> bfs_level(graph.vertex_count); 
    std::fill(bfs_level.begin(), bfs_level.end(), 
      std::numeric_limits<uint64_t>::max());

    std::vector<uint8_t> visited(graph.vertex_count);
    std::fill(visited.begin(), visited.end(), 0);
     
    visited[static_cast<size_t>(source)] = 1;
    bfs_level[static_cast<size_t>(source)] = 0;  

    bool finished = false;

    // BFS 
    do {
      finished = true;
      for (Vertex v = 0; v < graph.vertex_count; v++) {
        //std::cout << v << " : ";
        if (visited[static_cast<size_t>(v)] == 1) {
          uint64_t v_level = bfs_level[static_cast<size_t>(v)]; 

          for (Edge e = graph.vertices[v]; e < graph.vertices[v + 1]; e++) {
            Vertex v_nbr = graph.edges[e];
            //std::cout << v_nbr << ", ";
            
            uint64_t v_nbr_level = bfs_level[static_cast<size_t>(v_nbr)];
            
            if (visited[static_cast<size_t>(v_nbr)] == 0) {
              visited[static_cast<size_t>(v_nbr)] = 1;              
              finished = false;
            } 

            if (v_level < v_nbr_level) {
               bfs_level[static_cast<size_t>(v_nbr)] = v_level + 1;
               if (!finished) {
                 finished = false;		
               } 
            }  

          } // for 
        } // if  
        //std::cout << std::endl;
      } // for
    } while(!finished);

    // Test
    /*std::cout << "BFS Levels" << std::endl;  
    for (size_t v = 0; v < bfs_level.size(); v++) {
      std::cout << v << " " << bfs_level[v] << std::endl;	
    }
    std::cout << std::endl;*/ 
    // Test 
   
    // get the (or a) shortest path from source to target
    return find_path<Vertex, Edge, EdgeListTuple, Graph>(graph, source, target, 
      bfs_level[static_cast<size_t>(target)]);

  } 

};

}} // end namespace prunejuice::pattern 
