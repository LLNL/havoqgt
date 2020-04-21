#pragma once

namespace prunejuice { namespace pattern {

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename EdgeListTuple, typename EdgeListMap>
class basic_graph {

public:

  static constexpr size_t VERTEX_COUNT = 64;
  typedef std::bitset<VERTEX_COUNT> VertexBitSet;

  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> EdgeSet; // TODO: ?

  //basic_graph();
  //~basic_graph();
  
  //template <typename Vertex, typename Edge>
  static Edge edge_hash_vertex_bitset(Vertex s, Vertex t) {
    VertexBitSet vertex_bitset;
    vertex_bitset.set(static_cast<size_t>(s));
    vertex_bitset.set(static_cast<size_t>(t));
    return static_cast<Edge>(vertex_bitset.to_ullong());
  }

  template <typename EdgeListFilter>
  static EdgeListTuple filter_edge_list(EdgeListTuple& edge_list,
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

  template <typename EdgeSet> 
  static bool is_subset_edgeset(EdgeListMap& a, EdgeSet& b) {
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

  //template <typename EdgeSet>
  // Important: lvalue const references could bind to rvalues
  static EdgeSet edgelisttuple_to_edgeset(EdgeListTuple const& edge_list) {
    EdgeSet edgeset(0);

    size_t r = 1;
    for (auto edge : edge_list) {
      Vertex s = std::get<0>(edge);
      Vertex t = std::get<1>(edge);

      Edge edge_hash = edge_hash_vertex_bitset(s, t);

      auto find_edge_hash = edgeset.find(edge_hash);    
      if (find_edge_hash == edgeset.end()) {
        edgeset.insert({edge_hash, std::forward_as_tuple(s, t, r)});
        //update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered)); // TODO: ?
        r++; 
      }  
    } 

    assert(edge_list.size() == edgeset.size());   

    return edgeset;
  }    

  // no edgelist_unique
  void generate_graph(Vertex& vertex_count, VertexList& vertices,
    VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
    Vertex& max_vertex) {

    // TODO: imporve
    //typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> T;
    //T dummy_edge_list_unique;
    EdgeListMap dummy_edge_list_unique;
    dummy_edge_list_unique.insert({0, std::forward_as_tuple(0, 0)});

    generate_graph(vertex_count, vertices, vertex_degree, edges, edge_list,
      dummy_edge_list_unique, max_vertex);
  }

  /**
   * Generate CSR graph from a given edgelist
   */
  // TODO: this function assumes the input edgelist is undirected? 
  // Note: it also works for a directed edgelist
  void generate_graph(Vertex& vertex_count, VertexList& vertices,
    VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
    EdgeListMap& edge_list_unique, Vertex& max_vertex) {

    //Vertex max_vertex;
    //Edge edge_count = edge_list.size();

    // sort edges by source
    //std::cout << "Sorting edges by source vertex ..." << std::endl;
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

    // generate vertex list
    //std::cout << "Generating vertex list ..." << std::endl;
    //std::cout << "Creating CSR row pointers ..." << std::endl;
    vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree,
      max_vertex);
    //std::cout << "Size of vertex list: " << vertices.size() << std::endl;
    //std::cout << "Size of vertex degree list: " << vertex_degree.size() << std::endl;

    //sort targets in increasing order // TODO: move this inside generate_edge_list  
    //std::cout << "Sorting neighbors ..." << std::endl;
    {
    //#pragma omp parallel for // TODO: uncomment
    for (size_t v = 0; v < vertices.size() - 1; v++) {
      size_t start = vertices[v];
      size_t end = vertices[v + 1];
      std::stable_sort(edge_list.begin() + start,
                       edge_list.begin() + end,
        [](const std::tuple<Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
    } // for
    }

    // generate edge list
    //std::cout << "Generating edge list ..." << std::endl;
    generate_edge_list(edge_list, edges);
    //std::cout << "Size of edge list: " << edges.size() << std::endl;

    if (edge_list_unique.size() < 1) {
      //std::cout << "Generating unique edges ..." << std::endl;
      //generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
      //  (edge_list, edge_list_unique);
      generate_unique_edge_list(edge_list, edge_list_unique);  
      //std::cout << "Number of unique edges: " << edge_list_unique.size()
      //  << std::endl;
    }

    //  std::cout << "CSR Graph generation completed." << std::endl;
    //  std::cout << "Number of vertices: " << vertex_count << std::endl;
    //  std::cout << "Number of edges: " << edge_count << std::endl;
    //  std::cout << "Max vertex: " << max_vertex << std::endl;         
  }

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

  void generate_edge_list(EdgeListTuple& edge_list, EdgeList& edges) {
    for (size_t e = 0; e < edge_list.size(); e++) {
      edges.push_back(std::get<1>(edge_list[e]));
    } // for
  }

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
    } // for
  }

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
    
};

}} // end namespace prunejuice::utilities
