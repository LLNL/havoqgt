#pragma once

#include "basic_graph.hpp"
#include "template_constraint.hpp"

namespace prunejuice { namespace pattern {

template <typename Vertex, typename Edge, typename VertexData,
  typename VertexList, typename EdgeList, typename EdgeListTuple, 
  typename EdgeListMap, typename VertexDataList>
class template_graph : public basic_graph<Vertex, Edge, VertexList,
    EdgeList, EdgeListTuple, EdgeListMap> {

public:
 
  typedef basic_graph<Vertex, Edge, VertexList,
    EdgeList, EdgeListTuple, EdgeListMap> basic_graph_t;

  //static constexpr size_t VERTEX_COUNT = 64;
  
  typedef Vertex Vertex_t;
  typedef Edge Edge_t;  
  typedef VertexData VertexData_t;

  typedef EdgeListTuple EdgeListTuple_t;
  typedef EdgeListMap EdgeListMap_t;
  typedef EdgeList EdgeList_t;
  typedef VertexList VertexList_t;

  typedef std::vector<EdgeListTuple_t> EdgeListTupleVector_t;

  typedef std::unordered_map<Edge_t, std::tuple<Vertex_t, Vertex_t, size_t>> EdgeSet_t; // TODO: ?
  typedef std::vector<EdgeSet_t> EdgeSetVector_t; // TODO: ?
  
  typedef prunejuice::pattern::template_constraint<Vertex_t, Edge_t, EdgeSet_t,
    EdgeSetVector_t, EdgeListTupleVector_t> TemplateConstraint_t; // TODO: ?
  typedef std::vector<TemplateConstraint_t> TemplateConstraintVector_t; // TODO: ?

  typedef std::unordered_map<VertexData_t, std::unordered_map<Edge_t, 
    std::tuple<Vertex_t, Vertex_t>>> VetexDataVertexPairs_t; // TODO: ?

  template_graph() : is_null(true) {}

  //template_graph(EdgeListTuple _edgelist) : 
  //  edgelist(_edgelist), edgelist_unique(0), edges(0), vertices(0), 
  //  vertex_degree(0), template_constraints(0), is_null(false) {
  //}

  template_graph(EdgeListTuple _edgelist, Vertex _max_vertex = 0) : 
    edgelist(_edgelist), edgelist_unique(0), edges(0), vertices(0),
    vertex_degree(0), template_constraints(0), template_nonlocal_constraints(0), 
    max_vertex(_max_vertex), is_null(false) {
   
    basic_graph_t::generate_graph(vertex_count, vertices, vertex_degree, edges, 
      edgelist, edgelist_unique, max_vertex); 
 
    edge_count = static_cast<Edge>(edges.size());
    max_degree = static_cast<Edge>(0);
  }

  template_graph(EdgeListTuple _edgelist, VertexDataList _vertex_data, 
    Vertex _max_vertex = 0) : 
    edgelist(_edgelist), edgelist_unique(0), edges(0), vertices(0),
    vertex_degree(0), template_constraints(0), template_nonlocal_constraints(0),
    max_vertex(_max_vertex), vertex_data(_vertex_data), vertex_data_vertex_pairs(0),  
    is_null(false) {
   
    basic_graph_t::generate_graph(vertex_count, vertices, vertex_degree, edges, 
      edgelist, edgelist_unique, max_vertex); 
 
    edge_count = static_cast<Edge>(edges.size());
    max_degree = static_cast<Edge>(0); // TODO:

    generate_vertex_pairs_with_identical_metadata();
  }

  template_graph(EdgeListTuple _edgelist, EdgeListMap _edgelist_unique, 
    EdgeList _edges, VertexList _vertices, VertexList _vertex_degree) :
    edgelist(_edgelist), edgelist_unique(_edgelist_unique), edges(_edges),
    vertices(_vertices), vertex_degree(_vertex_degree), template_constraints(0), 
    template_nonlocal_constraints(0), is_null(false) {
  }

  template_graph(EdgeListTuple _edgelist, EdgeListMap _edgelist_unique,
    EdgeList _edges, Edge _edge_count, VertexList _vertices, 
    Vertex _vertex_count, Vertex _max_vertex, VertexList _vertex_degree) :
    edgelist(_edgelist), edgelist_unique(_edgelist_unique), edges(_edges),
    edge_count(_edge_count), vertices(_vertices), vertex_count(_vertex_count),
    max_vertex(_max_vertex), vertex_degree(_vertex_degree), 
    template_constraints(0), template_nonlocal_constraints(0), is_null(false) {
  }

#ifdef ENABLE_BLOCK
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
#endif
  ~template_graph() {}

  EdgeListTuple edgelist;
  EdgeListMap edgelist_unique;
  EdgeList edges;
  VertexList vertices;
  VertexList vertex_degree;
  VertexDataList vertex_data;   
  
  Edge edge_count;
  Vertex vertex_count;
  Vertex max_vertex;
  Edge max_degree;
  const bool is_null;

  VetexDataVertexPairs_t vertex_data_vertex_pairs;  

  // the list of template constraints
  std::vector<size_t> template_constraints; // TODO: ?
  TemplateConstraintVector_t template_nonlocal_constraints;

private:

  /**
   * Identify vertices with identical metadata      
   */
  void generate_vertex_pairs_with_identical_metadata() {
    assert(vertex_data.size() > 0);
    for (size_t v = 0; v < vertex_data.size(); v++) {
      auto find_vertex_data = vertex_data_vertex_pairs.find(vertex_data[v]);
      if (find_vertex_data == vertex_data_vertex_pairs.end()) {
        vertex_data_vertex_pairs.insert({vertex_data[v], std::unordered_map<Edge, std::tuple<Vertex, Vertex>>()});
      }  
    } // for

    for (size_t i = 0; i < vertex_data.size() - 1; i++) {
      Vertex u = static_cast<Vertex>(i);          
      for (size_t j = i + 1; j < vertex_data.size(); j++) {            
        if (i == j) { 
          continue;  
        } else if (i < j) {
          Vertex v = static_cast<Vertex>(j); 
          if (vertex_data[i] == vertex_data[j]) {

            std::bitset<64> vertex_pair_bitset; // TODO: write a function to generate edge_hash 
	    vertex_pair_bitset.set(i);	
            vertex_pair_bitset.set(j);
            Edge vertex_pair_hash = static_cast<Edge>(vertex_pair_bitset.to_ullong());;

            auto find_vertex_pair = vertex_data_vertex_pairs[vertex_data[i]].find(vertex_pair_hash);      
            if (find_vertex_pair == vertex_data_vertex_pairs[vertex_data[i]].end()) {               
              vertex_data_vertex_pairs[vertex_data[i]].insert({vertex_pair_hash, std::forward_as_tuple(u, v)}); 
            } 
          }        
        }    
      } // for 
    } // for
  } 
 
#ifdef ENABLE_BLOCK
private:

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
//      generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
//        (edge_list, edge_list_unique);
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
#endif
};

}} // end namespace prunejuice::pattern
