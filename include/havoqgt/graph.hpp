#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <vector>
#include <tuple>

#include <boost/algorithm/string.hpp>

template<typename Vertex, typename Edge, typename VertexData>
class graph {
  public:
    graph(std::string vertex_input_filename,  std::string edge_input_filename, 
      std::string vertex_data_input_filename, const bool _directed = true) :
      directed(_directed), 
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0) {
      std::cout << "Building CSR graph ... " << std::endl;
      
      std::cout << "Reading vertex list ... " << std::endl;
      vertex_count = read_vertex_list(vertex_input_filename, vertices,
        vertex_degree);

      std::cout << "Reading vertex data list ... " << std::endl;
      read_vertex_data_list(vertex_data_input_filename, vertex_data);

      std::cout << "Reading edge list ... " << std::endl;
      edge_count = read_edge_list(edge_input_filename, edges);

      std::cout << "Completed building graph." << std::endl;
      output_stat(); 
    }

    graph(std::string edge_input_filename, std::string vertex_input_filename,
      std::string vertex_data_input_filename, 
      std::string stat_input_filename,
      const bool _directed = true, const bool _mutable = false) :
      directed(_directed),
      vertex_count(0),
      edge_count(0),
      diameter(0),  
      vertices(0),
      vertex_degree(0),
      vertex_data(0), 
      edges(0), 
      edge_list(0) {
//      std::cout << "Building CSR graph ... " << std::endl;
      
//      std::cout << "Reading edge list ... " << std::endl;
      edge_count = read_edge_list(edge_input_filename);

      //std::cout << "Generating vertex list ..." << std::endl;
      vertex_count = generate_vertex_list() ;       

//      std::cout << "Reading vertex list ... " << std::endl;
//      vertex_count = read_vertex_list(vertex_input_filename);

//      std::cout << "Reading vertex data list ... " << std::endl;
      read_vertex_data_list(vertex_data_input_filename);

      read_stat(stat_input_filename);

//      std::cout << "Completed building graph." << std::endl;
//      output_stat(); 
    }

    ~graph() {
//      std::cout << "Disposing graph ... " << std::endl;
    }    

    const bool directed; 
    Vertex vertex_count;
    Edge edge_count;
    Edge diameter;  
    std::vector<Edge> vertices;
    std::vector<Edge> vertex_degree;
    std::vector<VertexData> vertex_data;
    std::vector<Vertex> edges;
    std::vector<std::tuple<Vertex, Vertex>> edge_list; 

  private:
    Vertex read_vertex_list(std::string vertex_input_filename) { 
      std::ifstream vertex_input_file(vertex_input_filename, std::ifstream::in);
      std::string line;
      while (std::getline(vertex_input_file, line)) {
        std::istringstream iss(line);
        Vertex v_source(0);
        Edge v_degree(0), v_offset(0);
        iss >> v_source >> v_degree >> v_offset;
        vertices.push_back(v_offset);
        vertex_degree.push_back(v_degree);
      }
      vertex_input_file.close();
      vertex_degree.erase(vertex_degree.end() - 1);
      return vertices.size() <= 0 ? 0 : vertices.size() - 1; // vertex count
    }

    void read_vertex_data_list(std::string vertex_data_input_filename) {
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

    Edge read_edge_list(std::string edge_input_filename) {
      std::ifstream edge_input_file(edge_input_filename, std::ifstream::in);
      std::string line;
      while(std::getline(edge_input_file, line)) {
        std::istringstream iss(line);
        Vertex s(0), t(0);
        iss >> s >> t;
        edges.push_back(t);
        edge_list.push_back(std::forward_as_tuple(s, t));
      }
      edge_input_file.close();
      return edges.size(); // edge count
    } 

    Vertex generate_vertex_list() {
      Vertex vertex_count = 0;
      Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
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
      //<< "\n";
      //vertex_file.close();   
      //vertex_data_file.close();

      return vertex_count;  
    }

    void read_stat(std::string stat_input_filename) {
      std::ifstream stat_input_file(stat_input_filename, 
        std::ifstream::in);
      std::string line;
      const char delim = ':';
      while(std::getline(stat_input_file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while(std::getline(iss, token, delim)) {
          tokens.push_back(token);
        }          
        //std::cout << tokens[0] << " " << tokens[1] << std::endl; 
        assert(tokens.size() > 1);
        boost::trim(tokens[0]);
        boost::trim(tokens[1]); 
        if (boost::iequals(tokens[0], "diameter")) {
          diameter = std::stoull(tokens[1]);         
        } 
      }
      stat_input_file.close();
    } 

    void output_stat() {
      std::cout << "Number of vertices: " << vertex_count << std::endl;
      std::cout << "Number of edges: " << edge_count << std::endl;
      std::cout << "Number of vertex data: " << vertex_data.size() << std::endl;
    }
};
