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
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "omp.h"

#include <assert.h>
#include <unistd.h>

#include <boost/json/serialize.hpp>

//#include <iostream>
//#include <string>
#include <metall/metall.hpp>
#include <metall/container/experimental/jgraph/jgraph.hpp>

using namespace metall::container::experimental;

using graph_type = jgraph::jgraph<metall::manager::allocator_type<std::byte>>;

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

template <typename Graph, typename Vertex>
bool bfs_level_synchronous_single_source(Graph* graph, 
  const Vertex& source_veretx, const Vertex& target_vertex) {

  size_t current_level = 0;   
  std::unordered_map<Vertex, size_t> bfs_level;
  bfs_level[source_veretx] = current_level;    

  bool finished = false;
  bool target_vertex_found = false; 

  do {

    finished = true; 

    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); 
      vitr != vend; ++vitr) {

      //std::cout << "Vertex value = " << vitr->value() << std::endl;
      auto& v_value = vitr->value();
      //std::cout << "Vertex value = " << v_value << std::endl;
      auto v_id = v_value.as_object()["id"].as_string();
      //std::cout << v_id << std::endl;
      //auto v_label = v_value.as_object()["type"].as_string();
      //auto v_label = v_value.as_object()["labels"].as_array()[0].as_string(); 
      //if (v_label == "Food") {
      //  std::cout << v_label << std::endl;    
      //}
      
      // TODO: do not use vitr->key(); use v_value.as_object()["id"]   

      //if (vitr->key() == source_veretx) {
        //std::cout << source_veretx << std::endl;
      //}
  
      //auto find_vertex = bfs_level.find(vitr->key());
      auto find_vertex = bfs_level.find(v_id);
      if (find_vertex == bfs_level.end()) {
        continue;  
      } else {
        //std::cout << vitr->key() << " " << bfs_level[vitr->key()] << std::endl;          
        //std::cout << v_id << " " << bfs_level[v_id] << std::endl; 
      }  

      //if (bfs_level[vitr->key()] == current_level) {
      if (bfs_level[v_id] == current_level) {

        //size_t nbr_count = 0;  

        for (auto eitr = graph->edges_begin(vitr->key()), 
          eend = graph->edges_end(vitr->key()); eitr != eend; ++eitr) {
          //std::cout << "Edge ID = " << eitr->key() << std::endl; 
          //std::cout << "Edge value = " << eitr->value() << std::endl;
          //nbr_count++;

          auto& e_value = eitr->value();
          //auto v_nbr = e_value.as_object()["start"].as_object()["id"].as_string();
          auto v_nbr = e_value.as_object()["end"].as_object()["id"].as_string();
          //std::cout << v_nbr << std::endl;

          auto find_v_nbr = bfs_level.find(v_nbr);
          if (find_v_nbr == bfs_level.end()) {
            bfs_level[v_nbr] = current_level + 1;
            //std::cout << bfs_level[v_nbr] << std::endl; 

            if (finished == true) {
              finished = false;
            }
          }

          if (v_nbr == target_vertex) {
            //std::cout << v_nbr << std::endl;
            target_vertex_found = true;
            break;
          }
  
        } // for  

        //std::cout << vitr->key() << " # neighbors: " << nbr_count << std::endl;

      }

      if (target_vertex_found) {
        finished = true; 
        break;
      } 

    } // for

    //std::cout << current_level << " " << target_vertex_found << " " << 
    //  finished << std::endl; 

    current_level++;    

  } while(!finished);

  if (target_vertex_found) {
    std::cout << target_vertex << " found at depth " << 
      bfs_level[target_vertex] << std::endl;
  } else {
    std::cout << target_vertex << " was not found at maximum depth " << 
      (current_level - 1) << std::endl;
  } 

  std::cout << "#visited vertices: " << bfs_level.size() << std::endl;

  //for (auto i : bfs_level) {
  //  std::cout << i.first << " " << i.second << std::endl; 
  //}

  return target_vertex_found; 
} // bfs_level_synchronous_single_source  

template <typename Graph, typename SearchResult>
bool vertex_keyword_search(Graph* graph, const std::string label, 
  const std::string keyword, SearchResult& vertex_search_result) {

  size_t current_level = 0;   
  //std::unordered_map<Vertex, size_t> bfs_level;
  //bfs_level[source_veretx] = current_level;    

  bool finished = false;
  bool target_vertex_found = false; 

  size_t tmp_count = 0;

  size_t match_count = 0; 

  std::string regex_string = keyword;  
  //std::string regex_string = "[a-z_][a-z_0-9]*\\.[a-z0-9]+";
  //std::string regex_string = "[a-z0-9_-]+";
  //std::string regex_string = "(pou)(.*)";
  //std::string regex_string = "DOZFQ-NTUHNPA"; // works  
  //std::regex regex_object(regex_string, std::regex_constants::icase);
  std::regex regex_object(regex_string, std::regex_constants::ECMAScript | std::regex_constants::icase); 

  //do {

    finished = true; 

    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); 
      vitr != vend; ++vitr) {

      //std::cout << "Vertex value = " << vitr->value() << std::endl;
      auto& v_value = vitr->value();
      //std::cout << "Vertex value = " << v_value << std::endl;
      //auto v_id = v_value.as_object()["id"].as_string();
      //std::cout << v_id << std::endl;
      //auto v_label = v_value.as_object()["type"].as_string();
      auto v_label = v_value.as_object()["labels"].as_array()[0].as_string(); 
      //if (v_label == "Food") {
      //  std::cout << v_label << std::endl;    
      //}
      if (v_label != label) {
        continue;
      } 

      auto v_id = v_value.as_object()["id"].as_string();
      
      std::string v_value_serialized = serialize(v_value);
      //std::cout << v_value_serialized << std::endl;

      bool match_found = std::regex_search(v_value_serialized, regex_object, std::regex_constants::match_any);
      //std::cout << match_found << std::endl;
      if (match_found) {
        //std::cout << v_id << std::endl;
        vertex_search_result[v_id] = v_value_serialized; 
        //match_count++;
      } 

      // TODO: do not use vitr->key(); use v_value.as_object()["id"]   

      //if (vitr->key() == source_veretx) {
        //std::cout << source_veretx << std::endl;
      //}
  
      /*//auto find_vertex = bfs_level.find(vitr->key());
      auto find_vertex = bfs_level.find(v_id);
      if (find_vertex == bfs_level.end()) {
        continue;  
      } else {
        //std::cout << vitr->key() << " " << bfs_level[vitr->key()] << std::endl;          
        //std::cout << v_id << " " << bfs_level[v_id] << std::endl; 
      }*/  

      /*//if (bfs_level[vitr->key()] == current_level) {
      if (bfs_level[v_id] == current_level) {

        //size_t nbr_count = 0;  

        for (auto eitr = graph->edges_begin(vitr->key()), 
          eend = graph->edges_end(vitr->key()); eitr != eend; ++eitr) {
          //std::cout << "Edge ID = " << eitr->key() << std::endl; 
          //std::cout << "Edge value = " << eitr->value() << std::endl;
          //nbr_count++;

          auto& e_value = eitr->value();
          //auto v_nbr = e_value.as_object()["start"].as_object()["id"].as_string();
          auto v_nbr = e_value.as_object()["end"].as_object()["id"].as_string();
          //std::cout << v_nbr << std::endl;

          auto find_v_nbr = bfs_level.find(v_nbr);
          if (find_v_nbr == bfs_level.end()) {
            bfs_level[v_nbr] = current_level + 1;
            //std::cout << bfs_level[v_nbr] << std::endl; 

            if (finished == true) {
              finished = false;
            }
          }

          if (v_nbr == target_vertex) {
            //std::cout << v_nbr << std::endl;
            target_vertex_found = true;
            break;
          }
  
        } // for  

        //std::cout << vitr->key() << " # neighbors: " << nbr_count << std::endl;

      }*/

      /*if (target_vertex_found) {
        finished = true; 
        break;
      }*/

      //tmp_count++;
      //if (tmp_count > 10) {
        //break;
      //}     

    } // for

    //std::cout << current_level << " " << target_vertex_found << " " << 
    //  finished << std::endl; 

    //current_level++;    

  //} while(!finished);

  /*if (target_vertex_found) {
    std::cout << target_vertex << " found at depth " << 
      bfs_level[target_vertex] << std::endl;
  } else {
    std::cout << target_vertex << " was not found at maximum depth " << 
      (current_level - 1) << std::endl;
  } 

  std::cout << "#visited vertices: " << bfs_level.size() << std::endl;

  //for (auto i : bfs_level) {
  //  std::cout << i.first << " " << i.second << std::endl; 
  //}*/

  //std::cout << "# matches found: " << match_count << std::endl; 

  return target_vertex_found; 
} // bfs_level_synchronous_single_source  

int main() {
#ifdef ENABLE_BLOCK
  // Assumes that the format of the input file is JSON Lines
  //const std::string input_json_file_name("/usr/workspace/reza2/metall/metalljson/spoke.json");//("../spoke.json");
  const std::string input_json_file_name("/p/lustre2/havoqgtu/Spoke/spoke-20201130.json");

  {
    std::cout << "--- Create ---" << std::endl;
    metall::manager manager(metall::create_only, "/usr/workspace/reza2/metall/jgraph_spoke_undir_obj");

    std::ifstream ifs(input_json_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open file: " << input_json_file_name << std::endl;
      std::abort();
    }

    std::string json_string;

    auto *graph = manager.construct<graph_type>(metall::unique_instance)(manager.get_allocator());

    // Parse each line of the input file one by one
    while (std::getline(ifs, json_string)) {

      try {

      // Parse the JSON string and allocate a JSON value object.
      // Pass a Metall allocator so that the contents of the object is allocated in Metall space.
      auto json_value = json::parse(json_string, graph->get_allocator());

      if (json_value.as_object()["type"].as_string() == "node") {
        const auto &vertex_id = json_value.as_object()["id"].as_string();
        graph->register_vertex(vertex_id);
        graph->vertex_value(vertex_id) = std::move(json_value);
      } else if (json_value.as_object()["type"].as_string() == "relationship") {
        const auto &src_id = json_value.as_object()["start"].as_object()["id"].as_string();
        const auto &dst_id = json_value.as_object()["end"].as_object()["id"].as_string();
        const auto &edge_id = json_value.as_object()["id"].as_string();
        graph->register_edge(src_id, dst_id, edge_id);
        graph->register_edge(dst_id, src_id, edge_id); // uncomment to create undirected a graph
        graph->edge_value(edge_id) = std::move(json_value);
      }

      } catch (const std::exception& e) {
        std::cout << " a standard exception was caught, with message '" << e.what() << std::endl;
      }

    } // while

    std::cout << "#of vertices: " << graph->num_vertices() << std::endl;
    std::cout << "#of edges: " << graph->num_edges() << std::endl;
  }
#endif

  {
    std::cout << "\n--- Open ---" << std::endl;
    //metall::manager manager(metall::open_read_only, "/usr/workspace/reza2/metall/jgraph_spoke_undir_obj");
    metall::manager manager(metall::open_read_only, "/dev/shm/jgraph_spoke_undir_obj_2");

    const auto *graph = manager.find<graph_type>(metall::unique_instance).first;

    std::cout << "#of vertices: " << graph->num_vertices() << std::endl;
    std::cout << "#of edges: " << graph->num_edges() << std::endl;

    /*std::cout << "<Vertices>" << std::endl;
    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); vitr != vend; ++vitr) {
      std::cout << "Vertex ID = " << vitr->key() << std::endl;
      std::cout << "Vertex value = " << vitr->value() << std::endl;
    }*/

    /*size_t tmp_count = 0;

    // Access vertex values and edge values using the iterators
    std::cout << "\n<Edges>" << std::endl;
    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); vitr != vend; ++vitr) {
      std::cout << "Vertex ID = " << vitr->key() << std::endl;
      std::cout << "Vertex value = " << vitr->value() << std::endl; 
      for (auto eitr = graph->edges_begin(vitr->key()), eend = graph->edges_end(vitr->key()); eitr != eend; ++eitr) {
        std::cout << "Edge ID = " << eitr->key() << std::endl;
        std::cout << "Edge value = " << eitr->value() << std::endl;

        //json_value.as_object()["end"].as_object()["id"].as_string();
        auto& e_value = eitr->value();
        auto v_nbr_id = e_value.as_object()["end"].as_object()["id"].as_string();  
        std::cout << v_nbr_id << std::endl; 

        tmp_count++;
      }
      if (tmp_count > 10) {
        break;
      }
    }
  
    return 0;*/
 
#ifdef ENABLE_BLOCK   
    size_t edge_count = 0;
 
    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); vitr != vend; ++vitr) {

      //std::cout << "Vertex value = " << vitr->value() << std::endl;
      auto& v_value = vitr->value();
      //std::cout << "Vertex value = " << v_value << std::endl;
      //auto v_label = v_value.as_object()["type"].as_string();
      auto v_label = v_value.as_object()["labels"].as_array()[0].as_string(); 
      if (v_label == "Food") {
        //std::cout << v_label << std::endl;    
      }
  
      for (auto eitr = graph->edges_begin(vitr->key()), eend = graph->edges_end(vitr->key()); eitr != eend; ++eitr) {
        edge_count++;
      }  

      //if (edge_count > 10) {
      //  break;
      //}    
    }

    std::cout << "#of undirected edges: " << edge_count << std::endl;
   
    return 0; 
#endif

    // graph algorithm

    using Vertex = std::string; //std::string_view; // wrong
    using Edge = std::string; //std::string_view;
    using VertexVertexMap = std::unordered_map<Vertex, Vertex>;
    using VertexStringMap = std::unordered_map<Vertex, std::string>;

    Vertex source_vertex = "2283541"; //"2268378";
    Vertex target_veretx = "2268378"; //"2283541";

    std::string vertex_label = "Compound"; 
    //std::string vertex_label = "Reaction";

    //std::string search_keyword = "poun";     
    std::string search_keyword = "carbon";
    //std::string search_keyword = "WCMBWPHB-BJPFZR";
    //std::string search_keyword = "2283541"; 

    //VertexVertexMap vertex_predecessor_map;   
    VertexStringMap vertex_search_result; 

    std::chrono::time_point<std::chrono::steady_clock> global_start_time =
      std::chrono::steady_clock::now();

    //bfs_level_synchronous_single_source(graph, source_vertex, target_veretx);
    vertex_keyword_search(graph, vertex_label, search_keyword,
      vertex_search_result);

    std::chrono::time_point<std::chrono::steady_clock> global_end_time =
      std::chrono::steady_clock::now();
    double global_elapsed_time = getElapsedTimeSecond(global_start_time,
      global_end_time);
    std::cout << "Time: " << global_elapsed_time << " seconds."
      << std::endl;

    std::cout << "# matches found: " << vertex_search_result.size() << std::endl;

  }

  return 0;
}
