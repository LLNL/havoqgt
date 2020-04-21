#pragma once

#include <filesystem> // C++17, -lstdc++fs is required
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>

namespace stdfs = std::filesystem;

namespace prunejuice { namespace pattern {

template <typename Graph>
class file_utilities {

public:
  
  // directories
  static inline constexpr std::string_view PATTERN_DIR = "/pattern"; // C++17

  static inline constexpr std::string_view OUTPUT_DIR = "/output";
  static inline constexpr std::string_view OUTPUT_EDGES_DIR = "/all_ranks_active_edges"; 
  static inline constexpr std::string_view OUTPUT_EDGES_COUNT_DIR = "/all_ranks_active_edges_count";
  static inline constexpr std::string_view OUTPUT_VERTICES_DIR = "/all_ranks_active_vertices";
  static inline constexpr std::string_view OUTPUT_VERTICES_COUNT_DIR = "/all_ranks_active_vertices_count";  
  static inline constexpr std::string_view OUTPUT_SUBGRAPHS_DIR = "/all_ranks_subgraphs"; 
  static inline constexpr std::string_view OUTPUT_MSG_COUNT_DIR = "/all_ranks_messages";  
  static inline constexpr std::string_view OUTPUT_VERTEX_DATA_DIR = "/all_ranks_vertex_data";

  static inline constexpr std::string_view PRUNED_GRAPH_DIR = "/pruned_graph";

  static inline constexpr std::string_view PRUNED_GRAPH_UNION_DIR_NAME = "pruned_graph_union";
  static inline constexpr std::string_view EDGELIST_DIR = "/edgelist";
  static inline constexpr std::string_view VERTEX_DATA_DIR = "/vertex_data";   
  static inline constexpr std::string_view GRAPH_DIR = "/graph";

  // files
  static inline constexpr std::string_view PATTERN_EDGE_FILE = "/pattern_edge";
  static inline constexpr std::string_view PATTERN_EDGE_DATA_FILE = "/pattern_edge_data";
  static inline constexpr std::string_view PATTERN_VERTEX_DATA_FILE = "/pattern_vertex_data";
  static inline constexpr std::string_view PATTERN_STAT_FILE = "/pattern_stat";
  static inline constexpr std::string_view PATTERN_LC_FILE = "/pattern_local_constraint";
  static inline constexpr std::string_view PATTERN_NLC_FILE = "/pattern_nonlocal_constraint"; 

  static void create_template_directories(std::string root_dir, std::string base_dir) {
    //std::cout << "Creating Template Directories ..." << std::endl;
   
    base_dir = "/" + base_dir; 

    std::string pattern_dir(PATTERN_DIR);

    std::string output_dir(OUTPUT_DIR);
    std::string edges_dir(OUTPUT_EDGES_DIR); 
    std::string edges_count_dir(OUTPUT_EDGES_COUNT_DIR);
    std::string vertices_dir(OUTPUT_VERTICES_DIR);
    std::string vertices_count_dir(OUTPUT_VERTICES_COUNT_DIR);
    std::string subgraphs_dir(OUTPUT_SUBGRAPHS_DIR);
    std::string msg_count_dir(OUTPUT_MSG_COUNT_DIR);
    std::string vertex_data_dir(OUTPUT_VERTEX_DATA_DIR);

    std::string pruned_graph_dir(PRUNED_GRAPH_DIR);      
    
    stdfs::create_directory(root_dir + base_dir);
    
    stdfs::create_directory(root_dir + base_dir + pattern_dir);

    stdfs::create_directory(root_dir + base_dir + output_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + edges_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + edges_count_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + vertices_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + vertices_count_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + subgraphs_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + msg_count_dir);
    stdfs::create_directory(root_dir + base_dir + output_dir + vertex_data_dir);

    stdfs::create_directory(root_dir + base_dir + pruned_graph_dir);

    // TODO: overload '+' for std::string_view

    /*stdfs::create_directory(root_dir + base_dir + PATTERN_DIR); 

    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + EDGES_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + EDGES_COUNT_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + VERTICES_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + VERTICES_COUNT_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + SUBGRAPHS_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + MSG_COUNT_DIR);
    stdfs::create_directory(root_dir + base_dir + OUTPUT_DIR + VERTEX_DATA_DIR); 

    stdfs::create_directory(root_dir + base_dir + PRUNED_GRAPH_DIR);*/

    //std::cout << "Done creating directories." << std::endl;
  }

  static void create_graph_directory(std::string root_dir, std::string base_dir) {
    //std::cout << "Creating Graph Directories ..." << std::endl;

    base_dir = "/" + base_dir;
    std::string edgelist_dir(EDGELIST_DIR);
    std::string vertex_data_dir(VERTEX_DATA_DIR);
    std::string graph_dir(GRAPH_DIR);

    stdfs::create_directory(root_dir + base_dir);
    stdfs::create_directory(root_dir + base_dir + edgelist_dir);
    stdfs::create_directory(root_dir + base_dir + vertex_data_dir);
    stdfs::create_directory(root_dir + base_dir + graph_dir);

    //std::cout << "Done creating directories." << std::endl;  
  }

  /*template <typename TemplateGraph>
  static void write_template_to_file(TemplateGraph& template_graph, 
    std::string root_dir, std::string base_dir) {

    base_dir = "/" + base_dir;
    std::string pattern_dir(PATTERN_DIR);        

    std::string pattern_filepath = root_dir + base_dir + pattern_dir;

    std::string pattern_edge_filename(PATTERN_EDGE_FILE);
    std::string pattern_edge_data_filename(PATTERN_EDGE_DATA_FILE);
    std::string pattern_vertex_data_filename(PATTERN_VERTEX_DATA_FILE);   
    std::string pattern_stat_filename(PATTERN_STAT_FILE);
    std::string pattern_nlc_filename(PATTERN_LC_FILE);
    std::string pattern_lc_filename(PATTERN_NLC_FILE);

    // write edge
    // write edge_data
    // write vertex_data 
    // write stat
    // write local_constraint
    // write nonlocal_constraint
  }*/

  file_utilities(std::string root_dir, std::string base_dir, 
    std::string file_dir, std::string filename) {

    base_dir = "/" + base_dir;
    output_filepath = root_dir + base_dir + file_dir + filename;  
    output_file = std::ofstream(output_filepath, std::ofstream::out);  
     
    //std::cout << output_filepath << std::endl; 
  }   

  file_utilities() {}

  ~file_utilities() {
    output_file.close();
  } 

  std::string output_filepath;
  std::ofstream output_file; 

}; 

}} // end namespace prunejuice::pattern  
