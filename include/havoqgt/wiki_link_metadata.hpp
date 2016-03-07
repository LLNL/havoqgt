#ifndef WIKI_LINK_METADATA_INCLUDED
#define WIKI_LINK_METADATA_INCLUDED

#include <havoqgt/wiki_parallel_edge_list_reader.hpp>
#include <iostream>
#include <chrono>

#define get(y, x) y[index_map[x]]

namespace havoqgt {
  
  using clock_t = std::chrono::high_resolution_clock; 
  using time_point_t  = std::chrono::time_point<clock_t>;
  
  int index_map[9] = {0, 1, 4, 5, 2, 3, 6, 7, 8};

  class wiki_link_metadata {
  public:
    sha1key link_from;
    sha1key link_to;
    time_point_t added_at;
    time_point_t deleted_at;
    sha1key added_by;
    sha1key deleted_by;
    int type;
    bool redirect;
    int ns;
    
    static std::unordered_map<sha1key, uint64_t, sha1hasher> sha1_label_map;

    wiki_link_metadata() : type(0) {} //required due to templatization
    
    void initialize(std::string line) {
      std::stringstream ss(line);
      std::string token;
      std::vector<std::string> tokens;
      while( std::getline( ss, token, ' ') ) {
	tokens.push_back( token );
      }

      link_from  = sha1key( get(tokens, 0) );
      link_to    = sha1key( get(tokens, 1) );
      added_at   = time_point_t(parse_as_seconds( get(tokens, 2) ));
      deleted_at = time_point_t(parse_as_seconds( get(tokens, 3) ));
      added_by   = sha1key( get(tokens, 4) );
      deleted_by = sha1key( get(tokens, 5) );
      type       = -(std::atoi( get(tokens, 6).c_str() ));
      redirect   = std::atoi( get(tokens, 7).c_str() );
      ns         = std::atoi( get(tokens, 8).c_str() );
    }
    
    bool is_recorded() {
      return type > 0;
    }

    void register_recorded() {
      type = -type;
    }

    uint64_t get_src_label() const { return wiki_link_metadata::sha1_label_map.find( link_from )->second; }
    uint64_t get_dest_label() const { return wiki_link_metadata::sha1_label_map.find( link_to )->second; }
    time_point_t start_time() const { return added_at; }
    time_point_t end_time() const  { return deleted_at; }
    /*    template<typename T>
    static void preprocess(T args) {
      std::string str = static_cast<std::string>(args);
      wiki_parallel_edge_list_reader.read_pagekeyset(str, sha1_label_map); 
      }*/

  private:
    std::chrono::seconds parse_as_seconds( std::string& sec) {
      return std::chrono::seconds( static_cast<uint64_t>( std::atol( sec.c_str() ) ) );
   }
  } __attribute__ ((packed)) ;


};

#endif
