#ifndef WIKI_LINK_METADATA_INCLUDED
#define WIKI_LINK_METADATA_INCLUDED

#include <havoqgt/wiki_parallel_edge_list_reader.hpp>
#include <havoqgt/detail/hash.hpp>
#include <iostream>
#include <chrono>

#define get_index(y, x) y[index_map[x]]

namespace havoqgt {
  
  using clock_t = std::chrono::high_resolution_clock; 
  using time_point_t  = std::chrono::time_point<clock_t>;
  
  int index_map[9] = {0, 1, 4, 5, 2, 3, 6, 7, 8};

  class wiki_link_metadata {
  public:
    uint64_t link_from;
    uint64_t link_to;

    uint64_t added_at;
    uint64_t deleted_at;

    sha1key added_by;
    sha1key deleted_by;

    int type;
    bool redirect;
    int ns;
    
    static std::unordered_map<sha1key, uint64_t, sha1hasher> sha1_label_map;

    wiki_link_metadata() : type(0) {} //required due to templatization
    
    void initialize(std::string line) {
      try {
      std::stringstream ss(line);
      std::string token;
      std::vector<std::string> tokens;
      while( std::getline( ss, token, ' ') ) {
	tokens.push_back( token );
      }
      if( tokens.size() != 9 ) {
	std::cout << "Bad line in the file. Line: " << line << std::endl;
	return;
      }
      link_from  = std::atol( get_index(tokens, 0).c_str()) - 1;
      link_to    = std::atol( get_index(tokens, 1).c_str()) - 1;

      added_at   = std::atol( get_index(tokens, 2).c_str());
      deleted_at = std::atol( get_index(tokens, 3).c_str());
      if(deleted_at == 0) deleted_at = std::numeric_limits<uint64_t>::max();
      added_by   = sha1key( get_index(tokens, 4) );
      deleted_by = sha1key( get_index(tokens, 5) );

      type       = -(std::atoi( get_index(tokens, 6).c_str() ));
      redirect   = std::atoi( get_index(tokens, 7).c_str() );
      ns         = std::atoi( get_index(tokens, 8).c_str() );

      } catch( std::exception& e) {
	std::cout << "exception caught in wiki_link_metadata : " << e.what() << std::endl;
      }
    }
    
    bool is_recorded() {
      return type > 0;
    }

    void register_recorded() {
      type = -type;
    }

    uint64_t get_src_label() const { return link_from; }
    uint64_t get_dest_label() const { return link_to; }

    uint64_t start_time() const { return added_at; }
    uint64_t end_time() const { return deleted_at; }

  private:
    std::chrono::seconds parse_as_seconds( std::string& sec) {
      return std::chrono::seconds( static_cast<uint64_t>( std::atol( sec.c_str() ) ) );
   }
  }__attribute__ ((packed)) ;

};

#endif
