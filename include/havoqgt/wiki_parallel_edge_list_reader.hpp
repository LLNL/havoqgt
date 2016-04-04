#ifndef WIKI_PARALLEL_EDGE_LIST_READER_INCLUDED
#define WIKI_PARALLEL_EDGE_LIST_READER_INCLUDED

#include <havoqgt/parallel_edge_list_reader.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/hex.hpp>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>



namespace havoqgt {

class sha1key {
public:
  sha1key() { 
    for(int i = 0; i < 5; i++) keys[i] = 0; 
  }
  
  explicit sha1key( std::string key) {
    boost::algorithm::unhex( key.begin(), key.end(), keys );
  }

  std::size_t hash() const {
    std::size_t seed = 0;
    for(size_t i = 0 ; i < 5; i++ )
      boost::hash_combine( seed, keys[i]);
    return seed;
  }

  int compare( const sha1key& key) const {
    for(size_t i = 0; i < 5; i++) {
      if( keys[i] != key.keys[i] ) return keys[i] > key.keys[i] ? 1 : -1;
    }
    return 0;
  }

  bool operator==(const sha1key& other) const {
    return ( compare(other) == 0 ? true : false );
  }

private:
  uint32_t keys[5]; //magic number 5
} __attribute__((packed));

  struct sha1hasher{
    typedef sha1key argument_type;
    typedef std::size_t result_type;

    result_type operator()( argument_type const& _sha1key) const {
      return _sha1key.hash();
    }
  };


class wiki_parallel_edge_list_reader: public parallel_edge_list_reader {

public:
  wiki_parallel_edge_list_reader(const std::vector<std::string>& filenames, bool undirected
				 , std::unordered_map<sha1key, uint64_t, sha1hasher>& _pagekeymap):
    parallel_edge_list_reader(), pagekeymap(_pagekeymap){
    initialize(filenames, false);
  }

  static void read_pagekeymap( std::string filename, std::unordered_map<sha1key, uint64_t, sha1hasher>& pagekeymap) {
    std::ifstream pagekeylist(filename, std::ifstream::in);
    std::string line;
    uint64_t id = 0; // Should start from??
    while(std::getline(pagekeylist, line)) {
      try{
	pagekeymap.insert( std::make_pair( sha1key(line), id++) );
      }catch( std::exception& e) {
	std::cout<< filename << " : " << line << std::endl;
	std::cout<< e.what() << std::endl;
	exit(1);
      }
    }
    pagekeylist.close();
  }

protected:
  void populate_edge_data(std::string& line, edge_type& edge) {
    int i = 0;
    std::string token;
    std::stringstream line_ss(line);
    while( std::getline( line_ss, token, ' ') ) {
      if( i == 0) {
	edge.first = std::atol(token.c_str()) - 1;
      } else  if(i == 1) {
	edge.second = std::atol(token.c_str()) - 1;
        break;
      }
      i++;
    }
  }

private:
  std::unordered_map<sha1key, uint64_t, sha1hasher>& pagekeymap;
};


}




#endif
