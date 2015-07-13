#ifndef HAVOQGT_PARALLEL_FLOW_EDGE_LIST_READER_INCLUDED
#define HAVOQGT_PARALLEL_FLOW_EDGE_LIST_READER_INCLUDED

#include <vector>
#include <fstream>
#include <deque>
#include <string>
#include <sstream>
#include <utility>
#include <stdint.h>

#include <havoqgt/environment.hpp>

namespace havoqgt {

class parallel_flow_edge_list_reader {

public:
  typedef uint64_t                       vertex_descriptor;
  typedef std::pair<uint64_t, uint64_t> edge_type;

  class input_iterator_type : public std::iterator<std::input_iterator_tag, edge_type, ptrdiff_t, const edge_type*, const edge_type&> {
  public:
    input_iterator_type(parallel_flow_edge_list_reader* ptr_reader, uint64_t count)
      : m_ptr_reader(ptr_reader)
      , m_count(count) {
      if(m_count == 0) {
	get_next();
	m_count = 0; //reset to zero;
      }
    }

    const edge_type& operator*() const { return m_current; }
    
    input_iterator_type& operator++() {
      get_next();
      return *this;
    }

    input_iterator_type operator++(int) {
      input_iterator_type __tmp = *this;
      get_next();
      return __tmp;
    }

    edge_type *operator->() {
      return &m_current;
    }

    bool is_equal(const input_iterator_type& _x) const {
      return m_count == (_x.m_count);
    }

    /// Return true if x and y are both end or not end, or x and y are the same.
    friend bool
    operator==(const input_iterator_type& x, const input_iterator_type& y) {
      return x.is_equal(y);
    }

    /// Return false if x and y are both end or not end, or x and y are the same.
    friend bool
    operator!=(const input_iterator_type& x, const input_iterator_type& y) {
      return !x.is_equal(y);
    }

  private:
    input_iterator_type();
    
    void get_next() {
      bool ret = m_ptr_reader->try_read_edge(m_current);
      ++m_count;
      assert( m_current.first <= m_ptr_reader->max_vertex_id());
      assert( m_current.second <= m_ptr_reader->max_vertex_id());
    }

    parallel_flow_edge_list_reader* m_ptr_reader;
    uint64_t m_count;
    edge_type m_current;
  };

  parallel_flow_edge_list_reader(const std::vector<std::string>& filenames, bool _undirected) : undirected(_undirected) {
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();

    m_local_edge_count = 0;
    m_global_max_vertex = 0;

    for(size_t i = 0; i < filenames.size(); ++i){
      if(i % mpi_size == mpi_rank) {
	m_local_filenames.push_back(filenames[i]);
      }
    }

    // First pass to calc max vertex and count edges
    open_files();
    std::cout << "files open" << std::endl;
    edge_type edge;
    uint64_t local_max_vertex = 0;
    while(try_read_edge(edge)) {
      ++m_local_edge_count;
      local_max_vertex = std::max( edge.first, local_max_vertex);
      local_max_vertex = std::max(edge.second, local_max_vertex);
    }
    m_global_max_vertex = mpi::mpi_all_reduce(local_max_vertex, std::greater<uint64_t>(), 
					      havoqgt_env()->world_comm().comm());
  }

  uint64_t max_vertex_id() {
    return m_global_max_vertex;
  }

  size_t size() {
    return m_local_edge_count;
  }

  input_iterator_type begin() {
    open_files();
    return input_iterator_type(this, 0);
  }

  input_iterator_type end() {
    return input_iterator_type(this, m_local_edge_count);
  }

protected:
  bool try_read_edge(edge_type& edge) {
    //static int count = 0;
    std::string line;
    while(!m_ptr_ifstreams.empty()) {
      if(std::getline(*(m_ptr_ifstreams.front()), line )) {
        std::stringstream ssline(line);
	std::string token;
	int i = 0;
	while(std::getline(ssline, token, ';')){
	  if(i == 1) edge.first = std::stoi(token);
	  else if(i == 5) {
	    edge.second = std::stoi(token);
	    break; // break here
	  }
	  i++;
	}
	/*if(count < 100){
	  count++;
	  std::cout << edge.first <<", "<<edge.second<< std::endl;
	  }*/
	  return true; 
      } else {
	delete m_ptr_ifstreams.front();
	m_ptr_ifstreams.pop_front();
      }
    }
    return false;
  }

  void open_files() {
    if(!m_ptr_ifstreams.empty()) {
      HAVOQGT_ERROR_MSG("m_ptr_ifstreams not empty. ");
    }
    for( auto itr = m_local_filenames.begin(); itr!=m_local_filenames.end(); ++itr){
      std::ifstream* ptr = new std::ifstream(*itr);
      if(ptr->good()) {
	m_ptr_ifstreams.push_back(ptr);
      } else {
	std::cerr << "Erro opening filename: " << *itr;
      }
    }
  }

  std::vector<std::string> m_local_filenames;
  std::deque<std::ifstream*> m_ptr_ifstreams;
  uint64_t m_local_edge_count;
  uint64_t m_global_max_vertex;
  bool undirected;
};

} // end namespace havoqgt 

#endif //  HAVOQGT_PARALLEL_FLOW_EDGE_LIST_READER_INCLUDED
