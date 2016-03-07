#ifndef INGEST_METADATA_LIST_INCLUDED
#define INGEST_METADATA_LIST_INCLUDED

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>

namespace havoqgt {

template<typename metadata_t>
class ingest_metadata_list {

  class metadata_input_iterator :
    public std::iterator< std::input_iterator_tag, metadata_t, ptrdiff_t, const metadata_t*, const metadata_t&> {

  public:
    metadata_input_iterator(ingest_metadata_list* metadata_list_reader, uint64_t count)
      :m_ptr_reader(metadata_list_reader)
      ,m_count(count){
      if(m_count == 0) {
	get_next();
	m_count = 0;
      }
    }

    const metadata_t& operator*() const {
      return m_current;
    }

    metadata_input_iterator& operator++() {
      get_next();
      return *this;
    }

    metadata_input_iterator operator++(int) {
      metadata_input_iterator __tmp = *this;
      get_next();
      return __tmp;
    }

    metadata_t *operator->() {
      return &m_current;
    }

    bool is_equal(const metadata_input_iterator& _x) const {
      return m_count == (_x.m_count);
    }

    friend bool
    operator ==(const metadata_input_iterator& x, const metadata_input_iterator& y)
    { return x.is_equal(y); }

    friend bool
    operator !=(const metadata_input_iterator& x, const metadata_input_iterator& y)
    { return !x.is_equal(y); }

  private:
    metadata_input_iterator() { }

    void get_next() {
      bool retr = m_ptr_reader->try_read_metadata(m_current);
      ++m_count;
    }
    
    ingest_metadata_list* m_ptr_reader;
    uint64_t m_count;
    metadata_t m_current;
  };
  //protected:
  //bool try_read_flow(flow, bool);
  
public:
  typedef ingest_metadata_list<metadata_t>::metadata_input_iterator metadata_iterator_type;

  ingest_metadata_list(const std::vector<std::string>& filenames, uint64_t local_edge_count = -1)
    : m_local_edge_count(local_edge_count) {  // surajpoudel: -1 to signify values not passed
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();
    
    for(size_t i = 0; i < filenames.size(); i++){
      if( i % mpi_size == mpi_rank){
	m_local_filenames.push_back(filenames[i]);
      }
    }
    
    if(local_edge_count == -1) { // dont have knowledge of the local_edge_count
      open_files();
      m_local_edge_count = 0;
      metadata_t _metadata;
      while( try_read_metadata(_metadata, false) ) {
	++m_local_edge_count;
      }	
      std::cout << "Total Local Edge Count " << m_local_edge_count << std::endl;
    }
  }
  
  metadata_input_iterator begin() {
    open_files();
    return metadata_input_iterator(this, 0);
  }

  metadata_input_iterator end() {
    return metadata_input_iterator(this, m_local_edge_count);
  }

  size_t size() {
    return m_local_edge_count;
  }

protected:
  void open_files() {
    if(!m_ptr_ifstreams.empty()) {
      std::cout << "m_ptr_ifstreams not empty" << std::endl;
      HAVOQGT_ERROR_MSG("m_ptr_ifstreams not empty");
    }
    for( auto itr=m_local_filenames.begin(); itr != m_local_filenames.end(); ++itr){
      std::ifstream* ptr = new std::ifstream(*itr);
      if(ptr->good()) {
	//std::cout << "Opening File for metadata  : " << (*itr) << std::endl;
	m_ptr_ifstreams.push_back(ptr);
	m_filenames.push_back(*itr);
      } else {
	std::cerr << "Error opening filename: " << *itr;
      }
    }
  }
    
  bool try_read_metadata(metadata_t &curr_metadata, bool parse = true) {
    if(m_ptr_ifstreams.empty()) return false;
    
    std::ifstream *ptr = m_ptr_ifstreams.front();
    std::string line;
    while(std::getline(*ptr, line)) {
      if(parse)
	curr_metadata.initialize(line);
      return true;
    }
    if(parse)
      std::cout << "Closing reading file " << m_filenames.front() << std::endl;
    ptr->close();
    m_ptr_ifstreams.pop_front();
    m_filenames.pop_front();
    delete(ptr);
    return try_read_metadata(curr_metadata, parse);
  }
  

  std::deque<std::ifstream*> m_ptr_ifstreams;
  std::deque<std::string> m_filenames;
  std::vector<std::string> m_local_filenames;
  uint64_t m_local_edge_count;
};

}; // end namespace havoqgt

#endif
