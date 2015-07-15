#ifndef HAVOQGT_INGEST_FLOW_EDGE_LIST_INCLUDED
#define HAVOQGT_INGEST_FLOW_EDGE_LIST_INCLUDED

#include <iostream>
#include <sstream>
#include <cstdint>
#include <fstream>
#include <vector>
#include <deque>
#include <havoqgt/nano_time.hpp>

#include <havoqgt/environment.hpp>

namespace havoqgt {

struct IPV4{
  u_short vals[4];
} __attribute__ ((packed));

typedef struct IPV4 ip_type;
  
class flow {
protected:
  ip_type src_ip;
  uint32_t src_label;
  u_short src_port;
  uint64_t src_port_label;
  
  ip_type dest_ip;
  uint32_t dest_label;
  u_short dest_port;
  uint64_t dest_port_label;
  
  nano_time start;
  nano_time end;
  uint32_t total_bytes;
  u_short protocol;
  uint32_t packet_counts;

  void tokenize(std::string line, char delim, std::vector<std::string>& tokens){
    std::stringstream ssline(line);
    std::string token;
    int i = 0;
    while(std::getline(ssline, token, delim)){
      tokens.at(i++) = token;
    }
  }

  uint32_t to_integer(std::string token) {
    uint32_t val = -1;
    try{
      val = std::stoi(token);
    } catch (const std::invalid_argument& ia){
      std::cout << "Throwing invalid argument on token '" << token << "'." << std::endl;
    }
    return val;
  }

  uint64_t to_long(std::string token) {
    uint64_t val = -1;
    try{
      val = static_cast<uint64_t>(std::stol(token));
    }catch(const std::invalid_argument& ia) {
      std::cout << "Throwing invalid argument on token'" << token <<"'." << std::endl;
    }
    return val;
  }

  void split_ip_token(std::string ip_token, ip_type& ip) {
    std::stringstream ssline(ip_token);
    std::string token;
    int i = 0;
    while(std::getline(ssline, token, '.')) {
      ip.vals[i++] = static_cast<u_short>(to_integer(token));
    }
  }

public:
  flow() : protocol(0x8000) { }

  void initialize(std::string line) { 
    std::vector<std::string> tokens(13);
    tokenize(line, ';', tokens);
        
    split_ip_token(tokens[0], src_ip);
    src_label = to_integer(tokens[1]);
    src_port = static_cast<u_short>(to_integer(tokens[2]));
    src_port_label = to_long(tokens[3]);
					   
    split_ip_token(tokens[4], dest_ip);
    dest_label = to_integer(tokens[5]);
    dest_port = static_cast<u_short>(to_integer(tokens[6]));
    dest_port_label = to_long(tokens[7]);

    start = nano_time(tokens[8]);
    end = nano_time(tokens[9]);
    total_bytes = to_integer(tokens[10]);
    protocol = static_cast<u_short>(to_integer(tokens[11])) | (0x8000);
    packet_counts = to_integer(tokens[12]);
  }

  bool is_recorded() {
    return ( (protocol & 0x8000) == 0 ) ? true : false;
  }

  void register_recorded() {
    protocol = protocol ^ 0x8000;
  }

  uint32_t get_src_label() const { return src_label; }
  uint32_t get_dest_label() const { return dest_label; }
} __attribute__ ((packed));

class ingest_flow_edge_list {


  class flow_input_iterator :
    public std::iterator< std::input_iterator_tag, flow, ptrdiff_t, const flow*, const flow&> {

  public:
    flow_input_iterator(ingest_flow_edge_list* flow_list_reader, uint64_t count)
      :m_ptr_reader(flow_list_reader)
      ,m_count(count){
      if(m_count == 0) {
	get_next();
	m_count = 0;
      }
    }

    const flow& operator*() const {
      return m_current;
    }

    flow_input_iterator& operator++() {
      get_next();
      return *this;
    }

    flow_input_iterator operator++(int) {
      flow_input_iterator __tmp = *this;
      get_next();
      return __tmp;
    }

    flow *operator->() {
      return &m_current;
    }

    bool is_equal(const flow_input_iterator& _x) const {
      return m_count == (_x.m_count);
    }

    friend bool
    operator ==(const flow_input_iterator& x, const flow_input_iterator& y)
    { return x.is_equal(y); }

    friend bool
    operator !=(const flow_input_iterator& x, const flow_input_iterator& y)
    { return !x.is_equal(y); }

  private:
    flow_input_iterator() { }

    void get_next() {
      bool retr = m_ptr_reader->try_read_flow(m_current);
      ++m_count;
    }
    
    ingest_flow_edge_list* m_ptr_reader;
    uint64_t m_count;
    flow m_current;
  };
  //protected:
  //bool try_read_flow(flow, bool);
  
public:
  ingest_flow_edge_list(const std::vector<std::string>& filenames, uint64_t local_edge_count = -1)
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
      flow _flow;
      while( try_read_flow(_flow, false) ) {
	++m_local_edge_count;
      }	
      std::cout << "Total Local Edge Count " << m_local_edge_count << std::endl;
    }
  }
  
  flow_input_iterator begin() {
    open_files();
    return flow_input_iterator(this, 0);
  }

  flow_input_iterator end() {
    return flow_input_iterator(this, m_local_edge_count);
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
	std::cout << "Opening File for metadata  : " << (*itr) << std::endl;
	m_ptr_ifstreams.push_back(ptr);
	m_filenames.push_back(*itr);
      } else {
	std::cerr << "Error opening filename: " << *itr;
      }
    }
  }
    
  bool try_read_flow(flow &curr_flow, bool parse = true) {
    if(m_ptr_ifstreams.empty()) return false;
    
    std::ifstream *ptr = m_ptr_ifstreams.front();
    std::string line;
    while(std::getline(*ptr, line)) {
      if(parse)
	curr_flow.initialize(line);
      return true;
    }
    std::cout << "Closing reading file " << m_filenames.front() << std::endl;
    ptr->close();
    m_ptr_ifstreams.pop_front();
    m_filenames.pop_front();
    delete(ptr);
    return !m_ptr_ifstreams.empty();
  }
  

  std::deque<std::ifstream*> m_ptr_ifstreams;
  std::deque<std::string> m_filenames;
  std::vector<std::string> m_local_filenames;
  uint64_t m_local_edge_count;
};

}; // end namespace havoqgt

#endif  //  HAVOQGT_INGEST_FLOW_EDGE_LIST_INCLUDED
