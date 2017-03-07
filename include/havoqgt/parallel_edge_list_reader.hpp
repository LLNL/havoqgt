/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef HAVOQGT_PARALLEL_EDGE_LIST_READER_INCLUDED
#define HAVOQGT_PARALLEL_EDGE_LIST_READER_INCLUDED

#include <vector>
#include <fstream>
#include <deque>
#include <string>
#include <sstream>
#include <utility>
#include <stdint.h>
#include <iostream>

#include <havoqgt/environment.hpp>

namespace havoqgt {

std::vector<std::string> split(const std::string& line, const char delim) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream iss(line);

  while(std::getline(iss, token, delim)) {
    tokens.push_back(token);
  }

  return tokens;
}

/// Parallel edge list reader
///
template <typename edge_data_type = uint8_t>
class parallel_edge_list_reader {

public:
  typedef uint64_t                      vertex_descriptor;
  typedef std::tuple<uint64_t, uint64_t, edge_data_type> value_type;
  typedef value_type edge_type;
  typedef edge_data_type edge_data_value_type; 

  ///
  /// InputIterator class for rmat_edge_generator

  class input_iterator_type : public std::iterator<std::input_iterator_tag, edge_type, ptrdiff_t, const edge_type*, const edge_type&> {

  public:
    input_iterator_type(parallel_edge_list_reader* ptr_reader, uint64_t count)
    	: m_ptr_reader(ptr_reader)
    	, m_count(count)
      , m_make_undirected(false) {
      if(m_count == 0) {
        get_next();
        m_count = 0; //reset to zero
      }
    }

    const edge_type& operator*() const { return m_current; }
    //const uint64_t* operator->() const { return &(operator*()); }
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

    ///  Return true if x and y are both end or not end, or x and y are the same.
    friend bool
    operator==(const input_iterator_type& x, const input_iterator_type& y)
    { return x.is_equal(y); }

    ///  Return false if x and y are both end or not end, or x and y are the same.
    friend bool
    operator!=(const input_iterator_type& x, const input_iterator_type& y)
    { return !x.is_equal(y); }

  private:
    input_iterator_type();

    void get_next() {
      if(m_ptr_reader->m_undirected && m_make_undirected) {
        std::swap(std::get<0>(m_current), std::get<1>(m_current));
        m_make_undirected = false;
      } else {
        bool ret = m_ptr_reader->try_read_edge(m_current);
        ++m_count;
        m_make_undirected = true;
      }
      assert(std::get<0>(m_current) <= m_ptr_reader->max_vertex_id());
      assert(std::get<1>(m_current) <= m_ptr_reader->max_vertex_id());           
    }

    parallel_edge_list_reader* m_ptr_reader;
    uint64_t m_count;
    edge_type m_current;
    bool m_make_undirected;
  };


  /// @todo Add undirected flag
  parallel_edge_list_reader(const std::vector< std::string >& filenames, bool undirected)
    : m_undirected(undirected) {
    int shm_rank  = havoqgt_env()->node_local_comm().rank();
    int shm_size  = havoqgt_env()->node_local_comm().size();
    int node_rank = havoqgt_env()->node_offset_comm().rank();
    int node_size = havoqgt_env()->node_offset_comm().size();
    m_local_edge_count = 0;
    m_global_max_vertex = 0;      
    m_has_edge_data = false;
    
    // identify filenames to be read by local rank
    for(size_t i=0; i<filenames.size(); ++i) {
      //if(i % mpi_size == mpi_rank) {
      if(i % node_size == node_rank && (i / node_size) % shm_size == shm_rank) {
        m_local_filenames.push_back(filenames[i]);
      }
    }

    size_t global_num_files = mpi::mpi_all_reduce(m_local_filenames.size(), std::plus<size_t>(), MPI_COMM_WORLD);
    if(havoqgt_env()->world_comm().rank() == 0) {
      std::cout << "Ingesting from " << global_num_files << " files." << std::endl;   

      // rank 0 reads the first input file and determines if edge data exists 
      // and broadcasts to all other ranks
      std::string line;         
      std::ifstream input_file(m_local_filenames[0], std::ifstream::in);            
      if(std::getline(input_file, line)) {
        auto tokens = split(line, ' ');
        if(tokens.size() > 2) {
          m_has_edge_data = true;
        }
      }       
      input_file.close();     
    }

    // broadcast if edge data exists
    mpi::mpi_bcast(m_has_edge_data, 0, MPI_COMM_WORLD);

    // First pass to calc max vertex and count edges.
    open_files();
    //std::cout << "files open" << std::endl;

    edge_type edge;
    uint64_t local_max_vertex = 0;
    while(try_read_edge(edge)) {
      ++m_local_edge_count;
      local_max_vertex = std::max(std::get<0>(edge), local_max_vertex);
      local_max_vertex = std::max(std::get<1>(edge), local_max_vertex);
    }
    m_global_max_vertex = mpi::mpi_all_reduce(local_max_vertex, std::greater<uint64_t>(), MPI_COMM_WORLD);
  }


  /// Returns the begin of the input iterator
  input_iterator_type begin() {
    // Reset and prepare to have readers
    open_files();
    return input_iterator_type(this, 0);
  }

  /// Returns the end of the input iterator
  input_iterator_type end() {
    return input_iterator_type(this, m_local_edge_count);
  }

  /// @todo implement
  uint64_t max_vertex_id() {
    return m_global_max_vertex;
  }

  size_t size() {
    return m_local_edge_count;
  }

  bool has_edge_data() {
    return m_has_edge_data; 
  }

protected:
  
  bool try_read_edge(edge_type& edge) {
    std::string line;
    uint64_t source;
    uint64_t target;
    edge_data_type weight;
    while(!m_ptr_ifstreams.empty()) {
      if(std::getline(*(m_ptr_ifstreams.front()), line)) {
        std::stringstream ssline(line); 
        if (m_has_edge_data) {  
          ssline >> source >> target >> weight;
          //std::cout << source << " " << target << " " << weight << std::endl;
        } else {
          ssline >> source >> target;
          weight = static_cast<edge_data_type>(1); // 1 - true
          //std::cout << source << " " << target << " " << weight << std::endl;
        }
        edge = std::forward_as_tuple(source, target, weight);
        return true;
      } else { //No remaining lines, close file.
        delete m_ptr_ifstreams.front();
        m_ptr_ifstreams.pop_front();
      }
    }
    return false;
  }
 
  void open_files() {
    if(!m_ptr_ifstreams.empty()) {
      HAVOQGT_ERROR_MSG("m_ptr_ifstreams not empty.");
    }
    for(auto itr=m_local_filenames.begin(); itr!=m_local_filenames.end(); ++itr) {
      std::ifstream* ptr = new std::ifstream(*itr);
      if(ptr->good()) {
        m_ptr_ifstreams.push_back(ptr);
      } else {
        std::cerr << "Error opening filename: " << *itr;
      }
    }
  }
  
  std::vector< std::string >     m_local_filenames;
  std::deque< std::ifstream* > m_ptr_ifstreams;
  uint64_t m_local_edge_count;
  uint64_t m_global_max_vertex;
  bool m_undirected;
  bool m_has_edge_data;
};

} //end namespace havoqgt

#endif //HAVOQGT_PARALLEL_EDGE_LIST_READER_INCLUDED
