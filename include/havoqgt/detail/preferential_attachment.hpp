// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <utility>
#include <boost/random.hpp>



namespace havoqgt { namespace detail {

template <typename DUMMY=uint64_t>
class preferential_attachment_helper {
public:
  typedef uint64_t                  vertex_descriptor;
  typedef typename std::pair<uint64_t,uint64_t> edge_type;

  preferential_attachment_helper(uint64_t k, uint64_t m, double beta, uint64_t rng_seed=5489)
    : m_rng(rng_seed) 
  {
    m_k = k;
    m_num_edges = m;
    m_koffset = m_k*(m_k+1)/2;
    m_ptr_mask =  ~(std::numeric_limits<vertex_descriptor>::max() >> 1);
    m_alpha = (beta / double(m_k) + double(1)) / (beta / double(m_k) + 2);
    //std::cout << "beta = " << beta << ", alpha = " << m_alpha << std::endl;
  }

  edge_type gen_edge(uint64_t _edge_index) {
    edge_type to_return;
    to_return.first = calc_source(_edge_index);
    if(_edge_index >= m_koffset)  {
      //
      // Generate random edge_list location based on beta model
      boost::random::uniform_01<boost::random::mt19937> rand_prob(m_rng);
      boost::random::uniform_int_distribution<uint64_t> uid(0,_edge_index-1);
      uint64_t rand = uid(m_rng) * 2;
      if(rand_prob() > m_alpha) {
        ++rand;
      } 
      if(rand % 2 == 0) { //this is a source vertex, we can calc!
        to_return.second = calc_source(rand/2);
      } else {
        uint64_t edge_rand = rand/2;
        if(edge_rand < m_koffset) {  //this is an early edge we can calc!
          to_return.second = calc_target(edge_rand);
        } else {
          to_return.second = make_pointer(edge_rand);
        }
      }
    } else {
      to_return.second = calc_target(_edge_index);
    }
    return to_return;
  }

  uint64_t calc_source(uint64_t i) {
    uint64_t to_return;
    if(i+1>m_koffset) {
      to_return = ((i-m_koffset)/m_k)+m_k+1;
    } else {
      to_return = (uint64_t) floor(double(-0.5f) + 
                  sqrt(double(0.25f)+double(2)*double(i))+1);
    }
    return to_return;
  }


  uint64_t calc_target(uint64_t i) {
    uint64_t to_return;
    if(i+1>m_koffset) {
      assert(false);
      to_return = i;
    } else {
      double tmp = double(-0.5f) + sqrt(double(0.25f)+double(2)*double(i))+1;
      to_return = (uint64_t) ((tmp - floor(tmp)) * floor(tmp));
    }
    return to_return;
  }

  vertex_descriptor make_pointer(vertex_descriptor i) {
    return i | m_ptr_mask;
  }

  bool is_pointer(vertex_descriptor i) {
    return (i & m_ptr_mask) > 0;
  }

  //dereferences point
  uint64_t value_of_pointer(vertex_descriptor i, uint64_t _num_partitions=1) {
    i = i & ~m_ptr_mask;
    uint64_t my_partition = i%_num_partitions;
    uint64_t my_partition_offset = i/_num_partitions;
    uint64_t edges_per_partition = m_num_edges / _num_partitions;
    return my_partition * edges_per_partition + my_partition_offset;
  }


private:
  boost::random::mt19937 m_rng;
  uint64_t               m_k;
  uint64_t               m_koffset;
  vertex_descriptor      m_ptr_mask;
  uint64_t               m_num_edges;
  double                 m_alpha;
};


}} //end namespace havoqgt::detail
