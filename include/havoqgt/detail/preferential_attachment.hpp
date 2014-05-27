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
