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

#ifndef HAVOQGT_MPI_CONSTRUCT_DYNAMICGRAPH_VEC_HPP_INCLUDED
#define HAVOQGT_MPI_CONSTRUCT_DYNAMICGRAPH_VEC_HPP_INCLUDED

#include <havoqgt/mpi.hpp>
//#include <havoqgt/distributed_edge_list.hpp>
#include <havoqgt/detail/iterator.hpp>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/map.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/offset_ptr.hpp>
#include <boost/interprocess/containers/set.hpp>

#include <stdint.h>
#include <utility>
#include <limits>

#include <havoqgt/robin_hood_hashing.hpp>

 namespace havoqgt {
  namespace mpi {


class IOInfo {
  public:
    IOInfo();
    void init();
    void reset_baseline();
    void get_status(int &r, int &w);
    void log_diff(bool final);

  private:
    int read_previous_mb_;
    int written_previous_mb_;
    int read_total_mb_;
    int written_total_mb_;
};


namespace bip = boost::interprocess;


#ifndef WITHOUT_DUPLICATE_INSERTION
  #define WITHOUT_DUPLICATE_INSERTION 1
#endif

#ifndef DEBUG_INSERTEDEDGES
  #define DEBUG_INSERTEDEDGES 0
  #warning DEBUG_INSERTEDEDGES is enabled.
#endif
#if DEBUG_INSERTEDEDGES == 1
  static const std::string kFnameDebugInsertedEdges = "/usr/localdisk/fusion/graph_out.debug_edges";
#endif

template <typename SegmentManager>
class construct_dynamicgraph {
public:

  template<typename T>
  using SegmentAllocator = bip::allocator<T, SegmentManager>;
  typedef bip::managed_mapped_file mapped_t;

  // class DegreeTable {
  //   public:
  //     DegreeTable() {
  //       max_vertex_id_ = -1;
  //       // num_verticies_ = 0;
  //       // num_edges_ = 0;
  //       capacity_ = 0;
  //     };

  //     ~DegreeTable() {
  //       if (capacity_ > 0)
  //         delete[] degree_table_;
  //     };

  //     inline uint64_t count_up(uint64_t vertex_id) {
  //       if (vertex_id > max_vertex_id_) {
  //         max_vertex_id_ = vertex_id;
  //         grow();
  //       }

  //       return ++degree_table_[vertex_id]; 
  //     };

  //     inline uint64_t count_down(uint64_t vertex_id) {
  //       if (vertex_id > max_vertex_id_) {
  //         return 0;
  //       }
  //       return --degree_table_[vertex_id];
  //     };

  //     inline uint64_t degree(uint64_t vertex_id){
  //       if (vertex_id > max_vertex_id_) {
  //         return 0;
  //       } else {
  //         return degree_table_[vertex_id];
  //       }
  //     };

  //   private:
  //     inline void grow() {
  //       uint64_t old_capacity = capacity_;
  //       uint64_t *old_table = degree_table_;

  //       capacity_ |= (capacity_ == 0);
  //       while (capacity_ < max_vertex_id_+1)
  //         capacity_ *= 2ULL;

  //       calloc();

  //       std::memcpy(degree_table_, old_table, old_capacity * sizeof(int64_t));
  //       delete[] old_table;
  //     }

  //     inline void calloc() {
  //       degree_table_ = new uint64_t[capacity_];
  //       std::memset(degree_table_, 0, capacity_ * sizeof(int64_t));
  //     }

  //     uint64_t *degree_table_;
  //     int64_t max_vertex_id_;
  //     // uint64_t num_verticies_;
  //     // uint64_t num_edges_;
  //     size_t capacity_;
  // };

  // ---------  Data Structures ------------ //
  typedef bip::vector<uint64_t, SegmentAllocator<uint64_t>> uint64_vector_t;
  typedef bip::vector<uint64_vector_t, SegmentAllocator<uint64_vector_t>> adjacency_matrix_vec_vec_t;
  
  typedef std::pair<const uint64_t, uint64_vector_t> map_value_vec_t;
  typedef boost::unordered_map<uint64_t, uint64_vector_t, boost::hash<uint64_t>, std::equal_to<uint64_t>, SegmentAllocator<map_value_vec_t>> adjacency_matrix_map_vec_t;

  typedef robin_hood_hash<int64_t, int64_t, SegmentManager> robin_hood_hashing_t;

  typedef boost::unordered_map<uint64_t, uint64_t> degree_map_t;

  // typedef std::pair<uint64_t, uint64_t> uint64_pair_t;
  // typedef bip::vector<uint64_pair_t, SegmentAllocator<uint64_pair_t>> vec_pair_t;
  // typedef bip::set<uint64_pair_t, SegmentAllocator<uint64_pair_t>> set_pair_t;


  enum DataStructureMode {
    kUseVecVecMatrix,
    kUseMapVecMatrix,
    kUseRobinHoodHash,
    kUseDegreeAwareModel
  };

  //--  Constructors -- //
  construct_dynamicgraph(mapped_t& asdf, SegmentAllocator<void>& seg_allocator, const DataStructureMode mode);

  //--  Deconstructors -- //
  ~construct_dynamicgraph();


  /// add edges
  template <typename Container>
  inline void add_edges_adjacency_matrix(Container& edges)
  {
    switch(data_structure_type_) {
      case kUseVecVecMatrix:
        add_edges_adjacency_matrix_vector_vector(edges);
        break;

      case kUseMapVecMatrix:
        add_edges_adjacency_matrix_map_vector(edges);
        break;

      case kUseRobinHoodHash:
        add_edges_robin_hood_hash(edges);
        break;

      case kUseDegreeAwareModel:
        add_edges_hybrid(edges);
        break;

      default:
        std::cerr << "Unknown data structure type" << std::endl;
        assert(false);
        exit(-1);
    }

  };

  void print_profile();


private:

  /// add edges vector-vector adjacency-matrix
  template <typename Container>
  void add_edges_adjacency_matrix_vector_vector(Container& edges);

  /// add edges unsorted_map-vector adjacency-matrix
  template <typename Container>
  void add_edges_adjacency_matrix_map_vector(Container& edges);

  /// add edges array by using robin hood hash
  template <typename Container>
  void add_edges_robin_hood_hash(Container& edges);

  /// --- TODO: This is a temporarily code ---
  template <typename EdgeType>
  void add_edges_adjacency_matrix_map_vector_core(const EdgeType& edge)
  {
#if DEBUG_INSERTEDEDGES == 1
    fout_debug_insertededges_ << edge.first << "\t" << edge.second << std::endl;
#endif
    auto value = adjacency_matrix_map_vec_->find(edge.first);
    if (value == adjacency_matrix_map_vec_->end()) { // new vertex
      uint64_vector_t vec(1, edge.second, seg_allocator_);
      adjacency_matrix_map_vec_->insert(map_value_vec_t(edge.first, vec));
    } else {
      uint64_vector_t& adjacency_list_vec = value->second;
#if WITHOUT_DUPLICATE_INSERTION == 1
      if (boost::find<uint64_vector_t>(adjacency_list_vec, edge.second) != adjacency_list_vec.end() )
        return;
#endif
      adjacency_list_vec.push_back(edge.second);
    }
  }
  /// --- TODO: This is a temporarily code ---
  template <typename EdgeType>
  inline void add_edges_robin_hood_hash_core(const EdgeType& edge)
  {
#if DEBUG_INSERTEDEDGES == 1
    fout_debug_insertededges_ << edge.first << "\t" << edge.second << std::endl;
#endif
#if WITHOUT_DUPLICATE_INSERTION == 1
    robin_hood_hashing_->insert_unique(edge.first, edge.second);
#else
    robin_hood_hashing_->insert(edge.first, edge.second);
#endif
  }

  template <typename Container>
  void add_edges_hybrid(Container& edges);

  inline void flush_pagecache() {
    asdf_.flush();
  }

  mapped_t& asdf_;
  const SegmentAllocator<void>& seg_allocator_;
  const DataStructureMode data_structure_type_;
  degree_map_t degree_map;

  adjacency_matrix_vec_vec_t *adjacency_matrix_vec_vec_;
  uint64_vector_t *init_vec;
  adjacency_matrix_map_vec_t *adjacency_matrix_map_vec_;
  robin_hood_hashing_t *robin_hood_hashing_;

  IOInfo *io_info_;
  double total_exectution_time_;

#if DEBUG_INSERTEDEDGES == 1
  std::ofstream fout_debug_insertededges_;
#endif

};



/// Frees the container of edges
template <typename Container>
void free_edge_container(Container &edges) {};


template<>
void free_edge_container<std::vector<std::pair<uint64_t, uint64_t> > >(std::vector<std::pair<uint64_t, uint64_t> > &edges){
  std::vector< std::pair<uint64_t, uint64_t> >empty(0);
  edges.swap(empty);
};


} // namespace mpi
} // namespace havoqgt

#include <havoqgt/impl/construct_dynamicgraph.ipp>

#endif //HAVOQGT_MPI_CONSTRUCT_DYNAMICGRAPH_VEC_HPP_INCLUDED
