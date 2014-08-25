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
#endif

#if DEBUG_INSERTEDEDGES == 1
  #warning DEBUG_INSERTEDEDGES is enabled.
  static const std::string kFnameDebugInsertedEdges = "/l/ssd/graph_out.debug_edges";
  //static const std::string kFnameDebugInsertedEdges = "/usr/localdisk/fusion/graph_out.debug_edges";
#endif


template <typename SegmentManager>
class construct_dynamicgraph {
public:

  template<typename T>
  using SegmentAllocator = bip::allocator<T, SegmentManager>;
  typedef bip::managed_mapped_file mapped_t;

  class bitmap_mgr {
  public:
    bitmap_mgr()
    {
      table_ = nullptr;
      max_id_ = 0;
      capacity_ = 0;
    }

    ~bitmap_mgr() {
      if (table_ != nullptr)
        delete[] table_;
    }

    inline void set(const uint64_t _id)
    {
      #define roundup(x,n) (((x)+((n)-1))&(~((n)-1)))
      if (_id >= max_id_) {
        max_id_ = roundup(_id, 64ULL);
        grow_table();
      }
      table_[_id/64ULL] |= (1ULL << (_id & 63ULL));
    }

    inline void unset(const uint64_t _id)
    {
      table_[_id/64ULL] ^= (1ULL << (_id & 63ULL));
    }

    inline bool get(const uint64_t _id)
    {
      if(static_cast<int64_t>(_id) >= max_id_) return false;
      return table_[_id/64ULL] & (1ULL << (_id & 63ULL));
    }

    inline size_t size() const
    {
      return capacity_;
    }

  private:

    void grow_table() {
      const uint64_t old_capacity = capacity_;
      const uint64_t *old_table = table_;

      capacity_ |= (capacity_ == 0);
      while (capacity_ < max_id_ / 64ULL)
        capacity_ *= 2ULL;
      max_id_ = capacity_*64ULL;

      calloc();

      std::memcpy(table_, old_table, old_capacity * sizeof(uint64_t));
      delete[] old_table;
    }

    inline void calloc() {
        table_ = new uint64_t[capacity_];
        std::memset(table_, 0, capacity_ * sizeof(uint64_t));
    }

    uint64_t *table_;
    uint64_t max_id_;
    size_t capacity_;

  };


  // ---------  Data Structures ------------ //
  typedef bip::vector<uint64_t, SegmentAllocator<uint64_t>> uint64_vector_t;
  typedef bip::vector<uint64_vector_t, SegmentAllocator<uint64_vector_t>> adjacency_matrix_vec_vec_t;
  
  typedef std::pair<const uint64_t, uint64_vector_t> map_value_vec_t;
  typedef boost::unordered_map<uint64_t, uint64_vector_t, boost::hash<uint64_t>, std::equal_to<uint64_t>, SegmentAllocator<map_value_vec_t>> adjacency_matrix_map_vec_t;

  typedef robin_hood_hash<int64_t, int64_t, SegmentManager> robin_hood_hashing_t;

  //typedef std::pair<const uint64_t, robin_hood_hashing_t> map_value_rbh_t;
  typedef robin_hood_hash<int64_t, robin_hood_hashing_t, SegmentManager> adjacency_matrix_rbh_rbh_t;
  //typedef boost::unordered_map<uint64_t, robin_hood_hashing_t, boost::hash<uint64_t>, std::equal_to<uint64_t>, SegmentAllocator<map_value_rbh_t>> adjacency_matrix_map_rbh_t;

  enum DataStructureMode {
    kUseVecVecMatrix,
    kUseMapVecMatrix,
    kUseRobinHoodHash,
    kUseDegreeAwareModel
  };


  ///  ------------------------------------------------------ ///
  ///              Constructor / Deconstructor
  ///  ------------------------------------------------------ ///
  //--  Constructor -- //
  construct_dynamicgraph(
    mapped_t& asdf,
    SegmentAllocator<void>& seg_allocator,
    const DataStructureMode mode,
    const uint64_t n);

  //--  Deconstructors -- //
  ~construct_dynamicgraph();


  ///  ------------------------------------------------------ ///
  ///              Public Member Functions
  ///  ------------------------------------------------------ ///
  /// add edges controller
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
        add_edges_degree_aware_rbh_first(edges);
        //add_edges_degree_aware_bitmap_first(edges);
        //add_edges_degree_aware_adj_first(edges);
        //add_edges_degree_aware_rbh_first_rbh_mtrx(edges);
        break;

      default:
        std::cerr << "Unknown data structure type" << std::endl;
        exit(-1);
    }

  };

  void print_profile();


private:

  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///

  /// add edges vector-vector adjacency-matrix
  template <typename Container>
  void add_edges_adjacency_matrix_vector_vector(Container& edges);

  /// add edges unsorted_map-vector adjacency-matrix
  template <typename Container>
  void add_edges_adjacency_matrix_map_vector(Container& edges);

  /// add edges array by using robin hood hash
  template <typename Container>
  void add_edges_robin_hood_hash(Container& edges);

  /// Add edges to robihn-hood-hash or adjacency-matrix depends on degree of souce vertex
  /// check robihn hood hashing first
  template <typename Container>
  void add_edges_degree_aware_rbh_first(Container& edges);

  /// Add edges to robihn-hood-hash or adjacency-matrix depends on degree of souce vertex
  /// check bitmap first
  template <typename Container>
  void add_edges_degree_aware_bitmap_first(Container& edges);

  /// Add edges to robihn-hood-hash or adjacency-matrix depends on degree of souce vertex
  /// check adjacency-matrxit first
  template <typename Container>
  void add_edges_degree_aware_adj_first(Container& edges);

  template <typename Container>
  void add_edges_degree_aware_rbh_first_rbh_mtrx(Container& edges);

  inline void flush_pagecache() {
    asdf_.flush();
  }


  ///  ------------------------------------------------------ ///
  ///              Private Member Variables
  ///  ------------------------------------------------------ ///
  mapped_t& asdf_;
  SegmentAllocator<void>& seg_allocator_;
  const DataStructureMode data_structure_type_;
  const uint64_t kLowDegreeThreshold;

  adjacency_matrix_vec_vec_t *adjacency_matrix_vec_vec_;
  uint64_vector_t *init_vec;
  adjacency_matrix_map_vec_t *adjacency_matrix_map_vec_;
  robin_hood_hashing_t *robin_hood_hashing_;
  bitmap_mgr *is_exist_bmp_;
  adjacency_matrix_rbh_rbh_t *adjacency_matrix_rbh_rbh_;

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
