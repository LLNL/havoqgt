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

#include <havoqgt/RHHAdjacencyMatrix.hpp>

#include <havoqgt/RHH/RHHAllocHolder.hpp>
#include <havoqgt/RHH/RHHMain.hpp>

namespace havoqgt {
namespace mpi {
    namespace bip = boost::interprocess;

#ifndef WITHOUT_DUPLICATE_INSERTION
  #define WITHOUT_DUPLICATE_INSERTION 1
#endif

#ifndef DEBUG_DUMPUPDATEREQUESTANDRESULTS
  #define DEBUG_DUMPUPDATEREQUESTANDRESULTS 0
#endif
#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
  #warning DEBUG_DUMPUPDATEREQUESTANDRESULTS is enabled.
    static const std::string kFnameDebugInsertedEdges = "/l/ssd/graph_out.debug_edges";
  //static const std::string kFnameDebugInsertedEdges = "/usr/localdisk/fusion/graph_out.debug_edges";
#endif

#ifndef DEBUG_DETAILPROFILE
  #define DEBUG_DETAILPROFILE 0
#endif

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


struct EdgeUpdateRequest
{
  EdgeUpdateRequest(){}

  EdgeUpdateRequest(std::pair<uint64_t, uint64_t> _edge, bool _is_delete)
  {
    edge = _edge;
    is_delete = _is_delete;
  }

  std::pair<uint64_t, uint64_t> edge;
  bool is_delete;
};

bool edgerequest_asc( const EdgeUpdateRequest& left, const EdgeUpdateRequest& right ) {
  return left.edge.first < right.edge.first;
}

template <typename SegmentManager>
    class construct_dynamicgraph {
    public:

  // ---------  Typedefs and Enums ------------ //

  template<typename T>
      using SegmentAllocator = bip::allocator<T, SegmentManager>;
      typedef bip::managed_mapped_file mapped_t;

      typedef bip::vector<uint64_t, SegmentAllocator<uint64_t>> uint64_vector_t;
      typedef bip::vector<uint64_vector_t, SegmentAllocator<uint64_vector_t>> adjacency_matrix_vec_vec_t;

      typedef std::pair<const uint64_t, uint64_vector_t> map_value_vec_t;
      typedef boost::unordered_map<uint64_t, uint64_vector_t, boost::hash<uint64_t>, std::equal_to<uint64_t>, SegmentAllocator<map_value_vec_t>> adjacency_matrix_map_vec_t;

      typedef robin_hood_hash<uint64_t, uint64_t, SegmentManager> rhh_single_array_t;

      typedef RHHAdjacencyMatrix<uint64_t, uint64_t, SegmentAllocator<void>> rhh_matrix_t;

      enum DataStructureMode {
        kUseVecVecMatrix,
        kUseMapVecMatrix,
        kUseRHHAsArray,
        kUseRHHAsMatrix,
        kUseHybridDegreeAwareModel
      };

  ///  ------------------------------------------------------ ///
  ///              Constructor / Destructor
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
      inline void add_edges_adjacency_matrix(Container req_itr, size_t length)
      {
        switch(data_structure_type_) {
          case kUseVecVecMatrix:
          add_edges_adjacency_matrix_vector_vector(req_itr, length);
          break;

          case kUseMapVecMatrix:
          add_edges_adjacency_matrix_map_vector(req_itr, length);
          break;

          case kUseRHHAsArray:
          add_edges_rhh_single_array(req_itr, length);
          break;

          case kUseRHHAsMatrix:
          add_edges_rhh_matrix(req_itr, length);
          break;

          case kUseHybridDegreeAwareModel:
          add_edges_degree_aware_hybrid(req_itr, length);
          break;

          default:
          std::cerr << "Unknown data structure type" << std::endl;
          exit(-1);
        }

      }
      void print_profile();
      void reset_profile();

    private:

      static const char kNoValue;

  ///  ------------------------------------------------------ ///
  ///              Private Member Functions
  ///  ------------------------------------------------------ ///

  /// add edges vector-vector adjacency-matrix
  template <typename Container>
      void add_edges_adjacency_matrix_vector_vector(Container req_itr, size_t length);

  /// add edges unsorted_map-vector adjacency-matrix
  template <typename Container>
      void add_edges_adjacency_matrix_map_vector(Container req_itr, size_t length);

  /// add edges array by using robin hood hash
  template <typename Container>
      void add_edges_rhh_single_array(Container req_itr, size_t length);

  template <typename Container>
      void add_edges_rhh_matrix(Container req_itr, size_t length);

  template <typename Container>
      void add_edges_degree_aware_hybrid(Container req_itr, size_t length);

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

      rhh_single_array_t *rhh_single_array;

      rhh_matrix_t *rhh_matrix_;

      RHH::AllocatorsHolder *alloc_holder;
      RHH::RHHMain<uint64_t, uint64_t> *hybrid_matrix;

  // IOInfo *io_info_;
      double total_exectution_time_;

#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
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
