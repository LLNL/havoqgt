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

#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED


/**
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */

namespace havoqgt {
namespace mpi {

  static const std::string kDeviceName = "md0";

  IOInfo::IOInfo()
  : read_total_mb_(0.0), written_total_mb_(0.0)
  {
    init();
  }

  void IOInfo::init() {
    read_total_mb_ = 0.0;
    written_total_mb_ = 0.0;
    get_status(read_previous_mb_, written_previous_mb_);
  }

  void IOInfo::reset_baseline() {
    get_status(read_previous_mb_, written_previous_mb_);
  }

  void IOInfo::get_status(int &r, int &w) {
    FILE *pipe;
    char str[250];
    std::string fname = "iostat -m | grep " + kDeviceName + " 2>&1";
    pipe = popen(fname.c_str(), "r" );

    float temp;
    fscanf(pipe, "%s", str);
    fscanf(pipe, "%f %f %f", &temp, &temp, &temp);
    fscanf(pipe, "%d %d \n", &r, &w);
    pclose(pipe);
  };

  void IOInfo::log_diff(bool final = false) {
    int read_current_mb, written_current_mb;
    get_status(read_current_mb, written_current_mb);
    
    int read    = read_current_mb    - read_previous_mb_;
    int written = written_current_mb - written_previous_mb_;

    read_total_mb_    += read;
    written_total_mb_ += written;

    std::cout << "MB Read:\t"     << read    << std::endl;
    std::cout << "MB Written:\t"  << written << std::endl;
    if (final) {
      std::cout << "Total MB Read:\t"    << read_total_mb_     << std::endl;
      std::cout << "Total MB Written:\t" << written_total_mb_  << std::endl;    
    }

    read_previous_mb_    = read_current_mb;
    written_previous_mb_ = written_current_mb;

  };



/**
 * Constructor
 *
 * @param seg_allocator       Reference to segment allocator
 */
template <typename SegementManager>
construct_dynamicgraph<SegementManager>::
construct_dynamicgraph (
  mapped_t& asdf,
  SegmentAllocator<void>& seg_allocator,
  const DataStructureMode mode,
  const uint64_t n )
  : asdf_(asdf)
  , seg_allocator_(seg_allocator)
  , data_structure_type_(mode)
  , kLowDegreeThreshold(n)
{
  
  switch(data_structure_type_) {
    case kUseVecVecMatrix:
      adjacency_matrix_vec_vec_ = new adjacency_matrix_vec_vec_t(seg_allocator);
      init_vec = new uint64_vector_t(seg_allocator);
      break;

    case kUseMapVecMatrix:
      adjacency_matrix_map_vec_ = new adjacency_matrix_map_vec_t(seg_allocator);
      break;

    case kUseRobinHoodHash:
      robin_hood_hashing_ = new robin_hood_hashing_t(seg_allocator);
      break;

    case kUseDegreeAwareModel:
      adjacency_matrix_map_vec_ = new adjacency_matrix_map_vec_t(seg_allocator);
      robin_hood_hashing_ = new robin_hood_hashing_t(seg_allocator);
      is_exist_bmp_ = new bitmap_mgr();
      break;

    default:
      std::cerr << "Unknown data structure type" << std::endl;
      assert(false);
      exit(-1);
  }

  io_info_ = new IOInfo();
  total_exectution_time_  = 0.0;

#if DEBUG_INSERTEDEDGES == 1
  fout_debug_insertededges_.open(kFnameDebugInsertedEdges+"_raw");
#endif

}

/**
 * Deconstructor
 */
template <typename SegementManager>
construct_dynamicgraph<SegementManager>::
~construct_dynamicgraph()
{
  if (adjacency_matrix_vec_vec_ != NULL) delete(adjacency_matrix_vec_vec_);
  if (init_vec != NULL) delete(init_vec);
  if (adjacency_matrix_map_vec_ != NULL) delete(adjacency_matrix_map_vec_);
  if (robin_hood_hashing_ != NULL) delete(robin_hood_hashing_);
  if (io_info_ != NULL) delete(io_info_);
}

/**
 * Add edges into vector-vector adjacency matrix using
 * boost:interprocess:vector with from and unsorted sequence of edges.
 *
 * @param edges               input edges to partition
*/
template <typename SegementManager>
template <typename Container>
void construct_dynamicgraph<SegementManager>::
add_edges_adjacency_matrix_vector_vector(Container& edges)
{
  // TODO: make initializer
  if (adjacency_matrix_vec_vec_->size() == 0) {
      adjacency_matrix_vec_vec_->resize(1, *init_vec);
  }

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;

    while (adjacency_matrix_vec_vec_->size() <= edge.first) {
      adjacency_matrix_vec_vec_->resize(adjacency_matrix_vec_vec_->size() * 2, *init_vec);
    }
    uint64_vector_t& adjacency_list_vec = adjacency_matrix_vec_vec_->at(edge.first);
#if WITHOUT_DUPLICATE_INSERTION == 1
    // add a edge without duplication
    if (boost::find<uint64_vector_t>(adjacency_list_vec, edge.second) == adjacency_list_vec.end()) {
      adjacency_list_vec.push_back(edge.second);
    }
#else
    adjacency_list_vec.push_back(edge.second);
#endif
  }
  flush_pagecache();
  double time_end = MPI_Wtime();

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;  
  total_exectution_time_ += time_end - time_start;

  io_info_->log_diff();

  //free_edge_container(edges);
}


template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_adjacency_matrix_map_vector(Container& edges)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;
    auto value = adjacency_matrix_map_vec_->find(edge.first);
    if (value == adjacency_matrix_map_vec_->end()) { // new vertex
      uint64_vector_t vec(1, edge.second, seg_allocator_);
      adjacency_matrix_map_vec_->insert(map_value_vec_t(edge.first, vec));
    } else {
      uint64_vector_t& adjacency_list_vec = value->second;

#if WITHOUT_DUPLICATE_INSERTION == 1
      // add a edge without duplication
      if (boost::find<uint64_vector_t>(adjacency_list_vec, edge.second) != adjacency_list_vec.end() )
        continue;
#endif

      adjacency_list_vec.push_back(edge.second);

    }
  }  
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_->log_diff();

}

template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_robin_hood_hash(Container& edges)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;
#if DEBUG_INSERTEDEDGES == 1
    fout_debug_insertededges_ << edge.first << "\t" << edge.second << std::endl;
#endif
#if WITHOUT_DUPLICATE_INSERTION == 1
    robin_hood_hashing_->insert_unique(edge.first, edge.second);
#else
    robin_hood_hashing_->insert(edge.first, edge.second);
#endif
  }
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_->log_diff();
  //robin_hood_hashing_->disp_elements();

}


template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_degree_aware_rbhs_first(Container& edges)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {

    const int64_t source_vtx = itr->first;
    const int64_t target_vtx = itr->second;

#if DEBUG_INSERTEDEDGES == 1
        fout_debug_insertededges_  << source_vtx << "\t" << target_vtx << std::endl;
#endif

    const uint64_t num_edges_rbhs = robin_hood_hashing_->count(source_vtx);

    if (num_edges_rbhs == 0) {
      // ---- These are 2 candidates ---- //
      // 1. new vertex
      // 2. non-low-degree edges

      auto value = adjacency_matrix_map_vec_->find(source_vtx);
      if (value == adjacency_matrix_map_vec_->end()) {
        // -- 1. new vertex -- //
        robin_hood_hashing_->insert(source_vtx, target_vtx);
      } else { 
        // -- 2. non-low-degree edges -- //
        uint64_vector_t& adjacency_list_vec = value->second;
#if WITHOUT_DUPLICATE_INSERTION == 1
        if (boost::find<uint64_vector_t>(adjacency_list_vec, target_vtx) != adjacency_list_vec.end() )
          continue; // Since this edge is duplicated, skip adding this edge.
#endif
        adjacency_list_vec.push_back(target_vtx);
      }

    } else if (num_edges_rbhs < kLowDegreeThreshold) {
      // -- Low-degree edge -- //
#if WITHOUT_DUPLICATE_INSERTION == 1
      robin_hood_hashing_->insert_unique(source_vtx, target_vtx);
#else
      robin_hood_hashing_->insert(source_vtx, target_vtx);
#endif

    } else {
      // -- degree exceed the low-degree-threshold -- //
      // Note: this implementation dose not care about order of edges insertion.
      uint64_vector_t adjacency_list_vec(1, target_vtx, seg_allocator_);

      // Move edges to the adjacency-matrix from the robin-hood-hashing //
      auto itr_rb = robin_hood_hashing_->find(source_vtx);
      while (itr_rb.is_valid_index()) {
         const int64_t trg_vtx = *itr_rb;
#if WITHOUT_DUPLICATE_INSERTION == 1
         // Since we have already avoided duplicated edges when we add edges into robin_hood_hashing array,
         // in this time, we only compare to the latest edge.
         if (trg_vtx == target_vtx) goto NEXT_MOVING; // a duplicated edge
#endif
        adjacency_list_vec.push_back(trg_vtx);

NEXT_MOVING:
        robin_hood_hashing_->erase(itr_rb); // Delete the edge from robin_hood_hashing array
        ++itr_rb;
      } // End of edges moving loop
      adjacency_matrix_map_vec_->insert(map_value_vec_t(source_vtx, adjacency_list_vec));      
    }

  } // End of a edges insertion step
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_->log_diff();
}



template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_degree_aware_bitmap_first(Container& edges)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {

    const int64_t source_vtx = itr->first;
    const int64_t target_vtx = itr->second;

#if DEBUG_INSERTEDEDGES
        fout_debug_insertededges_  << source_vtx << "\t" << target_vtx << std::endl;
#endif
    //std::cout << source_vtx << "\t" << target_vtx << std::endl;
    
    const bool is_exist = is_exist_bmp_->get(source_vtx);
    if (is_exist) {
      // -- 1. Non-new vertex -- //

      const uint64_t num_edges_rbhs = robin_hood_hashing_->count(source_vtx);
      if (num_edges_rbhs < kLowDegreeThreshold) {
        // 1-1. Low-degree edge or non-low-degree edge
        
        auto value = adjacency_matrix_map_vec_->find(source_vtx);
        if (value == adjacency_matrix_map_vec_->end()) {
          // -- 1-1-1. Low-degree edge -- //
#if WITHOUT_DUPLICATE_INSERTION == 1
          robin_hood_hashing_->insert_unique(source_vtx, target_vtx);
#else
           robin_hood_hashing_->insert(source_vtx, target_vtx);
#endif
        } else {  
          // -- 1-1-2. non-low-degree edges -- //
          uint64_vector_t& adjacency_list_vec = value->second;
  #if WITHOUT_DUPLICATE_INSERTION
          if (boost::find<uint64_vector_t>(adjacency_list_vec, target_vtx) != adjacency_list_vec.end() )
            continue; // Since this edge is duplicated, skip adding this edge.
  #endif
          adjacency_list_vec.push_back(target_vtx);
        }

      } else {
        // -- 1-2. degree exceed the low-degree-threshold -- //
        // Note: this implementation dose not care about order of edges insertion.
        uint64_vector_t adjacency_list_vec(1, target_vtx, seg_allocator_);

        // Move edges to the adjacency-matrix from the robin-hood-hashing //
        auto itr_rb = robin_hood_hashing_->find(source_vtx);
        while (itr_rb.is_valid_index()) {
           const int64_t trg_vtx = *itr_rb;
#if WITHOUT_DUPLICATE_INSERTION
          // Since we have already avoided duplicated edges when we add edges into robin_hood_hashing array,
          // in this time, we only compare to the latest edge.
          if (trg_vtx == target_vtx) goto NEXT_MOVING; // a duplicated edge
#endif
          adjacency_list_vec.push_back(trg_vtx);

NEXT_MOVING:
          robin_hood_hashing_->erase(itr_rb); // Delete the edge from robin_hood_hashing array
          ++itr_rb;
        } // End of edges moving loop
        adjacency_matrix_map_vec_->insert(map_value_vec_t(source_vtx, adjacency_list_vec));      
      }


    } else {
      // -- 2. new vertex -- //

      robin_hood_hashing_->insert(source_vtx, target_vtx);
      is_exist_bmp_->set(source_vtx);
    }

  } // End of a edges insertion step
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_->log_diff();
}



template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_degree_aware_adj_first(Container& edges)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const int64_t source_vtx = itr->first;
    const int64_t target_vtx = itr->second;

#if DEBUG_INSERTEDEDGES == 1
        fout_debug_insertededges_  << source_vtx << "\t" << target_vtx << std::endl;
#endif


    // Check whether we already have edges at adjacency-matrix
    auto value = adjacency_matrix_map_vec_->find(source_vtx);
    if (value == adjacency_matrix_map_vec_->end()) {
#if 0
      int64_t degree = ++degree_map[edge.first];
#else
      const uint64_t new_degree = robin_hood_hashing_->count(source_vtx) + 1;
#endif

      // If new_degree is kLowDegreeThreshold or less, add the edge into robin_hood_hashing array.
      // If new_degree is more than  kLowDegreeThreshold, move edges from robin-hood to adjacency-matrix
      if (new_degree <= kLowDegreeThreshold) {

#if WITHOUT_DUPLICATE_INSERTION == 1
        robin_hood_hashing_->insert_unique(source_vtx, target_vtx);
#else
        robin_hood_hashing_->insert(source_vtx, target_vtx);
#endif

      } else {
        // -- move edges to adhacency-matrix -- //

        // Since this is a first time to add source_vtx's edges into adjacency-matrix,
        // allocate source_vtx's adjacency-list-vector
        // Note: this implementation dose not care about order of edges insertion.
        uint64_vector_t adjacency_list_vec(1, target_vtx, seg_allocator_);
        

        auto itr_rb = robin_hood_hashing_->find(source_vtx);
        while (itr_rb.is_valid_index()) {
           const int64_t trg_vtx = *itr_rb;
           //std::cout << trg_vtx;
#if WITHOUT_DUPLICATE_INSERTION == 1
           // Since we have already avoided duplicated edges when we add edges into robin_hood_hashing array,
           // in this time, we only compare to the latest edge.
           if (trg_vtx == target_vtx) goto NEXT_MOVING; // a duplicated edge
#endif
          adjacency_list_vec.push_back(trg_vtx);

NEXT_MOVING:
          robin_hood_hashing_->erase(itr_rb); // Delete the edge from robin_hood_hashing array
          ++itr_rb;
        } // End of edges moving loop

        adjacency_matrix_map_vec_->insert(map_value_vec_t(source_vtx, adjacency_list_vec));
      }

    } else { // We already have edges at adjacency-matrix
      // If we already have edges at adjacency-matrix, 
      // degree of the edge is more than kDegreeThreshold.
      uint64_vector_t& adjacency_list_vec = value->second;
#if WITHOUT_DUPLICATE_INSERTION == 1
      if (boost::find<uint64_vector_t>(adjacency_list_vec, target_vtx) != adjacency_list_vec.end() )
        continue; // because this edge is duplicated, goto next edge.
#endif
      adjacency_list_vec.push_back(target_vtx);
    }

  } // End of a edges insertion step
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_->log_diff();
}




template <typename SegmentManager>
void construct_dynamicgraph<SegmentManager>::
print_profile()
{
  std::cout << "TIME: Total Execution time (sec.) =\t" << total_exectution_time_ << std::endl;
  io_info_->log_diff(true);

#if DEBUG_INSERTEDEDGES == 1
  std::ofstream tmp;
  tmp.open(kFnameDebugInsertedEdges+"_graph", std::ios::trunc);
  tmp.close();

  if (data_structure_type_ == kUseMapVecMatrix || data_structure_type_ == kUseDegreeAwareModel) {
    std::ofstream fout;
    fout.open(kFnameDebugInsertedEdges+"_graph", std::ios::out | std::ios::app);
    for (auto itr = adjacency_matrix_map_vec_->begin(); itr != adjacency_matrix_map_vec_->end(); ++itr) {
      uint64_vector_t& adjacency_list_vec = (*itr).second;
      for (auto itr2 = adjacency_list_vec.begin(); itr2 != adjacency_list_vec.end(); ++itr2) {
        fout << (*itr).first << "\t" << *itr2 << std::endl;
      }
    }
    fout.close();
  }
  if (data_structure_type_ == kUseRobinHoodHash || data_structure_type_ == kUseDegreeAwareModel) {
    robin_hood_hashing_->dump_elements(kFnameDebugInsertedEdges+"_graph");
  }
  fout_debug_insertededges_.close();
#endif
}


} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
