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


  template <typename SegementManager>
    const char construct_dynamicgraph<SegementManager>::kNoValue = 0;

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
    adjacency_matrix_vec_vec_ = new adjacency_matrix_vec_vec_t(seg_allocator_);
    init_vec = new uint64_vector_t(seg_allocator_);
    break;

    case kUseMapVecMatrix:
    adjacency_matrix_map_vec_ = new adjacency_matrix_map_vec_t(seg_allocator_);
    break;

    case kUseRobinHoodHash:
    robin_hood_hashing_ = new robin_hood_hashing_t(seg_allocator_);
    break;

    case kUseDegreeAwareModel:
    adjacency_matrix_map_vec_ = new adjacency_matrix_map_vec_t(seg_allocator_);
    robin_hood_hashing_ = new robin_hood_hashing_t(seg_allocator_);
    break;

    case kUseDegreeAwareModelRbhMtx:
    robin_hood_hashing_ = new robin_hood_hashing_t(seg_allocator_);
    adjacency_matrix_rbh_rbh_ = new adjacency_matrix_rbh_rbh_t(seg_allocator_);
    std::srand(1);
    std::cout << "Random seed is 1" << std::endl;
    // std::srand(std::time(0));
    // std::cout << "Random seed is time" << std::endl;
    break;

    default:
    std::cerr << "Unknown data structure type" << std::endl;
    assert(false);
    exit(-1);
  }

  io_info_ = new IOInfo();
  total_exectution_time_  = 0.0;

#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
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
  delete adjacency_matrix_vec_vec_;
  delete init_vec;
  delete adjacency_matrix_map_vec_;
  delete robin_hood_hashing_;
  delete io_info_;
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
 add_edges_adjacency_matrix_vector_vector(Container req_itr, size_t length)
 {
  // TODO: make initializer
  if (adjacency_matrix_vec_vec_->size() == 0) {
    adjacency_matrix_vec_vec_->resize(1, *init_vec);
  }

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (size_t k = 0; k < length; ++k, ++req_itr) {
    const auto &edge = req_itr->edge;

    while (adjacency_matrix_vec_vec_->size() <= edge.first) {
      adjacency_matrix_vec_vec_->resize(adjacency_matrix_vec_vec_->size() * 2, *init_vec);
    }
    uint64_vector_t& adjacency_list_vec = adjacency_matrix_vec_vec_->at(edge.first);
#if WITHOUT_DUPLICATE_INSERTION == 1
    /// add a edge without duplication
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

}

template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_adjacency_matrix_map_vector(Container req_itr, size_t length)
{

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (size_t k = 0; k < length; ++k, ++req_itr) {
    const auto &edge = req_itr->edge;
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
add_edges_robin_hood_hash(Container req_itr, size_t length)
{
  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (size_t k = 0; k < length; ++k, ++req_itr) {
    const auto& edge = req_itr->edge;
#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
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
}


template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_degree_aware_rbh_first(Container req_itr, size_t length)
{
  uint64_t count_new = 0;
  uint64_t count_low = 0;
  uint64_t count_non_low = 0;
  uint64_t count_move = 0;

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (size_t k = 0; k < length; ++k, ++req_itr) {

    const int64_t source_vtx = req_itr->edge.first;
    const int64_t target_vtx = req_itr->edge.second;

#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
    fout_debug_insertededges_  << source_vtx << "\t" << target_vtx << std::endl;
#endif

    const uint64_t num_edges_rbh = robin_hood_hashing_->count(source_vtx);

    if (num_edges_rbh == 0) {
      // ---- These are 2 candidates ---- //
      // 1. new vertex
      // 2. non-low-degree edges

      auto value = adjacency_matrix_map_vec_->find(source_vtx);
      if (value == adjacency_matrix_map_vec_->end()) {
        // -- 1. new vertex -- //
        robin_hood_hashing_->insert(source_vtx, target_vtx);
        ++count_new;
      } else { 
        // -- 2. non-low-degree edges -- //
        uint64_vector_t& adjacency_list_vec = value->second;
#if WITHOUT_DUPLICATE_INSERTION == 1
        if (boost::find<uint64_vector_t>(adjacency_list_vec, target_vtx) != adjacency_list_vec.end() )
          continue; // Since this edge is duplicated, skip adding this edge.
#endif
        adjacency_list_vec.push_back(target_vtx);
        ++count_non_low;
      }


    } else if (num_edges_rbh < kLowDegreeThreshold) {
      // -- Low-degree edge -- //
#if WITHOUT_DUPLICATE_INSERTION == 1
      robin_hood_hashing_->insert_unique(source_vtx, target_vtx);
#else
      robin_hood_hashing_->insert(source_vtx, target_vtx);
#endif
      ++count_low;

    } else {
      // -- degree exceed the low-degree-threshold -- //
      auto itr_map = adjacency_matrix_map_vec_->insert(map_value_vec_t(source_vtx, uint64_vector_t(seg_allocator_))).first;
      uint64_vector_t &adjacency_list_vec = itr_map->second;
      adjacency_list_vec.reserve(kLowDegreeThreshold+1);

      // - Move edges to the adjacency-matrix from the robin-hood-hashing - //
      auto itr_rb = robin_hood_hashing_->find(source_vtx);
      for (int64_t i = 0; i < kLowDegreeThreshold; ++i) {
        auto trg_vtx = itr_rb->value;
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
      adjacency_list_vec.push_back(target_vtx);
      ++count_move;
    }

  } // End of a edges insertion loop
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;
  std::cout << "Count: # new inserted edges =\t" << count_new << std::endl;
  std::cout << "Count: # low degree edges =\t" << count_low << std::endl;
  std::cout << "Count: # non-low degree edges =\t" << count_non_low << std::endl;
  std::cout << "Count: # moved vertices =\t" << count_move << std::endl;
  std::cout << "Count: # moved edges =\t" << count_move * kLowDegreeThreshold << std::endl;
  std::cout << "Allocated size for low-degree-edges =\t" << robin_hood_hashing_->allocated_size() << std::endl;
  io_info_->log_diff();
  total_exectution_time_ += time_end - time_start;
}


template <typename SegmentManager>
template <typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_degree_aware_rbh_first_rbh_mtrx(Container req_itr, size_t length)
{
  uint64_t count_new = 0;
  uint64_t count_low = 0;
  uint64_t count_non_low = 0;
  uint64_t count_move = 0;
  uint64_t count_delete = 0;

  io_info_->reset_baseline();
  double time_start = MPI_Wtime();
  for (size_t k = 0; k < length; ++k, ++req_itr) {

    const uint64_t source_vtx = req_itr->edge.first;
    const uint64_t target_vtx = req_itr->edge.second;

#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
    fout_debug_insertededges_ << source_vtx << "\t" << target_vtx << "\t0" << std::endl;
#endif

    if (req_itr->is_delete) {
      count_delete += delete_edge_from_rbh_rbh_mtrx(source_vtx, target_vtx);
#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
      fout_debug_insertededges_ << source_vtx << "\t" << target_vtx << "\t1" << std::endl;
#endif
      continue;
    }

    const uint64_t num_edges_rbh = robin_hood_hashing_->count(source_vtx);
    if (num_edges_rbh == 0) {
      auto itr_rbh_mtrx = adjacency_matrix_rbh_rbh_->find(source_vtx);
      if (!itr_rbh_mtrx.is_valid_index()) {
        // -- 1. new vertex -- //
        robin_hood_hashing_->insert(source_vtx, target_vtx);

        ++count_new;
      } else {
        // -- 2. Non-low-degree edges -- //
        auto& adjacency_list_rbh = itr_rbh_mtrx->value;


#if WITHOUT_DUPLICATE_INSERTION == 1
        count_non_low += adjacency_list_rbh.insert_unique(target_vtx, kNoValue);
#else
        adjacency_list_rbh.insert(target_vtx, kNoValue);
        ++count_non_low;
#endif
      }


    } else if (num_edges_rbh < kLowDegreeThreshold) {
      // -- 3. Low-degree edge -- //
#if WITHOUT_DUPLICATE_INSERTION == 1
      count_low += robin_hood_hashing_->insert_unique(source_vtx, target_vtx);
#else
      robin_hood_hashing_->insert(source_vtx, target_vtx);
      ++count_low;
#endif
    } else {
      // -- 4. degree exceed the low-degree-threshold -- //

      /// If a new edge is same as existed edge which are stored in low-degree data structure
      /// and low-degree threshold is 1, we don't make instance into adjacency-matrix.
      auto itr_rb = robin_hood_hashing_->find(source_vtx);
#if WITHOUT_DUPLICATE_INSERTION == 1
      if (itr_rb->value == target_vtx && num_edges_rbh == 1) continue;
#endif
      adjacency_matrix_rbh_rbh_->insert(source_vtx, robin_hood_hashing_novalue_t(seg_allocator_));
      robin_hood_hashing_novalue_t &adjacency_list_rbh = adjacency_matrix_rbh_rbh_->find(source_vtx)->value;
      adjacency_list_rbh.insert(target_vtx, kNoValue);
      ++count_move;

      // - Move edges to the adjacency-matrix from the robin-hood-hashing - //
      for (int64_t i = 0; i < kLowDegreeThreshold; ++i, ++itr_rb) {
       const uint64_t trg_vtx = itr_rb->value;
#if WITHOUT_DUPLICATE_INSERTION == 1
        if (trg_vtx == target_vtx) goto NEXT_MOVING;
#endif
        adjacency_list_rbh.insert(trg_vtx, kNoValue);
        ++count_move;

        NEXT_MOVING:
        robin_hood_hashing_->erase(itr_rb); // Delete the edge from robin_hood_hashing array
      } // End of the edges moving loop
    }

  } // End of a edges insertion loop
  flush_pagecache();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;
  std::cout << "Count: # new inserted edges =\t" << count_new << std::endl;
  std::cout << "Count: # low degree edges =\t" << count_low << std::endl;
  std::cout << "Count: # non-low degree edges =\t" << count_non_low << std::endl;
  std::cout << "Count: # moved edges =\t" << count_move << std::endl;
  std::cout << "Count: # deleted edges =\t" << count_delete << std::endl;
  std::cout << "Allocated size for low-degree-edges =\t" << robin_hood_hashing_->allocated_size() << std::endl;
  std::cout << "Allocated size for non-low-degree-edges =\t" << adjacency_matrix_rbh_rbh_->allocated_size() << std::endl;
  io_info_->log_diff();
  total_exectution_time_ += time_end - time_start;
}

template <typename SegmentManager>
bool construct_dynamicgraph<SegmentManager>::
delete_edge_from_rbh_rbh_mtrx(const int64_t source_vtx, const int64_t target_vtx)
{

  /// Delete from the robin hood hashing 1-D array
  if (robin_hood_hashing_->erase(source_vtx, target_vtx)) { return true; }
  
  /// If there is no 'source_vtx', then return false.
  auto itr_mtrx_rb = adjacency_matrix_rbh_rbh_->find(source_vtx);
  if (!itr_mtrx_rb.is_valid_index()) { return false; }

  /// Delete from the adjavency-list
  /// If the edge is not existed, then return false.
  auto& adjacency_list_rbh = itr_mtrx_rb->value;
  if (!adjacency_list_rbh.erase(target_vtx, kNoValue)) { return false; }
  
  /// If degree on the adjacency-list become equal to low-degree threshold,
  /// move edges from the adjacency-list to the 1-D hash array.
  uint64_t degree = adjacency_list_rbh.size();
  if (degree == kLowDegreeThreshold) {
    auto itr = adjacency_list_rbh.begin();
    for (int64_t i = 0; i < kLowDegreeThreshold; ++i, ++itr) {
      robin_hood_hashing_->insert(source_vtx, itr->key); /// Note: we don't need to call unique insert function.
      //adjacency_list_rbh.erase(itr);
    }
    /// delete a element from the adjacency-matrix 
    adjacency_matrix_rbh_rbh_->erase(itr_mtrx_rb);
  }

  return true;
}


template <typename SegmentManager>
void construct_dynamicgraph<SegmentManager>::
print_profile()
{
  std::cout << "TIME: Total Execution time (sec.) =\t" << total_exectution_time_ << std::endl;
  io_info_->log_diff(true);

  std::cout << "WITHOUT_DUPLICATE_INSERTION is " << WITHOUT_DUPLICATE_INSERTION << std::endl;
  std::cout << "DEBUG_DUMPUPDATEREQUESTANDRESULTS is " << DEBUG_DUMPUPDATEREQUESTANDRESULTS << std::endl;
#if DEBUG_DUMPUPDATEREQUESTANDRESULTS
  std::cout << "kFnameDebugInsertedEdges is " << kFnameDebugInsertedEdges << std::endl;
#endif
  if (data_structure_type_ == kUseRobinHoodHash ||
      data_structure_type_ == kUseDegreeAwareModel ||
      data_structure_type_ == kUseDegreeAwareModelRbhMtx) {
    std::cout << "USE_SEPARATE_HASH_ARRAY is " << USE_SEPARATE_HASH_ARRAY << std::endl;
    std::cout << "USE_TOMBSTONE is " << USE_TOMBSTONE << std::endl;
    std::cout << "# elements in Robin-Hood-Hashing = " << robin_hood_hashing_->size() << std::endl;
  }

#if DEBUG_DUMPUPDATEREQUESTANDRESULTS == 1
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
  if (data_structure_type_ == kUseRobinHoodHash || data_structure_type_ == kUseDegreeAwareModel || data_structure_type_ == kUseDegreeAwareModelRbhMtx) {
    std::cout << "# elements in Robin-Hood-Hashing = " << robin_hood_hashing_->size() << std::endl;
    robin_hood_hashing_->dump_elements(kFnameDebugInsertedEdges+"_graph");
    //robin_hood_hashing_->disp_elements();
  }
  if (data_structure_type_ == kUseDegreeAwareModelRbhMtx) {
    std::ofstream fout;
    fout.open(kFnameDebugInsertedEdges+"_graph", std::ios::out | std::ios::app);
    for(auto itr = adjacency_matrix_rbh_rbh_->begin(), itr_end=adjacency_matrix_rbh_rbh_->end();
        itr!=itr_end;
        ++itr) 
    {
          auto& list = itr->value;
          for (auto itr_lst = list.begin(), itr_lst_end = list.end();
            itr_lst != itr_lst_end;
            ++itr_lst) 
          {
              fout << itr->key << "\t" << itr_lst->key << std::endl;
          }
    }
    fout.close();
    //robin_hood_hashing_->disp_elements();
  }

  fout_debug_insertededges_.close();

#endif
}

template <typename SegmentManager>
void construct_dynamicgraph<SegmentManager>::
reset_profile()
{
  total_exectution_time_ = 0.0;
}

} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
