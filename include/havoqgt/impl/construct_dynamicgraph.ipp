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

  void IOInfo::get_status(int &r, int &w) {
    FILE *pipe;
    char str[250];
    pipe = popen("iostat -m | grep md0 2>&1", "r" );

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
 * Constructor : sets secment allocator
 *
 * @param seg_allocator       Reference to segment allocator
 */
template <typename SegementManager>
construct_dynamicgraph<SegementManager>::
construct_dynamicgraph(const SegmentAllocator<void>& seg_allocator, const int mode)
  : adjacency_matrix_vec_vec_(seg_allocator)
  , adjacency_matrix_map_vec_(seg_allocator)
  , data_structure_type_(mode)
  , io_info_() 
{
  total_exectution_time_  = 0.0;
}



template <typename SegementManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph<SegementManager>::
add_edges_adjacency_matrix(
  ManagedMappedFile& asdf,
  const SegmentAllocator<void>& seg_allocator,
  Container& edges)
{
  if (data_structure_type_ == kUseVecVecMatrix) {
    add_edges_adjacency_matrix_vector_vector(asdf, seg_allocator, edges);
  } else if (data_structure_type_ == kUseMapVecMatrix) {
    add_edges_adjacency_matrix_map_vector(asdf, seg_allocator, edges);    
  } else {
    std::cerr << "Unknown data structure type" << std::endl;
    exit(-1);
  }
}


/**
 * Add edges into vector-vector adjacency matrix using
 * boost:interprocess:vector with from and unsorted sequence of edges.
 *
 * @param edges               input edges to partition
 * @param max_vertex          Max vertex ID
*/
template <typename SegementManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph<SegementManager>::
add_edges_adjacency_matrix_vector_vector(
  ManagedMappedFile& asdf,
  const SegmentAllocator<void>& seg_allocator,
  Container& edges)
{
  // XXX: make initializer
  uint64_vector_t init_vec(seg_allocator);
  if (adjacency_matrix_vec_vec_.size() == 0)
      adjacency_matrix_vec_vec_.resize(1, init_vec);

  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;

    while (adjacency_matrix_vec_vec_.size() <= edge.first) {
      adjacency_matrix_vec_vec_.resize(adjacency_matrix_vec_vec_.size() * 2, init_vec);
    }
    
    adjacency_matrix_vec_vec_[edge.first].push_back(edge.second);
  } 
  asdf.flush();
  double time_end = MPI_Wtime();

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;  
  total_exectution_time_ += time_end - time_start;

  io_info_.log_diff();

  //free_edge_container(edges);
}


template <typename SegmentManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph<SegmentManager>::
add_edges_adjacency_matrix_map_vector(
  ManagedMappedFile& asdf, 
  const SegmentAllocator<void>& seg_allocator, 
  Container& edges)
{

  uint64_vector_t init_vec(seg_allocator);

  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;
    auto value = adjacency_matrix_map_vec_.find(edge.first);
    if (value == adjacency_matrix_map_vec_.end()) {
      uint64_vector_t vec(1, edge.second, seg_allocator);
      adjacency_matrix_map_vec_.insert(map_value_t(edge.first, vec));
    } else {
      uint64_vector_t *adjacency_list_vec = &value->second;
      adjacency_list_vec->push_back(edge.second);
    }
  }  
  asdf.flush();
  double time_end = MPI_Wtime();  

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;

  total_exectution_time_ += time_end - time_start;
  io_info_.log_diff();

}

template <typename SegmentManager>
void construct_dynamicgraph<SegmentManager>::
print_profile()
{
  std::cout << "TIME: Total Execution time (sec.) =\t" << total_exectution_time_ << std::endl;
  io_info_.log_diff(true);
}


template <typename SegementManager>
const int construct_dynamicgraph<SegementManager>::kUseVecVecMatrix = 1;
template <typename SegementManager>
const int construct_dynamicgraph<SegementManager>::kUseMapVecMatrix = 2;


} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
