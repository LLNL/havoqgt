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

#ifndef HAVOQGT_MPI_GEN_PREFERENTIAL_ATTACHMENT_EDGE_LIST_HPP_INCLUDED
#define HAVOQGT_MPI_GEN_PREFERENTIAL_ATTACHMENT_EDGE_LIST_HPP_INCLUDED

#include <havoqgt/mpi.hpp>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/detail/preferential_attachment.hpp>
#include <boost/random.hpp>
#include <algorithm>

namespace havoqgt { namespace mpi {


template <typename EdgeType>
void gen_preferential_attachment_edge_list(std::vector<EdgeType>& local_edges,
                                           uint64_t in_base_seed,
                                           uint64_t in_node_scale,
                                           uint64_t in_edge_scale,
                                           double in_beta,
                                           double in_prob_rewire,
                                           MPI_Comm in_comm)
{
  namespace ad = havoqgt::detail;
  int mpi_rank, mpi_size;
  CHK_MPI( MPI_Comm_rank(in_comm, &mpi_rank) );
  CHK_MPI( MPI_Comm_size(in_comm, &mpi_size) );

  //
  // Calc global number of nodes and edges
  uint64_t global_num_nodes = uint64_t(1) << in_node_scale;
  uint64_t global_num_edges = uint64_t(1) << in_edge_scale;
  uint64_t k = global_num_edges / global_num_nodes;
  uint64_t edges_per_rank = global_num_edges / mpi_size;

  ad::preferential_attachment_helper<uint64_t> pa(k, global_num_edges, in_beta,
                                             in_base_seed*mpi_rank + mpi_rank);

  //
  // Generate inital pa edges.  Includes unresolved links
  double init_time_start = MPI_Wtime();
  local_edges.resize(edges_per_rank);
  for(uint64_t i=0; i<edges_per_rank; ++i) {
    uint64_t edge_index = uint64_t(mpi_rank) + i*uint64_t(mpi_size);
    local_edges[i] = pa.gen_edge(edge_index);
  }
  double init_time_end = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Initial rng for edges took " << init_time_end - init_time_start
              << " seconds. " << std::endl;
  }

  //
  // We loop until all edges are resolved
  uint64_t iteration_count(1), count_missing(0), count_total_mising(0);
  do {
    //
    // Count number still missing
    for(uint64_t i=0; i<edges_per_rank; ++i) {
      if(pa.is_pointer(local_edges[i].second)) {
        ++count_total_mising;
      }
    }

    count_missing = 0;
    double time_start = MPI_Wtime();
    //
    // Generate vectors to exchange
    std::vector<std::vector<uint64_t> > tmp_send_to_rank(mpi_size);
    std::vector<std::vector<uint64_t> > tmp_mark_send_to_rank(mpi_size);
    uint64_t count_to_send(0);
    for(uint64_t i=0; i<local_edges.size(); ++i) {
      if(pa.is_pointer(local_edges[i].second)) {
        int owner = pa.value_of_pointer(local_edges[i].second) % mpi_size;
        tmp_mark_send_to_rank[owner].push_back(i);
        tmp_send_to_rank[owner].push_back(pa.value_of_pointer(local_edges[i].second));
        ++count_to_send;
      }
    }
    //
    // Combine send/recv vectors
    std::vector<uint64_t> to_send(count_to_send); to_send.reserve(1);
    std::vector<uint64_t> mark_to_send(count_to_send);
    std::vector<int>      sendcnts(mpi_size,0); sendcnts.reserve(1);
    to_send.clear(); mark_to_send.clear();
    for(int i = 0; i < mpi_size; ++i) {
      for(size_t j = 0; j < tmp_send_to_rank[i].size(); ++j) {
        to_send.push_back(tmp_send_to_rank[i][j]);
        mark_to_send.push_back(tmp_mark_send_to_rank[i][j]);
      }
      sendcnts[i] = tmp_send_to_rank[i].size();
    }
    //
    // Release memory from temp vectors
    {
      std::vector<std::vector<uint64_t> > empty;
      tmp_send_to_rank.swap(empty);
    }
    {
      std::vector<std::vector<uint64_t> > empty;
      tmp_mark_send_to_rank.swap(empty);
    }
    //
    // Exchange vectors
    std::vector<uint64_t> to_recv; // not sized
    std::vector<int>      recvcnts; // not sized
    havoqgt::mpi::mpi_all_to_all(to_send, sendcnts, to_recv, recvcnts, MPI_COMM_WORLD);
    //
    // Look up pointers, pointer jump!
    for(size_t i=0; i<to_recv.size(); ++i) {
      uint64_t local_index = to_recv[i]/mpi_size;
      to_recv[i] = local_edges[local_index].second;
    }
    //
    // Return exchange vectors
    havoqgt::mpi::mpi_all_to_all(to_recv, recvcnts, to_send, sendcnts, MPI_COMM_WORLD);
    //
    // Return new pointers to marked location
    for(size_t i=0; i<to_send.size(); ++i) {
        local_edges[ mark_to_send[i] ].second = to_send[i];;
    }
    //
    // Count number still missing
    for(uint64_t i=0; i<edges_per_rank; ++i) {
      if(pa.is_pointer(local_edges[i].second)) {
        ++count_missing;
      }
    }
    count_missing = havoqgt::mpi::mpi_all_reduce(count_missing, std::plus<uint64_t>(), MPI_COMM_WORLD);

    //
    // Iteration ouput
    double time_end = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "Iteration " << iteration_count << " took "
                << time_end - time_start << " seconds.   "
                << count_missing << " still missing. " << std::endl;
    }
    ++iteration_count;
  } while(count_missing > 0);
  //
  // Output total nuumber of global exchanges
  count_total_mising = havoqgt::mpi::mpi_all_reduce(count_total_mising, std::plus<uint64_t>(), MPI_COMM_WORLD);
  if(mpi_rank == 0) {
    std::cout << "Total missing (how much data globally exchanged) = " << count_total_mising << std::endl;
  }


  //
  // Randomly rewire
  if(in_prob_rewire > double(0)) {
    boost::mt19937 rng(in_base_seed + uint64_t(mpi_rank) * 3ULL);
    boost::random::uniform_int_distribution<> rand_node(0,global_num_nodes-1);
    boost::random::uniform_01<> rand_real;
    for(size_t i=0; i<local_edges.size(); ++i) {
      if(rand_real(rng) < in_prob_rewire) {
        EdgeType rand_edge(rand_node(rng), rand_node(rng));
        local_edges[i] = rand_edge;
      }
    }
  }

  //
  // Scramble edges
  for(size_t i=0; i<local_edges.size(); ++i) {
    //TODO:  This needs a mod because we never exactly make the correct
    //number of vertices.   We need a better solution!
    local_edges[i].first  %= global_num_nodes;
    local_edges[i].second %= global_num_nodes;
    local_edges[i].first  = ad::hash_nbits(local_edges[i].first,  in_node_scale);
    local_edges[i].second = ad::hash_nbits(local_edges[i].second, in_node_scale);
  }

  // //
  // // Symmetrizing b/c PA is always undirected.
  // uint64_t old_size = local_edges.size();
  // local_edges.reserve(old_size * 2);
  // for(uint64_t i=0; i<old_size; ++i) {
  //   EdgeType edge;
  //   edge.first = local_edges[i].second;
  //   edge.second = local_edges[i].first;
  //   local_edges.push_back( edge );
  // }

  // //
  // // Shuffle edges because we need to mix the hubs in for future partitioning
  // std::random_shuffle(local_edges.begin(), local_edges.end());
}

}} //end namespace havoqgt::mpi


#endif //end HAVOQGT_MPI_GEN_PREFERENTIAL_ATTACHMENT_EDGE_LIST_HPP_INCLUDED

