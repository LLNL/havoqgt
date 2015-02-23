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

// /usr/local/tools/mvapich2-gnu-1.9/bin/mpicxx ../havoqgt/scripts/catalyst.llnl.gov/sort_webgraph.cpp -std=c++11 -O3
// export MPICH_CXX=/opt/rh/devtoolset-1.1/root/usr/bin/g++

#include <sched.h>
#include <vector>
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <utility>
#include <string>

#define CHK_MPI(a) { if (a != MPI_SUCCESS) {\
                      char* error_string = NULL; \
                      int len = 0; \
                      MPI_Error_string(a, error_string, &len); \
                      std::cerr << __FILE__ << ", line " << __LINE__  \
                           <<" MPI ERROR = " << error_string << std::endl; \
                           exit(-1); \
                     } }

inline MPI_Datatype mpi_typeof(char) {return MPI_CHAR;}
inline MPI_Datatype mpi_typeof(signed short) {return MPI_SHORT;}
inline MPI_Datatype mpi_typeof(signed int) {return MPI_INT;}
inline MPI_Datatype mpi_typeof(signed long) {return MPI_LONG;}
inline MPI_Datatype mpi_typeof(unsigned char) {return MPI_UNSIGNED_CHAR;}
inline MPI_Datatype mpi_typeof(unsigned short) {return MPI_UNSIGNED_SHORT;}
inline MPI_Datatype mpi_typeof(unsigned) {return MPI_UNSIGNED;}
inline MPI_Datatype mpi_typeof(unsigned long) {return MPI_UNSIGNED_LONG;}
inline MPI_Datatype mpi_typeof(unsigned long long) {return MPI_UNSIGNED_LONG_LONG; }
inline MPI_Datatype mpi_typeof(signed long long) {return MPI_LONG_LONG_INT;}
inline MPI_Datatype mpi_typeof(double) {return MPI_DOUBLE;}
inline MPI_Datatype mpi_typeof(long double) {return MPI_LONG_DOUBLE;}
inline MPI_Datatype mpi_typeof(std::pair<int,int>) {return MPI_2INT;}
inline MPI_Datatype mpi_typeof(std::pair<float,int>) {return MPI_FLOAT_INT;}
inline MPI_Datatype mpi_typeof(std::pair<double,int>) {return MPI_DOUBLE_INT;}
inline MPI_Datatype mpi_typeof(std::pair<long double,int>) {return MPI_LONG_DOUBLE_INT;}
inline MPI_Datatype mpi_typeof(std::pair<short,int>) {return MPI_SHORT_INT;}

///
/// All to All exchange of std::vector<T>
///
/// \param[in]  in_p_vec input vector
/// \param[out] out_p_vec output vector
/// \param[in]  mpi_comm MPI communicator
///
/// \todo Add ability to optimize for sparse sends.  Put MPI_Request
///       objects into a resizing vecor to be skipped for 0 byte
///       send/recvs
template <typename T>
void mpi_all_to_all(std::vector< std::vector<T> >& in_p_vec,
                    std::vector< std::vector<T> >& out_p_vec,
                    MPI_Comm mpi_comm) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );
  assert( mpi_size == (int) in_p_vec.size() );

  std::vector<size_t> per_rank_send_counts(mpi_size);
  std::vector<size_t> per_rank_recv_counts(mpi_size);

  for(int i=0; i<mpi_size; ++i) {
    per_rank_send_counts[i] = in_p_vec[i].size();
  }

  CHK_MPI( MPI_Alltoall( &(per_rank_send_counts[0]), 1, mpi_typeof(size_t()),
                         &(per_rank_recv_counts[0]), 1, mpi_typeof(size_t()),
                         mpi_comm ) );
  //Allocate recv memory
  out_p_vec.resize(mpi_size);
  for(int i=0; i<mpi_size; ++i) {
    out_p_vec[i].resize(per_rank_recv_counts[i]);
  }

/*
 *  //Aggregisive method, good for small-med p?
 *  std::vector<MPI_Request> vec_req;
 *  for(int i=0; i<mpi_size; ++i) {
 *    int send_to_rank   = (mpi_rank + i) % mpi_size;
 *    int recv_from_rank = (mpi_rank - i + mpi_size) % mpi_size;
 *    //Post Irecvs
 *    if(out_p_vec[recv_from_rank].size() > 0) {
 *      MPI_Request req;
 *      CHK_MPI( MPI_Irecv( &(out_p_vec[recv_from_rank][0]),
 *                          (out_p_vec[recv_from_rank].size() * sizeof(T)),
 *                          MPI_BYTE,
 *                          recv_from_rank, 0, mpi_comm, &req ) );
 *      vec_req.push_back(req);
 *    }
 *
 *    //Post Isends
 *    if(in_p_vec[send_to_rank].size() > 0) {
 *      MPI_Request req;
 *      CHK_MPI( MPI_Isend( &(in_p_vec[send_to_rank][0]),
 *                          (in_p_vec[send_to_rank].size() * sizeof(T)),
 *                          MPI_BYTE,
 *                          send_to_rank, 0, mpi_comm, &req) );
 *      vec_req.push_back(req);
 *    }
 *  }
 *  CHK_MPI( MPI_Waitall(vec_req.size(), &(vec_req[0]), MPI_STATUSES_IGNORE) );
 */

  //Basic method -- good for large p?
  //For each rank, in parallel do:
  for(int i=0; i<mpi_size; ++i) {
    MPI_Request request[2];
    int send_to_rank   = (mpi_rank + i) % mpi_size;
    int recv_from_rank = (mpi_rank - i + mpi_size) % mpi_size;
    CHK_MPI( MPI_Isend( &(in_p_vec[send_to_rank][0]),
                        (in_p_vec[send_to_rank].size() * sizeof(T)),
                        MPI_BYTE,
                        send_to_rank, 0, mpi_comm, &(request[0]) ) );
    CHK_MPI( MPI_Irecv( &(out_p_vec[recv_from_rank][0]),
                        (out_p_vec[recv_from_rank].size() * sizeof(T)),
                        MPI_BYTE,
                        recv_from_rank, 0, mpi_comm, &(request[1]) ) );
    CHK_MPI( MPI_Waitall(2, request, MPI_STATUSES_IGNORE) );
  }
}

  using edge_type = std::pair<uint64_t, uint64_t>;
  using edges_vec_type = std::vector<edge_type>;
  using mpi_buf_vec_type = std::vector<edges_vec_type>;

int main (int argc, char** argv)
{
  int mpi_rank(0), mpi_size(0);
  double t_start;

  CHK_MPI(MPI_Init(&argc, &argv));
  {

  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
  CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

  if (mpi_rank == 0) {

    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
    // havoqgt::get_environment().print();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  mpi_buf_vec_type send_buf_vec(mpi_size);
  mpi_buf_vec_type recev_buf_vec(mpi_size);

  for (int i = 0; i < mpi_size; ++i) {
    send_buf_vec[i] = edges_vec_type();
    recev_buf_vec[i] = edges_vec_type();
  }

  std::string source_file_name("/p/lscratchf/iwabuchi/WebDataCommons/2012/edges/sorted/part-r-00");
  if (mpi_rank < 10) {
    source_file_name += "0";
  }
  if (mpi_rank < 100) {
    source_file_name += "0";
  }
  source_file_name += std::to_string(mpi_rank);
  std::cout << "rank: " << mpi_rank << " open " << source_file_name << std::endl;

  std::ifstream source_file;
  source_file.open(source_file_name);
  assert( source_file.is_open() );

   if (mpi_rank == 0) {
    t_start = MPI_Wtime();
    std::cout << "Loading data from sourcefiles...\n";
   }
  size_t count = 0;
  while (source_file.good()) {
    uint64_t src_vtx, dst_vtx;
    source_file >> src_vtx >> dst_vtx;
    uint64_t target_rank = rand() % mpi_size;
    send_buf_vec[target_rank].push_back(std::make_pair(src_vtx, dst_vtx));
    ++count;
  }
  source_file.close();
  std::cout << "rank: " << mpi_rank << " will sent " << count << " edges" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
   if (mpi_rank == 0) {
    std::cout << "done. " << MPI_Wtime() - t_start << " sec\n";
   }


  if (mpi_rank == 0) {
    t_start = MPI_Wtime();
    std::cout << "All to all start\n";
  }
  mpi_all_to_all(send_buf_vec, recev_buf_vec, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "All to all done:" << MPI_Wtime() - t_start << " sec\n";
  }


  std::string result_file_name("/p/lscratchf/iwabuchi/WebDataCommons/2012/edges/randomized/part-r-00");
  if (mpi_rank < 10) {
    result_file_name += "0";
  }
  if (mpi_rank < 100) {
    result_file_name += "0";
  }
  result_file_name += std::to_string(mpi_rank);
  std::cout << "rank: " << mpi_rank << " open " << result_file_name << std::endl;

  std::ofstream result_file;
  result_file.open(result_file_name);
  assert( result_file.is_open() );

  edges_vec_type randomized_edges_vec;
  for (auto it_src_proc = recev_buf_vec.begin(), it_src_proc_end = recev_buf_vec.end(); it_src_proc != it_src_proc_end; ++it_src_proc) {
    edges_vec_type edges = *it_src_proc;
    for (auto it_edges = edges.begin(), it_edges_end = edges.end(); it_edges != it_edges_end; ++it_edges) {
      randomized_edges_vec.push_back(*it_edges);
    }
  }
  std::random_shuffle(randomized_edges_vec.begin(), randomized_edges_vec.end());
  std::cout << "rank: " << mpi_rank << " received " << randomized_edges_vec.size() << " edges" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

   if (mpi_rank == 0) {
    t_start = MPI_Wtime();
    std::cout << "Storing data to files...\n";
   }
  for (auto it = randomized_edges_vec.begin(), it_end = randomized_edges_vec.end(); it != it_end; ++it) {
    result_file << it->first << "\t" << it->second << std::endl;
  }
  result_file.close();
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "done. " << MPI_Wtime() - t_start << " sec\n";
  }

  } //END Main MPI
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));

  CHK_MPI(MPI_Finalize());
  if (mpi_rank == 0) {
    std::cout << "FIN." << std::endl;
  }
  return 0;
}

