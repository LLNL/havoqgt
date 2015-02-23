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

#ifndef _HAVOQGT_MPI_HPP_
#define _HAVOQGT_MPI_HPP_

#include <havoqgt/detail/null_ostream.hpp>
#include <sched.h>
#include <vector>
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdint.h>


#define CHK_MPI(a) { if (a != MPI_SUCCESS) {\
                      char* error_string = NULL; \
                      int len = 0; \
                      MPI_Error_string(a, error_string, &len); \
                      std::cerr << __FILE__ << ", line " << __LINE__  \
                           <<" MPI ERROR = " << error_string << std::endl; \
                           exit(-1); \
                     } }


namespace havoqgt { namespace mpi {

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

template <typename T>
MPI_Op mpi_typeof(std::greater<T>) { return MPI_MAX; }

template <typename T>
MPI_Op mpi_typeof(std::less<T>) { return MPI_MIN; }

template <typename T>
MPI_Op mpi_typeof(std::plus<T>) { return MPI_SUM; }

template <typename T>
MPI_Op mpi_typeof(std::multiplies<T>) { return MPI_PROD; }

template <typename T>
MPI_Op mpi_typeof(std::logical_and<T>) { return MPI_LAND; }

template <typename T>
MPI_Op mpi_typeof(std::logical_or<T>) { return MPI_LOR; }

class mpi_communicator{
public:
  mpi_communicator(int argc, char** argv)
    : init(true) {
    check_mpi( MPI_Init(&argc, &argv) );
    m_comm = MPI_COMM_WORLD;
  }
  mpi_communicator(MPI_Comm comm)
    : m_comm(MPI_COMM_WORLD),
      init(false) { }

  ~mpi_communicator() {
    if(init) {
      check_mpi( MPI_Finalize() );
    }
  }

private:
  /// Not Copyable
  mpi_communicator(const mpi_communicator&);

  /// Checks MPI return codes
  void check_mpi(int ret) {
    if(ret != MPI_SUCCESS) {
      char estring[MPI_MAX_ERROR_STRING];
      int len(0);
      MPI_Error_string(ret, estring, &len);
      if(m_rank==0) {
        std::cerr << "MPI ERROR = " << estring << std::endl;
      }
      MPI_Abort(m_comm, -1);
    }
  }

  MPI_Comm m_comm;
  uint32_t m_size;
  uint32_t m_rank;
  bool init;
};

void mpi_yield_barrier(MPI_Comm mpi_comm) {
  MPI_Request request;
  CHK_MPI(MPI_Ibarrier(mpi_comm, &request));

  for(; ;) {
    MPI_Status status;
    int is_done = false;
    MPI_Test(&request, &is_done, &status);
    if (is_done) {
      return;
    } else {
      sched_yield();
    }
  }
}


//no std:: equivalent for MPI_BAND, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MAXLOC, MPI_MINLOC

template <typename T, typename Op>
T mpi_all_reduce(T in_d, Op in_op, MPI_Comm mpi_comm) {
  T to_return;
  CHK_MPI(
    MPI_Allreduce( &in_d, &to_return, 1, mpi_typeof(T()), mpi_typeof(in_op),
      mpi_comm)
  );
  return to_return;
}

template <typename Vec, typename Op>
void mpi_all_reduce_inplace(Vec &vec, Op in_op, MPI_Comm mpi_comm) {
  CHK_MPI(
    MPI_Allreduce(MPI_IN_PLACE, &(vec[0]), vec.size(),
      mpi_typeof(typename Vec::value_type()), mpi_typeof(in_op), mpi_comm)
  );
}

template <typename T, typename Op>
void mpi_all_reduce(std::vector<T>& in_vec, std::vector<T>& out_vec, Op in_op,
    MPI_Comm mpi_comm) {
  out_vec.resize(in_vec.size());
  CHK_MPI(
    MPI_Allreduce( &(in_vec[0]), &(out_vec[0]), in_vec.size(),
      mpi_typeof(in_vec[0]), mpi_typeof(in_op), mpi_comm)
  );
}


template<typename T, typename Partitioner>
class owner_sort {
public:
  owner_sort(Partitioner _owner):owner(_owner) { }

  inline bool operator()(const T& a, const T& b) const {
    return owner(a) < owner(b);
  }

private:
  Partitioner owner;
};

template <typename T>
void mpi_all_to_all(std::vector<T>& in_vec, std::vector<int>& in_sendcnts,
                    std::vector<T>& out_vec, std::vector<int>& out_recvcnts,
                    MPI_Comm mpi_comm) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );

  void* send_vec  =  in_vec.size() > 0 ? (void*) &in_vec[0] : NULL;
  int* send_cnts =  in_sendcnts.size() > 0 ? &in_sendcnts[0] : NULL;

  //
  // Exchange send counts
  out_recvcnts.resize(mpi_size);
  CHK_MPI( MPI_Alltoall( (void*) send_cnts, sizeof(int), MPI_BYTE,
                         (void*) &(out_recvcnts[0]), sizeof(int), MPI_BYTE,
                         mpi_comm ) );

  //
  // Calc send & recv disps
  std::vector<int> sdispls(mpi_size,0), rdispls(mpi_size,0);
  std::partial_sum(in_sendcnts.begin(), in_sendcnts.end(), sdispls.begin());
  for(size_t i=0; i<sdispls.size(); ++i) {
    sdispls[i] -= in_sendcnts[i]; //set to 0 offset
  }
  std::partial_sum(out_recvcnts.begin(), out_recvcnts.end(), rdispls.begin());
  for(size_t i=0; i<rdispls.size(); ++i) {
    rdispls[i] -= out_recvcnts[i]; //set to 0 offset
  }

  out_vec.resize(std::accumulate(out_recvcnts.begin(), out_recvcnts.end(), 0));

  int* send_displs =  sdispls.size() > 0 ? &sdispls[0] : NULL;
  CHK_MPI( MPI_Alltoallv(send_vec, send_cnts, send_displs,
                  mpi_typeof(T()), &(out_vec[0]), &(out_recvcnts[0]),
                  &(rdispls[0]), mpi_typeof(T()), mpi_comm) );

}


template <typename T, typename Partitioner>
void mpi_all_to_all(std::vector<T>& inout_vec, std::vector<T>& temp_vec,
      Partitioner &owner, MPI_Comm mpi_comm) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );

  std::vector<int> send_counts(mpi_size,0);
  std::vector<int> send_disps(mpi_size,0);

  std::vector<int> recv_counts(mpi_size,0);
  std::vector<int> recv_disps(mpi_size,0);

  //sort send vector by owner
  //std::sort(inout_vec.begin(), inout_vec.end(), owner_sort<T,Partitioner>(owner));

  //calc send counts
  for(size_t i=0; i<inout_vec.size(); ++i) {
    send_counts[owner(inout_vec[i])] += sizeof(T); //sizeof(t) lets us use PODS
    //if(i>0) {
    //  assert(owner(inout_vec[i-1]) <= owner(inout_vec[i]));
    //}
  }

  //cacl send disps
  std::partial_sum(send_counts.begin(), send_counts.end(), send_disps.begin());
  for(size_t i=0; i<send_disps.size(); ++i) {
    send_disps[i] -= send_counts[i]; //set to 0 offset
  }

  temp_vec.resize(inout_vec.size());
  for(size_t i=0; i<inout_vec.size(); ++i) {
    std::vector<int> temp_arrange(mpi_size,0);
    int dest_rank = owner(inout_vec[i]);
    size_t dest_offset = send_disps[dest_rank] + temp_arrange[dest_rank];
    temp_arrange[dest_rank] += sizeof(T);
    dest_offset /= sizeof(T);
    temp_vec[dest_offset] = inout_vec[i];
  }


  //exchange send counts
  CHK_MPI( MPI_Alltoall( (void*) &(send_counts[0]), sizeof(int), MPI_BYTE,
                         (void*) &(recv_counts[0]), sizeof(int), MPI_BYTE,
                         mpi_comm ) );

  //Allocate recv vector
  int total_recv = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
  inout_vec.resize(total_recv/sizeof(T));

  //cacl recv disps
  std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_disps.begin());
  for(size_t i=0; i<recv_disps.size(); ++i) {
    recv_disps[i] -= recv_counts[i]; //set to 0 offset
  }

  //perform actual a3lloallv
  CHK_MPI( MPI_Alltoallv( (void*) &(temp_vec[0]), &(send_counts[0]), &(send_disps[0]), MPI_BYTE,
                          (void*) &(inout_vec[0]),&(recv_counts[0]), &(recv_disps[0]), MPI_BYTE,
                          mpi_comm) );
}

template <typename T, typename Partitioner>
void mpi_all_to_all_better(std::vector<T>& in_vec, std::vector<T>& out_vec,
      Partitioner &owner, MPI_Comm mpi_comm) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );

  std::vector<int> send_counts(mpi_size,0);
  std::vector<int> send_disps(mpi_size,0);
  std::vector<int> recv_counts(mpi_size,0);
  std::vector<int> recv_disps(mpi_size,0);

  //sort send vector by owner
  //std::sort(in_vec.begin(), in_vec.end(), owner_sort<T,Partitioner>(owner));

  //calc send counts
  for(size_t i=0; i<in_vec.size(); ++i) {
    send_counts[owner(in_vec[i], true)] += sizeof(T); //sizeof(t) lets us use PODS
    //if(i>0) {
    //  assert(owner(inout_vec[i-1]) <= owner(inout_vec[i]));
    //
  }


  //cacl send disps
  std::partial_sum(send_counts.begin(), send_counts.end(), send_disps.begin());
  for(size_t i=0; i<send_disps.size(); ++i) {
    send_disps[i] -= send_counts[i]; //set to 0 offset
  }


  { //rearrange instead of sorting
    std::vector<T> order_vec(in_vec.size());
    std::vector<int> temp_arrange(mpi_size,0);
    for(size_t i=0; i<in_vec.size(); ++i) {
      int dest_rank = owner(in_vec[i], false);
      assert (dest_rank >=0 && dest_rank < mpi_size);

      size_t dest_offset = send_disps[dest_rank] + temp_arrange[dest_rank];
      temp_arrange[dest_rank] += sizeof(T);
      dest_offset /= sizeof(T);
      order_vec[dest_offset] = in_vec[i];
    }
    in_vec.swap(order_vec);
  }


  //exchange send counts
  CHK_MPI( MPI_Alltoall( (void*) &(send_counts[0]), sizeof(int), MPI_BYTE,
                         (void*) &(recv_counts[0]), sizeof(int), MPI_BYTE,
                         mpi_comm ) );


  //Allocate recv vector
  int total_recv = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
  out_vec.resize(total_recv/sizeof(T));

  //cacl recv disps
  std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_disps.begin());
  for(size_t i=0; i<recv_disps.size(); ++i) {
    recv_disps[i] -= recv_counts[i]; //set to 0 offset
  }

  //perform actual alltoallv
  void* send_ptr = in_vec.empty()  ? NULL :(void*) &(in_vec[0]);
  void* recv_ptr = out_vec.empty() ? NULL :(void*) &(out_vec[0]);
  CHK_MPI( MPI_Alltoallv( send_ptr, &(send_counts[0]), &(send_disps[0]), MPI_BYTE,
                          recv_ptr, &(recv_counts[0]), &(recv_disps[0]), MPI_BYTE,
                          mpi_comm) );
}


//This code is very similar to that of:
//http://trac.mcs.anl.gov/projects/mpich2/browser/mpich2/trunk/src/mpi/coll/alltoall.c
//
// It is re-implemted here because <cough> BG/P doesn't support mpi-2.2 which
// added MPI_IN_PLACE for Alltoall.
template <typename T>
void mpi_all_to_all_in_place(std::vector<T>& in_out_vec, size_t count,
                             MPI_Comm mpi_comm) {
  const int HAVOQGT_TAG=42;
  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );
  assert( in_out_vec.size() == count * mpi_size );

  for(int i=0; i<mpi_size; ++i) {
    //start inner loop at i to avoid re-exchanging data
    for(int j=i; j<mpi_size; ++j) {
      if( mpi_rank == i) {
        MPI_Status status;
        CHK_MPI( MPI_Sendrecv_replace(&(in_out_vec[j*count]), count,
                                      mpi_typeof(T()),
                                      j, HAVOQGT_TAG,
                                      j, HAVOQGT_TAG,
                                      mpi_comm, &status) );

      } else if(mpi_rank == j) {
        MPI_Status status;
        CHK_MPI( MPI_Sendrecv_replace(&(in_out_vec[i*count]), count,
                                      mpi_typeof(T()),
                                      i, HAVOQGT_TAG,
                                      i, HAVOQGT_TAG,
                                      mpi_comm, &status) );


      }
    }
  }
}

/// TODO:  Add tests
template <typename T>
void mpi_all_gather(T _t, std::vector<T>& out_p_vec, MPI_Comm mpi_comm) {
  int mpi_size(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );

  out_p_vec.resize(mpi_size);
  CHK_MPI( MPI_Allgather(&_t, sizeof(_t), MPI_BYTE,
                         &(out_p_vec[0]),  sizeof(_t), MPI_BYTE, mpi_comm) );
}

/// TODO:  Add tests, especially with non mpi types, POD
template <typename T>
void mpi_all_gather(std::vector<T>& in_send,
                    std::vector<T>& out_recv_gather, MPI_Comm mpi_comm) {

  int mpi_size(0), mpi_rank(0);
  CHK_MPI( MPI_Comm_size( mpi_comm, &mpi_size) );
  CHK_MPI( MPI_Comm_rank( mpi_comm, &mpi_rank) );

  int my_size = in_send.size();
  std::vector<int> recv_counts(mpi_size,0);
  std::vector<int> recv_disps(mpi_size,0);

  //Gather recv counts for all ranks
  mpi_all_gather(my_size, recv_counts, mpi_comm);

  //Allocate recv vector
  int total_recv = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
  if(total_recv > 0) {
    out_recv_gather.resize(total_recv);

    //cacl recv disps
    std::partial_sum(recv_counts.begin(), recv_counts.end(), recv_disps.begin());
    for(size_t i=0; i<recv_disps.size(); ++i) {
      recv_disps[i] -= recv_counts[i]; //set to 0 offset
    }

    void* send_buff = in_send.size() == 0 ? NULL : &(in_send[0]);
    CHK_MPI( MPI_Allgatherv( send_buff, my_size, mpi_typeof(T()),
                             &(out_recv_gather[0]), &(recv_counts[0]), &(recv_disps[0]),
                             mpi_typeof(T()), mpi_comm));
  }
}


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

template <typename T>
void mpi_bcast(T& data, int root, MPI_Comm comm)
{
  CHK_MPI( MPI_Bcast(&data, sizeof(data), MPI_BYTE, root, comm));
}


inline std::ostream& cout_rank0() {
  int mpi_rank(0);
  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );

  if(mpi_rank == 0) {
    return std::cout;
  }
  return havoqgt::detail::get_null_ostream();
}

inline std::ostream& cout_rank0_barrier() {
  CHK_MPI( MPI_Barrier( MPI_COMM_WORLD ));
  int mpi_rank(0);
  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );

  if(mpi_rank == 0) {
    return std::cout;
  }
  return havoqgt::detail::get_null_ostream();
}


inline int mpi_comm_rank() {
  int rank;
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
  return rank;
}

inline int mpi_comm_size() {
  int size;
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &size) );
  return size;
}


} //end namespace mpi

///
/// Start "new" communicator here
///

class communicator {
public:
  communicator(MPI_Comm in_comm)
    : m_comm( in_comm )
  {
    CHK_MPI( MPI_Comm_rank(m_comm, &m_rank) );
    CHK_MPI( MPI_Comm_size(m_comm, &m_size) );
  }
  
  communicator() {}
  
  MPI_Comm comm() const { return m_comm; }
  int      size() const { return m_size; }
  int      rank() const { return m_rank; }
  
  void barrier() const {
    MPI_Barrier(m_comm);
  }
private:
  MPI_Comm m_comm;
  int      m_size;
  int      m_rank;
    
};


} // end namespace havoqgt 

#endif //_HAVOQGT_MPI_HPP_

