// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <deque>
#include <tuple>
#include <vector>
#include <ygm/mpi.hpp>

namespace ygm {

template <typename MSG>
class comm_exchanger {
  using count_pair = std::pair<size_t, size_t>;

 public:
  comm_exchanger(MPI_Comm comm, int tag) : m_tag(tag), m_local_count(0) {
    CHK_MPI(MPI_Comm_dup(comm, &m_comm));
    CHK_MPI(MPI_Comm_size(m_comm, &m_comm_size));
    CHK_MPI(MPI_Comm_rank(m_comm, &m_comm_rank));
    m_vec_send.resize(m_comm_size);
    init_recv_counts();
  }

  ~comm_exchanger() { cancel_recv_counts(); }

  void queue(int rank, const MSG &msg) {
    // if (m_local_count == 0) { init_recv_counts(); }
    ++m_local_count;
    // m_local_count += sizeof(MSG);
    m_vec_send[rank].push_back(msg);
  }

  std::size_t queue_bytes(int rank, const MSG &msg) {
    // if (m_local_count == 0) { init_recv_counts(); }
    ++m_local_count;
    // m_local_count += sizeof(MSG);
    m_vec_send[rank].push_back(msg);
    // return sizeof(MSG);
    return 1;
  }

  // RETURNS:  Total exchange count of communicator
  // WARNING:  Could optimzie not to send to self....   However, more general
  // this way....
  template <typename RecvHandlerFunc>
  uint64_t exchange(RecvHandlerFunc recv_func, uint64_t extracount = 0) {
    std::deque<std::tuple<MPI_Request, MSG *, size_t>>
                             q_req_irecv_data;  // req, buff, recv_size
    std::vector<MPI_Request> vec_req_isend_data;

    uint64_t to_return = 0;
    // Send counts
    for (int i = 0; i < m_comm_size; ++i) {
        int rank = (i + m_comm_rank)%m_comm_size;
      count_pair to_send;
      to_send.first = m_vec_send[rank].size();
      to_send.second =
          m_local_count + extracount;  // FIXME, should count this better
      CHK_MPI(MPI_Send(&to_send, sizeof(count_pair), MPI_BYTE, rank, 2 * m_tag,
                       m_comm));
    }

    // Big do/while loop adds a bit of async.  Recvs can start while sends are
    // still in progress.
    do {
      // Wait for all counts to come in, post recvs.
      /*while*/if(!m_req_irecv_counts.empty()) {
        MPI_Status status;
        auto       req_pair = m_req_irecv_counts.front();
        int        flag;
        CHK_MPI(MPI_Test(&(req_pair.first), &flag, &status));
        if (flag) {
          int    recv_counts_vec_pos = req_pair.second;
          int    recv_rank           = status.MPI_SOURCE;
          size_t recv_size =
              m_vec_recv_counts[recv_counts_vec_pos].first * sizeof(MSG);
          to_return += m_vec_recv_counts[recv_counts_vec_pos]
                           .second;  // add up global exc count
          if (recv_size > 0) {
            MSG *       buff = (MSG *)malloc(recv_size);
            if(buff == NULL) {
              std::cerr << "comm_exchanger:: unable to malloc" << std::endl << std::flush;  exit(-1);
            }
            MPI_Request req;
            CHK_MPI(MPI_Irecv((void *)buff, recv_size, MPI_BYTE, recv_rank,
                              2 * m_tag + 1, m_comm, &req));
            q_req_irecv_data.push_back(std::make_tuple(req, buff, recv_size));
          }
          /// If I received a rank's count, they should be close to accepting my
          /// messages
          size_t send_size = m_vec_send[recv_rank].size() * sizeof(MSG);
          if (send_size > 0) {
            MPI_Request req;
            CHK_MPI(MPI_Isend((void *)&(m_vec_send[recv_rank][0]), send_size,
                              MPI_BYTE, recv_rank, 2 * m_tag + 1, m_comm,
                              &req));
            vec_req_isend_data.push_back(req);
          }
          m_req_irecv_counts.pop_front();
        } else {
          //break;
        }
      }

      // Wait for all recvs to come in
      // WARNING:  this might be better as a if/then, posting the recvs is
      // actually higher priority...
      while (!q_req_irecv_data.empty()) {
        auto req_tuple = q_req_irecv_data.front();
        q_req_irecv_data.pop_front();
        int flag;
        CHK_MPI(MPI_Test(&(std::get<0>(req_tuple)), &flag, MPI_STATUS_IGNORE));
        if (flag) {
          MSG *  recvbuff  = std::get<1>(req_tuple);
          size_t recv_size = std::get<2>(req_tuple) / sizeof(MSG);
          for (size_t i = 0; i < recv_size; ++i) {
            recv_func(recvbuff[i]);
          }
          free(recvbuff);
        } else {
          q_req_irecv_data.push_back(req_tuple);
          break;
        }
      }
    } while (!m_req_irecv_counts.empty() || !q_req_irecv_data.empty());

    // Wait on all my sends
    if (!vec_req_isend_data.empty()) {
      CHK_MPI(MPI_Waitall(vec_req_isend_data.size(), &(vec_req_isend_data[0]),
                          MPI_STATUS_IGNORE));
    }

    // clear send queues
    for (size_t i = 0; i < m_vec_send.size(); ++i) {
      m_vec_send[i].clear();
    }
    m_local_count = 0;

    m_tag += 2;
    init_recv_counts();

    return to_return;
  }

 private:
  void init_recv_counts() {
    // Only inits once per exchange
    if (m_req_irecv_counts.size() != m_comm_size) {
      m_vec_recv_counts.clear();
      m_vec_recv_counts.resize(m_comm_size, std::make_pair(0, 0));
      for (int i = 0; i < m_comm_size; ++i) {
        MPI_Request req;
        CHK_MPI(MPI_Irecv((void *)&(m_vec_recv_counts[i]), sizeof(count_pair),
                          MPI_BYTE, MPI_ANY_SOURCE, 2 * m_tag, m_comm, &req));
        m_req_irecv_counts.push_back(std::make_pair(req, i));
      }
    }
  }

  void cancel_recv_counts() {
    for (int i = 0; i < m_req_irecv_counts.size(); ++i) {
      CHK_MPI(MPI_Cancel(&(m_req_irecv_counts[i].first)));
    }
  }

  // Basic Comm Data
  MPI_Comm m_comm;
  int      m_comm_rank;
  int      m_comm_size;
  int      m_tag;

  // send queue
  std::vector<std::vector<MSG>> m_vec_send;
  size_t                        m_local_count;

  // exchange data
  std::vector<count_pair> m_vec_recv_counts;
  std::deque<std::pair<MPI_Request, int>> m_req_irecv_counts;  // req, rank
};                                                             // namespace ygm
}  // namespace ygm