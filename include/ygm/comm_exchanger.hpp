#pragma once

#include <deque>
#include <ygm/mpi.hpp>
#include <tuple>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>

namespace ygm {

template <typename MSG>
class atav_comm_exchanger {
 public:
  atav_comm_exchanger(MPI_Comm comm, uint64_t dummy) {
    CHK_MPI(MPI_Comm_dup(comm, &m_comm));
    CHK_MPI(MPI_Comm_size(m_comm, &m_comm_size));
    CHK_MPI(MPI_Comm_rank(m_comm, &m_comm_rank));
    // std::cout << "m_comm_size = " << m_comm_size << std::endl;
    m_send_dest_counts.resize(m_comm_size, 0);
  }

  ~atav_comm_exchanger() { CHK_MPI(MPI_Comm_free(&m_comm)); }

  void queue(int rank, const MSG &msg) {
    m_send_queue.push_back(msg);
    m_send_queue_dests.push_back(rank);
    m_send_dest_counts[rank] += sizeof(MSG);
  }

  std::size_t queue_bytes(int rank, const MSG &msg) {
    queue(rank, msg);
    return 1;
  }

  template <typename RecvHandlerFunc>
  uint64_t exchange(RecvHandlerFunc recv_func, uint64_t send_cnt = 0) {
    uint64_t local_count = m_send_queue.size();
    uint64_t comm_count(0);
    CHK_MPI(MPI_Allreduce((void *)&local_count, (void *)&comm_count, 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_SUM, m_comm));
    // if (m_comm_rank == 0) {
    //   std::cout << "comm_count = " << comm_count << std::endl;
    // }
    if (comm_count > 0) {
      std::vector<int> send_disps(m_comm_size, 0);
      std::vector<int> recv_counts(m_comm_size, 0);
      std::vector<int> recv_disps(m_comm_size, 0);

      // cacl send disps
      std::partial_sum(m_send_dest_counts.begin(), m_send_dest_counts.end(),
                       send_disps.begin());
      for (size_t i = 0; i < send_disps.size(); ++i) {
        send_disps[i] -= m_send_dest_counts[i];  // set to 0 offset
      }

      // std::vector<MSG> send_vec(m_send_queue.size());
      MSG *send_ptr = NULL;
      if (m_send_queue.size() > 0) {
        send_ptr = (MSG *)malloc(m_send_queue.size() * sizeof(MSG));
        if (send_ptr == NULL) {
          std::cerr << "Unable to malloc send_ptr" << std::endl;
          exit(-1);
        }
      }
      {  // rearrange instead of sorting
        std::vector<int> temp_arrange(m_comm_size, 0);
        for (size_t i = 0; i < m_send_queue.size(); ++i) {
          int dest_rank = m_send_queue_dests[i];
          assert(dest_rank >= 0 && dest_rank < m_comm_size);
          size_t dest_offset = send_disps[dest_rank] + temp_arrange[dest_rank];
          temp_arrange[dest_rank] += sizeof(MSG);
          dest_offset /= sizeof(MSG);
          send_ptr[dest_offset] = m_send_queue[i];
        }
      }

      // exchange send counts
      CHK_MPI(MPI_Alltoall(m_send_dest_counts.data(), sizeof(int), MPI_BYTE,
                           recv_counts.data(), sizeof(int), MPI_BYTE, m_comm));

      // Allocate recv vector
      int total_recv =
          std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
      // std::vector<MSG> recv_vec(total_recv / sizeof(MSG));
      MSG *recv_ptr = NULL;
      if (total_recv > 0) {
        recv_ptr = (MSG *)malloc(total_recv);
        if (recv_ptr == NULL) {
          std::cerr << "Unable to malloc recv_ptr" << std::endl;
          exit(-1);
        }
      }
      // std::cout << "recv_vec.size() = " << recv_vec.size() << std::endl;

      // cacl recv disps
      std::partial_sum(recv_counts.begin(), recv_counts.end(),
                       recv_disps.begin());
      for (size_t i = 0; i < recv_disps.size(); ++i) {
        recv_disps[i] -= recv_counts[i];  // set to 0 offset
      }

      // perform actual alltoallv
      // void *send_ptr = send_vec.empty() ? NULL : send_vec.data();
      // void *recv_ptr = recv_vec.empty() ? NULL : recv_vec.data();
      CHK_MPI(MPI_Alltoallv(
          send_ptr, &(m_send_dest_counts[0]), &(send_disps[0]), MPI_BYTE,
          recv_ptr, &(recv_counts[0]), &(recv_disps[0]), MPI_BYTE, m_comm));
      if (send_ptr) free(send_ptr);

      for (size_t i = 0; i < total_recv / sizeof(MSG); ++i) {
        recv_func(recv_ptr[i]);
      }
      if (recv_ptr) free(recv_ptr);

      // reset everything
      m_send_queue.clear();
      m_send_queue_dests.clear();
      for (int &c : m_send_dest_counts) { c = 0; }
    }
    uint64_t comm_mailbox_count;
    CHK_MPI(MPI_Allreduce((void *)&send_cnt, (void *)&comm_mailbox_count, 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_SUM, m_comm));
    return comm_mailbox_count;
  }

 private:
  MPI_Comm m_comm;
  int m_comm_size;
  int m_comm_rank;
  std::vector<MSG> m_send_queue;
  std::vector<int> m_send_queue_dests;
  std::vector<int> m_send_dest_counts;
};

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
};


}  // namespace ygm
