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
  enum tag { partial = 0, complete = 1, empty = 2 };

 public:
  comm_exchanger(MPI_Comm comm, size_t send_size) : m_send_size(send_size) {
    CHK_MPI(MPI_Comm_dup(comm, &m_comm));
    CHK_MPI(MPI_Comm_dup(comm, &m_comm_next));
    CHK_MPI(MPI_Comm_size(m_comm, &m_comm_size));
    CHK_MPI(MPI_Comm_rank(m_comm, &m_comm_rank));
    // if (m_comm_rank == 0) {
    //   std::cout << "m_send_size = " << m_send_size << std::endl;
    // }
    m_vec_send.resize(m_comm_size);
    for (auto &vecs : m_vec_send) { vecs.reserve(m_send_size); }
    init_recvs();
    create_send_count_order();
  }

  ~comm_exchanger() {
    cancel_recvs();
    if (!m_req_irecv_buff.empty() || !m_req_isend_buff.empty()) {
      std::cerr << "ERROR ~comm_exchanger shutdown" << std::endl;
    }
    CHK_MPI(MPI_Comm_free(&m_comm));
    CHK_MPI(MPI_Comm_free(&m_comm_next));
  }

  void queue(int rank, const MSG &msg) {
    m_local_empty = false;
    m_vec_send[rank].push_back(msg);
    if (m_vec_send[rank].size() == m_send_size) {
      // std::cout << whoami() << " -- sending a partial chunk" << std::endl;
      post_isend(rank, tag::partial);
    }
  }

  std::size_t queue_bytes(int rank, const MSG &msg) {
    queue(rank, msg);
    return 1;
  }

  void post_isend(int rank, tag send_tag) {
    std::vector<MSG> *data = allocate_new_send_buffer();
    data->swap(m_vec_send[rank]);
    MPI_Request req;

    void *send_ptr = data->empty() ? NULL : data->data();
    CHK_MPI(MPI_Isend(send_ptr, data->size() * sizeof(MSG), MPI_BYTE, rank,
                      send_tag, m_comm, &req));
    m_req_isend_buff.push_back(std::make_pair(req, data));
  }

  std::vector<MSG> *allocate_new_send_buffer() {
    if (!m_req_isend_buff.empty()) {
      MPI_Request req = m_req_isend_buff.front().first;
      std::vector<MSG> *pvec = m_req_isend_buff.front().second;
      int flag;
      CHK_MPI(MPI_Test(&req, &flag, MPI_STATUS_IGNORE));
      if (flag) {  // send has finished, can recycle buffer
        m_req_isend_buff.pop_front();
        pvec->clear();
        return pvec;
      }
    }
    // If here, recycling old send failed.
    std::vector<MSG> *to_return = new std::vector<MSG>();
    to_return->reserve(m_send_size);
    return to_return;
  }

  // RETURNS:  Total exchange count of communicator
  // WARNING:  Could optimzie not to send to self....   However, more general
  // this way....
  template <typename RecvHandlerFunc>
  uint64_t exchange(RecvHandlerFunc recv_func, uint64_t extracount = 0) {
    bool could_be_empty = m_local_empty && (extracount == 0);
    tag send_tag;
    if (could_be_empty) {
      send_tag = tag::empty;
    } else {
      send_tag = tag::complete;
    }

    // Post all my isends
    for (int i = 0; i < m_comm_size; ++i) {
      post_isend(m_send_count_order[i], send_tag);
    }

    size_t count_not_empty(0);
    while (!m_req_irecv_buff.empty()) {
      MPI_Request req = m_req_irecv_buff.front().first;
      void *buffer = m_req_irecv_buff.front().second;
      MPI_Status status;
      CHK_MPI(MPI_Wait(&req, &status));
      int count;
      CHK_MPI(MPI_Get_count(&status, MPI_BYTE, &count));
      int recvtag = status.MPI_TAG;  // this might not be portable!

      // recv msgs
      MSG *recvbuff = (MSG *)buffer;
      size_t recv_size = count / sizeof(MSG);
      for (size_t i = 0; i < recv_size; ++i) { recv_func(recvbuff[i]); }

      if (recvtag == tag::partial) {
        // std::cout << whoami() << " -- received a partial chunk" <<
        // std::endl;
        post_irecv(buffer);
      } else if (recvtag == tag::complete) {
        count_not_empty++;
        free(buffer);
      } else if (recvtag == tag::empty) {
        free(buffer);
      } else {
        std::cerr << "Unknown tag!" << std::endl;
        exit(-1);
      }
      m_req_irecv_buff.pop_front();
    }

    // Wait on all my sends
    while (!m_req_isend_buff.empty()) {
      MPI_Request req = m_req_isend_buff.front().first;
      std::vector<MSG> *pvec = m_req_isend_buff.front().second;
      CHK_MPI(MPI_Wait(&req, MPI_STATUS_IGNORE));
      delete pvec;
      m_req_isend_buff.pop_front();
    }

    // clear send queues
    m_local_empty = true;

    init_recvs();
    // std::cout << "count_not_empty = " << count_not_empty << std::endl;
    return count_not_empty;
  }

 private:
  void init_recvs() {
    std::swap(m_comm, m_comm_next);
    assert(m_req_irecv_buff.empty());
    for (int i = 0; i < m_comm_size; ++i) { post_irecv(); }
  }

  void post_irecv(void *buffer = NULL) {
    if (buffer == NULL) {
      buffer = malloc(m_send_size * sizeof(MSG));
      if (buffer == NULL) {
        std::cerr << "comm_exchanger:  unable to malloc" << std::endl;
        exit(-1);
      }
    }
    MPI_Request req;
    CHK_MPI(MPI_Irecv(buffer, m_send_size * sizeof(MSG), MPI_BYTE,
                      MPI_ANY_SOURCE, MPI_ANY_TAG, m_comm, &req));
    m_req_irecv_buff.push_back(std::make_pair(req, buffer));
  }

  void cancel_recvs() {
    for (int i = 0; i < m_req_irecv_buff.size(); ++i) {
      CHK_MPI(MPI_Cancel(&(m_req_irecv_buff[i].first)));
      free(m_req_irecv_buff[i].second);
    }
    m_req_irecv_buff.clear();
  }

  void create_send_count_order() {
    std::mt19937 gen(m_comm_rank);

    m_send_count_order.resize(m_comm_size);
    for (int i = 0; i < m_comm_size; ++i) { m_send_count_order[i] = i; }

    std::shuffle(m_send_count_order.begin(), m_send_count_order.end(), gen);
  }

  // Basic Comm Data
  MPI_Comm m_comm;
  MPI_Comm m_comm_next;
  int m_comm_rank;
  int m_comm_size;
  size_t m_send_size;

  // count send order
  std::vector<int> m_send_count_order;

  // send queue
  std::vector<std::vector<MSG>> m_vec_send;
  std::deque<std::pair<MPI_Request, std::vector<MSG> *>> m_req_isend_buff;

  bool m_local_empty = true;

  // recv data
  std::deque<std::pair<MPI_Request, void *>> m_req_irecv_buff;  // req, rank
};                                                              // namespace ygm
}  // namespace ygm
