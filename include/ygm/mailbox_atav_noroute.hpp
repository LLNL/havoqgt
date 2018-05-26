#pragma once 
#include <ygm/mailbox_base.hpp>

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <ygm/mpi.hpp>
#include <numeric>
#include <vector>

using std::vector;

template <typename Data, typename RecvHandlerFunc>
class mailbox_atav_noroute : private mailbox_base<Data, RecvHandlerFunc> {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t dest : 30;
    Data     data;
  };//__attribute__((packed));

 public:
  mailbox_atav_noroute(RecvHandlerFunc recv_func, size_t batch_size)
      : mailbox_base<Data, RecvHandlerFunc>(recv_func, batch_size) {
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &m_local_size));
    m_send_buf.reserve(batch_size);
    m_to_send.resize(batch_size);
    m_recv_vec.resize(batch_size);
    m_send_counts.resize(m_mpi_size);
    m_send_disps.resize(m_mpi_size);
    m_recv_counts.resize(m_mpi_size);
    m_recv_disps.resize(m_mpi_size);
  }

  void send(uint32_t dest, Data data) {
    //std::cout << "Mailbox -- send()" << std::endl;
    start_timer();
    assert(dest < m_mpi_size);
    if (dest == m_mpi_rank) {
      m_recv_func(false, data);
      end_timer(m_time_ingest);
    } else {
      push_to_buf(m_send_buf, m_send_counts, m_send_count, dest,
                  message{0, 0, dest, data});
      // m_send_counts[dest]++;
      // m_send_buf.push_back(message{0, 0, dest, data});
      end_timer(m_time_ingest);
      if (m_send_count >= m_batch_size) {
        do_exchange();
      }
    }
  }

  void send_bcast(Data data) {
    //std::cout << "Mailbox -- send_bcast()" << std::endl;
    for (uint32_t i = 0; i < m_mpi_size; i++) {
      if (i == m_mpi_rank) continue;
      start_timer();
      push_to_buf(m_send_buf, m_send_counts, m_send_count, i,
                  message{1, 0, i, data});
      // m_send_counts[i]++;
      // m_send_buf.push_back(message{1, 0, i, data});
      end_timer(m_time_ingest);
      if (m_send_count >= m_batch_size) {
        do_exchange();
      }
    }
  }

  using mailbox_base<Data, RecvHandlerFunc>::global_empty;

  void print_status(std::string p) {
    std::string path = p;
    path += "noroute";
    path += "_N";
    path += std::to_string(m_mpi_size / m_local_size);
    path += "_";
    path += std::to_string(m_mpi_rank);
    path += ".csv";
    std::ofstream outf(path);
    outf << m_time_ingest << ", " << m_time_cpu << ", " << m_time_wait << ", "
         << m_time_comm << ", " << m_total_sent << ", " << m_total_recv
         << std::endl;
  }

 private:
  void cleanup() {
    start_timer();
    std::fill(m_send_counts.begin(), m_send_counts.end(), 0);
    m_send_buf.clear();
    m_send_count = 0;
    end_timer(m_time_cpu);
  }

  bool do_exchange() {
    ++m_count_exchanges;
    uint64_t local_count = m_send_count;
    uint64_t global_count(0);
    start_timer();
    CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    end_timer(m_time_wait);
    start_timer();
    CHK_MPI(MPI_Allreduce(&local_count, &global_count, 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
    end_timer(m_time_comm);

    if (global_count == 0) {
      assert(m_send_count == 0);
      return true;
    }
    start_timer();

    calc_disps(m_send_counts, m_send_disps);

    //
    // Sort send buf implicitly
    std::vector<int> send_counts(m_mpi_size, 0);
    for (const message &msg : m_send_buf) {
      write_to_vec(m_to_send, m_send_disps, send_counts, msg.dest, msg);
      // int addr = m_send_disps[msg.dest] + send_counts[msg.dest];
      // m_to_send[addr] = msg;
      // send_counts[msg.dest]++;
    }

    //
    // Scale up send counters/offsets
    // for (size_t i = 0; i < m_send_counts.size(); ++i) {
    //   m_send_counts[i] *= sizeof(message);
    //   m_send_disps[i] *= sizeof(message);
    // }
    scale_addresses(m_send_counts, m_send_disps, sizeof(message));

    end_timer(m_time_cpu);
    CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    end_timer(m_time_wait);
    start_timer();
    //
    // Exchange send/recv counts
    CHK_MPI(MPI_Alltoall((void *)&(m_send_counts[0]), sizeof(int), MPI_BYTE,
                         (void *)&(m_recv_counts[0]), sizeof(int), MPI_BYTE,
                         MPI_COMM_WORLD));
    end_timer(m_time_comm);
    start_timer();

    //
    // Allocate recv vector
    // int total_recv = std::accumulate(m_recv_counts.begin(),
    //                                  m_recv_counts.end(), 0);
    // std::vector<message> m_recv_vec(total_recv / sizeof(message));
    // m_max_alloc = std::max(m_max_alloc, uint64_t(m_recv_vec.size()));
    alloc_recv_vec(m_recv_vec, m_recv_counts, sizeof(message));
    calc_disps(m_recv_counts, m_recv_disps);

    // perform actual alltoallv
    void *send_ptr = m_to_send.empty() ? NULL : (void *)&(m_to_send[0]);
    void *recv_ptr = m_recv_vec.empty() ? NULL : (void *)&(m_recv_vec[0]);
    end_timer(m_time_cpu);
    CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    end_timer(m_time_wait);
    start_timer();
    CHK_MPI(MPI_Alltoallv(send_ptr, &(m_send_counts[0]), &(m_send_disps[0]),
                          MPI_BYTE, recv_ptr, &(m_recv_counts[0]),
                          &(m_recv_disps[0]), MPI_BYTE, MPI_COMM_WORLD));
    end_timer(m_time_comm);
    start_timer();

    //
    // Recv messages
    for (const message &msg : m_recv_vec) {
      m_recv_func(msg.bcast, msg.data);
    }
    m_total_sent += m_send_count;
    m_total_recv += m_recv_vec.size();
    end_timer(m_time_cpu);
    cleanup();
    return false;
  }

  using mailbox_base<Data, RecvHandlerFunc>::calc_disps;
  using mailbox_base<Data, RecvHandlerFunc>::start_timer;
  using mailbox_base<Data, RecvHandlerFunc>::end_timer;
  using mailbox_base<Data, RecvHandlerFunc>::push_to_buf;
  using mailbox_base<Data, RecvHandlerFunc>::scale_addresses;
  using mailbox_base<Data, RecvHandlerFunc>::write_to_vec;
  using mailbox_base<Data, RecvHandlerFunc>::alloc_recv_vec;

 private:
  using mailbox_base<Data, RecvHandlerFunc>::m_recv_func;
  using mailbox_base<Data, RecvHandlerFunc>::m_batch_size;
  using mailbox_base<Data, RecvHandlerFunc>::m_max_alloc;
  using mailbox_base<Data, RecvHandlerFunc>::m_mpi_size;
  using mailbox_base<Data, RecvHandlerFunc>::m_mpi_rank;
  using mailbox_base<Data, RecvHandlerFunc>::m_count_exchanges;

  vector<message>                              m_send_buf;
  vector<message, default_init_alloc<message>> m_to_send;
  vector<message, default_init_alloc<message>> m_recv_vec;
  vector<int>                                  m_send_counts;
  vector<int>                                  m_send_disps;
  vector<int>                                  m_recv_counts;
  vector<int>                                  m_recv_disps;
  size_t                                       m_send_count = 0;
  int                                          m_local_size;
  double                                       m_time_ingest = 0.0f;
  double                                       m_time_comm = 0.0f;
  double                                       m_time_wait = 0.0f;
  double                                       m_time_cpu = 0.0f;
  uint32_t                                     m_total_sent = 0;
  uint32_t                                     m_total_recv = 0;

 public:
  value_measurements get_measurements() {
    value_measurements vm;
    vm.ingest = m_time_ingest;
    vm.cpu = m_time_cpu;
    vm.comm = m_time_comm;
    vm.wait = m_time_wait;
    return vm;
  }
};