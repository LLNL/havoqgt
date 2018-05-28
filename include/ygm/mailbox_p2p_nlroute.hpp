#include <comm_exchanger.hpp>
#include <mpi.hpp>

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

using std::vector;

template <typename Data, typename RecvHandlerFunc>
class mailbox_p2p_nlroute {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t local : 6;
    uint32_t node : 24;  // Supports addressing <= 16777216 nodes w/ <= 64 cores
    Data     data;
  } __attribute__((packed));

 public:
  mailbox_p2p_nlroute(RecvHandlerFunc recv_func, size_t batch_size,
                      MPI_Comm local_comm, MPI_Comm remote_comm)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        m_local_comm(local_comm),
        m_remote_comm(remote_comm),
        m_local_exchanger(local_comm, 1),
        m_remote_exchanger(remote_comm, 2) {  // should change tag eventually
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size));
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
    CHK_MPI(MPI_Comm_size(m_local_comm, &m_local_size));
    CHK_MPI(MPI_Comm_rank(m_local_comm, &m_local_rank));
    CHK_MPI(MPI_Comm_size(m_remote_comm, &m_remote_size));
    CHK_MPI(MPI_Comm_rank(m_remote_comm, &m_remote_rank));
  }

  void send(uint32_t dest, Data data) {
    if (dest == m_mpi_rank) {
      m_recv_func(false, data);
    } else {
      uint32_t local = dest % m_local_size;
      uint32_t node = dest / m_local_size;
      if (local == m_local_rank) {
        m_remote_exchanger.queue(node, message{0, 0, local, node, data});
      } else {
        m_local_exchanger.queue(local, message{0, 0, local, node, data});
      }
      if (++m_send_count >= m_batch_size) {
        do_exchange();
      }
    }
  }

  void send_bcast(Data data) {
    for (uint32_t i = 0; i < m_local_size; i++) {
      if (i == m_local_rank) {
        for (uint32_t j = 0; j < m_remote_size; j++) {
          if (j == m_remote_rank) continue;
          m_remote_exchanger.queue(
              j, message{1, 0, uint32_t(m_local_rank), j, data});
          if (++m_send_count >= m_batch_size) do_exchange();
        }
        continue;
      }
      m_local_exchanger.queue(i,
                              message{1, 0, i, uint32_t(m_remote_rank), data});
      if (++m_send_count >= m_batch_size) {
        do_exchange();
      }
    }
  }

  bool global_empty() { return do_exchange(); }

 private:
  uint64_t do_exchange() {
    m_count_exchanges++;
    m_total_sent += m_send_count;
    m_send_count = 0;
    uint64_t total(0);
    total += m_local_exchanger.exchange([&](const message &msg) {
      if (msg.local == m_local_rank && msg.node == m_remote_rank) {
        // we are the destination
        if (msg.bcast) {
          for (uint32_t i = 0; i < m_remote_size; i++) {
            if (i == m_remote_rank) continue;
            m_remote_exchanger.queue(
                i, message{1, 0, uint32_t(m_local_rank), i, msg.data});
          }
        }
        // retire the message
        m_recv_func(msg.bcast, msg.data);
      } else {
        // forwarding with remote exchange
        m_remote_exchanger.queue(msg.node, msg);
      }
    });
    total += m_remote_exchanger.exchange(
        [&](const message &msg) { m_recv_func(msg.bcast, msg.data); });
    return total;
  }

 private:
  comm_exchanger<message> m_local_exchanger;
  comm_exchanger<message> m_remote_exchanger;
  RecvHandlerFunc         m_recv_func;
  size_t                  m_send_count = 0;
  size_t                  m_batch_size;
  uint64_t                m_max_alloc;
  uint64_t                m_count_exchanges = 0;
  uint32_t                m_total_sent = 0;
  int                     m_mpi_size;
  int                     m_mpi_rank;
  MPI_Comm                m_local_comm;
  int                     m_local_size;
  int                     m_local_rank;
  MPI_Comm                m_remote_comm;
  int                     m_remote_size;
  int                     m_remote_rank;
  // int                  m_local_size;
  // uint32_t m_total_recv = 0;

  //  public:
  //   value_measurements get_measurements() {
  //     value_measurements vm;
  //     vm.ingest = m_time_ingest;
  //     vm.cpu = m_time_cpu;
  //     vm.comm = m_time_comm;
  //     vm.wait = m_time_wait;
  //     return vm;
  //   }
};