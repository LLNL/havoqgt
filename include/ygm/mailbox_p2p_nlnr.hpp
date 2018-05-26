#include <ygm/comm_exchanger.hpp>
#include <ygm/mpi.hpp>

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

using std::vector;

template <typename Data, typename RecvHandlerFunc>
class mailbox_p2p_nlnr {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t local : 6;
    uint32_t node : 24;  // Supports addressing <= 16777216 nodes w/ <= 64 cores
    Data     data;
  };

 public:
  mailbox_p2p_nlnr(RecvHandlerFunc recv_func, size_t batch_size,
                   MPI_Comm local_comm, MPI_Comm remote_comm)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        m_local_comm(local_comm),
        m_remote_comm(remote_comm),
        m_first_local_exchanger(local_comm, 1),  // should change tag eventually
        m_remote_exchanger(remote_comm, 2),
        m_second_local_exchanger(local_comm, 3) {
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size));
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
    CHK_MPI(MPI_Comm_size(m_local_comm, &m_local_size));
    CHK_MPI(MPI_Comm_rank(m_local_comm, &m_local_rank));
    CHK_MPI(MPI_Comm_size(m_remote_comm, &m_remote_size));
    CHK_MPI(MPI_Comm_rank(m_remote_comm, &m_remote_rank));
    m_node_rank = m_mpi_rank / m_local_size;
    m_layer_rank = m_node_rank % m_local_size;
    m_parity = m_remote_rank & 1;
  }

    ~mailbox_p2p_nlnr() {
    if(m_mpi_rank == 0) {
      std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    }
  }

  void send(uint32_t dest, Data data) {
    if (dest == m_mpi_rank) {
      m_recv_func(false, data);
    } else {
      uint32_t local = dest % m_local_size;
      uint32_t node = dest / m_local_size;
      if (node == m_node_rank) {
        // true local
        m_second_local_exchanger.queue(local, message{0, 0, local, node, data});
        return;
      }
      uint32_t layer = node % m_local_size;
      uint32_t remote = node / m_local_size;
      if (m_layer_rank != layer) {
        remote *= 2;
        if (layer > m_layer_rank) ++remote;
      }
      if (m_local_rank == layer) {
        // remote
        m_remote_exchanger.queue(remote, message{0, 0, local, remote, data});
      } else {
        // forwarding local
        m_first_local_exchanger.queue(layer,
                                      message{0, 0, local, remote, data});
      }
      if (++m_send_count >= m_batch_size) do_exchange();
    }
  }

  // Verify that node/core address is not needed.
  void send_bcast(Data data) {
    for (uint32_t i = 0; i < m_local_size; i++) {
      if (i == m_local_rank) {
        for (uint32_t j = 0; j < m_remote_size; j++) {
          if (j == m_remote_rank) continue;
          if ((m_local_rank != m_layer_rank && (m_parity == (j & 1)))) continue;
          m_remote_exchanger.queue(j, message{1, 0, 0, 0, data});
          if (++m_send_count >= m_batch_size) do_exchange();
        }
        continue;
      }
      m_first_local_exchanger.queue(i, message{1, 0, 0, 0, data});
      if (++m_send_count >= m_batch_size) do_exchange();
    }
  }

  bool global_empty() { return do_exchange() == 0; }

 private:
  uint64_t do_exchange() {
    m_count_exchanges++;
    m_total_sent += m_send_count;
    m_send_count = 0;
    uint64_t total(0);
    // first local exchange
    total += m_first_local_exchanger.exchange([&](const message &msg) {
      if (msg.bcast) {
        for (uint32_t i = 0; i < m_remote_size; i++) {
          if (i == m_remote_rank) continue;
          if ((m_local_rank != m_layer_rank && (m_parity == (i & 1)))) continue;
          m_remote_exchanger.queue(i, msg);
        }
        m_recv_func(msg.bcast, msg.data);
      } else {
        // !!! right now, msg.node holds the remote id, not the node id !!!
        m_remote_exchanger.queue(msg.node, msg);
      }
    });
    // remote exchange
    total += m_remote_exchanger.exchange([&](const message &msg) {
      if (msg.bcast) {
        for (uint32_t i = 0; i < m_local_size; i++) {
          if (i == m_local_rank) continue;
          m_second_local_exchanger.queue(i, msg);
        }
        m_recv_func(msg.bcast, msg.data);
      } else {
        if (msg.local == m_local_rank) {
          m_recv_func(msg.bcast, msg.data);
        } else {
          // right now not modifying the msg's fields. Will probably need to do
          // this in a general implementation.
          m_second_local_exchanger.queue(msg.local, msg);
          // m_second_local_exchanger.queue(msg.local,
          //                                msg{0, 0, msg.local, m_node_rank,
          //                                msg.data});
        }
      }
    });
    // second local exchange
    total += m_second_local_exchanger.exchange(
        [&](const message &msg) { m_recv_func(msg.bcast, msg.data); });
    return total;
  }

 private:
  comm_exchanger<message> m_first_local_exchanger;
  comm_exchanger<message> m_remote_exchanger;
  comm_exchanger<message> m_second_local_exchanger;
  RecvHandlerFunc         m_recv_func;
  size_t                  m_send_count = 0;
  size_t                  m_batch_size;
  uint64_t                m_max_alloc;
  uint64_t                m_count_exchanges = 0;
  uint32_t                m_total_sent = 0;
  int                     m_mpi_size;
  int                     m_mpi_rank;
  int                     m_node_rank;
  int                     m_layer_rank;
  int                     m_parity;
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