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

namespace ygm {
template <typename Data, typename RecvHandlerFunc>
class mailbox_p2p_nrroute {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t local : 6;
    uint32_t node : 24;  // Supports addressing <= 16777216 nodes w/ <= 64 cores
    Data     data;
  };  //__attribute__((packed));

 public:
  mailbox_p2p_nrroute(RecvHandlerFunc recv_func, size_t batch_size)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        // should change tag eventually
        m_local_exchanger(comm_nl().mpi_comm(), 1),
        m_remote_exchanger(comm_nr().mpi_comm(), 2) {}

  ~mailbox_p2p_nrroute() {
    wait_empty();
    if (comm_world().rank() == 0) {
      std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    }
  }

  void send(uint32_t dest, Data data) {
    if (dest == comm_world().rank()) {
      m_recv_func(false, data);
    } else {
      uint32_t local = dest % comm_nl().size();
      uint32_t node  = dest / comm_nl().size();
      if (node == comm_nr().rank()) {
        m_local_exchanger.queue(local, message{0, 0, local, node, data});
      } else {
        m_remote_exchanger.queue(node, message{0, 0, local, node, data});
      }
      if (++m_send_count >= m_batch_size) {
        do_exchange();
      }
    }
  }

  void send_bcast(Data data) {
    for (uint32_t i = 0; i < comm_nr().size(); i++) {
      if (i == comm_nr().rank()) continue;
      m_remote_exchanger.queue(i, message{1, 0, 0, 0, data});
      ++m_send_count;
    }
    for (uint32_t j = 0; j < comm_nl().size(); j++) {
      if (j == comm_nl().rank()) continue;
      m_local_exchanger.queue(j, message{1, 0, 0, 0, data});
      ++m_send_count;
    }
    if (m_send_count >= m_batch_size) do_exchange();
    // bcast to self
    // m_recv_func(true,data);
  }

  bool global_empty() { return do_exchange() == 0; }

  void wait_empty() {
    do {
    } while (!global_empty());
  }

 private:
  /// WARNING, this count return is kinda flaky, not good...
  uint64_t do_exchange() {
    m_count_exchanges++;
    m_total_sent += m_send_count;
    uint64_t total = m_remote_exchanger.exchange(
        [&](const message &msg) {
          if (msg.bcast) {
            for (uint32_t i = 0; i < comm_nl().size(); i++) {
              if (i == comm_nl().rank())
                m_recv_func(msg.bcast, msg.data);
              else
                m_local_exchanger.queue(i, msg);
            }
          } else if (msg.local == comm_nl().rank() &&
                     msg.node == comm_nr().rank()) {
            // we are the destination
            m_recv_func(msg.bcast, msg.data);
          } else {
            // forwarding with local exchange
            m_local_exchanger.queue(msg.local, msg);
          }
        },
        m_send_count);
    total += m_local_exchanger.exchange(
        [&](const message &msg) { m_recv_func(msg.bcast, msg.data); }, total);
    m_send_count = 0;
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
  uint32_t                m_total_sent      = 0;

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
}  // namespace ygm
