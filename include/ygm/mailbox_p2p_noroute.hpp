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
class mailbox_p2p_noroute {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t dest : 30;
    Data     data;
  };  //__attribute__((packed));

 public:
  mailbox_p2p_noroute(RecvHandlerFunc recv_func, size_t batch_size,
                      MPI_Comm local_comm, MPI_Comm remote_comm)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        m_exchanger(MPI_COMM_WORLD, 1) {  // should change tag eventually
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size));
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
  }

  ~mailbox_p2p_noroute() {
    if (m_mpi_rank == 0) {
      std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    }
  }

  void send(uint32_t dest, Data data) {
    if (dest == m_mpi_rank) {
      m_recv_func(false, data);
    } else {
      m_exchanger.queue(dest, message{0, 0, dest, data});
      if (++m_send_count >= m_batch_size) { do_exchange(); }
    }
  }

  void send_bcast(Data data) {
    for (uint32_t i = 0; i < m_mpi_size; i++) {
      if (i == m_mpi_rank) continue;
      m_exchanger.queue(i, message{1, 0, i, data});
      ++m_send_count;
    }
    if (m_send_count >= m_batch_size) { do_exchange(); }
    // bcast to self
    // m_recv_func(true, data);
  }

  bool global_empty() { return do_exchange() == 0; }

 private:
  uint64_t do_exchange() {
    m_count_exchanges++;
    m_total_sent += m_send_count;
    m_send_count = 0;
    return m_exchanger.exchange(
        [&](const message &msg) { m_recv_func(msg.bcast, msg.data); });
  }

 private:
  comm_exchanger<message> m_exchanger;
  RecvHandlerFunc         m_recv_func;
  size_t                  m_send_count = 0;
  size_t                  m_batch_size;
  uint64_t                m_max_alloc;
  uint64_t                m_count_exchanges = 0;
  uint32_t                m_total_sent = 0;
  int                     m_mpi_size;
  int                     m_mpi_rank;
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
}  // namespace ygm
