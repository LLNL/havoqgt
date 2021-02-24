// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

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
template <typename Data, typename RecvHandlerFunc,
          template <typename> class EXCHANGER = comm_exchanger>
class mailbox_p2p_nrroute {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t local : 6;
    uint32_t node : 24;  // Supports addressing <= 16777216 nodes w/ <= 64 cores
    Data     data;
    template <class Archive>
    void save(Archive& ar) const {
      uint32_t dummy;
      encode(dummy);
      ar(dummy, data);
    }

    template <class Archive>
    void load(Archive& ar) {
      uint32_t dummy;
      ar(dummy, data);
      decode(dummy);
    }

    void encode(uint32_t& dummy) const {
      dummy = (node << 8) + (local << 2) + (interrupt << 1) + bcast;
    }

    void decode(uint32_t& input) {
      bcast     = input & 1;
      interrupt = (input >> 1) & 1;
      local     = (input >> 2) & 63;
      node      = (input >> 8);
    }
  };  //__attribute__((packed));

 public:
  mailbox_p2p_nrroute(RecvHandlerFunc recv_func, const size_t batch_size,
                      const int tag = 1)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        m_local_exchanger(comm_nl().mpi_comm(), tag),
        m_remote_exchanger(comm_nr().mpi_comm(), tag + 1),
        m_num_nodes(comm_world().size() / comm_nl().size()),
        m_per_node_batch_size(2 * m_batch_size / m_num_nodes),
        m_last_dest_node(0) {
    m_per_node_send_count.resize(m_num_nodes);
  }

  ~mailbox_p2p_nrroute() {
    wait_empty();
    // if (comm_world().rank() == 0) {
    //   std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    // }
    // std::cout << whoami() << " m_local_send = " << m_local_send
    //          << ", m_local_bcast = " << m_local_bcast << std::endl;
  }

  void send(uint32_t dest, const Data& data) {
    ++m_local_send;
    if (dest == comm_world().rank()) {
      m_recv_func(this, false, data);
    } else {
      if (in_exchange) {
        m_send_queue.push_back(std::make_pair(dest, data));
      } else {
        do_send(dest, data);
        if (m_send_count >= m_batch_size ||
            m_per_node_send_count[m_last_dest_node] >= m_per_node_batch_size) {
          do_exchange();
          std::fill(m_per_node_send_count.begin(), m_per_node_send_count.end(),
                    0);
        }
      }
    }
  }

  void send_bcast(const Data& data) {
    ++m_local_bcast;
    if (in_exchange) {
      m_bcast_queue.push_back(data);
    } else {
      do_send_bcast(data);
      if (m_send_count >= m_batch_size) do_exchange();
    }
  }

  bool global_empty() { return do_exchange() == 0; }

  void wait_empty() {
    do {
    } while (!global_empty());
  }

 private:
  void do_send(uint32_t dest, const Data& data) {
    assert(dest != comm_world().rank());
    uint32_t local = dest % comm_nl().size();
    uint32_t node  = dest / comm_nl().size();
    if (node == comm_nr().rank()) {
      m_send_count += m_local_exchanger.queue_bytes(
          local, message{0, 0, local, node, data});
      // m_local_exchanger.queue(local, message{0, 0, local, node, data});
    } else {
      m_send_count += m_remote_exchanger.queue_bytes(
          node, message{0, 0, local, node, data});
      // m_remote_exchanger.queue(node, message{0, 0, local, node, data});
    }
    // ++m_send_count;
    // if (++m_send_count >= m_batch_size) {
    //   do_exchange();
    // }
  }

  void do_send_bcast(const Data& data) {
    for (uint32_t i = 0; i < comm_nr().size(); i++) {
      if (i == comm_nr().rank()) continue;
      m_send_count +=
          m_remote_exchanger.queue_bytes(i, message{1, 0, 0, 0, data});
      // m_remote_exchanger.queue(i, message{1, 0, 0, 0, data});
      // ++m_send_count;
    }
    for (uint32_t j = 0; j < comm_nl().size(); j++) {
      if (j == comm_nl().rank()) continue;
      m_send_count +=
          m_local_exchanger.queue_bytes(j, message{1, 0, 0, 0, data});
      // m_local_exchanger.queue(j, message{1, 0, 0, 0, data});
      // ++m_send_count;
    }
    // if (m_send_count >= m_batch_size) do_exchange();
    // bcast to self
    // m_recv_func(true,data);
  }

  /// WARNING, this count return is kinda flaky, not good...
  uint64_t do_exchange() {
    in_exchange = true;
    m_count_exchanges++;
    m_total_sent += m_send_count;
    uint64_t total = m_remote_exchanger.exchange(
        [&](const message& msg) {
          if (msg.bcast) {
            for (uint32_t i = 0; i < comm_nl().size(); i++) {
              if (i == comm_nl().rank())
                m_recv_func(this, msg.bcast, msg.data);
              else
                m_local_exchanger.queue(i, msg);
            }
          } else if (msg.local == comm_nl().rank() &&
                     msg.node == comm_nr().rank()) {
            // we are the destination
            m_recv_func(this, msg.bcast, msg.data);
          } else {
            // forwarding with local exchange
            m_local_exchanger.queue(msg.local, msg);
          }
        },
        m_send_count);
    total += m_local_exchanger.exchange(
        [&](const message& msg) { m_recv_func(this, msg.bcast, msg.data); },
        total);
    m_send_count = 0;

    //
    // push out recursive queued
    auto queue_size = m_send_queue.size() + m_bcast_queue.size();
    // if (queue_size > 0) {
    //  std::cout << whoami()
    //            << " recursive m_send_queue.size() = " << m_send_queue.size()
    //            << " m_bcast_queue.size() = " << m_bcast_queue.size()
    //            << " exchange total = " << total << std::endl;
    //}
    for (const auto& p : m_send_queue) {
      do_send(p.first, p.second);
    }
    for (const auto& d : m_bcast_queue) {
      do_send_bcast(d);
    }
    m_send_queue.clear();
    m_bcast_queue.clear();
    in_exchange = false;
    // std::cout << whoami() << " - ending exchange, total = " << total
    //          << std::endl;
    return total;
  }

 private:
  // comm_exchanger<message> m_local_exchanger;
  // comm_exchanger<message> m_remote_exchanger;
  EXCHANGER<message> m_local_exchanger;
  EXCHANGER<message> m_remote_exchanger;
  RecvHandlerFunc    m_recv_func;
  size_t             m_send_count = 0;
  size_t             m_batch_size;
  uint64_t           m_max_alloc;
  uint64_t           m_count_exchanges = 0;
  uint32_t           m_total_sent      = 0;
  uint64_t           m_num_nodes;
  size_t             m_per_node_batch_size;
  uint64_t           m_last_dest_node;

  uint64_t m_local_send  = 0;
  uint64_t m_local_bcast = 0;

  bool in_exchange = false;
  std::vector<std::pair<uint32_t, Data>> m_send_queue;
  std::vector<Data>     m_bcast_queue;
  std::vector<uint64_t> m_per_node_send_count;

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