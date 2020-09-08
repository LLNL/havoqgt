// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <ygm/mpi.hpp>
#include <ygm/comm_exchanger.hpp>

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

using std::vector;

namespace ygm {
template < typename Data, typename RecvHandlerFunc,
         template <typename> class EXCHANGER = comm_exchanger >
class mailbox_p2p_noroute {
  struct message {
    uint32_t bcast : 1;
    uint32_t interrupt : 1;
    uint32_t dest : 30;
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
      dummy = (dest << 2) + (interrupt << 1) + bcast;
    }

    void decode(uint32_t& input) {
      bcast = input & 1;
      interrupt = (input >> 1) & 1;
      dest = (input >> 2);
    }
  };  //__attribute__((packed));

 public:
  mailbox_p2p_noroute(RecvHandlerFunc recv_func, const size_t batch_size,
		      const int tag = 1)
      : m_recv_func(recv_func),
        m_batch_size(batch_size),
        m_max_alloc(0),
        m_exchanger(comm_world().mpi_comm(), tag) {}

  ~mailbox_p2p_noroute() {
    wait_empty();
    //if (comm_world().rank() == 0) {
    //  std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    //}
  }

  void send(uint32_t dest, const Data& data) {
    if (dest == comm_world().rank()) {
      m_recv_func(this, false, data);
    } else {
      if (in_exchange) {
        m_send_queue.push_back(std::make_pair(dest, data));
      } else {
        do_send(dest, data);
        if (m_send_count >= m_batch_size) do_exchange();
      }
      /*
      m_exchanger.queue(dest, message{0, 0, dest, data});
      if (++m_send_count >= m_batch_size) { do_exchange(); }
      */
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
    /*
    for (uint32_t i = 0; i < comm_world().size(); i++) {
      if (i == comm_world().rank()) continue;
      m_exchanger.queue(i, message{1, 0, i, data});
      ++m_send_count;
    }
    if (m_send_count >= m_batch_size) { do_exchange(); }
    */
    // bcast to self
    // m_recv_func(true, data);
  }

  bool global_empty() { return do_exchange() == 0; }

  void wait_empty() {
    do {
    } while (!global_empty());
  }

 private:
  void do_send(uint32_t dest, const Data& data) {
    //m_exchanger.queue(dest, message{0, 0, dest, data});
    //++m_send_count;
    m_send_count += m_exchanger.queue_bytes(
        dest, message{0, 0, dest, data});
  }

  void do_send_bcast(const Data& data) {
    for (uint32_t i = 0; i < comm_world().size(); i++) {
      if (i == comm_world().rank()) continue;
      m_send_count +=
          m_exchanger.queue_bytes(i, message{1, 0, i, data});
      //m_exchanger.queue(i, message{1, 0, i, data});
      //++m_send_count;
    }
    //if (m_send_count >= m_batch_size) { do_exchange(); }
  }

  uint64_t do_exchange() {
    in_exchange = true;
    m_count_exchanges++;
    m_total_sent += m_send_count;
    uint64_t total(0);
    total =  m_exchanger.exchange(
        [&](const message &msg) { m_recv_func(this, msg.bcast, msg.data); });
    m_send_count = 0;
    //
    // push out recursive queued
    auto queue_size = m_send_queue.size() + m_bcast_queue.size();
    //if (queue_size > 0) {
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
  //comm_exchanger<message> m_exchanger;
  EXCHANGER<message>      m_exchanger;
  RecvHandlerFunc         m_recv_func;
  size_t                  m_send_count = 0;
  size_t                  m_batch_size;
  uint64_t                m_max_alloc;
  uint64_t                m_count_exchanges = 0;
  uint32_t                m_total_sent = 0;

  uint64_t m_local_send  = 0;
  uint64_t m_local_bcast = 0;

  bool                                   in_exchange = false;
  std::vector<std::pair<uint32_t, Data>> m_send_queue;
  std::vector<Data>                      m_bcast_queue;

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
