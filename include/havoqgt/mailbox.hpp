/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef HAVOQGT_MPI_MAILBOX_ROUTED_HPP_INCLUDED
#define HAVOQGT_MPI_MAILBOX_ROUTED_HPP_INCLUDED

#include <havoqgt/mpi.hpp>
#include <havoqgt/environment.hpp>
#include <vector>
#include <list>
#include <limits>
#include <deque>
#include <boost/unordered_map.hpp>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <boost/tuple/tuple.hpp>

namespace havoqgt { namespace mpi {





template <typename TMsg>
class mailbox_routed {

  // class routed_msg_type {
  // public:
  //   routed_msg_type(uint32_t _dest, const TMsg& _msg)
  //     : msg(_msg) { }
  //   TMsg     msg;
  //   //uint32_t m_dest;
  //   uint32_t dest() const { return msg.vertex.owner();/*m_dest;*/ }
  //   bool     is_tree_op() const { return msg.vertex.is_bcast() ||
  //                                        msg.vertex.is_parent_op(); }
  // };

  struct routed_msg_type {
    uint64_t dest      : 30;
    uint64_t bcast     : 1;
    uint64_t intercept : 1;
    TMsg        msg;
  };

  class twod_router{
  public:
    //const static int procs_per_node = 16;
    int procs_per_node;
    twod_router() { }
    twod_router(uint32_t rank, uint32_t size) {
      procs_per_node = std::min(size,uint32_t(24));
      uint64_t rank_to_print = 2;
      //if(rank == rank_to_print)
      //std::cout << "Rank " << rank_to_print << "'s bcast_targets: ";
      my_node_base_rank = rank - (rank % procs_per_node);;
      for(uint32_t i=0; i<procs_per_node; ++i) {
        uint64_t target = my_node_base_rank + i;
        m_bcast_targets.push_back(target);
        //if(rank==rank_to_print)
        //std::cout << target << " ";
      }
      //if(rank==rank_to_print)
      //  std::cout << std::endl << "Rank " << rank_to_print << "'s bcast_proxies: ";
      uint64_t node_offset = rank % procs_per_node;
      for(uint32_t i=0; i<size / procs_per_node; ++i) {
        uint64_t proxy = (i * procs_per_node) + node_offset;
        m_bcast_proxies.push_back(proxy);
        //if(rank==rank_to_print)
        //  std::cout << proxy << " ";
      }
      //if(rank == rank_to_print)
      //  std::cout << std::endl;
    }
    uint32_t proxy_rank(uint32_t dest) {
      uint64_t dest_offset = dest % procs_per_node;
      return my_node_base_rank + dest_offset;
    }
    const std::vector<uint32_t>& bcast_targets() const { return m_bcast_targets; }
    const std::vector<uint32_t>& bcast_proxies() const { return m_bcast_proxies; }
  private:
    std::vector<uint32_t> m_bcast_proxies;
    std::vector<uint32_t> m_bcast_targets;
    uint64_t my_node_base_rank;
  };


  class msg_buffer {
  public:
    msg_buffer():m_size(0),m_ptr(NULL) { }
    msg_buffer(void* _ptr):m_size(0),m_ptr(_ptr) { }

    size_t push_back(const routed_msg_type& _msg) {
      static_cast<routed_msg_type*>(m_ptr)[m_size] = _msg;
      return ++m_size;
    }

    routed_msg_type& operator[](size_t i) { return static_cast<routed_msg_type*>(m_ptr)[i]; }

    size_t size_in_bytes() { return m_size * sizeof(routed_msg_type); }

    size_t size() const { return m_size; }
    bool empty() const { return m_size == 0; }
    void clear() { m_size = 0; m_ptr = NULL; }
    bool is_init() const {return m_ptr != NULL; }
    void* get_ptr() const {return m_ptr; }


  private:
    size_t m_size;
    void* m_ptr;
  };
public:
  typedef TMsg message_type;

  mailbox_routed(int _mpi_tag):
              m_mpi_comm(MPI_COMM_WORLD),
              m_mpi_tag(_mpi_tag) {
    m_pending_partial_buffers = 0;
    m_num_pending_isend = 0;

    m_last_recv_count = 0;

    //statistics
    m_mpi_send_counter  = 0;
    m_tree_send_counter = 0;
    m_route_counter     = 0;
    m_send_counter      = 0;
    m_recv_counter      = 0;

    m_receiving = false;


    CHK_MPI(MPI_Comm_rank(m_mpi_comm, &m_mpi_rank) );
    CHK_MPI(MPI_Comm_size(m_mpi_comm, &m_mpi_size) );

    //Allocate buffer slots for each rank.
    //This does not acutally allocate the buffer's memory
    m_buffer_per_rank.resize(m_mpi_size);
    m_list_isend_request_per_rank.resize(m_mpi_size);
    m_pending_iterator_per_rank.resize(m_mpi_size,m_list_pending.end());


    for(size_t i=0; i<get_environment().mailbox_num_irecv(); ++i) {
      void* irecv_buff = NULL;
      int ret = posix_memalign(&irecv_buff, 32,
                     get_environment().mailbox_aggregation() * sizeof(routed_msg_type));
      if(ret !=0) {
        perror("posix_memalign-irecv"); exit(-1);
      }
      post_new_irecv(irecv_buff);
    }
    m_2d_comm = twod_router(m_mpi_rank, m_mpi_size);

    m_tree_parent = m_mpi_rank / 2;
    m_tree_child1 = (m_mpi_rank * 2) + 1;
    m_tree_child2 = (m_mpi_rank * 2) + 2;

    if(m_mpi_rank ==0) {
      std::cout << "sizeof(message_type) = " << sizeof(message_type) << std::endl;
    }

  }

  ~mailbox_routed() {
    assert(m_pending_partial_buffers == 0);
    while(m_num_pending_isend > 0) {
      cleanup_pending_isend_requests(true);
    }
    while(!m_list_irecv_request.empty()) {
      CHK_MPI( MPI_Cancel( &(m_list_irecv_request.front().first) ) );
      free(m_list_irecv_request.front().second);
      m_list_irecv_request.pop_front();
    }
    for(size_t i=0; i<m_vec_free_buffers.size(); ++i) {
      free(m_vec_free_buffers[i]);
    }

    if(get_environment().mailbox_print_stats()) {
      uint64_t g_mpi_send_counter  = mpi_all_reduce(m_mpi_send_counter, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t g_tree_send_counter = mpi_all_reduce(m_tree_send_counter, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t g_route_counter     = mpi_all_reduce(m_route_counter, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t g_send_counter      = mpi_all_reduce(m_send_counter, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t g_recv_counter      = mpi_all_reduce(m_recv_counter, std::plus<uint64_t>(), MPI_COMM_WORLD);
      if(m_mpi_rank == 0) {
        std::cout << "******************  Mailbox Statistics ********************" << std::endl;
        std::cout << "routed message size = " << sizeof(TMsg) << std::endl;
        std::cout << "mpi_send_counter  = " << g_mpi_send_counter  << std::endl;
        std::cout << "tree_send_counter = " << g_tree_send_counter << std::endl;
        std::cout << "route_counter     = " << g_route_counter     << std::endl;
        std::cout << "send_counter      = " << g_send_counter      << std::endl;
        std::cout << "recv_counter      = " << g_recv_counter      << std::endl;
        std::cout << "Average send size = " << double(g_send_counter + g_route_counter) / double(g_mpi_send_counter) << std::endl;
        std::cout << "***********************************************************" << std::endl;
      }
    }

    CHK_MPI( MPI_Barrier( m_mpi_comm ) );
  }

  /*void send_tree_fast(int raw_dest, const TMsg& _raw_msg) {
    ++m_tree_send_counter;
    routed_msg_type msg(raw_dest, _raw_msg);
    if(!m_buffer_per_rank[raw_dest].is_init()) {
      m_buffer_per_rank[raw_dest] = allocate_msg_buffer();
      ++m_pending_partial_buffers;
      m_list_pending.push_back(raw_dest);
      m_pending_iterator_per_rank[raw_dest] = --m_list_pending.end();
    }
    size_t size = m_buffer_per_rank[raw_dest].push_back(msg);

    if(size >= get_environment().mailbox_tree_aggregation() ) {
      post_isend(raw_dest);
      check_for_starvation();
    }
  }

  void send_tree_parent(const TMsg& _raw_msg) {
    send_tree_fast(m_tree_parent, _raw_msg);
  }

  void send_tree_children(const TMsg& _raw_msg) {
    if(m_tree_child1 < m_mpi_size) {
      send_tree_fast(m_tree_child1, _raw_msg);
      if(m_tree_child2 < m_mpi_size) {
        send_tree_fast(m_tree_child2, _raw_msg);
      }
    }
  }
  */

  void route_fast_path(uint32_t dest, const routed_msg_type& _msg) {
    ++m_route_counter;
    if(!m_buffer_per_rank[dest].is_init()) {
      m_buffer_per_rank[dest] = allocate_msg_buffer();
      ++m_pending_partial_buffers;
      m_list_pending.push_back(dest);
      m_pending_iterator_per_rank[dest] = --m_list_pending.end();
    }
    size_t size = m_buffer_per_rank[dest].push_back(_msg);

    if(size == get_environment().mailbox_aggregation() ) {
      post_isend(dest);
      check_for_starvation();
    }
  }

  template <typename OutputIterator>
  void bcast(TMsg _raw_msg, OutputIterator _oitr) {
    _raw_msg.vertex.set_bcast(true);
    routed_msg_type msg;
    msg.bcast = true;
    msg.msg = _raw_msg;
    msg.dest = m_mpi_size+1;
    for(uint32_t i=0; i<m_2d_comm.bcast_proxies().size(); ++i) {
      uint32_t proxy_rank = m_2d_comm.bcast_proxies()[i];
      if(proxy_rank == uint32_t(m_mpi_rank)) {
        bcast_to_targets(msg);
      } else {
        route_fast_path(proxy_rank, msg);
      }
    }
    if(m_receiving) return;   //prevent receive recursion @todo, make this a function called wait
    do {
      cleanup_pending_isend_requests();
      receive(_oitr);
    } while(m_num_pending_isend > get_environment().mailbox_num_isend());
  }

  void bcast_to_targets(routed_msg_type msg) {
    assert(msg.bcast);
    for(uint32_t i=0; i<m_2d_comm.bcast_targets().size(); ++i) {
      uint32_t target = m_2d_comm.bcast_targets()[i];
      msg.dest= target;
      route_fast_path(target, msg);
    }
  }




  void check_for_starvation() {
    /*if(!m_list_pending.empty()) {
      size_t to_check = m_list_pending.front();
      size_t size = m_send_ticket_per_rank.size();
      //size_t to_check = ++m_fair_checker % size;
      if(m_send_counter > m_send_ticket_per_rank[to_check] + size) {
        post_isend(to_check);
      }
    }*/
  }

  template <typename OutputIterator>
  void send(int raw_dest, const TMsg& _raw_msg, OutputIterator _oitr, bool _intercept) {
    ++m_send_counter;
    assert(raw_dest >= 0 && raw_dest < m_mpi_size);
    //++m_last_recv_count;
    assert(raw_dest != m_mpi_rank); // just dont send to self!
    //routed_msg_type msg(raw_dest, _raw_msg);
    routed_msg_type msg; //@todo fixme
    msg.msg  = _raw_msg;
    msg.dest = raw_dest;
    msg.intercept = _intercept;
    msg.bcast = false;
    int dest = m_2d_comm.proxy_rank(raw_dest);
    if(dest == m_mpi_rank) dest = raw_dest;
    if(!m_buffer_per_rank[dest].is_init()) {
      m_buffer_per_rank[dest] = allocate_msg_buffer();
      ++m_pending_partial_buffers;
      m_list_pending.push_back(dest);
      m_pending_iterator_per_rank[dest] = --m_list_pending.end();
    }
    size_t size = m_buffer_per_rank[dest].push_back(msg);

    if(size == get_environment().mailbox_aggregation() ) {
      post_isend(dest);
      check_for_starvation();
      if(m_receiving) return;   //prevent receive recursion
      do {
        //cleanup_pending_isend_requests_index(dest);
        cleanup_pending_isend_requests();
        receive(_oitr);
      //} while(m_list_isend_request_per_rank[dest].size() > 1);
      } while(m_num_pending_isend > get_environment().mailbox_num_isend());
    }
    //if(m_last_recv_count > 128) {
    //  receive(_oitr);
    //}
  }

  template <typename OutputIterator > //was vector but also want list
  void receive(OutputIterator _oitr, bool aggregsive=false) {
    m_receiving = true;
    m_last_recv_count = 0;
    int flag(0);
    //do {
      if(!m_list_irecv_request.empty()) {
        std::pair<MPI_Request, void* > pair_req = m_list_irecv_request.front();
        m_list_irecv_request.pop_front();
        MPI_Request* request_ptr = &(pair_req.first);
        MPI_Status status;
        CHK_MPI( MPI_Test( request_ptr, &flag, &status) );
        if(flag) {
          routed_msg_type* recv_ptr = static_cast<routed_msg_type*> (
                           pair_req.second );
          int count(0);
          CHK_MPI( MPI_Get_count(&status, MPI_BYTE, &count) );
          for(size_t i=0; i<count/sizeof(routed_msg_type); ++i) {
            if(recv_ptr[i].dest == uint32_t(m_mpi_rank) /*|| recv_ptr[i].is_tree_op()*/) {
              *_oitr = recv_ptr[i].msg;//.msg;
              ++_oitr;
              ++m_recv_counter;
            } else if(recv_ptr[i].bcast) {
              bcast_to_targets(recv_ptr[i]);
            } else if(recv_ptr[i].intercept) {
              if( _oitr.intercept(recv_ptr[i].msg) ) {
                route_fast_path(recv_ptr[i].dest, recv_ptr[i]);
              }
            } else {
              route_fast_path(recv_ptr[i].dest, recv_ptr[i]);
            }
          }
          post_new_irecv(recv_ptr);
        } else {
          m_list_irecv_request.push_front(pair_req);
        }
      }
    //} while(flag);
    m_receiving = false;
  }


  bool is_idle() {
    cleanup_pending_isend_requests();
    return m_pending_partial_buffers == 0 && m_num_pending_isend == 0;
  }

  void flush_buffers_if_idle() {
    if(!m_list_pending.empty()) {
      size_t index = m_list_pending.front();
      if(m_num_pending_isend < 1 && !m_buffer_per_rank[index].empty()) {
        post_isend(index);
      }
    }
  }

  int comm_rank() const { return m_mpi_rank; }
  int comm_size() const { return m_mpi_size; }


private:
  msg_buffer allocate_msg_buffer() {
    if(m_vec_free_buffers.empty()) {
      void* buff = NULL;
      int ret = posix_memalign(&buff, 32,
                      get_environment().mailbox_aggregation() * sizeof(routed_msg_type));
      if(ret !=0) {
        perror("posix_memalign"); exit(-1);
      }
      m_vec_free_buffers.push_back(buff);
    }
    msg_buffer to_return(m_vec_free_buffers.back());
    m_vec_free_buffers.pop_back();
    return to_return;
  }

  void free_msg_buffer(void* _ptr) {
    m_vec_free_buffers.push_back(_ptr);
  }

  void post_isend(int index) {
    if(m_buffer_per_rank[index].empty()) return;
    m_mpi_send_counter++;
    int dest = index;
    bool was_first_pending = false;
    if(m_pending_iterator_per_rank[dest] != m_list_pending.end()) {
      if(m_pending_iterator_per_rank[dest] == m_list_pending.begin()) {
        was_first_pending = true;
      }
      m_list_pending.erase(m_pending_iterator_per_rank[dest]);
      m_pending_iterator_per_rank[dest] = m_list_pending.end();
    }

    m_list_isends.push_back(index);
    boost::tuple<MPI_Request, void*,std::list<size_t>::iterator> isend_req_tuple;
    MPI_Request* request_ptr = &(isend_req_tuple.get<0>());
    isend_req_tuple.get<1>() = m_buffer_per_rank[index].get_ptr();
    isend_req_tuple.get<2>() = --m_list_isends.end();
    void* buffer_ptr = m_buffer_per_rank[index].get_ptr();
    int size_in_bytes = m_buffer_per_rank[index].size_in_bytes();

    CHK_MPI( MPI_Issend( buffer_ptr, size_in_bytes, MPI_BYTE, dest,
                        m_mpi_tag, m_mpi_comm, request_ptr) );

    --m_pending_partial_buffers;
    m_buffer_per_rank[index].clear();
    //int flag(0);
    //CHK_MPI( MPI_Test( request_ptr, &flag, MPI_STATUS_IGNORE) );
    //if(!flag) {
    m_list_isend_request_per_rank[index].push_back(isend_req_tuple);
    ++m_num_pending_isend;
    //}
    cleanup_pending_isend_requests();
    if(!was_first_pending/* && m_tree_parent != index && m_tree_child1 != index && m_tree_child2 != index*/) {
      post_isend(m_list_pending.front());
    }
  }

  bool cleanup_pending_isend_requests_index(size_t index) {
    bool to_return = false;
    if(m_list_isend_request_per_rank[index].empty()) return true;
    while(!m_list_isend_request_per_rank[index].empty()) {
      int flag(0);
      MPI_Request* request_ptr = &(m_list_isend_request_per_rank[index].front().get<0>());
      CHK_MPI( MPI_Test( request_ptr, &flag, MPI_STATUS_IGNORE) );
      if(flag) {
        free_msg_buffer(m_list_isend_request_per_rank[index].front().get<1>());
        m_list_isends.erase(m_list_isend_request_per_rank[index].front().get<2>());
        m_list_isend_request_per_rank[index].pop_front();
        --m_num_pending_isend;
        to_return = true;
      } else {
        break;
      }
    }
    return to_return;
  }

  /// If the number of current pending isends is less than user desired,
  /// just check the next fair bit.   If the number is greater than
  /// desired, aggregsively check all bits.
  void cleanup_pending_isend_requests(bool force_aggressive = false) {
    while(!m_list_isends.empty()) {
      bool found = cleanup_pending_isend_requests_index(m_list_isends.front());
      if(!found) break;
    }
  }

  void post_new_irecv(void* _buff) {
    std::pair<MPI_Request, void*> irecv_req;
    irecv_req.second = _buff;
    MPI_Request* request_ptr = &(irecv_req.first);
    int num_bytes = get_environment().mailbox_aggregation() * sizeof(routed_msg_type);
    CHK_MPI( MPI_Irecv( _buff, num_bytes, MPI_BYTE, MPI_ANY_SOURCE,
                        m_mpi_tag, m_mpi_comm, request_ptr) );
    m_list_irecv_request.push_back(irecv_req);
  }


  /// MPI configuration
  MPI_Comm m_mpi_comm;
  int m_mpi_tag;
  int m_mpi_rank;
  int m_mpi_size;

  size_t m_pending_partial_buffers;
  size_t m_num_pending_isend;
  size_t                  m_last_recv_count;

  std::vector<void*> m_vec_free_buffers;
  std::vector<msg_buffer> m_buffer_per_rank;
  //boost::unordered_map<uint32_t, msg_buffer> m_buffer_per_rank;

  std::vector< std::list< boost::tuple<MPI_Request, void*, std::list<size_t>::iterator > > > m_list_isend_request_per_rank;
  //boost::unordered_map<uint64_t, std::list< boost::tuple<MPI_Request, void*, std::list<size_t>::iterator > > > m_list_isend_request_per_rank;
  std::list <size_t> m_list_isends;
  std::list < std::pair<MPI_Request, void*> > m_list_irecv_request;

  twod_router m_2d_comm;

  std::vector< std::list<size_t>::iterator > m_pending_iterator_per_rank;
  //boost::unordered_map<uint64_t, std::list<size_t>::iterator > m_pending_iterator_per_rank;
  std::list <size_t> m_list_pending;

  int m_tree_parent;
  int m_tree_child1;
  int m_tree_child2;

  bool m_receiving;

  //Statistics
  uint64_t                m_mpi_send_counter;
  uint64_t                m_tree_send_counter;
  uint64_t                m_route_counter;
  uint64_t                m_send_counter;
  uint64_t                m_recv_counter;

};


}} //namespace havoqgt { namespace mpi {



#endif //HAVOQGT_MPI_MAILBOX_ROUTED_HPP_INCLUDED
