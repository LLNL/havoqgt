// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MAILBOX_HPP_INCLUDED
#define HAVOQGT_MAILBOX_HPP_INCLUDED


#include <havoqgt/mpi.hpp>


#include <vector>
#include <deque>
#include <list>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>


namespace havoqgt {


template <typename T>
class mailbox {
  const char* msm_fname = "havoqgt_mailbox_msm";
  template<class OPT> using offset_ptr = boost::interprocess::offset_ptr<OPT>;
  using managed_shared_memory = boost::interprocess::managed_shared_memory;

  struct msg_wrapper;
  struct msg_bundle_shm;
  struct msg_bundle_mpi;
  struct shm_exchange;
public:
  mailbox(int _mpi_tag) {
    m_mpi_tag = _mpi_tag;
    m_mpi_comm = comm_nr().mpi_comm();
    m_mpi_rank = comm_nr().rank();
    m_mpi_size = comm_nr().size();
    m_shm_rank = comm_nl().rank();
    m_shm_size = comm_nl().size();
    m_world_rank = comm_world().rank();
    m_world_size = comm_world().size();

    init_environment_config();

    init_managed_shared_memory();

    init_alloc_bundle_shared();

    init_alloc_bundle_mpi();

    // if(m_world_rank == 0) {
    //   std::cout << "Number of nodes = " << m_mpi_size << std::endl;
    //   std::cout << "Ranks per node = " << m_shm_size << std::endl;
    //   std::cout << "Sizeof(shm_exchange) = " << sizeof(shm_exchange) << std::endl;
    //   std::cout << "Sizeof(msg_wrapper) = " << sizeof(msg_wrapper) << std::endl;
    // }
    //lock_time = 0;
    m_send_recv_balance = 0;
    m_isend_count = 0;
    m_isend_bytes = 0;
  }

  ~mailbox() {
    // if(m_world_rank == 0) {
    // std::cout << "ISend Count = " << m_isend_count << std::endl;
    // std::cout << "ISend Bytes = " << m_isend_bytes << std::endl;
    // double ave = (m_isend_count) ? (double(m_isend_bytes) / double(m_isend_count)) : double(0);
    // std::cout << "Ave size = " << ave << std::endl;
    // std::cout << "m_shm_transfer_count = " << m_shm_transfer_count << std::endl;
    // std::cout << "m_shm_transfer_bytes = " << m_shm_transfer_bytes << std::endl;
    // std::cout << "m_free_shm_list.size() = " << m_free_shm_list.size() << std::endl;
    // }

    while(!m_shm_request_list.empty() || !m_isend_request_list.empty()) {
      cleanup_pending_isend_requests();
    }

    while(!m_irecv_request_list.empty()) {
      CHK_MPI( MPI_Cancel( &(m_irecv_request_list.front().first) ) );
      free(m_irecv_request_list.front().second);
      m_irecv_request_list.pop_front();
    }

    m_my_exchange->~shm_exchange(); //frees mutexes
    delete m_pmsm;        //frees managed_shared_memory
    comm_nl().barrier();
    boost::interprocess::shared_memory_object::remove(msm_fname); //removes shared memory file
    for(auto itr = m_free_mpi_list.begin(); itr != m_free_mpi_list.end(); ++itr) {
      free(*itr);

    }
    comm_world().barrier();
  }

private:
  void init_environment_config() {
    if(s_num_isend == 0) {
      s_num_isend = havoqgt_getenv("HAVOQGT_MAILBOX_NUM_ISEND", uint32_t(4));
      s_num_irecv = havoqgt_getenv("HAVOQGT_MAILBOX_NUM_IRECV", uint32_t(4));
      uint32_t mpi_bytes = havoqgt_getenv("HAVOQGT_MAILBOX_MPI_SIZE", uint32_t(4096));
      msg_bundle_mpi::set_capacity_bytes(mpi_bytes);
      uint32_t shm_bytes = havoqgt_getenv("HAVOQGT_MAILBOX_SHM_SIZE", uint32_t(4096));
      msg_bundle_shm::set_capacity_bytes(shm_bytes);
      s_b_print = havoqgt_getenv("HAVOQGT_MAILBOX_PRINT_STATS", bool(true));
      s_b_route_on_dest = havoqgt_getenv("HAVOQGT_MAILBOX_ROUTE_ON_DEST", bool(false));
      assert(s_num_isend > 0);

      //Quick print
      if(m_world_rank == 0) {
        std::cout << "HAVOQGT_MAILBOX_NUM_ISEND     = " << s_num_isend << std::endl;
        std::cout << "HAVOQGT_MAILBOX_NUM_IRECV     = " << s_num_irecv << std::endl;
        std::cout << "HAVOQGT_MAILBOX_MPI_SIZE      = " << mpi_bytes   << std::endl;
        std::cout << "HAVOQGT_MAILBOX_SHM_SIZE      = " << shm_bytes   << std::endl;
        std::cout << "HAVOQGT_MAILBOX_PRINT_STATS   = " << s_b_print << std::endl;
        std::cout << "HAVOQGT_MAILBOX_ROUTE_ON_DEST = " << s_b_route_on_dest << std::endl;

        std::cout << "msg_bundle_mpi::capacity = " << msg_bundle_mpi::capacity << std::endl;
        std::cout << "msg_bundle_shm::capacity = " << msg_bundle_shm::capacity << std::endl;
      }
    }
  }

  void init_managed_shared_memory() {
    //
    //Open Shared Mem segment
    size_t shm_file_size = msg_bundle_shm::padded_size * shm_exchange::capacity * m_shm_size * 4;

    m_shm_rank = comm_nl().rank();
    if(m_shm_rank == 0) {
      boost::interprocess::shared_memory_object::remove(msm_fname);
      m_pmsm = new managed_shared_memory( boost::interprocess::create_only, msm_fname, shm_file_size);
      m_poffset_exchange = m_pmsm->construct<offset_ptr<shm_exchange>>("exchange")[m_shm_size]();
    }
    comm_nl().barrier();
    if(m_shm_rank != 0) {
      m_pmsm = new managed_shared_memory( boost::interprocess::open_only, msm_fname );
      std::pair<offset_ptr<shm_exchange>*, managed_shared_memory::size_type> res;
      res = m_pmsm->find<offset_ptr<shm_exchange>>("exchange");
      assert(res.second == comm_nl().size());
      m_poffset_exchange = res.first;
    }
    comm_nl().barrier();
    void* aligned = m_pmsm->allocate_aligned(sizeof(shm_exchange), 4096);
    m_my_exchange = new (aligned) shm_exchange(m_shm_rank);
    m_poffset_exchange[m_shm_rank] = m_my_exchange;
    m_shm_transfer_count=0;
    m_shm_transfer_bytes=0;
    comm_nl().barrier();
    m_pp_exchange.resize(m_shm_size, nullptr);
    for(size_t i=0; i<m_shm_size; ++i) {
      m_pp_exchange[i] = m_poffset_exchange[i].get();
    }
    comm_nl().barrier();
  }

  void init_alloc_bundle_shared() {
    //
    // Allocate large chunk to help NUMA page pacement.
    // Local rank always touches pages before sending to other ranks.
    size_t num_to_preallocate = m_shm_size * 2;
    char* chunk = (char*) m_pmsm->allocate( msg_bundle_shm::padded_size * num_to_preallocate);
    for(size_t i=0; i<num_to_preallocate; ++i) {
      void* addr = chunk + i*msg_bundle_shm::padded_size;
      msg_bundle_shm* ptr = new (addr) msg_bundle_shm();
      assert(ptr->size == 0);
      m_free_shm_list.push_back(ptr);
    }
    //
    // Allocate per_rank data
    m_bundle_per_shm_rank.resize(m_shm_size, nullptr);
    m_pending_iterator_per_shm_rank.resize(m_shm_size, m_shm_pending_list.end());
  }

  void init_alloc_bundle_mpi() {
    size_t num_to_alloc = s_num_isend + s_num_irecv + m_mpi_size;
    m_free_mpi_list.reserve( num_to_alloc * 2);
    for(size_t i=0; i<num_to_alloc; ++i) {
      void* addr = nullptr;
      int ret = posix_memalign(&addr, 4096, msg_bundle_mpi::padded_size);
      if(ret != 0) {
        HAVOQGT_ERROR_MSG("posix_memalign failed");
      }
      m_free_mpi_list.push_back( new (addr) msg_bundle_mpi() );
    }
    //
    // Allocate per rank data
    m_bundle_per_mpi_rank.resize(m_mpi_size, nullptr);
    m_pending_iterator_per_mpi_rank.resize(m_mpi_size, m_mpi_pending_list.end());

    for(size_t i=0; i<s_num_irecv; ++i) {
      post_new_irecv(get_free_mpi_bundle());
    }
  }

  msg_bundle_mpi* get_free_mpi_bundle() {
    msg_bundle_mpi* to_return;
    if(m_free_mpi_list.empty()) {
      void* addr = nullptr;
      int ret = posix_memalign(&addr, 4096, msg_bundle_mpi::padded_size);
      if(ret != 0) {
        HAVOQGT_ERROR_MSG("posix_memalign failed");
      }
      to_return = new (addr) msg_bundle_mpi();
    } else {
      to_return = m_free_mpi_list.back();
      m_free_mpi_list.pop_back();
    }
    to_return->size = 0;
    return to_return;
  }

  void free_mpi_bundle(msg_bundle_mpi* tofree) {
    tofree->size = 0;
    m_free_mpi_list.push_back(tofree);
  }

  msg_bundle_shm* get_free_shm_bundle() {
    msg_bundle_shm* to_return;
    if(m_free_shm_list.empty()) {
      void* addr = (void*) m_pmsm->allocate( msg_bundle_shm::padded_size );
      to_return = new (addr) msg_bundle_shm();
    } else {
      to_return = m_free_shm_list.back();
      m_free_shm_list.pop_back();
    }
    to_return->size = 0;
    return to_return;
  }

  void free_shm_bundle(msg_bundle_shm* tofree) {
    tofree->size = 0;
    m_free_shm_list.push_back(tofree);
  }

public:

  size_t comm_size() const { return m_world_size; }
  size_t comm_rank() const { return m_world_rank; }

  template <typename OutputIterator>
  void bcast(T raw_msg, OutputIterator oitr) {
    raw_msg.vertex.set_bcast(true);
    msg_wrapper wrapped;
    wrapped.intercept = 0;
    wrapped.msg =  raw_msg;

    wrapped.bcast = 1;
    for(size_t i=0; i<m_shm_size; ++i) {
      if(i != m_shm_rank) {
        route_shm(i, wrapped, oitr);
      }
    }

    wrapped.bcast = 2;
    for(size_t i=0; i<m_mpi_size; ++i) {
      if(s_b_route_on_dest) {
        if(i%m_shm_size == m_shm_rank && i != m_mpi_rank) {
          route_mpi(i, wrapped, oitr);
        }
      } else {
        if(i != m_mpi_rank) {
          route_mpi(i, wrapped, oitr);
        }
      }
    }

    *oitr = raw_msg;
    ++oitr;
    hold_receive(oitr);
  }

  template <typename OutputIterator>
  void route_bcast(msg_wrapper wrapped, OutputIterator oitr) {
    size_t orig_bcast = wrapped.bcast;
    if(s_b_route_on_dest) {
      if(orig_bcast == 1) {
        wrapped.bcast = 2;
        for(size_t i=0; i<m_mpi_size; ++i) {
          if(i%m_shm_size == m_shm_rank && i != m_mpi_rank) {
            route_mpi(i, wrapped, oitr);
          }
        }
      } else if(orig_bcast == 2) {
        wrapped.bcast = 3;
        for(size_t i=0; i<m_shm_size; ++i) {
          if(i != m_shm_rank) {
            route_shm(i, wrapped, oitr);
          }
        }
      } else {
        if(orig_bcast != 3) HAVOQGT_ERROR_MSG("orig_bcast != 3");
      }
    } else {
      if(orig_bcast == 1) {
        wrapped.bcast = 2;
        for(size_t i=0; i<m_mpi_size; ++i) {
          if(i != m_mpi_rank) {
            route_mpi(i, wrapped, oitr);
          }
        }
      } else {
        if(orig_bcast != 2) HAVOQGT_ERROR_MSG("orig_bcast != 2");
      }
    }
    *oitr = wrapped.msg;
    ++oitr;
  }

  template <typename OutputIterator>
  void send(size_t world_dest, const T& raw_msg, OutputIterator oitr, bool intercept) {
    //std::cout << whoami() << " send to " << world_dest << std::endl;
    ++m_send_recv_balance;
    if(world_dest == m_world_rank) {
      *oitr = raw_msg;
      ++oitr;
    }

    //assume block partitioning
    msg_wrapper wrapped;
    wrapped.dest_node = world_dest / m_shm_size;
    wrapped.dest_core = world_dest % m_shm_size;
    wrapped.bcast = 0;
    wrapped.msg = raw_msg;
    wrapped.intercept = intercept;

    bool b_sent = false;

    if(s_b_route_on_dest)  {
      if(wrapped.dest_node == m_mpi_rank) {
        assert(wrapped.dest_core != m_shm_rank);
        b_sent = route_shm(wrapped.dest_core, wrapped, oitr);
      } else {
        size_t first_core_route = wrapped.dest_node % m_shm_size;
        if(first_core_route != m_shm_rank) {
          b_sent = route_shm(first_core_route, wrapped, oitr);
        } else {
          b_sent = route_mpi(wrapped.dest_node, wrapped, oitr);
        }
      }
    } else {
      if(wrapped.dest_core != m_shm_rank) {
        b_sent = route_shm(wrapped.dest_core, wrapped, oitr);
      } else if (wrapped.dest_node != m_mpi_rank) {
        b_sent = route_mpi(wrapped.dest_node, wrapped, oitr);
      } else {
        HAVOQGT_ERROR_MSG("Logic Error");
      }
    }

    /*if(b_sent) {
      do{
        receive(oitr);
        cleanup_pending_isend_requests();
      } while(shm_transfer_slots_full() || isend_slots_full());
    }*/
    hold_receive(oitr);
  }

  template <typename OutputIterator>
  void hold_receive(OutputIterator oitr) {
    while(shm_transfer_slots_full() || isend_slots_full()) {
      receive(oitr);
      cleanup_pending_isend_requests();
    }
  }

  template <typename OutputIterator>
  void receive(OutputIterator oitr) {
    receive_mpi(oitr);
    receive_shm(oitr);
  }

  bool is_idle() {
    cleanup_pending_isend_requests();
    return m_mpi_pending_list.empty() && m_isend_request_list.empty() && m_shm_pending_list.empty()
           && !m_my_exchange->probe() && shm_count_transfers_pending() == 0;
  }

  void flush_buffers_if_idle() {
    //static double last_time = 0;
    //double new_time = MPI_Wtime();
    //if(new_time - last_time < 1e-6) return;
    static bool flopper = false;
    cleanup_pending_isend_requests();
    if(m_isend_request_list.empty() && shm_count_transfers_pending() == 0 && !m_my_exchange->probe()) {
      flopper = !flopper;
      if(flopper) {
        if(!m_shm_pending_list.empty()) {
          try_transfer_shm(m_shm_pending_list.front());
          //last_time = new_time;
        }
      } else {
        if(!m_mpi_pending_list.empty()) {
          post_isend(m_mpi_pending_list.front());
          //last_time = new_time;
        }
      }
    }
/*

    if(m_isend_request_list.empty() && !m_shm_pending_list.empty() && shm_count_transfers_pending() == 0 && !m_my_exchange->probe()) {
      bool sent = try_transfer_shm(m_shm_pending_list.front());
      //if(sent) return;
    }

    cleanup_pending_isend_requests();

    if(m_isend_request_list.empty() && !m_mpi_pending_list.empty()) {
      if(!m_mpi_pending_list.empty()) {
        post_isend(m_mpi_pending_list.front());
      }
    }*/
  }

  uint64_t shm_count_transfers_pending() { return m_shm_request_list.size(); }

private:

  template <typename OutputIterator>
  bool route_shm(size_t shm_rank, msg_wrapper& wrapped, OutputIterator oitr) {
    assert(shm_rank != m_shm_rank);
    //std::cout << whoami() << " route_shm to rank " << shm_rank << std::endl;
    bool to_return = false;
    if(m_bundle_per_shm_rank[shm_rank] == nullptr) {
      bool print = true;
      msg_bundle_shm* bundle = get_free_shm_bundle();
      //std::cout << whoami() << " route_shm size = " << bundle->size << ", source_core = " << bundle->get_source_core() << std::endl;
      m_bundle_per_shm_rank[shm_rank] = bundle;
      assert(m_bundle_per_shm_rank[shm_rank]->size == 0);
      assert(m_pending_iterator_per_shm_rank[shm_rank] == m_shm_pending_list.end());
      m_shm_pending_list.push_back(shm_rank);
      m_pending_iterator_per_shm_rank[shm_rank] = --(m_shm_pending_list.end());
    }
    size_t size = m_bundle_per_shm_rank[shm_rank]->push_back(wrapped);

    if(size == msg_bundle_shm::capacity) {
      to_return = true;
      /*if(!fast_path && shm_rank != m_shm_pending_list.front()) {
         while(!try_transfer_shm( m_shm_pending_list.front())) {
          if(fast_path) {
            HAVOQGT_ERROR_MSG("route_shm try_transfer_shm failed");
          } else {
            std::cout << "Add me up 1" << std::endl;
            receive_shm(oitr);
            receive_mpi(oitr);
          }
        }
      }*/
      //double stuck_time = 0.0;
      while(!try_transfer_shm(shm_rank)) {
        //if(stuck_time == 0.0) stuck_time = MPI_Wtime();
        //receive_shm(oitr);
        //receive_mpi(oitr);
      }
      /*if(stuck_time > 0.0) {
        stuck_time = MPI_Wtime() - stuck_time;
        if(stuck_time > 0.5) {
          std::cout << whoami() << " is stuck in route_shm, try_transfer_shm for "<< stuck_time << " seconds"  << std::endl;
        }
      }*/
    }
    return to_return;
  }
  bool shm_transfer_slots_full() {
    uint64_t pending = shm_count_transfers_pending();
    return pending > s_num_isend;
  }

  bool isend_slots_full() {
    return m_isend_request_list.size() > s_num_isend;
  }

  template <typename OutputIterator>
  bool route_mpi(size_t mpi_rank, msg_wrapper& wrapped, OutputIterator oitr) {
    assert(mpi_rank != m_mpi_rank);
    //std::cout << whoami() << " route_mpi sending to " << mpi_rank << std::endl;
    bool to_return = false;
    if(m_bundle_per_mpi_rank[mpi_rank] == nullptr) {
      m_bundle_per_mpi_rank[mpi_rank] = get_free_mpi_bundle();
      assert(m_bundle_per_mpi_rank[mpi_rank]->size == 0);
      assert(m_pending_iterator_per_mpi_rank[mpi_rank] == m_mpi_pending_list.end());
      m_mpi_pending_list.push_back(mpi_rank);
      m_pending_iterator_per_mpi_rank[mpi_rank] = --(m_mpi_pending_list.end());
    }
    size_t size = m_bundle_per_mpi_rank[mpi_rank]->push_back(wrapped);

    if(size == msg_bundle_mpi::capacity) {
      to_return = true;
      /*if(!fast_path && mpi_rank != m_mpi_pending_list.front()) {
        post_isend(m_mpi_pending_list.front()); //fair sending
      }*/
      post_isend(mpi_rank);
    }
    return to_return;
  }

  bool try_transfer_shm(size_t rank) {
    assert(*(m_pending_iterator_per_shm_rank[rank]) == rank);
    msg_bundle_shm* to_transfer = m_bundle_per_shm_rank[rank];
    assert(to_transfer->size > 0);
    if(m_pp_exchange[rank]->try_send(to_transfer)) {
      m_bundle_per_shm_rank[rank] = nullptr;
      m_shm_pending_list.erase(m_pending_iterator_per_shm_rank[rank]);
      m_pending_iterator_per_shm_rank[rank] = m_shm_pending_list.end();
      m_shm_transfer_count++;
      m_shm_transfer_bytes+=to_transfer->size * sizeof(msg_wrapper);
      m_shm_request_list.push_back(to_transfer);
      return true;
    }
    //std::cout << whoami() << " transfer_shm FAILED dest_core = " << rank << std::endl;
    return false;
  }

  void post_isend(size_t rank) {
    assert(*(m_pending_iterator_per_mpi_rank[rank]) == rank);
    assert(m_pending_iterator_per_mpi_rank[rank] != m_mpi_pending_list.end());
    msg_bundle_mpi* to_transfer = m_bundle_per_mpi_rank[rank];
    m_bundle_per_mpi_rank[rank] = nullptr;
    std::pair<MPI_Request, msg_bundle_mpi*> req_pair;
    MPI_Request* request_ptr=&req_pair.first;
    req_pair.second = to_transfer;
    void* buffer_ptr = (void*) to_transfer;
    size_t size_in_bytes = to_transfer->message_size();
    //std::cout << whoami() << " posting ISend to node " << rank << std::endl;
    //std::cout << "Inspecting contents, byte_size = " << size_in_bytes << " size = " << to_transfer->size << ", first hop_count = " << to_transfer->data[0].msg.m_hop_count << std::endl;
    size_t world_dest;
    if(s_b_route_on_dest) {
      world_dest = (rank * m_shm_size)+(m_mpi_rank % m_shm_size);
    } else {
      world_dest = (rank * m_shm_size) + m_shm_rank;
    }
    if(true) {//m_isend_count % 2 == 0) {
      //CHK_MPI( MPI_Issend(buffer_ptr, size_in_bytes, MPI_BYTE, rank, m_mpi_tag, m_mpi_comm, request_ptr) );
      CHK_MPI( MPI_Issend(buffer_ptr, size_in_bytes, MPI_BYTE, world_dest, m_mpi_tag, MPI_COMM_WORLD, request_ptr) );
    } else {
      //CHK_MPI( MPI_Isend(buffer_ptr, size_in_bytes, MPI_BYTE, rank, m_mpi_tag, m_mpi_comm, request_ptr) );
      CHK_MPI( MPI_Isend(buffer_ptr, size_in_bytes, MPI_BYTE, world_dest, m_mpi_tag, MPI_COMM_WORLD, request_ptr) );
    }
    m_isend_count++;
    m_isend_bytes+=size_in_bytes;

    m_isend_request_list.push_back(req_pair);
    m_mpi_pending_list.erase(m_pending_iterator_per_mpi_rank[rank]);
    m_pending_iterator_per_mpi_rank[rank] = m_mpi_pending_list.end();
    //std::cout << whoami() << " Inspecting m_mpi_pending_list:  size = " << m_mpi_pending_list.size()  << std::endl;
    /*for(auto itr = m_mpi_pending_list.begin(); itr != m_mpi_pending_list.end(); ++itr) {
      //std::cout << *itr << std::endl;
    }*/
    //std::cout << "Inspecting m_isend_request_list: size = " << m_isend_request_list.size() << std::endl;
  }

  template <typename OutputIterator>
  void receive_mpi(OutputIterator oitr) {
    m_send_recv_balance = 0;
    int flag(0);
    //do {
    MPI_Request* request_ptr = &(m_irecv_request_list.front().first);
    CHK_MPI( MPI_Test( request_ptr, &flag, MPI_STATUS_IGNORE) );
    if(flag) {
      msg_bundle_mpi* ptr = m_irecv_request_list.front().second;
      for(size_t i=0; i<ptr->size; ++i) {
        if(ptr->data[i].bcast > 0) {
            route_bcast(ptr->data[i], oitr);
        } else  if(s_b_route_on_dest) {
          assert(ptr->data[i].dest_node == m_mpi_rank);
          if(ptr->data[i].dest_core == m_shm_rank) {
            *oitr = ptr->data[i].msg;
            ++oitr;
          } else {
            if(ptr->data[i].intercept) {
              if(!oitr.intercept(ptr->data[i].msg)) { continue; }
            }
            route_shm(ptr->data[i].dest_core, ptr->data[i], oitr);
          }
        } else {
          assert(ptr->data[i].dest_core == m_shm_rank);
          assert(ptr->data[i].dest_node == m_mpi_rank);
          *oitr = ptr->data[i].msg;
          ++oitr;
        }
      }
      m_irecv_request_list.pop_front();
      ptr->size = 0;
      post_new_irecv(ptr);
    }
    //} while(flag);
  }

  template <typename OutputIterator>
  void receive_shm(OutputIterator oitr) {
    m_send_recv_balance = 0;
    std::vector<msg_bundle_shm*> to_recv;
    m_my_exchange->try_recv(to_recv);
    for(auto itr = to_recv.begin(); itr != to_recv.end(); ++itr) {


#if 0
      msg_bundle_shm* recvptr = *itr;
      std::vector<msg_wrapper> copy(recvptr->data, recvptr->data+recvptr->size);
      m_pp_exchange[recvptr->source_core]->free(recvptr);


    //if(recvptr != nullptr) {
      for(int i=0; i</*recvptr->size*/copy.size(); ++i) {
        if(s_b_route_on_dest) {
          if(/*recvptr->data*/copy[i].dest_node == m_mpi_rank && /*recvptr->data*/copy[i].dest_core == m_shm_rank) {
            *oitr = /*recvptr->data*/copy[i].msg;
            ++oitr;
          } else {
            if(copy[i].dest_node == m_mpi_rank) {
              std::cout << "logic problem!!!" << std::endl;
            }
            route_mpi(/*recvptr->data*/copy[i].dest_node, /*recvptr->data*/copy[i], oitr);
          }
        } else {
          assert(/*recvptr->data*/copy[i].dest_core == m_shm_rank);
          if(/*recvptr->data*/copy[i].dest_node == m_mpi_rank) {
            *oitr = /*recvptr->data*/copy[i].msg;
            ++oitr;
          } else {
            route_mpi(/*recvptr->data*/copy[i].dest_node, /*recvptr->data*/copy[i], oitr);
          }
        }
      }
      //m_pp_exchange[recvptr->source_core]->free(recvptr);
    //}

#else
      msg_bundle_shm* recvptr = *itr;
      assert(recvptr != nullptr);
      assert(recvptr->size > 0);
      for(int i=0; i<recvptr->size; ++i) {
        if(recvptr->data[i].bcast > 0) {
          route_bcast(recvptr->data[i], oitr);
        } else if(s_b_route_on_dest) {
          if(recvptr->data[i].dest_node == m_mpi_rank/* && recvptr->data[i].dest_core == m_shm_rank*/) {
            if(recvptr->data[i].dest_core != m_shm_rank) {
              //std::cout << whoami() << "recvptr->data[i].dest_core != m_shm_rank: dest = " << recvptr->data[i].msg.dest() << ", dest_node = " << recvptr->data[i].dest_node << ", dest_core = " << recvptr->data[i].dest_core << ", hop length = " << recvptr->data[i].msg.m_hop_count << std::endl;
            }
            assert(recvptr->data[i].dest_core == m_shm_rank);
            *oitr = recvptr->data[i].msg;
            ++oitr;
          } else {
            if(recvptr->data[i].intercept) {
              if(!oitr.intercept(recvptr->data[i].msg)) { continue; }
            }
            route_mpi(recvptr->data[i].dest_node, recvptr->data[i], oitr);
          }
        } else {
          assert(recvptr->data[i].dest_core == m_shm_rank);
          if(recvptr->data[i].dest_node == m_mpi_rank) {
            *oitr = recvptr->data[i].msg;
            ++oitr;
          } else {
            if(recvptr->data[i].intercept) {
              if(!oitr.intercept(recvptr->data[i].msg)) { continue; }
            }
            route_mpi(recvptr->data[i].dest_node, recvptr->data[i], oitr);
          }
        }
      }
      recvptr->mark_read();
 #endif
    }
  }



  void cleanup_pending_isend_requests() {
    while(!m_shm_request_list.empty()) {
      if(m_shm_request_list.front()->is_read()) {
        free_shm_bundle(m_shm_request_list.front());
        m_shm_request_list.pop_front();
      } else {
        break;
      }
    }

    while(!m_isend_request_list.empty()) {
      int flag(0);
      MPI_Request* request_ptr = &(m_isend_request_list.front().first);
      CHK_MPI( MPI_Test( request_ptr, &flag, MPI_STATUS_IGNORE) );
      if(flag) {
        free_mpi_bundle(m_isend_request_list.front().second);
        m_isend_request_list.pop_front();
      } else {
        break;
      }
    }
  }

  void post_new_irecv(msg_bundle_mpi* buff) {
    std::pair<MPI_Request, msg_bundle_mpi*> irecv_req;
    irecv_req.second = buff;
    MPI_Request* request_ptr = &irecv_req.first;
    //CHK_MPI( MPI_Irecv( (void*) buff, msg_bundle_mpi::padded_size, MPI_BYTE, MPI_ANY_SOURCE, m_mpi_tag, m_mpi_comm, request_ptr) );
    CHK_MPI( MPI_Irecv( (void*) buff, msg_bundle_mpi::padded_size, MPI_BYTE, MPI_ANY_SOURCE, m_mpi_tag, MPI_COMM_WORLD, request_ptr) );
    m_irecv_request_list.push_back(irecv_req);
  }


private:

  int m_world_rank;
  int m_world_size;

  //
  // MPI Related
  int m_mpi_tag;
  MPI_Comm m_mpi_comm;
  int m_mpi_size;
  int m_mpi_rank;
  std::vector<msg_bundle_mpi*>                        m_bundle_per_mpi_rank;
  std::vector<msg_bundle_mpi*>                        m_free_mpi_list;
  std::vector< std::list<size_t>::iterator >          m_pending_iterator_per_mpi_rank;
  std::list<size_t>                                   m_mpi_pending_list;
  std::deque< std::pair<MPI_Request, msg_bundle_mpi*>> m_irecv_request_list;
  std::deque< std::pair<MPI_Request, msg_bundle_mpi*>> m_isend_request_list;


  //
  // Shared Mem Related
  managed_shared_memory*    m_pmsm;
  offset_ptr<shm_exchange>* m_poffset_exchange;
  std::vector<shm_exchange*> m_pp_exchange;
  shm_exchange*             m_my_exchange;
  int                       m_shm_rank;
  int                       m_shm_size;
  std::vector<msg_bundle_shm*>                        m_bundle_per_shm_rank;
  std::vector< std::list<size_t>::iterator >          m_pending_iterator_per_shm_rank;
  std::list<size_t>                                   m_shm_pending_list;
  uint64_t                  m_shm_transfer_count;
  uint64_t                  m_shm_transfer_bytes;
  std::vector<msg_bundle_shm*>                        m_free_shm_list;
  std::deque< msg_bundle_shm* >                       m_shm_request_list;

  //
  // Satic Configs
  static uint32_t s_num_isend;
  static uint32_t s_num_irecv;
  static bool     s_b_print;
  static bool s_b_route_on_dest;

  //
  //
  uint64_t m_send_recv_balance;
  uint64_t m_isend_count;
  uint64_t m_isend_bytes;


};

template <typename T>
uint32_t mailbox<T>::s_num_isend = 0;
template <typename T>
uint32_t mailbox<T>::s_num_irecv = 0;
template <typename T>
bool     mailbox<T>::s_b_print = 0;
template <typename T>
bool     mailbox<T>::s_b_route_on_dest = 0;

template<typename T>
struct mailbox<T>::shm_exchange {
  using mutex = boost::interprocess::interprocess_mutex;
  using scoped_lock = boost::interprocess::scoped_lock<mutex>;
  static const uint32_t capacity = 2048;
  volatile uint64_t recv_end;
  char pad2[65-sizeof(recv_end)];
  volatile uint64_t recv_beg;
  uint64_t source_core;
  mutex recv_mutex;
  offset_ptr<msg_bundle_shm> recv_list[capacity];
  shm_exchange(uint32_t core) {
    recv_end = 0;
    recv_beg = 0;
    source_core = core;
    for(size_t i=0; i<capacity; ++i) {
      recv_list[i] = offset_ptr<msg_bundle_shm>();
    }
  }

#if 0
  bool try_send(msg_bundle_shm* to_send) {
    assert(to_send->get_source_core() == comm_nl().rank());
    //std::cout << whoami() << " sending to " << source_core << ", size = " << to_send->size << std::endl;
    //double lock_start = MPI_Wtime();
    scoped_lock lock(recv_mutex, boost::interprocess::try_to_lock);
    bool print = false;
    while(!lock) {
      if(print) {
        print = false;
        std::cout << whoami() << " lock failed" << std::endl;
      }
      lock.try_lock();
      //sched_yield();
    }
    //lock_time += MPI_Wtime() - lock_start;
    if(/*size_recv < capacity*/ /*recv_end >= recv_beg &&*/ recv_end - recv_beg < capacity) {
      recv_list[recv_end++ % capacity] = to_send;
      return true;
    }
    //std::cout << whoami() << ": Try_send Failed, sending_to_core = " << source_core << ", recv_end = " << recv_end << ", recv_beg = " << recv_beg << std::endl;
    return false;
  }
#else
  bool try_send(msg_bundle_shm* to_send) {
    __sync_synchronize();
    if(recv_end - recv_beg >= capacity-100) {
      return false;
    }
    uint64_t pos = __sync_fetch_and_add(&recv_end, 1);
    pos %= capacity;
    if(recv_list[pos] == offset_ptr<msg_bundle_shm>()) {
      recv_list[pos] = to_send;
      __sync_synchronize();
      return true;
    }
    std::cout << "recv_end = " << recv_end << ", recv_beg = " << recv_beg << std::endl;
    HAVOQGT_ERROR_MSG("Failed to try_send");
  }

#endif
  bool probe() {
    __sync_synchronize();
    return recv_end > recv_beg;
  }
#if 0
  void try_recv(std::vector<msg_bundle_shm*>& to_recv) {
    //double lock_start = MPI_Wtime();
    to_recv.clear();
    __sync_synchronize();
    if(recv_end == recv_beg) return;
    //to_recv.reserve(recv_end - recv_beg);
    scoped_lock lock(recv_mutex, boost::interprocess::try_to_lock);
    while(!lock) {
      lock.try_lock();
    }
    //lock_time += MPI_Wtime() - lock_start;
    while(/*size_recv > 0*/recv_end > recv_beg) {
      /*//std::cout << whoami() << " recv size = " << size_recv << std::endl;
      to_recv.resize(size_recv);
      for(size_t i=0; i<size_recv; ++i) {
        to_recv[i] = recv_list[i].get();
        //std::cout << "recving from: " << recv_list[i]->get_source_core() << ", size = " << recv_list[i]->size <<", size = " << recv_list[i]->size <<  std::endl;
        //std::cout << "recving from: " << to_recv[i]->get_source_core() << ", size = " << to_recv[i]->size <<", size = " << to_recv[i]->size <<  std::endl;
      }
      size_recv = 0;
      */
      //return recv_list[/*--size_recv*/recv_beg++ % capacity].get();
      to_recv.push_back(recv_list[recv_beg++ % capacity].get());
    }
    //return nullptr;
  }
#else
  void try_recv(std::vector<msg_bundle_shm*>& to_recv) {
    to_recv.clear();
    __sync_synchronize();
    to_recv.reserve(recv_end - recv_beg);
    /*while*/if(recv_end > recv_beg) { // data flow prob???
      uint64_t pos = (recv_beg) % capacity;
      if(recv_list[pos] != offset_ptr<msg_bundle_shm>()) {
        to_recv.push_back(recv_list[pos].get());
        recv_list[pos] = offset_ptr<msg_bundle_shm>();
        ++recv_beg;
        __sync_synchronize();
      } else {
        //break; put back in for while
      }
    }
  }
#endif
};


template<typename T>
struct mailbox<T>::msg_wrapper {
  uint64_t dest_node : 16;
  uint64_t dest_core : 8;
  uint64_t bcast     : 2;
  uint64_t intercept : 1;
  T        msg;
};

template<typename T>
struct mailbox<T>::msg_bundle_shm {
  // Data Members
  volatile uint32_t size;
  msg_wrapper data[0];
  msg_bundle_shm()
    : size(0)
  {
  }
  // Static Members & Functions
  static uint32_t capacity;
  static uint32_t padded_size;
  static void set_capacity_bytes(uint32_t bytes) {
    msg_bundle_shm::padded_size = bytes;
    msg_bundle_shm::capacity = (padded_size - sizeof(msg_bundle_mpi)) / sizeof(msg_wrapper);
  }
  size_t push_back(const msg_wrapper& _d) {
    data[size] = _d;
    return ++size;
  }

  void mark_read() {
    __atomic_store_n(&size, 0, __ATOMIC_SEQ_CST);
    /*assert(size > 0);
    __sync_synchronize();
    size = 0;
    __sync_synchronize();*/
  }

  bool is_read() {
    return __atomic_load_n(&size,__ATOMIC_SEQ_CST) == 0;
    /*__sync_synchronize();
    bool to_return = (size == 0);
    __sync_synchronize();
    return to_return;*/

  }

};

template <typename T>
uint32_t mailbox<T>::msg_bundle_shm::capacity = 0;

template <typename T>
uint32_t mailbox<T>::msg_bundle_shm::padded_size = 0;

template<typename T>
struct mailbox<T>::msg_bundle_mpi {
  uint32_t size;
  msg_wrapper data[0];
  msg_bundle_mpi()
    : size(0)
  { }
  // Static Members & Functions
  static uint32_t capacity;
  static uint32_t padded_size;
  static void set_capacity_bytes(uint32_t bytes) {
    msg_bundle_mpi::padded_size = bytes;
    msg_bundle_mpi::capacity = (padded_size - sizeof(msg_bundle_shm)) / sizeof(msg_wrapper); //matches msg_bundle_shm's
  }
  size_t message_size() { return sizeof(msg_wrapper[size]) + sizeof(msg_bundle_mpi); }
  size_t push_back(const msg_wrapper& _d) {
    data[size] = _d;
    return ++size;
  }
};

template <typename T>
uint32_t mailbox<T>::msg_bundle_mpi::capacity = 0;

template <typename T>
uint32_t mailbox<T>::msg_bundle_mpi::padded_size = 0;

} //namespace havoqgt

#endif //HAVOQGT_MAILBOX_HPP_INCLUDED

