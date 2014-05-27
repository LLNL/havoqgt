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


#ifndef HAVOQGT_OMP_MAILBOX_HPP_INCLUDED
#define HAVOQGT_OMP_MAILBOX_HPP_INCLUDED

#include <omp.hpp>
#include <vector>
#include <limits>
#include <list>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>

/**
 * @todo Make aggregation buffer size an environment parameter
 * @todo Remove std::list and manage order by intrusive list in bucket structure
 */

template<typename T>
class page_aligned_allocator {
public : 
    //    typedefs
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

public : 
    //    convert an page_aligned_allocator<T> to page_aligned_allocator<U>
    template<typename U>
    struct rebind {
        typedef page_aligned_allocator<U> other;
    };

public : 
    inline explicit page_aligned_allocator() {}
    inline ~page_aligned_allocator() {}
    inline explicit page_aligned_allocator(page_aligned_allocator const&) {}
    template<typename U>
    inline explicit page_aligned_allocator(page_aligned_allocator<U> const&) {}

    //    address
    inline pointer address(reference r) { return &r; }
    inline const_pointer address(const_reference r) { return &r; }

    //    memory allocation
    inline pointer allocate(size_type cnt, void* = 0 ) { 
        void* mem = 0;
        if (posix_memalign(&mem, 4096, sizeof(T) * cnt) != 0)
        {
            throw std::bad_alloc(); // or something
        }
        memset(mem,0,1);
        return (pointer) mem;
    }
    inline void deallocate(pointer p, size_type) { 
      free(p);
    }

    //    size
    inline size_type max_size() const { 
        return std::numeric_limits<size_type>::max() / sizeof(T);
 }

    //    construction/destruction
    inline void construct(pointer p, const T& t) { new(p) T(t); }
    inline void destroy(pointer p) { p->~T(); }

    inline bool operator==(page_aligned_allocator const&) { return true; }
    inline bool operator!=(page_aligned_allocator const& a) { return !operator==(a); }
};    //    end of class page_aligned_allocator 

namespace havoqgt { namespace omp {


const int mailbox_buffer_bytes = 4*1024;

namespace detail {
  template <typename T>
  class mailbox_tp {
  public:
    class bucket {
    public:
      T* ptr;
      uint16_t size;
      uint16_t owner;
      uint16_t next;
      uint16_t prev;
      void push_back(const T& data) {
        assert(size < mailbox_buffer_bytes / sizeof(T));
        ptr[size] = data;
        ++size;
      }
      bool full() const { return size == mailbox_buffer_bytes/sizeof(T); }
      void clear() { size = 0; next = -1; prev = -1; }
    };
  public:
    mailbox_tp() 
    {
      assert(omp_in_parallel() == 1);
      m_send_buckets.reserve( omp_get_num_threads() );
      m_pending_position.resize( omp_get_num_threads(), m_pending_list.end());
      m_to_receive_buckets.reserve(omp_get_num_threads() * 2);
      m_received_buckets.reserve( omp_get_num_threads() * 2);
      m_free_buckets.reserve( omp_get_num_threads() * 2);
      for(size_t i=0; i<omp_get_num_threads(); ++i) {
        void* data;
        posix_memalign( &data, 4096, mailbox_buffer_bytes );
        bzero(data,mailbox_buffer_bytes);
        bucket b;
        b.ptr = (T*) data;
        b.size = 0;
        b.owner = omp_get_thread_num();
        m_send_buckets.push_back(b);
      }
      for(size_t i=0; i<omp_get_num_threads()*2; ++i) {
        void* data;
        posix_memalign( &data, 4096, mailbox_buffer_bytes );
        bzero(data,mailbox_buffer_bytes);
        bucket b;
        b.ptr = (T*) data;
        b.size = 0;
        b.owner = omp_get_thread_num();
        m_free_buckets.push_back(b);
      }
      omp_init_lock(&m_receive_lock);
      omp_init_lock(&m_free_lock);
    }

    ~mailbox_tp() {
      omp_destroy_lock(&m_receive_lock);
      omp_destroy_lock(&m_free_lock);
      assert(m_to_receive_buckets.empty());
      assert(m_received_buckets.empty());
      for(size_t i=0; i<m_send_buckets.size(); ++i) {
        free(m_send_buckets[i].ptr);
      }
      for(size_t i=0; i<m_free_buckets.size(); ++i) {
        free(m_free_buckets[i].ptr);
      }
    }

    template <typename Mailbox>
    void send(size_t dest, const T& data, Mailbox* ptr_mailbox) {
      if(m_send_buckets[dest].size == 0) {
        m_pending_list.push_back(dest);
        m_pending_position[dest] = --(m_pending_list.end());
      }
      m_send_buckets[dest].push_back(data);
      if(m_send_buckets[dest].full()) { 
        send_bucket(dest, ptr_mailbox);
      }
    }

    template <typename Mailbox>
    void idle_send(Mailbox* ptr_mailbox ) {
      if(!m_pending_list.empty()) {
        size_t first_dest = m_pending_list.front();
        send_bucket(first_dest, ptr_mailbox);
      }
    }

    bool is_idle() { return m_pending_list.empty(); }

    template <typename Mailbox>
    void send_bucket(size_t dest, Mailbox* ptr_mailbox) {
        assert(m_pending_position[dest] != m_pending_list.end());
        assert(*m_pending_position[dest] == dest);
        bool first = dest == m_pending_list.front();
        m_pending_list.erase(m_pending_position[dest]);
        m_pending_position[dest] = m_pending_list.end();
        // send the full bucket
        ptr_mailbox->move_bucket(dest,m_send_buckets[dest]);

        // grab free bucket for sending
        bool found_free = false;
        do {
          receive(ptr_mailbox);  //Should we receive w/ ever bucket send?
          omp_set_lock(&m_free_lock);
          if(!m_free_buckets.empty()) {
            found_free = true;
            m_send_buckets[dest] = m_free_buckets.back();
            m_free_buckets.pop_back();
          } 
          omp_unset_lock(&m_free_lock);
          //if(!found_free) {
          //  receive(ptr_mailbox);
          //}
        } while(!found_free);

        if(!first) { //then send first
          size_t first_dest = m_pending_list.front();
          send_bucket(first_dest, ptr_mailbox);
        }
    }

    template <typename Mailbox>
    void receive(Mailbox* ptr_mailbox) {
      omp_set_lock(&m_receive_lock);
      m_received_buckets.swap(m_to_receive_buckets);
      omp_unset_lock(&m_receive_lock);
      //process m_received_buckets
      for(size_t i=0; i<m_received_buckets.size(); ++i) {
        for(size_t j=0; j<m_received_buckets[i].size; ++j) {
          // do something with m_received_buckets[i][j]
          ptr_mailbox->mailbox_receive(m_received_buckets[i].ptr[j]);
        }
        ptr_mailbox->free_bucket(m_received_buckets[i]);
      }
      m_received_buckets.clear();
    }

    //this is where an another thread pushes buckets m_to_receive_buckets
    void receive_bucket(const bucket& b) {
      omp_set_lock(&m_receive_lock);
      m_to_receive_buckets.push_back(b);
      omp_unset_lock(&m_receive_lock);
    }

    void free_bucket(bucket& b) {
      b.size=0;
      omp_set_lock(&m_free_lock);
      m_free_buckets.push_back(b);
      omp_unset_lock(&m_free_lock);
    }
  private:
    std::vector< bucket , page_aligned_allocator<bucket> > m_send_buckets;
    std::list<size_t> m_pending_list;
    std::vector< std::list<size_t>::iterator , page_aligned_allocator<std::list<size_t>::iterator> > m_pending_position;
    //uint64_t m_fair_sender __attribute__ ((aligned (64)));


    // Lock manages following three vectors
    omp_lock_t m_receive_lock __attribute__ ((aligned (64)));
    std::vector< bucket, page_aligned_allocator<bucket> > m_to_receive_buckets __attribute__ ((aligned (64)));
    std::vector< bucket, page_aligned_allocator<bucket> > m_received_buckets __attribute__ ((aligned (64)));
    omp_lock_t m_free_lock __attribute__ ((aligned (64)));
    std::vector< bucket, page_aligned_allocator<bucket> > m_free_buckets __attribute__ ((aligned (64)));
  };
}

template <typename T, typename H>
class mailbox {
friend class detail::mailbox_tp<T>;
public:
  mailbox(H* ptr_h) 
    : m_vec_tp(num_threads())
    , m_receive_handler(ptr_h)
  {
    //std::cout << "sizeof(detail::mailbox_tp<T>) = " << sizeof(detail::mailbox_tp<T>) << std::endl;
    assert(omp_in_parallel() == 0);
    #pragma omp parallel
    {
      void* data;
      posix_memalign(&data, 4096, sizeof(detail::mailbox_tp<T>));
      m_vec_tp[thread_num()] = new (data) detail::mailbox_tp<T>();
    }
  }

  ~mailbox() {
    assert(omp_in_parallel() == 0);
    #pragma omp parallel
    {
      m_vec_tp[thread_num()]->~mailbox_tp();
      free(m_vec_tp[thread_num()]);
    }
  }

  void send(size_t dest, const T& data) 
  {
    m_vec_tp[thread_num()]->send(dest, data, this);
  }

  void receive() {
    m_vec_tp[thread_num()]->receive(this);
  }

  void idle_receive() {
    m_vec_tp[thread_num()]->idle_send(this);
    m_vec_tp[thread_num()]->receive(this);
  }

  bool is_idle() {
    return m_vec_tp[thread_num()]->is_idle();
  }
protected:
  void move_bucket(size_t dest, typename detail::mailbox_tp<T>::bucket& b) {
    m_vec_tp[dest]->receive_bucket(b);
  }

  void free_bucket(typename detail::mailbox_tp<T>::bucket& b) {
    m_vec_tp[b.owner]->free_bucket(b);
  }

  void mailbox_receive(const T& data) {
    m_receive_handler->mailbox_receive(data);
  }

private:
  std::vector< detail::mailbox_tp<T>* > m_vec_tp;
  H* m_receive_handler;
};

} }



#endif
