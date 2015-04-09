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

#include <havoqgt/environment.hpp>
#include <havoqgt/mpi.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/mailbox.hpp>
#include <havoqgt/termination_detection.hpp>
#include <havoqgt/detail/hash.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <assert.h>



#include <random>
#include <iostream>

using namespace havoqgt;

class msg_hopper {
public:
  msg_hopper(uint32_t seed) {
    m_hop_count = 1;
    m_hop_state = seed;
    do {
     uint32_t old_hop_state = m_hop_state;
     m_hop_state = detail::hash32(old_hop_state);
     if(old_hop_state == m_hop_state) {
       m_hop_state++;
     }
     m_dest = m_hop_state % s_mpi_size;
    } while(m_dest ==  havoqgt_env()->world_comm().rank());  
  }

  msg_hopper() {}


  uint32_t dest() const { return m_dest; }
  uint32_t get_bcast() const { return false; }
  void set_bcast(uint32_t bcast) { }
  void set_dest(uint32_t dest) { }
  bool is_intercept() const { return false; }
  void set_intercept(bool) { }

  bool finished() const { return m_hop_count >= s_hop_limit; }

  msg_hopper next() const { 
    msg_hopper to_return;
    to_return.m_hop_count = m_hop_count + 1;
    to_return.m_hop_state = m_hop_state;
    do {
      uint32_t old_hop_state = to_return.m_hop_state;
      to_return.m_hop_state = detail::hash32(old_hop_state);
      if(old_hop_state == to_return.m_hop_state) {
        to_return.m_hop_state++;
      }
      to_return.m_dest = to_return.m_hop_state % s_mpi_size;
    } while (to_return.m_dest ==  havoqgt_env()->world_comm().rank());
    return to_return;
  }

  uint32_t m_dest;
  uint32_t m_hop_count;
  uint32_t m_hop_state;
  uint32_t pad;
  static uint32_t s_hop_limit;
  static uint32_t s_mpi_size;
};

uint32_t msg_hopper::s_hop_limit;
uint32_t msg_hopper::s_mpi_size;


class receive_iterator
      : public std::iterator<std::output_iterator_tag, void, void, void, void> {
public:
  receive_iterator(mpi::termination_detection<uint64_t>* _td, std::vector<msg_hopper>* _p) 
    :  ptd(_td)
    ,  ppending(_p) { }

  template <typename T>
  receive_iterator& operator=(const T& value) {
    if(value.finished()) {
      ptd->inc_completed();
      //std::cout << "Hop Finished!!!" << std::endl;
    } else {
      //std::cout << "Just received hop " << value.m_hop_count << ", rank = " << havoqgt_env()->world_comm().rank() << std::endl;
      T next = value.next();
      ppending->push_back(next);
    }
  } 

  template <typename T>
  bool intercept(const T& value) { return true; }

  receive_iterator& operator*() { return *this; }

  receive_iterator& operator++() { return *this; }

  receive_iterator operator++(int) { return *this; }

  mpi::termination_detection<uint64_t>* ptd;
  std::vector<msg_hopper>* ppending;
};



int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    havoqgt::get_environment();
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();
    if(mpi_rank == 0) {
      havoqgt::get_environment().print();
    }

    if(argc != 3) {
      std::cerr << "Usage: <num hoppers> <hop length>" << std::endl;
      exit(-1);
    }

    uint64_t count = boost::lexical_cast<uint64_t>(argv[1]);
    msg_hopper::s_hop_limit = boost::lexical_cast<uint32_t>(argv[2]);
    msg_hopper::s_mpi_size = mpi_size;
    if(mpi_rank == 0) {
    std::cout << "Hop Count = " << count << std::endl;
    std::cout << "Hop Limit = " << msg_hopper::s_hop_limit << std::endl;
    }

    mpi::mailbox_routed<msg_hopper> mailbox(1);
    mpi::termination_detection<uint64_t> td(MPI_COMM_WORLD);
    std::vector<msg_hopper> pending;

    uint64_t count_per_rank = count / uint64_t(mpi_size);
    boost::counting_iterator<uint64_t> seq_beg(uint64_t(mpi_rank) * count_per_rank);
    boost::counting_iterator<uint64_t> seq_end(uint64_t(mpi_rank+1) * count_per_rank);

    havoqgt_env()->world_comm().barrier();   
    double time_start = MPI_Wtime();

    do {
      do {
        if(seq_beg != seq_end) {
          msg_hopper msg(*seq_beg);
          ++seq_beg;
          mailbox.send(msg.dest(), msg, receive_iterator(&td, &pending), false);
          td.inc_queued();
        }
        while(!pending.empty()) {
          msg_hopper msg;
          msg = pending.back();
          pending.pop_back();
          mailbox.send(msg.dest(), msg, receive_iterator(&td, &pending), false);
        }
      } while(!pending.empty() || seq_beg != seq_end);
      mailbox.receive(receive_iterator(&td, &pending));
      mailbox.flush_buffers_if_idle();
    } while(!pending.empty() || !mailbox.is_idle() || !td.test_for_termination());

    havoqgt_env()->world_comm().barrier();   
    double time_end = MPI_Wtime();

    if(mpi_rank == 0) {
      double time = time_end - time_start;
      uint64_t hops = uint64_t(msg_hopper::s_hop_limit) * uint64_t(count);
      std::cout << "Ellapsed time = " << time << std::endl;
      std::cout << "Total hops = " << hops << std::endl;
      std::cout << "Rate = " << double(hops) / time << " hops/sec" << std::endl;
    }

  }
  havoqgt::havoqgt_finalize();

  return 0;
}
