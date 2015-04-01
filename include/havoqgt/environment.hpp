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

#ifndef HAVOQGT_ENVIRONMENT_HPP_INCLUDED
#define HAVOQGT_ENVIRONMENT_HPP_INCLUDED

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cstdlib>
#include <boost/lexical_cast.hpp>

#include <havoqgt/error.hpp>
#include <havoqgt/mpi.hpp>

namespace havoqgt {
/** 
 * @todo delete this 
 * @todo add environment variable documentation 
 */
class old_environment {
public:
  old_environment() {
    m_mailbox_num_irecv   = get_env_var<uint32_t>("HAVOQGT_MAILBOX_NUM_IRECV", 8);
    m_mailbox_num_isend   = get_env_var<uint32_t>("HAVOQGT_MAILBOX_NUM_ISEND", 8);
    m_mailbox_aggregation = get_env_var<uint32_t>("HAVOQGT_MAILBOX_AGGREGATION", 1024);
    m_mailbox_tree_aggregation = get_env_var<uint32_t>("HAVOQGT_MAILBOX_TREE_AGGREGATION", 64);
    m_mailbox_print_stats = get_env_var<bool>    ("HAVOQGT_MAILBOX_PRINT_STATS", false);
  }

  uint32_t mailbox_num_irecv()   const { return m_mailbox_num_irecv; }
  uint32_t mailbox_num_isend()   const { return m_mailbox_num_isend; }
  uint32_t mailbox_aggregation() const { return m_mailbox_aggregation; }
  uint32_t mailbox_tree_aggregation() const { return m_mailbox_tree_aggregation; }
  bool     mailbox_print_stats() const { return m_mailbox_print_stats; }

  template <typename T>
  inline T get_env_var(const char* key, T default_val) const;

  void print() const;

private:
  uint32_t  m_mailbox_num_irecv;
  uint32_t  m_mailbox_num_isend;
  uint32_t  m_mailbox_aggregation;
  uint32_t  m_mailbox_tree_aggregation;
  bool      m_mailbox_print_stats;
};

inline void
old_environment::print() const {
  std::cout << "HAVOQGT_MAILBOX_NUM_IRECV        "<< " = " << m_mailbox_num_irecv << std::endl;
  std::cout << "HAVOQGT_MAILBOX_NUM_ISEND        "<< " = " << m_mailbox_num_isend << std::endl;
  std::cout << "HAVOQGT_MAILBOX_AGGREGATION      "<< " = " << m_mailbox_aggregation << std::endl;
  std::cout << "HAVOQGT_MAILBOX_TREE_AGGREGATION "<< " = " << m_mailbox_tree_aggregation << std::endl;
  std::cout << "HAVOQGT_MAILBOX_PRINT_STATS      "<< " = " << m_mailbox_print_stats << std::endl;
}

template <typename T>
inline 
T
old_environment::get_env_var(const char* key, T default_val) const {
  char* val = std::getenv( key );
  if(val != NULL) {
    try {
      default_val = boost::lexical_cast<T>(val);
    } catch (...) {
      std::stringstream err;
      err << "havoqgt::environment -- Unable to parse environment variable: "
          << key << "=" << val << std::endl;
      throw std::runtime_error(err.str());
    }
  }

  return default_val;
}

inline 
old_environment& 
get_environment() {
  static old_environment env;
  return env;
}


///   Start 'new' environmet here

class environment {
public:
  environment(int *argc, char*** argv) {
    // Init MPI
    MPI_Init(argc, argv);
    
    // Create communicators
    m_world_comm = communicator(MPI_COMM_WORLD);
    MPI_Comm tmp_node_local_comm  = split_node_local_comm(MPI_COMM_WORLD);
    m_node_local_comm = communicator(tmp_node_local_comm);
    MPI_Comm tmp_node_offset_comm = split_node_offset_comm(MPI_COMM_WORLD, tmp_node_local_comm);
    m_node_offset_comm = communicator(tmp_node_offset_comm);
  }
  
  ~environment() {
    m_world_comm.barrier();
    MPI_Finalize();
  }
  
  const communicator& world_comm()       const { return m_world_comm; }
  const communicator& node_local_comm()  const { return m_node_local_comm; }
  const communicator& node_offset_comm() const { return m_node_offset_comm; }

  std::string whoami() const {
    std::stringstream sstr;
    sstr << "(R" <<  world_comm().rank() 
         << ",N" << node_offset_comm().rank() 
         << ",C" << node_local_comm().rank() << ")";
    return sstr.str();
  }
  
private:
  
  MPI_Comm split_node_local_comm(MPI_Comm world_comm)
  {
    MPI_Comm to_return;
    int rank;
    MPI_Comm_rank(world_comm, &rank);

    int color = 0;
    char* cnodeid   = NULL;//getenv("SLURM_NODEID");
    char chostname[256];
    gethostname(chostname, 256);
    if(cnodeid) {
      color = boost::lexical_cast<int>(cnodeid);
    } else if(chostname) {
      std::string shostname(chostname);
      std::hash<std::string> hasher;
      color = hasher(shostname);
      if(color < 0) color *= -1;
    }

    int key = rank;
    int ret = MPI_Comm_split(world_comm, color, key, &to_return);

    int node_local_rank, node_local_size;
    MPI_Comm_rank(to_return, &node_local_rank);
    MPI_Comm_size(to_return, &node_local_size);

    if(rank < node_local_size)
    {
      if(rank != node_local_rank) {
        std::cerr << "ERROR:   Only blocked task mapping is supported" << std::endl; exit(-1);
      }
    }



    return to_return;
  }

  MPI_Comm split_node_offset_comm(MPI_Comm world_comm, MPI_Comm node_local_comm)
  {
    MPI_Comm to_return;
    int world_rank, node_local_rank;
    MPI_Comm_rank(world_comm, &world_rank);
    MPI_Comm_rank(node_local_comm, &node_local_rank);

    int color = node_local_rank;
    int key   = world_rank;

    int ret = MPI_Comm_split(world_comm, color, key, &to_return);

    return to_return;
  }

  bool is_node_mapping_round_robin(MPI_Comm world_comm, MPI_Comm node_local_comm)
  {

  }

  bool is_node_mapping_block(MPI_Comm world_comm, MPI_Comm node_local_comm)
  {
    int world_rank, world_size, node_local_rank, node_local_size;
    MPI_Comm_rank(world_comm, &world_rank);
    MPI_Comm_size(world_comm, &world_size);
    MPI_Comm_rank(node_local_comm, &node_local_rank);
    MPI_Comm_size(node_local_comm, &node_local_size);

  }

  communicator   m_world_comm;
  communicator   m_node_local_comm;
  communicator   m_node_offset_comm;
};


namespace detail {
  inline environment*& priv_havoqgt_env() {
    static environment* penv;
    return penv;
  }
} 

inline environment* havoqgt_env() {
  return detail::priv_havoqgt_env();
}


inline void havoqgt_init( int *argc, char ***argv ) {
  detail::priv_havoqgt_env() = new environment(argc, argv);
}

inline void havoqgt_finalize() {
  delete detail::priv_havoqgt_env();
}


template <typename T>
inline 
T
havoqgt_getenv(const char* key, T default_val) {
  char* val = std::getenv( key );
  if(val != NULL) {
    try {
      default_val = boost::lexical_cast<T>(val);
    } catch (...) {
      std::stringstream err;
      err << "havoqgt::environment -- Unable to parse environment variable: "
          << key << "=" << val << std::endl;
      throw std::runtime_error(err.str());
    }
  }

  return default_val;
}



} //namespace havoqgt

#endif
