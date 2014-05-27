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

/** @page page1 Environment Variables
 *  @todo add documentation
 * List of Environment variables
 *    - HAVOQGT_VERBOSE           -- adsf
 *    - HAVOQGT_MAILBOX_NUM_IRECV -- adsf
 */


namespace havoqgt {

class environment {
public:
  environment() {
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
environment::print() const {
  std::cout << "HAVOQGT_MAILBOX_NUM_IRECV        "<< " = " << m_mailbox_num_irecv << std::endl;
  std::cout << "HAVOQGT_MAILBOX_NUM_ISEND        "<< " = " << m_mailbox_num_isend << std::endl;
  std::cout << "HAVOQGT_MAILBOX_AGGREGATION      "<< " = " << m_mailbox_aggregation << std::endl;
  std::cout << "HAVOQGT_MAILBOX_TREE_AGGREGATION "<< " = " << m_mailbox_tree_aggregation << std::endl;
  std::cout << "HAVOQGT_MAILBOX_PRINT_STATS      "<< " = " << m_mailbox_print_stats << std::endl;
}

template <typename T>
inline 
T
environment::get_env_var(const char* key, T default_val) const {
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
environment& 
get_environment() {
  static environment env;
  return env;
}


} //namespace havoqgt

#endif