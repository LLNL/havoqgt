
/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Steven Feldman <feldman12@llnl.gov>.
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

#ifndef __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__
#define __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__

namespace havoqgt {
namespace mpi {

class local_source_id {
 public:
  explicit local_source_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(std::pair<T, T> i) const { return i.first / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class local_dest_id {
 public:
  explicit local_dest_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(std::pair<T, T> i) const { return i.second / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class get_local_id {
 public:
  explicit get_local_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(T i) const { return i / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class owner_source_id {
 public:
  explicit owner_source_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(std::pair<T, T> i) const { return i.first % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class owner_dest_id {
 public:
  explicit owner_dest_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(std::pair<T, T> i) const { return i.second % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class get_owner_id {
 public:
  explicit get_owner_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(T i) const { return i % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

}  // namespace mpi
}  // namespace havoqgt

#endif  // __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__
