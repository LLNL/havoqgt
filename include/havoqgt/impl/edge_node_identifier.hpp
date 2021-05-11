// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__
#define __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__

namespace havoqgt {

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

}  // namespace havoqgt

#endif  // __HAVOQGT_IMP_EDGE_NODE_IDENTIFIER_HPP__
