// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fcntl.h>
#include <mpi.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iostream>
#include <vector>

#define CHK_MPI(a)                                               \
  {                                                              \
    if (a != MPI_SUCCESS) {                                      \
      char *error_string = NULL;                                 \
      int   len          = 0;                                    \
      MPI_Error_string(a, error_string, &len);                   \
      std::cerr << __FILE__ << ", line " << __LINE__             \
                << " MPI ERROR = " << error_string << std::endl; \
      exit(-1);                                                  \
    }                                                            \
  }

namespace ygm {

inline MPI_Datatype mpi_typeof(char) { return MPI_CHAR; }
inline MPI_Datatype mpi_typeof(signed short) { return MPI_SHORT; }
inline MPI_Datatype mpi_typeof(signed int) { return MPI_INT; }
inline MPI_Datatype mpi_typeof(signed long) { return MPI_LONG; }
inline MPI_Datatype mpi_typeof(unsigned char) { return MPI_UNSIGNED_CHAR; }
inline MPI_Datatype mpi_typeof(unsigned short) { return MPI_UNSIGNED_SHORT; }
inline MPI_Datatype mpi_typeof(unsigned) { return MPI_UNSIGNED; }
inline MPI_Datatype mpi_typeof(unsigned long) { return MPI_UNSIGNED_LONG; }
inline MPI_Datatype mpi_typeof(unsigned long long) {
  return MPI_UNSIGNED_LONG_LONG;
}
inline MPI_Datatype mpi_typeof(signed long long) { return MPI_LONG_LONG_INT; }
inline MPI_Datatype mpi_typeof(double) { return MPI_DOUBLE; }
inline MPI_Datatype mpi_typeof(long double) { return MPI_LONG_DOUBLE; }
inline MPI_Datatype mpi_typeof(std::pair<int, int>) { return MPI_2INT; }
inline MPI_Datatype mpi_typeof(std::pair<float, int>) { return MPI_FLOAT_INT; }
inline MPI_Datatype mpi_typeof(std::pair<double, int>) {
  return MPI_DOUBLE_INT;
}
inline MPI_Datatype mpi_typeof(std::pair<long double, int>) {
  return MPI_LONG_DOUBLE_INT;
}
inline MPI_Datatype mpi_typeof(std::pair<short, int>) { return MPI_SHORT_INT; }

class comm;
const comm &       comm_world();
const comm &       comm_nl();
const comm &       comm_nr();
const comm &       comm_nlnr();
inline std::string whoami();

void init(int *argc, char ***argv);

class comm {
 public:
  enum class split { nl, nr, nlnr };

  comm(MPI_Comm comm = MPI_COMM_WORLD);

  comm(split _split);

  ~comm();

  int      size() const { return m_size; }
  int      rank() const { return m_rank; }
  MPI_Comm mpi_comm() const { return m_comm; }

  void barrier() const { chk_ret(MPI_Barrier(m_comm)); }

  template <typename T>
  T all_reduce(T in_d, MPI_Op in_op) const {
    T to_return;
    chk_ret(
        MPI_Allreduce(&in_d, &to_return, 1, mpi_typeof(T()), in_op, m_comm));
    return to_return;
  }

  template <typename T>
  void broadcast(T &data, int root) const {
    chk_ret(MPI_Bcast(&data, sizeof(T), MPI_BYTE, root, m_comm));
  }

  template <typename T>
  void broadcast(std::vector<T> &vec_data, int root) const {
    size_t size(0), capacity(0);
    if (rank() == root) {
      size     = vec_data.size();
      capacity = vec_data.capacity();
    } else {
      vec_data.clear();
    }
    broadcast(size, root);
    broadcast(capacity, root);
    if (rank() != root) {
      vec_data.reserve(capacity);
      vec_data.resize(size);
    }
    chk_ret(
        MPI_Bcast(&(vec_data[0]), sizeof(T) * size, MPI_BYTE, root, m_comm));
  }

 private:
  /// Not Copyable
  comm();
  comm(const comm &) = delete;             // non construction-copyable
  comm &operator=(const comm &) = delete;  // non copyable

  /// Checks MPI return codes

  void chk_ret(int ret, const char *loc = NULL) const;

  void init_nlnr_comm();
  void init_nl_comm();
  void init_nr_comm();
  void init_rank_size();

  MPI_Comm m_comm;
  int      m_size;
  int      m_rank;
  bool     free_comm_in_destructor;
};

void abort(int errorcode = -1) { MPI_Abort(MPI_COMM_WORLD, errorcode); }

}  // namespace ygm

#include <ygm/impl/mpi.ipp>
