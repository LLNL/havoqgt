// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_METALL_DISTRIBUTED_DB_HPP_INCLUDED
#define HAVOQGT_METALL_DISTRIBUTED_DB_HPP_INCLUDED

#include <metall/utility/metall_mpi_adaptor.hpp>

namespace havoqgt {

class db_create {};
class db_open {};
class db_open_read_only {};

class metall_distributed_db {
 private:
  using impl_t = metall::utility::metall_mpi_adaptor;

 public:
  using manager_type = metall::utility::metall_mpi_adaptor::manager_type;
  template <typename T = std::byte>
  using allocator = manager_type::allocator_type<T>;

  metall_distributed_db(db_open,
                        const std::string &data_store_dir,
                        const MPI_Comm &comm = MPI_COMM_WORLD)
      : m_impl(metall::open_only, data_store_dir, comm) {}

  metall_distributed_db(db_open_read_only,
                        const std::string &data_store_dir,
                        const MPI_Comm &comm = MPI_COMM_WORLD)
      : m_impl(metall::open_read_only, data_store_dir, comm) {}

  metall_distributed_db(db_create,
                        const std::string &data_store_dir,
                        const MPI_Comm &comm = MPI_COMM_WORLD)
      : m_impl(metall::create_only, data_store_dir, comm) {}

  metall_distributed_db(db_create,
                        const std::string &data_store_dir,
                        const double gbyte_per_rank,
                        const MPI_Comm &comm = MPI_COMM_WORLD)
      : m_impl(metall::create_only, data_store_dir, uint64_t(gbyte_per_rank * 1024 * 1024) * 1024ULL, comm) {}

  static bool transfer(const std::string& src_base_fname,
                       const std::string& dest_base_fname,
                       const MPI_Comm &comm = MPI_COMM_WORLD) {
    return metall::utility::metall_mpi_adaptor::copy(src_base_fname.c_str(), dest_base_fname.c_str(), comm);
  }

  bool snapshot(const char *destination_dir_path) {
    return m_impl.snapshot(destination_dir_path);
  }

  static bool remove(const std::string &dir_path,
                     const MPI_Comm &comm = MPI_COMM_WORLD) {
    return impl_t::remove(dir_path.c_str(), comm);
  }

  static int partitions(const std::string &data_store_dir, const MPI_Comm &comm = MPI_COMM_WORLD) {
    return impl_t::partitions(data_store_dir.c_str(), comm);
  }

  manager_type *get_manager() { return &(m_impl.get_local_manager()); }

  const manager_type *get_manager() const {
    return &(m_impl.get_local_manager());
  }

  template <typename T = std::byte>
  allocator<T> get_allocator() {
    return m_impl.get_local_manager().get_allocator<T>();
  }

 private:
  impl_t m_impl;
};

}  // namespace havoqgt
#endif  // HAVOQGT_METALL_DISTRIBUTED_DB_HPP_INCLUDED
