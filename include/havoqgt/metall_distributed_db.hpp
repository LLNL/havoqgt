/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */
#ifndef HAVOQGT_METALL_DISTRIBUTED_DB_HPP_INCLUDED
#define HAVOQGT_METALL_DISTRIBUTED_DB_HPP_INCLUDED

#include <havoqgt/distributed_db.hpp>
#include <metall_utility/metall_mpi_adaptor.hpp>

namespace havoqgt {

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
                        const std::size_t capacity,
                        const MPI_Comm &comm = MPI_COMM_WORLD)
      : m_impl(metall::create_only, data_store_dir, capacity, comm) {}

  static bool transfer(const char *const src_base_fname,
                       const char *const dest_base_fname,
                       const MPI_Comm &comm = MPI_COMM_WORLD) {
    return metall::utility::metall_mpi_adaptor::copy(src_base_fname, dest_base_fname, comm);
  }

  bool snapshot(const char *destination_dir_path) {
    return m_impl.snapshot(destination_dir_path);
  }

  static bool remove(const char *dir_path,
                     const MPI_Comm &comm = MPI_COMM_WORLD) {
    return impl_t::remove(dir_path, comm);
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
