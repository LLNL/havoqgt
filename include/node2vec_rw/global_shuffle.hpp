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

#ifndef HAVOQGT_GLOBAL_SHUFFLE_HPP
#define HAVOQGT_GLOBAL_SHUFFLE_HPP

#include <sched.h>
#include <vector>
#include <mpi.h>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <utility>
#include <string>
#include <random>
#include <sstream>
#include <iomanip>      // std::setfill, std::setw
#include <boost/container/vector.hpp>
#include <havoqgt/mpi.hpp>

namespace node2vec_rw {

/// \brief Very naive implementation of mpi all to all for 64 bit data size
template <typename input_iterator_type, typename back_insert_iterator_type>
void mpi_large_all_to_all(input_iterator_type input_begin, input_iterator_type input_end,
                          std::function<int(const typename std::iterator_traits<input_iterator_type>::value_type &)> destination,
                          MPI_Comm mpi_comm,
                          back_insert_iterator_type output_begin) {
  using value_type = typename std::iterator_traits<input_iterator_type>::value_type;

  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_size(mpi_comm, &mpi_size));
  CHK_MPI(MPI_Comm_rank(mpi_comm, &mpi_rank));

  std::vector<std::vector<value_type>> send_buf(mpi_size);
  std::vector<std::vector<value_type>> recv_buf;
  const int max_count = (int)((1ULL << 31ULL) / mpi_size / (int)std::log2(mpi_size));

  auto input_itr = input_begin;
  while (true) {
    // --- Copy data to the send buffer --- //
    while (true) {
      if (input_itr == input_end) break;

      const int dest_rank = destination(*input_itr);
      send_buf[dest_rank].push_back(*input_itr);
      ++input_itr;

      if (send_buf[dest_rank].size() >= max_count) break;
    }

    havoqgt::mpi_all_to_all(send_buf, recv_buf, mpi_comm);
    for (const auto &single_buf : recv_buf) {
      std::transform(single_buf.begin(), single_buf.end(), output_begin,
                     [](value_type value) -> value_type { return value; });
    }
    send_buf.clear();
    recv_buf.clear();

    // --- Global termination check --- //
    char local_finished_flag = (input_itr == input_end);
    const char global_finished_flag = havoqgt::mpi_all_reduce(local_finished_flag,
                                                              std::logical_and<char>{},
                                                              mpi_comm);
    if (global_finished_flag) break;
  }
}

} // namespace node2vec_rw

#endif //HAVOQGT_GLOBAL_SHUFFLE_HPP
