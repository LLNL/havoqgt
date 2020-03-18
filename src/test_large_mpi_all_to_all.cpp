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

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <iterator>
#include <havoqgt/mpi.hpp>
#include <node2vec_rw/global_shuffle.hpp>

int main(int argc, char **argv) {

  havoqgt::init(&argc, &argv);
  {
    int mpi_size(0), mpi_rank(0);
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

    const std::size_t num_global_elements = (1ULL << 20ULL);

    std::mt19937 rnd_gen(123);

    std::vector<std::pair<int, uint64_t>> send_buf;
    for (std::size_t i = 0; i < num_global_elements; ++i) {
      if (rnd_gen() % mpi_size == mpi_rank) {
        const int dest_rank = i % mpi_size;
        const int local_index = i / mpi_size;
        send_buf.emplace_back(std::make_pair(dest_rank, local_index));
      }
    }

    auto partitioner = [](const std::pair<int, uint64_t> &data) -> int {
      return data.first;
    };

    std::vector<std::pair<int, uint64_t>> recv_buf;
    node2vec_rw::mpi_large_all_to_all(send_buf.begin(),
                                      send_buf.end(),
                                      partitioner,
                                      MPI_COMM_WORLD,
                                      std::back_inserter(recv_buf));

    std::sort(recv_buf.begin(), recv_buf.end(),
              [](const std::pair<int, uint64_t> &lhd, const std::pair<int, uint64_t> &rhd) {
                return lhd.second < rhd.second;
              });

    uint64_t local_count = 0;
    for (int r = 0; r < mpi_size; ++r) {
      if (r == mpi_rank) {
        for (std::size_t i = 0; i < recv_buf.size(); ++i) {
          if (recv_buf[i].first != mpi_rank || recv_buf[i].second != i) {
            std::ostringstream oss;
            oss << "Test failed at mpi_rank = " << mpi_rank << ", i = " << i << " : "
                << recv_buf[i].first << " " << recv_buf[i].second;
            std::cerr << oss.str() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
          }
        }
        local_count = recv_buf.size();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    const uint64_t global_count = havoqgt::mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
    if (global_count != num_global_elements) {
      havoqgt::cout_rank0() << "Test failed (global count does not match): "
                            << global_count << " != " << num_global_elements << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    havoqgt::cout_rank0() << "Passed the test" << std::endl;
  }

  return 0;
}