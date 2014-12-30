/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Re-written by Steven Feldman <feldman12@llnl.gov>.
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
#ifndef __HAVOQGT_IMP_LOG_STEP_HPP__
#define __HAVOQGT_IMP_LOG_STEP_HPP__

#include <havoqgt/environment.hpp>

namespace havoqgt {
namespace mpi {

class LogStep {
 public:
  LogStep(const std::string &stp_str, MPI_Comm mpi_comm, int mpi_rank)
  : stp_str_(stp_str)
  , mpi_comm_(mpi_comm)
  , mpi_rank_(mpi_rank) {
    MPI_Barrier(mpi_comm_);
    if (mpi_rank_ == 0) {
      time_ = MPI_Wtime();
      std::cout << "Starting:  " << stp_str_ << std::endl << std::flush;
    }
    if (mpi_rank_ % havoqgt_env()->node_local_comm().size()  == 0) {
      get_io_stat_info(mb_read_, mb_written_);
      std::cout << "\t[" << mpi_rank_ << "] Dirty Pages: " << get_dirty_pages()
        << "kb." << std::endl << std::flush;
    }
    MPI_Barrier(mpi_comm_);
  }

  ~LogStep() {
    MPI_Barrier(mpi_comm_);

    if (mpi_rank_ == 0) {
      time_ = MPI_Wtime() - time_;

      std::cout << "Finished: " << stp_str_ <<  " in " << time_ << " seconds."
        << std::endl << std::flush;
    }
    if (mpi_rank_ % havoqgt_env()->node_local_comm().size() == 0) {
      int read = -1;
      int written = -1;
      get_io_stat_info(read, written);

      std::cout
        << "\t[" << mpi_rank_ << "]" << std::endl
        << "\tSpace used: " << get_disk_utilization() << "gB." << std::endl
        << "\tDirty Pages: " << get_dirty_pages() << "kB." << std::endl
        << "\tMB Read: " << (read - mb_read_) << std::endl
        << "\tMB Written: " << (written - mb_written_) << std::endl
        << std::flush;
    }
    MPI_Barrier(mpi_comm_);
  }

 private:
  double time_;
  int mb_read_;
  int mb_written_;

  std::string stp_str_;
  MPI_Comm mpi_comm_;
  int mpi_rank_;
};


}  // namespace mpi
}  // namespace havoqgt

#endif  // __HAVOQGT_IMP_LOG_STEP_HPP__
