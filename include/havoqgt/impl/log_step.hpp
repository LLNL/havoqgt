// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT
#ifndef __HAVOQGT_IMP_LOG_STEP_HPP__
#define __HAVOQGT_IMP_LOG_STEP_HPP__

namespace havoqgt {

class LogStep {
 public:
  LogStep(const std::string &stp_str, MPI_Comm mpi_comm, int mpi_rank)
      : stp_str_(stp_str), mpi_comm_(mpi_comm), mpi_rank_(mpi_rank) {
    MPI_Barrier(mpi_comm_);
    if (mpi_rank_ == 0) {
      time_ = MPI_Wtime();
      std::cout << "Starting:  " << stp_str_ << std::endl << std::flush;
    }
    MPI_Barrier(mpi_comm_);
  }

  ~LogStep() {
    MPI_Barrier(mpi_comm_);

    if (mpi_rank_ == 0) {
      time_ = MPI_Wtime() - time_;

      std::cout << "Finished: " << stp_str_ << " in " << time_ << " seconds."
                << std::endl
                << std::flush;
    }
    MPI_Barrier(mpi_comm_);
  }

 private:
  double time_;
  int    mb_read_;
  int    mb_written_;

  std::string stp_str_;
  MPI_Comm    mpi_comm_;
  int         mpi_rank_;
};

}  // namespace havoqgt

#endif  // __HAVOQGT_IMP_LOG_STEP_HPP__
