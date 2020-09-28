// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <iostream>
#include <string>

#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Init(&argc, &argv));
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::string fname_input;
  std::string fname_output;

  int pos = 1;
  if (argc == 3) {
    fname_input  = argv[pos++];
    fname_output = argv[pos++];

  } else {
    if (mpi_rank == 0) {
      std::cout << "Parameter Error: <input_file> <output_file>" << std::endl;
    }
    return -1;
  }

  if (mpi_rank == 0) {
    std::cout << "Transferring Data:";
    std::cout << "Source = " << fname_input << std::endl;
    std::cout << "Dest = " << fname_output << std::endl;
  }

  if (!havoqgt::distributed_db::transfer(fname_input, fname_output)) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) std::cout << "Transfer Graph Fin." << std::endl;

  MPI_Finalize();

  return 0;
};
