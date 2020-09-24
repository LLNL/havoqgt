// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <iostream>
#include <string>

#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>

enum operation_mode { invalid, copy, rm, stats };

void usage() {
  if (havoqgt::comm_world().rank() == 0) {
    std::cerr << "Usage: [one of operation modes: 'c', 'r', or 's'] [directory "
                 "path 1] "
                 "[directory path 2 (required for copy mode)]\n"
              << " -c            - copy mode\n"
              << " -r            - remove mode\n"
              << " -s            - stats mode\n"
              << " -h            - print help and exit\n";
  }
}

void parse_cmd_line(int argc, char** argv, operation_mode& op_mode,
                    std::string& dir_path1, std::string& dir_path2) {
  if (havoqgt::comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  op_mode = operation_mode::invalid;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "crsh ")) != -1) {
    switch (c) {
      case 'c':
        op_mode = operation_mode::copy;
        break;

      case 'r':
        op_mode = operation_mode::rm;
        break;

      case 's':
        op_mode = operation_mode::stats;
        break;

      case 'h':
        prn_help = true;
        break;

      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }

  if (prn_help) {
    usage();
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int index = optind;
  if (index < argc) {
    dir_path1.assign(argv[index++]);
  }
  if (index < argc) {
    dir_path2.assign(argv[index]);
  }

  if (op_mode == operation_mode::invalid) {
    if (havoqgt::comm_world().rank() == 0)
      std::cerr << "Invalid operation mode" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (op_mode == operation_mode::copy &&
      (dir_path1.empty() || dir_path2.empty())) {
    if (havoqgt::comm_world().rank() == 0)
      std::cerr << "Invalid arguments for copy mode" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (op_mode == operation_mode::rm && dir_path1.empty()) {
    if (havoqgt::comm_world().rank() == 0)
      std::cerr << "Invalid arguments for remove mode" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  if (op_mode == operation_mode::stats && dir_path1.empty()) {
    if (havoqgt::comm_world().rank() == 0)
      std::cerr << "Invalid arguments for stats mode" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Init(&argc, &argv));
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  operation_mode op_mode;
  std::string    dir_path1;
  std::string    dir_path2;

  parse_cmd_line(argc, argv, op_mode, dir_path1, dir_path2);

  // Check if the correct number of MPI procs is used.
  if (mpi_rank == 0) {
    const int required_mpi_procs =
        havoqgt::distributed_db::partitions(dir_path1, MPI_COMM_WORLD);
    if (mpi_size != required_mpi_procs) {
      std::cerr << "Invalid #of MPI processes. Must be " << required_mpi_procs
                << ", not " << mpi_size << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (op_mode == operation_mode::copy) {
    if (mpi_rank == 0) {
      havoqgt::cout_rank0() << "Copy Data Store:";
      havoqgt::cout_rank0() << " Source = " << dir_path1 << std::endl;
      havoqgt::cout_rank0() << " Dest = " << dir_path2 << std::endl;
    }

    if (!havoqgt::distributed_db::transfer(dir_path1, dir_path2)) {
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  } else if (op_mode == operation_mode::rm) {
    if (mpi_rank == 0) {
      havoqgt::cout_rank0() << "Remove Data Store: " << dir_path1 << std::endl;
    }
    if (!havoqgt::distributed_db::remove(dir_path1, MPI_COMM_WORLD)) {
      if (mpi_rank == 0) {
        std::cerr << "Failed to remove data store" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
  } else if (op_mode == operation_mode::stats) {
    // Do some other status check, such as consistency?
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
};
