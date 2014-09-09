#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <fstream>      // std::ifstream
#include <sstream>      // std::stringstream
#include <iostream>     // std::cout, std::ios
#include <havoqgt/mpi.hpp>

#include <sys/mman.h>

const bool verbose = false;

const int processes_per_node = 24;
const int node_partions = 8;

void execute_command(int backup, std::string fname_input,
  std::string fname_output);

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Init(&argc, &argv));
  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
  CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

  int backup;
  std::string fname_input;
  std::string fname_output;

  int pos = 1;
  if (argc == 3) {
    fname_input = argv[pos++];
    fname_output = argv[pos++];

  } else {
    if (mpi_rank == 0) {
      std::cout << "Parameter Error: <input_file> <output_file>" << std::endl;
    }
    return -1;
  }

  if (mpi_rank == 0) {
    std::cout << "Restoring File:";
    std::cout << "Source = '" << fname_input.c_str() << std::endl;
    std::cout << "Dest = '" << fname_output.c_str() << std::endl;
  }

  fname_input += "_" + std::to_string(mpi_rank);
 {
    std::ifstream fin(fname_input);
    bool is_good = fin.good();
    fin.close();
    if (!is_good) {
      std::cout << "[" << mpi_rank << "] File not found: " << fname_input
        << "." << std::endl;
      return -1;
    }
  }

  bool hit = false;
  for (int i = 0; i < node_partions; i++) {
    if (mpi_rank == 0) {
      const int lower_bound = (i)*(processes_per_node/node_partions);
      const int upper_bound = (i+1)*(processes_per_node/node_partions)-1;

      std::cout << "***(" << i << "/" << node_partions << ") Processes: "
        << lower_bound << " - " << upper_bound
        << " executing on each node." << std::endl;
    }
    if ((mpi_rank % processes_per_node) % node_partions == i){
      if (hit) {
        std::cout << "[" << mpi_rank
          << "] error: attempting to backup/restor twice." << std::endl;
        exit(-1);
      } else {
        hit = true;
        if (verbose)
          std::cout << "[" << mpi_rank << "]";
        execute_command(backup, fname_input, fname_output);
      }

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0)
    std::cout << "Transfer Graph Fin." << std::endl;
  MPI_Finalize();
  return 0;

};

void execute_command(int backup, std::string fname_input,
  std::string fname_output) {

  std::string tar_str;
  if (verbose)
    std::cout << "Restoring File:";

  tar_str = "cp " + fname_input + " " + fname_output;

  if (verbose)
    std::cout << "Source = '" << fname_input.c_str() << "' "
      << "Destination = '" << fname_output.c_str() << "'" << std::endl
      << "Command = '" << tar_str.c_str() << "'" << std::endl;

  system(tar_str.c_str());
  std::cout << std::flush;

};  // execute_command;
