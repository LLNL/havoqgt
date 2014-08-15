#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <fstream>      // std::ifstream
#include <sstream>      // std::stringstream
#include <iostream>     // std::cout, std::ios
#include <havoqgt/mpi.hpp>

#include <sys/mman.h>

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  CHK_MPI(MPI_Init(&argc, &argv));
  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
  CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

  std::string fname_input, fname_output;

  int pos = 1;
  if (argc == 3) {
    fname_input = argv[pos++];
    fname_input += "_" + std::to_string(mpi_rank);
    fname_output = argv[pos++];

    std::cout << "[" << mpi_rank
      << "] Copying File From '" << fname_input.c_str()
      << "' to '" << fname_output.c_str() << "'." << std::endl;
  } else {
    if (mpi_rank == 0) {
      std::cout << "Parameter Error: <input_file> <output_directory>:"
        << std::endl;
    }
    return -1;
  }


  {
    std::ifstream fin(fname_input);
    if (!fin.good()) {
      fin.close();
      std::cout << "[" << mpi_rank << "] File not found: " << fname_input
        << "." << std::endl;
      return -1;
    }
    fin.close();
  }


  unsigned found = fname_input.rfind("/");
  if (found == std::string::npos) {
    found = 0;
  }
  else {
    found++;
  }

  std::string out_filename = fname_output+"/"+fname_input.substr(found);
  std::cout << "Out File: " << out_filename << std::endl;
  std::ifstream fin(out_filename);
  if (fin.good()) {
    time_t timer;
    time(&timer);
    fin.close();

    std::string moved_out_file = out_filename + "_" + std::to_string(timer);

    std::cout << "[" << mpi_rank << "] File already exist, moving it to  '"
      << moved_out_file.c_str() << "." << std::endl;

    std::string mv_str = "mv " + out_filename + " " + moved_out_file;
    system(mv_str.c_str());
    return -1;
  } else {
    fin.close();
  }

  std::string mv_str = "mv " + fname_input + " " + out_filename;
  system(mv_str.c_str());

  std::cout << "Fin." << std::endl;

}  // main
