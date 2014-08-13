#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <sstream>      // std::stringstream
#include <iostream>     // std::cout, std::ios
#include <sys/mman.h>

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  CHK_MPI(MPI_Init(&argc, &argv));
  {
    CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
    CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

    std::string fname_input, fname_output;

    if (argc < 3) {
      fname_input = argv[pos++] + "_" + mpi_rank;
      fname_output = argv[pos++];

      std::cout << "[" << mpi_rank
        << "] Copying File From " << fname_input.c_str() <<
        << " to " << fname_output.c_str() << "." << std::endl;
    } else {
      if (mpi_rank == 0) {
        std::cout << "Transfering Files <input_file> <output_directory>:"
          << std::endl;
      }
      return -1;
    }

    ifstream fin(fname_input);
    if (!fin.good()) {
      fin.close();
      std::cout << "[" << mpi_rank << "] File not found: " << fname_input
        << "." << std::endl;
      return -1;
    }
    fin.close();

    unsigned found = fname_output.rfind("/");
    if (found == std::string::npos || )
      found = 0;
    else
      found++;

    unsigned found = str.rfind("/");
    if (found == std::string::npos)
      found = 0;
    else
      found++;

    filename = str.substr(found);

    ifstream fin(fname_output);
    if (!fin.good()) {
      time_t timer;
      time(&timer);
      fin.close();
      std::cout << "[" << mpi_rank << "] File already exist, moving it: " << fname_output
        << "." << std::endl;
      std::stringstream mv_str = "mv " << fname_input << " " << fname_output;
      system(mv_str);
      return -1;
    }
    fin.close();

    std::stringstream mv_str = "mv " << fname_input << " " << fname_output;
    system(mv_str);

  } // CHK_MPI
}  // main
