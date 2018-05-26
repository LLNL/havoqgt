#pragma once

#include <fcntl.h>
#include <mpi.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iostream>

// #define CHK_MPI(a)                                               \
//   {                                                              \
//     if (a != MPI_SUCCESS) {                                      \
//       char *error_string = NULL;                                 \
//       int len = 0;                                               \
//       MPI_Error_string(a, error_string, &len);                   \
//       std::cerr << __FILE__ << ", line " << __LINE__             \
//                 << " MPI ERROR = " << error_string << std::endl; \
//       exit(-1);                                                  \
//     }                                                            \
//   }

/// @warning  This algorithm is O(P) because its general purpose.   Could be
///           Specialized to improve performance.
///
inline MPI_Comm build_node_local_comm() {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  MPI_Comm local_comm;

  int min_shm_rank = -1;
  int shm_com_size = 0;
  const char *shm_name = "build_node_local_comm";
  shm_unlink(shm_name);
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
  //
  // This executes rank-by-rank sequentially.
  {
    // Blocks except for rank 0
    if (mpi_rank > 0) {
      CHK_MPI(MPI_Recv(NULL, 0, MPI_BYTE, mpi_rank - 1, 1, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE));
      // std::cout << "Rank " << mpi_rank << " Ready" << std::endl;
    }
    int shm_fd = -1;
    const int shm_size = 4096;
    bool this_rank_created = false;

    shm_fd = shm_open(shm_name, O_RDWR, 0666);
    if (shm_fd == -1) {
      // std::cout << "Rank " << mpi_rank << ": Failed to open;  attempting
      // to create" << std::endl; Try to create
      shm_fd = shm_open(shm_name, O_CREAT | O_RDWR, 0666);
      if (shm_fd == -1) {
        std::cerr << "Failed to open & create -- Total Failure" << std::endl;
        exit(-1);
      }
      this_rank_created = true;
      if (ftruncate(shm_fd, shm_size) == -1) {
        std::cerr << "Created, but ftruncate failed" << std::endl;
      }
    }

    // By now, should be ready & correct size
    void *ptr =
        mmap(0, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED) {
      printf("Map failed\n");
      exit(-1);
    }

    std::pair<int, int> *p_min_rank_size = (std::pair<int, int> *)ptr;
    if (this_rank_created) {
      p_min_rank_size->first = mpi_rank;
      p_min_rank_size->second = 1;
    } else {
      p_min_rank_size->first = std::min(p_min_rank_size->first, mpi_rank);
      p_min_rank_size->second++;
    }

    //
    // Notifies rank+1 of completion
    if (mpi_rank < mpi_size - 1) {
      // std::cout << "Rank " << mpi_rank << " COMPLETE" << std::endl;
      CHK_MPI(MPI_Send(NULL, 0, MPI_BYTE, mpi_rank + 1, 1, MPI_COMM_WORLD));
    }

    // All ranks have completed.  Each pulls their node loacal's data
    CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    min_shm_rank = p_min_rank_size->first;
    shm_com_size = p_min_rank_size->second;

    //
    // Close shared segment
    if (munmap(ptr, shm_size) != 0) {
      std::cerr << "Rank " << mpi_rank << "munmap failed" << std::endl;
    }
    if (close(shm_fd) != 0) {
      std::cerr << "Rank " << mpi_rank << "close failed" << std::endl;
    }
    CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    if (this_rank_created && shm_unlink(shm_name) != 0) {
      std::cerr << "Rank " << mpi_rank << "shm_unlink failed" << std::endl;
    }
    // std::cout << "Rank " << mpi_rank << ", min_shm_rank = " << min_shm_rank
    //           << ", shm_com_size = " << shm_com_size << std::endl;
  }

  ///@todo
  ///  Test if round robin or blocked
  ///  test if equal number of ranks per node

  int local_color = min_shm_rank;
  // std::cout << "Rank = " << mpi_rank << " local_color = " << local_color
  //           << std::endl;
  int key = mpi_rank;
  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, local_color, key, &local_comm));

  return local_comm;
}

inline MPI_Comm build_node_remote_comm(MPI_Comm local_comm) {
  int mpi_rank(0), remote_color(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_rank(local_comm, &remote_color));

  MPI_Comm remote_comm;
  int key = mpi_rank;
  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, remote_color, key, &remote_comm));

  return remote_comm;
}

inline MPI_Comm build_node_bremote_comm(MPI_Comm local_comm) {
  int mpi_rank(0), local_offset(0), local_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_rank(local_comm, &local_offset));
  CHK_MPI(MPI_Comm_size(local_comm, &local_size));

  int node_offset = (mpi_rank / local_size) % local_size;

  int color = std::max(local_offset, node_offset) * local_size +
              std::min(local_offset, node_offset);
  int key = mpi_rank;
  MPI_Comm bremote_comm;

  // std::cout << "Rank = " << mpi_rank << " color = " << color
  //           << std::endl
  //           << "\t" << "node_offset = " << node_offset
  //           << ", node_offset % local_size = " << node_offset % local_size
  //           << ", local_offset = " << local_offset
  //           << std::endl;

  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, color, key, &bremote_comm));

  return bremote_comm;
}