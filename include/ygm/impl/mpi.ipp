// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once
#include <iostream>
#include <sstream>

namespace ygm {
const comm& comm_world() {
  static comm world(MPI_COMM_WORLD);
  return world;
}
const comm& comm_nl() {
  static comm nl(comm::split::nl);
  return nl;
}
const comm& comm_nr() {
  static comm nr(comm::split::nr);
  return nr;
}
const comm& comm_nlnr() {
  static comm nlnr(comm::split::nlnr);
  return nlnr;
}

inline comm::comm(MPI_Comm comm)
    : m_comm(comm), free_comm_in_destructor(false) {
  init_rank_size();
}

inline comm::comm(split _split) : free_comm_in_destructor(true) {
  if (_split == split::nl) {
    init_nl_comm();
  } else if (_split == split::nr) {
    init_nr_comm();
  } else if (_split == split::nlnr) {
    init_nlnr_comm();
  } else {
    std::cerr << "comm::comm(split) ERROR" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  init_rank_size();
}

inline comm::~comm() {
  if (free_comm_in_destructor) {
    int is_finalized;
    MPI_Finalized(&is_finalized);
    if (!is_finalized) {
      chk_ret(MPI_Comm_free(&m_comm));
    }
  }
}

inline void comm::init_rank_size() {
  chk_ret(MPI_Comm_size(m_comm, &m_size));
  chk_ret(MPI_Comm_rank(m_comm, &m_rank));
}

inline void comm::chk_ret(int ret, const char* loc) const {
  if (ret != MPI_SUCCESS) {
    char estring[MPI_MAX_ERROR_STRING];
    int  len(0);
    MPI_Error_string(ret, estring, &len);
    if (m_rank == 0) {
      std::cerr << "MPI ERROR = " << estring;
      if (loc) {
        std::cerr << ", LOCATION = " << loc;
      }
      std::cerr << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}

inline void comm::init_nl_comm() {
  /// @warning  This algorithm is O(P) because its general purpose.   Could be
  ///           Specialized to improve performance.
  ///

  int world_rank = comm_world().rank();
  int world_size = comm_world().size();

  MPI_Comm local_comm;

  int         min_shm_rank = -1;
  int         shm_com_size = 0;
  const char* shm_name     = "init_nl_comm";
  shm_unlink(shm_name);
  comm_world().barrier();
  //
  // This executes rank-by-rank sequentially.
  {
    // Blocks except for rank 0
    if (world_rank > 0) {
      CHK_MPI(MPI_Recv(NULL, 0, MPI_BYTE, world_rank - 1, 1, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE));
      // std::cout << "Rank " << world_rank << " Ready" << std::endl;
    }
    int       shm_fd            = -1;
    const int shm_size          = 4096;
    bool      this_rank_created = false;

    shm_fd = shm_open(shm_name, O_RDWR, 0666);
    if (shm_fd == -1) {
      // std::cout << "Rank " << world_rank << ": Failed to open;  attempting
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
    void* ptr =
        mmap(0, shm_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED) {
      printf("Map failed\n");
      exit(-1);
    }

    std::pair<int, int>* p_min_rank_size = (std::pair<int, int>*)ptr;
    if (this_rank_created) {
      p_min_rank_size->first  = world_rank;
      p_min_rank_size->second = 1;
    } else {
      p_min_rank_size->first = std::min(p_min_rank_size->first, world_rank);
      p_min_rank_size->second++;
    }

    //
    // Notifies rank+1 of completion
    if (world_rank < world_size - 1) {
      // std::cout << "Rank " << world_rank << " COMPLETE" << std::endl;
      CHK_MPI(MPI_Send(NULL, 0, MPI_BYTE, world_rank + 1, 1, MPI_COMM_WORLD));
    }

    // All ranks have completed.  Each pulls their node loacal's data
    comm_world().barrier();
    min_shm_rank = p_min_rank_size->first;
    shm_com_size = p_min_rank_size->second;

    //
    // Close shared segment
    if (munmap(ptr, shm_size) != 0) {
      std::cerr << "Rank " << world_rank << "munmap failed" << std::endl;
    }
    if (close(shm_fd) != 0) {
      std::cerr << "Rank " << world_rank << "close failed" << std::endl;
    }
    comm_world().barrier();
    if (this_rank_created && shm_unlink(shm_name) != 0) {
      std::cerr << "Rank " << world_rank << "shm_unlink failed" << std::endl;
    }
    // std::cout << "Rank " << world_rank << ", min_shm_rank = " << min_shm_rank
    //           << ", shm_com_size = " << shm_com_size << std::endl;
  }
  ///@todo
  ///  Test if round robin or blocked
  ///  test if equal number of ranks per node

  int local_color = min_shm_rank;
  // std::cout << "Rank = " << world_rank << " local_color = " << local_color
  //           << std::endl;
  int key = world_rank;
  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, local_color, key, &m_comm));
}

inline void comm::init_nr_comm() {
  int world_rank   = comm_world().rank();
  int remote_color = comm_nl().rank();

  int key = world_rank;
  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, remote_color, key, &m_comm));
}

inline void comm::init_nlnr_comm() {
  int world_rank   = comm_world().rank();
  int local_offset = comm_nl().rank();
  int local_size   = comm_nl().size();

  int node_offset = (world_rank / local_size) % local_size;

  int color = std::max(local_offset, node_offset) * local_size +
              std::min(local_offset, node_offset);
  int key = world_rank;

  CHK_MPI(MPI_Comm_split(MPI_COMM_WORLD, color, key, &m_comm));
  // std::cout << "Rank = " << world_rank << " color = " << color
  //           << std::endl
  //           << "\t" << "node_offset = " << node_offset
  //           << ", node_offset % local_size = " << node_offset % local_size
  //           << ", local_offset = " << local_offset
  //           << std::endl;
}

inline std::ostream& cout0() {
  static std::ostringstream dummy;
  dummy.clear();
  if (comm_world().rank() == 0) {
    return std::cout;
  }
  return dummy;
}

template <typename T>
inline void do_once(T func) {
  if (comm_world().rank() == 0) {
    func();
  }
}

namespace detail {
class init_final {
 public:
  init_final(int* argc, char*** argv) {
    int is_initialized;
    MPI_Initialized(&is_initialized);
    if (!is_initialized) {
      MPI_Init(argc, argv);
    }
  }
  ~init_final() {
    int is_finalized;
    MPI_Finalized(&is_finalized);
    if (!is_finalized) {
      MPI_Finalize();
    }
  }
};
}

inline void init(int* argc, char*** argv) {
  static detail::init_final i_f(argc, argv);
  comm_world();
  comm_nl();
  comm_nr();
  comm_nlnr();

  int rank            = comm_world().rank();
  int node_local_size = comm_nl().size();
  int node_local_rank = comm_nl().rank();

  if (rank < node_local_size) {
    if (rank != node_local_rank) {
      std::cerr << "ERROR:   Only blocked task mapping is supported"
                << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

inline std::string whoami() {
  std::stringstream sstr;
  sstr << "(R" << comm_world().rank() << ",N" << comm_nr().rank() << ",C"
       << comm_nl().rank() << ")";
  return sstr.str();
}

}  // namespace ygm