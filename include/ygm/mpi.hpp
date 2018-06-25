#pragma once

#include <fcntl.h>
#include <mpi.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iostream>

#define CHK_MPI(a)                                               \
  {                                                              \
    if (a != MPI_SUCCESS) {                                      \
      char *error_string = NULL;                                 \
      int   len          = 0;                                    \
      MPI_Error_string(a, error_string, &len);                   \
      std::cerr << __FILE__ << ", line " << __LINE__             \
                << " MPI ERROR = " << error_string << std::endl; \
      exit(-1);                                                  \
    }                                                            \
  }

namespace ygm {

class comm;
const comm &comm_world();
const comm &comm_nl();
const comm &comm_nr();
const comm &comm_nlnr();
inline std::string whoami();

void init(int *argc, char ***argv);

class comm {
 public:
  enum class split { nl, nr, nlnr };

  comm(MPI_Comm comm = MPI_COMM_WORLD);

  comm(split _split);

  ~comm() {
    if (free_comm_in_destructor) {
      chk_ret(MPI_Comm_free(&m_comm));
    }
  }

  int      size() const { return m_size; }
  int      rank() const { return m_rank; }
  MPI_Comm mpi_comm() const { return m_comm; }

  void barrier() const { chk_ret(MPI_Barrier(m_comm)); }

 private:
  /// Not Copyable
  comm();
  comm(const comm &) = delete;             // non construction-copyable
  comm &operator=(const comm &) = delete;  // non copyable

  /// Checks MPI return codes

  void chk_ret(int ret, const char *loc = NULL) const;

  void init_nlnr_comm();
  void init_nl_comm();
  void init_nr_comm();
  void init_rank_size();

  MPI_Comm m_comm;
  int      m_size;
  int      m_rank;
  bool     free_comm_in_destructor;
};

}  // namespace ygm

#include <ygm/impl/mpi.ipp>