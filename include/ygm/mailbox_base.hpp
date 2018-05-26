#pragma once

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <ygm/default_init_alloc.hpp>
#include <fstream>
#include <iostream>
#include <ygm/mpi.hpp>
#include <numeric>
#include <vector>

using std::vector;

struct value_measurements {
  double ingest;
  double comm;
  double local_comm;
  double first_local_comm;
  double last_local_comm;
  double remote_comm;
  double wait;
  double local_wait;
  double first_local_wait;
  double last_local_wait;
  double remote_wait;
  double cpu;
  double local_cpu;
  double first_local_cpu;
  double last_local_cpu;
  double remote_cpu;
  uint64_t local_sent;
  uint64_t first_local_sent;
  uint64_t last_local_sent;
  uint64_t remote_sent;
  uint64_t local_recv;
  uint64_t first_local_recv;
  uint64_t last_local_recv;
  uint64_t remote_recv;
};

template <typename Data, typename RecvHandlerFunc>
class mailbox_base {
 public:
  mailbox_base(RecvHandlerFunc recv_func, size_t batch_size)
      : m_recv_func(recv_func), m_max_alloc(0), m_batch_size(batch_size) {
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size));
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank));
  }

  void set_communicators(MPI_Comm local_comm, MPI_Comm remote_comm) {
    m_local_comm = local_comm;
    m_remote_comm = remote_comm;
    CHK_MPI(MPI_Comm_size(m_local_comm, &m_local_size));
    CHK_MPI(MPI_Comm_rank(m_local_comm, &m_local_rank));
    CHK_MPI(MPI_Comm_size(m_remote_comm, &m_remote_size));
    CHK_MPI(MPI_Comm_rank(m_remote_comm, &m_remote_rank));
  }

  virtual void send(uint32_t dest, Data data) = 0;

  bool global_empty() { return do_exchange(); }

  ~mailbox_base() {
    uint64_t global_max_alloc;
    CHK_MPI(MPI_Allreduce(&m_max_alloc, &global_max_alloc, 1,
                          MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD));
    if (m_mpi_rank == 0) {
      std::cout << "m_max_alloc = " << global_max_alloc << std::endl;
      std::cout << "m_count_exchanges = " << m_count_exchanges << std::endl;
    }
  }

  virtual void print_status(std::string p) = 0;

  virtual void cleanup() = 0;

  virtual bool do_exchange() = 0;

  virtual value_measurements get_measurements() = 0;

  template <typename M, typename T>
  inline void push_to_buf(vector<M> &buf, vector<int> &counts, size_t &counter,
                          T &offset, const M &msg) {
    counts[offset]++;
    buf.push_back(msg);
    counter++;
  }

  template <typename M, typename T>
  inline void write_to_vec(vector<M, default_init_alloc<M>> &vec,
                           vector<int> &disps, vector<int> &counter, T &offset,
                           const M &msg) {
    vec[disps[offset] + counter[offset]] = msg;
    counter[offset]++;
  }

  inline void scale_addresses(vector<int> &counts, vector<int> &disps,
                              size_t scale) {
    for (size_t i = 0; i < counts.size(); ++i) {
      counts[i] *= scale;
      disps[i] *= scale;
    }
  }

  template <typename M>
  inline void alloc_recv_vec(vector<M, default_init_alloc<M>> &recv_vec,
                             vector<int> &recv_counts, size_t scale) {
    int total = std::accumulate(recv_counts.begin(), recv_counts.end(), 0);
    recv_vec.resize(total / scale);
    m_max_alloc = std::max(m_max_alloc, uint64_t(recv_vec.size()));
  }

  template <typename T>
  inline void calc_disps(const vector<T> &counts, vector<T> &disps) {
    // cacl send disps
    std::partial_sum(counts.begin(), counts.end(), disps.begin());
    for (size_t i = 0; i < disps.size(); ++i) {
      disps[i] -= counts[i];  // set to 0 offset
    }
  }

  inline void start_timer() { m_timer = MPI_Wtime(); }

  inline void end_timer(double &timer) { timer += MPI_Wtime() - m_timer; }

  RecvHandlerFunc m_recv_func;
  size_t m_batch_size;
  uint64_t m_max_alloc;
  int m_mpi_size;
  int m_mpi_rank;
  double m_timer = 0.0f;
  uint64_t m_count_exchanges = 0;

  MPI_Comm m_local_comm;
  int m_local_size;
  int m_local_rank;
  MPI_Comm m_remote_comm;
  int m_remote_size;
  int m_remote_rank;
};