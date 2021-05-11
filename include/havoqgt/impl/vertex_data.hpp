// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename Allocator>
template <typename T, typename AllocatorOther >
class delegate_partitioned_graph<Allocator>::vertex_data {
 public:
  typedef T value_type;
  vertex_data() {}

  T&       operator[] (const vertex_locator& locator)
  {
    if(locator.is_delegate()) {
      assert(locator.local_id() < m_delegate_data.size());
      return m_delegate_data[locator.local_id()];
    }
    assert(locator.local_id() < m_owned_vert_data.size());
    return m_owned_vert_data[locator.local_id()];
  }
  
  
  const T& operator[] (const vertex_locator& locator) const 
  {
    if(locator.is_delegate()) {
      assert(locator.local_id() < m_delegate_data.size());
      return m_delegate_data[locator.local_id()];
    }
    assert(locator.local_id() < m_owned_vert_data.size());
    return m_owned_vert_data[locator.local_id()];
  }

  void reset(const T& r) {
    for(size_t i=0; i<m_owned_vert_data.size(); ++i) {
      m_owned_vert_data[i] = r;
    }
    for(size_t i=0; i<m_delegate_data.size(); ++i) {
      m_delegate_data[i] = r;
    }
  }

  void clear() {
    for(size_t i=0; i<m_owned_vert_data.size(); ++i) {
      m_owned_vert_data[i].clear();
    }
    for(size_t i=0; i<m_delegate_data.size(); ++i) {
      m_delegate_data[i].clear();
    }
  }

  void all_reduce() {
    std::vector<T> tmp_in(m_delegate_data.begin(), m_delegate_data.end());
    std::vector<T> tmp_out(tmp_in.size(), 0);
    mpi_all_reduce(tmp_in, tmp_out, std::plus<T>(), MPI_COMM_WORLD);
    std::copy(tmp_out.begin(), tmp_out.end(), m_delegate_data.begin());
  }

  void all_max_reduce() {
    std::vector<T> tmp_in(m_delegate_data.begin(), m_delegate_data.end());
    std::vector<T> tmp_out(tmp_in.size(), 0);
    mpi_all_reduce(tmp_in, tmp_out, std::greater<T>(), MPI_COMM_WORLD);
    std::copy(tmp_out.begin(), tmp_out.end(), m_delegate_data.begin());
  }

  void all_min_reduce() {
    std::vector<T> tmp_in(m_delegate_data.begin(), m_delegate_data.end());
    std::vector<T> tmp_out(tmp_in.size(), 0);
    mpi_all_reduce(tmp_in, tmp_out, std::less<T>(), MPI_COMM_WORLD);
    std::copy(tmp_out.begin(), tmp_out.end(), m_delegate_data.begin());
  }
    
  vertex_data(const delegate_partitioned_graph& dpg, AllocatorOther allocate = AllocatorOther() )
    : m_owned_vert_data(allocate)
    , m_delegate_data(allocate) {
    m_owned_vert_data.resize(dpg.m_owned_info.size());
    m_delegate_data.resize(dpg.m_delegate_info.size());
    }
      
  T local_accumulate() const {
    T to_return = T();
    to_return = std::accumulate(m_owned_vert_data.begin(), m_owned_vert_data.end(), to_return);
    to_return = std::accumulate(m_delegate_data.begin(), m_delegate_data.end(), to_return);
    return to_return;
  }
  
  T global_accumulate() const {
    T local = local_accumulate();
    return mpi_all_reduce(local,std::plus<T>(), MPI_COMM_WORLD);
  }

 private:
  bip::vector<T, other_allocator<AllocatorOther, T>> m_owned_vert_data;
  bip::vector<T, other_allocator<AllocatorOther, T>> m_delegate_data;
};

}  // namespace havoqgt
#endif  // HAVOQGT_MPI_IMPL_VERTEX_DATA_HPP_
