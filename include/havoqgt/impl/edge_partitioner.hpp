// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
#define __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
#include <map>
#include <deque>

namespace havoqgt {

/**
 * @class delegate_partitioned_graph
 * @details Put details here for class
 */

class source_partitioner {
 public:
  explicit source_partitioner(int p):m_mpi_size(p) { }
  int operator()(uint64_t i) const { return i % m_mpi_size; }

 private:
  int m_mpi_size;
};

template <typename edge_type>
class edge_source_partitioner {
 public:
  explicit edge_source_partitioner(int p):m_mpi_size(p) { }

  int operator()(edge_type i, bool is_counting) const {
    return std::get<0>(i) % m_mpi_size;
  }

 private:
  int m_mpi_size;
};

class edge_target_partitioner {
 public:
  explicit edge_target_partitioner(int p):m_mpi_size(p) { }
  int operator()(std::pair<uint64_t, uint64_t> i) const {
    return i.second % m_mpi_size;
  }

 private:
  int m_mpi_size;
};

/**
 * This class is used to determine where to send a high edge.
 * If the edge's destination is owned by another node, then the edge is sent to
 * that node. Otherwise, it is sent to a node based on the transfer_info.
 *
 * Transfer_info is a map of delgate ids (not vertex ids) to a dequeue.
 *   i.e. delgate_id = m_map_delegate_locator(vertex_id)
 * Each dequeue contains one or more OverflowSendInfo object which contains
 *   to_send_id, to_send_count, temp_to_send_count. which determines who will
 *   recieve the extra edges.
 */
 typedef struct OverflowSendInfo{
  OverflowSendInfo(int sid, int count)
    : to_send_id(sid)
    , to_send_count(count)
    , temp_to_send_count(0) {}

  int to_send_id;
  int32_t to_send_count;
  int32_t temp_to_send_count;
}OverflowSendInfo;

template <typename edge_type>
class high_edge_partitioner {
 public:
  explicit high_edge_partitioner(int s, int r,
    std::map<uint64_t, std::deque<OverflowSendInfo>> *transfer_info)
    : m_mpi_size(s)
    , m_mpi_rank(r)
    , m_transfer_info(transfer_info)
    /*, m_dof(dof)*/ { }

  /**
   * Determines where to send an edge.
   * If the edge's destination is owned by another node, then the edge is sent to
   * that node. Otherwise, it is sent to a node based on the transfer_info.
   *
   * @param  is_counting determines how to adjust the send_count variavles
   * @return the node to send the passed edge to.
   */
  int operator()(edge_type i, bool is_counting = true) {

    int dest = int(std::get<1>(i) % m_mpi_size);
    if (dest == m_mpi_rank) {
      // If the current node is the destination, then determine the destination
      // by examing the transfer_info object
      const uint64_t delegate_id = std::get<0>(i);
      if (m_transfer_info->count(delegate_id) == 0) {
        return m_mpi_rank;
      }

      assert(m_transfer_info->at(delegate_id).size() > 0);

      if (is_counting) {
        // If it is counting then use the temp_to_send_count, which is reset
        // the next time this called with an edge with the same delegate_id and
        // is_counting is set to false.
        for (size_t j = 0; j < m_transfer_info->at(delegate_id).size(); j++) {
          if (m_transfer_info->at(delegate_id)[j].temp_to_send_count <
                m_transfer_info->at(delegate_id)[j].to_send_count) {
            (m_transfer_info->at(delegate_id)[j].temp_to_send_count)++;

            // If this has an edge to send, return the send_id
            return m_transfer_info->at(delegate_id)[j].to_send_id;
          }
        }
        return m_mpi_rank;
      } else {
        // Not counting, so update the edge counts
        dest = m_transfer_info->at(delegate_id).front().to_send_id;
        //m_dof->send_of_delegate(delegate_id, dest);

        const int to_send_count =
          m_transfer_info->at(delegate_id).front().to_send_count--;
        assert(m_transfer_info->at(delegate_id).front().to_send_count >= 0);

        // Cleanup, if no more edges for this destination then remove it.
        // If no more edges for this delaget remove it!
        if (m_transfer_info->at(delegate_id).front().to_send_count == 0) {
          m_transfer_info->at(delegate_id).pop_front();
          if (m_transfer_info->at(delegate_id).size() == 0) {
            m_transfer_info->erase(delegate_id);
          }
        } else {
          // Otherwise reset the temp variable
          m_transfer_info->at(delegate_id).front().temp_to_send_count = 0;
        }
      }  // else not counting

    } else {  // if dest == rank
      // m_dof->send_delegate(delegate_id, dest);
    }

    assert(dest >= 0);
    assert(dest != m_mpi_rank);
    return dest;
  }  // operator()

 private:
  const int m_mpi_size;
  const int m_mpi_rank;
  std::map<uint64_t, std::deque<OverflowSendInfo>> * m_transfer_info;
};  // class high_edge_partitioner


class dest_pair_partitioner {
 public:
  template<typename T>
  int operator()(std::pair<int, T> i) const { return i.first; }
};



}  // namespace havoqgt
   //
#endif  // __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
