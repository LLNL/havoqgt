
/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Steven Feldman <feldman12@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
#define __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
#include <map>
#include <deque>

namespace havoqgt {
namespace mpi {

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



}  // namespace mpi
}  // namespace havoqgt
   //
#endif  // __HAVOQGT_IMP_EDGE_PARTITIONER_HPP__
