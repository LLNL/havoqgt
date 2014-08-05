#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED

/**
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */

#include <sys/mman.h>

#ifdef DEBUG
 #warning Debug MACRO is enabled.
#endif

#define FLUSH_BY_DIRTY_CHECK 1
// #define FLUSH_ROUND_ROBIN 1

namespace havoqgt {
namespace mpi {


class IOInfo {
 public:
  IOInfo() {
    init(mb_read_, mb_written_);
  }

  void init(int &r, int &w) {
    FILE *pipe;
    char str[250];
    pipe = popen("iostat -m | grep md0 2>&1", "r" );

    float temp;
    fscanf(pipe, "%s", str);
    fscanf(pipe, "%f %f %f", &temp, &temp, &temp);
    fscanf(pipe, "%d %d \n", &r, &w);
    pclose(pipe);
  };

  void log_diff(bool final = false) {
    int read = -1;
    int written = -1;
    init(read, written);

    if (final) {
      std::cout << "Total MB Read: " << (read - mb_read_) <<
                "\nTotal MB Written: " << (written - mb_written_) << std::endl << std::flush;
    } else {
      std::cout << "\tMB Read: " << (read - mb_read_) <<
                "\n\tMB Written: " << (written - mb_written_) << std::endl << std::flush;
    }
  };

 private:
  int mb_read_;
  int mb_written_;
};

class DebugOverFlow {
 public:
  DebugOverFlow(int mpi_size, int mpi_rank, int num_delgates)
    : m_mpi_size(mpi_size)
    , m_mpi_rank(mpi_rank)
    , m_num_delgates(num_delgates) {
      sent_to.resize(num_delgates, std::vector<int>(mpi_size, 0));
      recieved_from.resize(num_delgates, std::vector<int>(mpi_size, 0));
      sent_to_of.resize(num_delgates, std::vector<int>(mpi_size, 0));
      recieved_from_of.resize(num_delgates, std::vector<int>(mpi_size, 0));
    }

  void send_delegate(uint64_t delegate_id, int dest) {
    assert(dest < m_mpi_size);
    assert(dest >= 0);
    assert(delegate_id < m_num_delgates);
    assert(delegate_id >= 0);

    sent_to[delegate_id][dest]++;
  }

  void recv_delegate(uint64_t delegate_id, int source) {
    assert(source < m_mpi_size);
    assert(source >= 0);
    assert(delegate_id < m_num_delgates);
    assert(delegate_id >= 0);

    recieved_from[delegate_id][source]++;
  }

  void send_of_delegate(uint64_t delegate_id, int dest) {
    assert(dest < m_mpi_size);
    assert(dest >= 0);
    assert(delegate_id < m_num_delgates);
    assert(delegate_id >= 0);

    sent_to_of[delegate_id][dest]++;
  }

  void recv_of_delegate(uint64_t delegate_id, int source) {
    assert(source < m_mpi_size);
    assert(source >= 0);
    assert(delegate_id < m_num_delgates);
    assert(delegate_id >= 0);

    recieved_from_of[delegate_id][source]++;
  }

  std::string log() {
    std::string temp = "[" + std::to_string(m_mpi_rank) + "]\n";

    temp += "\tExchanged Delegates:";
    for (size_t i = 0; i < m_num_delgates; i++) {
      temp += "\t\t" + std::to_string(i) + ": [";
      for (size_t j =0; j < m_mpi_size; j++) {
        temp += std::to_string(sent_to[i][j]) + ", ";
      }
      temp += "\n";
    }

    temp += "\tRecieved Delegates:";
    for (size_t i = 0; i < m_num_delgates; i++) {
      temp += "\t\t" + std::to_string(i) + ": [";
      for (size_t j =0; j < m_mpi_size; j++) {
        temp += std::to_string(recieved_from[i][j]) + ", ";
      }
      temp += "\n";
    }

    temp += "\t Sent overflow delegates:";
    for (size_t i = 0; i < m_num_delgates; i++) {
      temp += "\t\t" + std::to_string(i) + ": [";
      for (size_t j =0; j < m_mpi_size; j++) {
        temp += std::to_string(sent_to_of[i][j]) + ", ";
      }
      temp += "\n";
    }

    temp += "\t Recv overflow delegates:";
    for (size_t i = 0; i < m_num_delgates; i++) {
      temp += "\t\t" + std::to_string(i) + ": [";
      for (size_t j =0; j < m_mpi_size; j++) {
        temp += std::to_string(recieved_from_of[i][j]) + ", ";
      }
      temp += "\n";
    }

    return temp;
  }

  const int m_mpi_rank, m_mpi_size, m_num_delgates;
  std::vector< std::vector<int> > sent_to, recieved_from;
  std::vector< std::vector<int> > sent_to_of, recieved_from_of;

};

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

class edge_source_partitioner {
 public:
  explicit edge_source_partitioner(int p):m_mpi_size(p) { }
  int operator()(std::pair<uint64_t, uint64_t> i, bool is_counting) const {
    return i.first % m_mpi_size;
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
class high_edge_partitioner {
 public:
  explicit high_edge_partitioner(int s, int r,
    std::map<uint64_t, std::deque<OverflowSendInfo>> *transfer_info
    /*, DebugOverFlow *dof */)
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
  int operator()(std::pair<uint64_t, uint64_t> i, bool is_counting = true) {


    int dest = int(i.second % m_mpi_size);
    if (dest == m_mpi_rank) {
      // If the current node is the destination, then determine the destination
      // by examing the transfer_info object
      const uint64_t delegate_id = i.first;
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
  //DebugOverFlow *m_dof;
};  // class high_edge_partitioner


class dest_pair_partitioner {
 public:
  template<typename T>
  int operator()(std::pair<int, T> i) const { return i.first; }
};

class local_source_id {
 public:
  explicit local_source_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(std::pair<T, T> i) const { return i.first / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class local_dest_id {
 public:
  explicit local_dest_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(std::pair<T, T> i) const { return i.second / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class get_local_id {
 public:
  explicit get_local_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(T i) const { return i / m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};



class owner_source_id {
 public:
  explicit owner_source_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(std::pair<T, T> i) const { return i.first % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class owner_dest_id {
 public:
  explicit owner_dest_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(std::pair<T, T> i) const { return i.second % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class get_owner_id {
 public:
  explicit get_owner_id(int p):m_mpi_size(p) {}
  template<typename T>
  int operator()(T i) const { return i % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

/**
 * This function iterates (1) through the edges and calculates the following:
 *
 *   m_local_incoming_count: For each vertedx that is assigned to this node it
 *   determines the number of incoming edges
 *
 *   high_vertex_count: the number of high edges generated b
 *
 *  @param mpi_comm: MPI communication
 *  @param unsorted_itr: The edge generator iteratior
 *  @param unsorted_itr_end: The end of the edge generator iterator
 *  @param global_hubs: A map of all hub vertex. This map is filled in this
 *  function
 *  @param delegate_degree_threshold: The mininum number of edges a high degree
 *  vertex has incomming.
 *  @return n/a
 */
template <typename SegmentManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegmentManager>::
count_edge_degrees(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hubs,
                 uint64_t delegate_degree_threshold) {
  double time_start = MPI_Wtime();
  using boost::container::map;
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  IOInfo *io_info_temp = new IOInfo();
  if (mpi_rank == 0) {
    std::cout << "Starting:  count_edge_degrees" << std::endl << std::flush;
  }


  uint64_t high_vertex_count(0);

  // Loop until no processor is producing edges
  while (!detail::global_iterator_range_empty(unsorted_itr,
        unsorted_itr_end, mpi_comm)) {
    std::vector<
      boost::container::map< int, std::pair<uint64_t, uint64_t> >
    > maps_to_send(m_mpi_size);
    int maps_to_send_element_count = 0;

    // Generate Enough information to send
    for (size_t i = 0; i < edge_chunk_size && unsorted_itr != unsorted_itr_end; i++) {
      // Update this vertex's outgoing edge count (first member of the pair)
      uint64_t local_id = local_source_id(m_mpi_size)(*unsorted_itr);
      int owner    = owner_source_id(m_mpi_size)(*unsorted_itr);
      if (owner == m_mpi_rank) {
        m_local_outgoing_count[local_id]++;
      } else {
        maps_to_send.at(owner)[local_id].first++;
      }

      // Update the vertex's incoming edge count (second member of the pair)
      local_id = local_dest_id(m_mpi_size)(*unsorted_itr);
      owner    = owner_dest_id(m_mpi_size)(*unsorted_itr);
      if (owner == m_mpi_rank) {
        m_local_incoming_count[local_id]++;
        if (m_local_incoming_count[local_id] == delegate_degree_threshold) {
          high_vertex_count++;
        }
      } else {
        int c = maps_to_send.at(owner)[local_id].second++;
        if (c == 0) {
          maps_to_send_element_count++;
        }
      }

      unsorted_itr++;
    }  // for until threshold is reached


    //Send Vertex degree count information to other nodes.
    send_vertex_info(mpi_comm, high_vertex_count, delegate_degree_threshold,
        maps_to_send, maps_to_send_element_count);

  }  // while more edges


  // Now, the m_local_incoming_count contains the total incoming and outgoing
  // edges for each vertex owned by this node.
  // Using this information we identify the hubs.
  std::vector<uint64_t> temp_hubs;
  temp_hubs.reserve(high_vertex_count);
  for (size_t i = 0; i < m_local_incoming_count.size(); i++) {
    // const uint64_t outgoing = m_local_outgoing_count[i];
    const uint64_t incoming = m_local_incoming_count[i];

    if (incoming >= delegate_degree_threshold) {
      const uint64_t global_id = (i * m_mpi_size) + m_mpi_rank;
      assert(global_id != 0);
      temp_hubs.push_back(global_id);
    }
  }

  assert(temp_hubs.size() == high_vertex_count);

  // Gather the hub liss and add them to the map,
  std::vector<uint64_t> vec_global_hubs;
  mpi_all_gather(temp_hubs, vec_global_hubs, mpi_comm);
  // Insert gathered global hubs to set
  global_hubs.insert(vec_global_hubs.begin(), vec_global_hubs.end());

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    const double total_time = time_end-time_start;
    std::cout << "count_edge_degrees time = " << total_time << std::endl << std::flush;
    io_info_temp->log_diff();
  }
}  // count_edge_degrees

/**
 * This function is used to send/recv information about vertexes during the
 * count_edge_degrees function.
 *
 * @param mpi_comm: the mpi_communication group
 * @param high_vertex_count: tracks the number of high vertices, used to
 *  determine how much space to allocate for the delegate degree info later on.
 * @param delegate_degree_threshold: The mininum number of edges a high degree
 *  vertex has incomming.
 * @param maps_to_send: a vector of maps of vertex ids to pairs of incoming and
 *  outgoing edge counts.
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
send_vertex_info(MPI_Comm mpi_comm, uint64_t& high_vertex_count,
  uint64_t delegate_degree_threshold, std::vector<
  boost::container::map< int, std::pair<uint64_t, uint64_t> >  >&
  maps_to_send, int maps_to_send_element_count) {


  int to_send_pos = 0;
  std::vector<uint64_t> to_send(maps_to_send_element_count*3, 0);
  std::vector<int> to_send_count(m_mpi_size, 0);

  assert(maps_to_send.size() == m_mpi_size);
  for (size_t i = 0; i < maps_to_send.size(); i++) {
    for (auto itr = maps_to_send[i].begin(); itr != maps_to_send[i].end(); itr++) {
      assert(to_send_pos < to_send.size());
      std::pair<int, std::pair<uint64_t, uint64_t>> triple = (*itr);
      to_send[to_send_pos++] = uint64_t(triple.first);
      to_send[to_send_pos++] = triple.second.first;
      to_send[to_send_pos++] = triple.second.second;
    }
    to_send_count[i] = maps_to_send[i].size()*3;
  }

  std::vector<uint64_t> to_recv;
  std::vector<int> out_recvcnts;

  mpi_all_to_all(to_send, to_send_count,to_recv, out_recvcnts, mpi_comm);

  for (size_t k = 0; k < to_recv.size(); ) {
    const uint64_t local_id = to_recv[k++];
    const uint64_t source_count = to_recv[k++];
    const uint64_t dest_count = to_recv[k++];
    assert(local_id < m_local_incoming_count.size());

    // If its not currently a high vertex but by adding this it becomes one
    // then increment high_vertex_count
    if (m_local_incoming_count[local_id] < delegate_degree_threshold
      && m_local_incoming_count[local_id] + dest_count >=
      delegate_degree_threshold) {

      high_vertex_count++;
    }
    m_local_outgoing_count[local_id] += source_count;
    m_local_incoming_count[local_id] += dest_count;
  }  // for each recieved element.

}  // send_vertex_info

/**
 * This function initlizes the transfer_info object by determining which nodes
 * needs extra edges and which nodes have extra edges. This information is
 * exchanged so that nodes that need more edges will know how many extra edges
 * it will recieve.
 *
 * @param mpi_comm: the mpi communication group
 * @param &owned_high_count: tracks the number of high edges owned by this node
 * it is updated when we give to or recieve edges from another node.
 * @param owned_low_count: The number of low edges owened by this node
 * @param transfer_info: used to track to whome and how many edges are given
 * to another node
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
calculate_overflow(MPI_Comm mpi_comm, uint64_t &owned_high_count,
    const uint64_t owned_low_count,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info) {
  double time_start = MPI_Wtime();
  IOInfo *io_info_temp = new IOInfo();
  if (m_mpi_rank == 0) {
    std::cout << "Starting:  calculate_overflow" << std::endl << std::flush;
  }
  //
  // Get the number of high and low edges for each node.
  //
  std::vector<uint64_t> low_count_per_rank, high_count_per_rank;
  mpi_all_gather(uint64_t(owned_low_count), low_count_per_rank, mpi_comm);
  mpi_all_gather(owned_high_count, high_count_per_rank, mpi_comm);

  // Determine the total of edges accorss all nodes.
  const uint64_t owned_total_edges = owned_high_count + owned_low_count;
  uint64_t global_edge_count = mpi_all_reduce(owned_total_edges,
      std::plus<uint64_t>(), mpi_comm);

  // Determine the desired number of edges at each node.
  const uint64_t target_edges_per_rank = global_edge_count / m_mpi_size;

  //
  // Compure the edge count exchange
  //
  int heavy_idx(0), light_idx(0);
  for(; heavy_idx < m_mpi_size && light_idx < m_mpi_size; ++heavy_idx) {

    while(low_count_per_rank[heavy_idx] + high_count_per_rank[heavy_idx]
            > target_edges_per_rank) {
      // while heavy_idx has edges to give
      const int64_t total_edges_low_idx = low_count_per_rank[light_idx]
                                  + high_count_per_rank[light_idx];
      if(total_edges_low_idx < target_edges_per_rank) {
        // if the low_idx needs edges

        if(high_count_per_rank[heavy_idx] == 0) {
          // If the heavy_idx has no more edges to give then break the while loop
          // causing the heavy_idx to be incremented.
          break;
        }

        // Determine the most that can be given.
        uint64_t max_to_offload = std::min(high_count_per_rank[heavy_idx],
            high_count_per_rank[heavy_idx] + low_count_per_rank[heavy_idx] -
            target_edges_per_rank);
        // Determine the most that can be recived
        uint64_t max_to_receive = target_edges_per_rank -
            high_count_per_rank[light_idx] - low_count_per_rank[light_idx];
        // Determine the most that can be moved
        uint64_t to_move = std::min(max_to_offload, max_to_receive);

        assert(to_move != 0);
        // Update the local count variables
        high_count_per_rank[heavy_idx]-=to_move;
        high_count_per_rank[light_idx]+=to_move;

        assert(heavy_idx != light_idx);
        if (heavy_idx == m_mpi_rank) { // This node is sending
          std::vector<uint64_t> send_list;
          // Generate a list of [delegate_ids, edge_counts] to send
          generate_send_list(send_list, to_move, light_idx, transfer_info);

          // Send the information
          int64_t send_len = send_list.size();
          CHK_MPI(MPI_Send(&send_len, 1, mpi_typeof(send_len), light_idx, 0, mpi_comm));
          CHK_MPI(MPI_Send(send_list.data(), send_len, mpi_typeof(to_move),
              light_idx, 0, mpi_comm));

          // Adjust the owned_high_count info.
          owned_high_count -= to_move;
        } else if (light_idx == m_mpi_rank) {  // This node is reciving
          MPI_Status status;
          int64_t recv_length;

          // Recieve the information
          CHK_MPI(MPI_Recv(&recv_length, 1, mpi_typeof(recv_length), heavy_idx,
              0, mpi_comm, &status));
          std::vector<uint64_t> recv_list(recv_length);
          CHK_MPI(MPI_Recv(recv_list.data(), recv_length, mpi_typeof(to_move),
               heavy_idx, 0, mpi_comm, &status));

          // Update my delagate edge counts.
          uint64_t sanity_count = 0;
          for (int i = 0; i < recv_length;) {
            const uint64_t vert_id = recv_list[i++];
            const uint64_t count = recv_list[i++];
            m_delegate_info[vert_id] += count;
            sanity_count += count;
          }

          // Adjust the owned_high_count info.
          owned_high_count += to_move;
          assert(sanity_count == to_move);
        } // else this node is not involved.
        MPI_Barrier(mpi_comm);
      } else {
        ++light_idx;
        if (light_idx == m_mpi_size) {
          break;
        }
      } // else
    }  // While
  }  // For

#ifdef DEBUG
  const uint64_t owned_total_edges2 = owned_high_count + owned_low_count;
  uint64_t sanity_global_edge_count = mpi_all_reduce(owned_total_edges2,
      std::plus<uint64_t>(), mpi_comm);
  assert(sanity_global_edge_count == global_edge_count);
  assert(owned_high_count == high_count_per_rank[m_mpi_rank]);
  uint64_t high_count_2 = 0;
  for (size_t i = 0; i < m_delegate_info.size()-1; i++) {
    high_count_2 += m_delegate_info[i];
  }

  assert(owned_high_count == high_count_2);

  std::vector<uint64_t> low_count_per_rank2, high_count_per_rank2;
  mpi_all_gather(uint64_t(owned_low_count), low_count_per_rank2, mpi_comm);
  mpi_all_gather(owned_high_count, high_count_per_rank2, mpi_comm);

  for (size_t i = 0; i < m_mpi_size; i++) {
    assert(low_count_per_rank2[i] == low_count_per_rank[i]);
    assert(high_count_per_rank2[i] == high_count_per_rank[i]);
  }

#endif

  double time_end = MPI_Wtime();
  if (m_mpi_rank == 0) {
    const double total_time = time_end-time_start;
    std::cout << "calculate_overflow time = " << total_time << std::endl << std::flush;
    io_info_temp->log_diff();
  }
}  // calculate overflow

/**
 * This function, in a non optimized manner, generates a list of vertexs and
 * counts to send to a node that needs more delegate edges.
 *
 * Edges are sent in the parition_high_degree function, this only tells the node
 * how many edges to expect for each delegate vertex.
 *
 * It can be improved by optimizing how selection is done.
 *
 * @param send_list: the location to store the [dest, count] pairs
 * @param num_send: the number of edges to needed to send
 * @param send_id: the location to send the edges. (used by transfer_info)
 * @transfer_info: A map of delegate vertex to a dequeue of [dest,count] pairs
 * used to track which nodes will be recieving extra high edges.
 *
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
generate_send_list(std::vector<uint64_t> &send_list, uint64_t num_send,
    int send_id,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info ) {
  // Tracking variables used to determine how much space to allocate.
  uint64_t send_count = 0;
  uint64_t ver_count = 0;
  for (uint64_t i = 0; i < m_delegate_info.size()-1 && send_count <  num_send;) {
    if (m_delegate_info[i] != 0) {
      send_count += m_delegate_info[i];
      ver_count++;
    }
    i++;
  }

  // Initilze the send_list
  send_list.reserve(ver_count * 2);
  send_count = 0;
  for (uint64_t i = 0; i < m_delegate_info.size()-1 && send_count <  num_send;) {
    if (m_delegate_info[i] != 0) {  // if there are edges to give
      uint64_t edges_to_give = m_delegate_info[i];

      if (send_count + edges_to_give > num_send) {  // reduce edges if it will
        edges_to_give = num_send - send_count;  // go over the num to send.
      }

      m_delegate_info[i] -= edges_to_give;  // update this node's edge count
      send_count += edges_to_give;  // update the send_count

      send_list.push_back(i);  // add the information to the send_list
      send_list.push_back(edges_to_give);


      //Add this informatio to the transfer_info map
      if (transfer_info.count(i) == 0 ) {
        transfer_info[i] = std::deque<OverflowSendInfo>();
      }
      transfer_info[i].push_back(OverflowSendInfo(send_id, edges_to_give));
    }  // if there are edges to give for this delegate
    i++;
  }
  assert(num_send == send_count);
}  // generate_send_list

/**
 * This function allocates and initlizes several of data members. It is called
 * after the count_edge_degrees function, which determined the size of these
 * data memebers.
 *
 * @param global_hubs: the set of hub vertices
 * @param delegate_degree_threshold: The edge limit when a vertex becomes a hub.
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
initialize_edge_storage(boost::unordered_set<uint64_t>& global_hubs,
  uint64_t delegate_degree_threshold) {

  //
  // Setup and Compute Low Edge CSR information
  //

  // Allocate the index into the low edge csr.
  // +2: because the m_max_vertex is indexible and the last position must hold
  // the number of low edges.
  //m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0)); moved to
  //constructor


  // Initilize the m_owned_info, bny iterating through owned vertexes and
  //  if it is now a hub, then it incremenets the edge count by the number of
  //  outgoing edges.
  uint64_t edge_count = 0;
  for (uint64_t vert_id = 0; vert_id < m_owned_info.size(); vert_id++) {
    const uint64_t outgoing = m_local_outgoing_count[vert_id];
    const uint64_t incoming = m_local_incoming_count[vert_id];

    m_owned_info[vert_id] = vert_info(false, 0, edge_count);

    if (incoming < delegate_degree_threshold) {
      edge_count += outgoing;
    } else {
      #ifdef DEBUG
        const uint64_t global_id = (vert_id * m_mpi_size) + m_mpi_rank;
        assert(global_id != 0);
        if (global_id < m_max_vertex) {
          // IF vert_id == size-1 then the above will be true
          // And this assert will hit incorrectly
          assert(global_hubs.count(global_id) != 0);
        }
      #endif
    }
  }  // for over m_owned_info

  // Allocate the low edge csr to accommdate the number of edges
  // This will be filled by the partion_low_edge function
  for (int i = 0; i < m_mpi_size; i++) {
    if (i == m_mpi_rank) {
      m_owned_targets.resize(edge_count, vertex_locator());
      flush_advise_vector_dont_need(m_owned_targets);
    }
    MPI_Barrier(m_mpi_comm);
  }


  //
  // Setup and Compute Hub information
  //
  std::vector<uint64_t> vec_sorted_hubs(global_hubs.begin(), global_hubs.end());
  std::sort(vec_sorted_hubs.begin(), vec_sorted_hubs.end());

  // Allocates and initilize the delegate (AKA hub) vertex infromation
  m_delegate_degree.resize(vec_sorted_hubs.size(), 0);
  flush_advise_vector_dont_need(m_delegate_degree);
  m_delegate_label.resize(vec_sorted_hubs.size());

  // Loop over the hub vertexes, initilizing the delegate_degree tracking
  // structures
  for(size_t i=0; i<vec_sorted_hubs.size(); ++i) {
    uint64_t t_local_id = vec_sorted_hubs[i] / uint64_t(m_mpi_size);
    int t_owner = uint32_t(vec_sorted_hubs[i] % uint32_t(m_mpi_size));
    vertex_locator new_ver_loc(true, i, t_owner);

    m_map_delegate_locator[vec_sorted_hubs[i]] = new_ver_loc;
    m_delegate_label[i] = vec_sorted_hubs[i];

    //
    // Tag owned delegates
    //
    if (t_owner == m_mpi_rank) {
      m_owned_info[t_local_id].is_delegate = 1;
      m_owned_info[t_local_id].delegate_id = i;
    }
  }  // for over vec_sorted_hubs
  // Allocate space for the delegate csr index.
  // This is initlized during the paritioning of the low edges and then adjusted
  // in initialize_delegate_target
  m_delegate_info.resize(m_map_delegate_locator.size()+1, 0);
  flush_advise_vector_dont_need(m_delegate_info);
  flush_advise_vector_dont_need(m_delegate_label);
}

/**
 * This function initlizes the member variables that are used to hold the high
 * degree edges.
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
initialize_delegate_target(int64_t edges_high_count) {
  // Currently, m_delegate_info holds the count of high degree edges assigned
  // to this node for each vertex.
  // Below converts it into an index into the m_delegate_targets array
  assert( m_delegate_info[ m_delegate_info.size()-1] == 0);
  int64_t edge_count = 0;
  for (size_t i=0; i < m_delegate_info.size(); i++) {
    uint64_t num_edges = m_delegate_info[i];
    m_delegate_info[i] = edge_count;
    edge_count += num_edges;
    assert(edge_count <= edges_high_count);
  }

  // Allocate space for storing high degree edges
  // This will be filled in the partion_high_degree function
  for (int i = 0; i < m_mpi_size; i++) {
    if (i == m_mpi_rank) {
      m_delegate_targets.resize(edges_high_count);
      flush_vector(m_delegate_targets);
      flush_vector(m_delegate_info);
      assert(edges_high_count == edge_count);
    }
    MPI_Barrier(m_mpi_comm);
  }


}  // initialize_delegate_target


/**
 * This function iterates (2) through the edges and sends the low degree edges
 * to the nodes that own them.
 *
 * At the same time it tracks the number of outgoing edges for each delegate
 * vertex and exchanges that information with the other nodes.
 *
 */
template <typename SegmentManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegmentManager>::
partition_low_degree_count_high(MPI_Comm mpi_comm,
                 InputIterator orgi_unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 uint64_t delegate_degree_threshold,
                 uint64_t &edges_high_count) {

  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);

  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

  IOInfo *io_info_temp = new IOInfo();
  if (mpi_rank == 0) {
    std::cout << "Starting:  partition_low_degree" << std::endl << std::flush;
  }

  // Temp Vector for storing offsets


  // Used to store high_edge count
  std::vector<uint64_t> tmp_high_count_per_rank(mpi_size,0);

  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;
  double last_loop_time = MPI_Wtime();


  for (int node_turn = 0; node_turn < node_partitions; node_turn++) {

    InputIterator unsorted_itr = orgi_unsorted_itr;
    if (m_mpi_rank == 0) {
      std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
        << " Node Turn: " << node_turn <<  " ): "<<
        " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;
    }
    if ( (m_mpi_rank % processes_per_node) % node_partitions == node_turn) {
      std::cout << "\t " << m_mpi_rank << " recieving edges (" <<
        node_turn << " % " << processes_per_node << " % " << node_partitions << " = "
        << ((m_mpi_rank % processes_per_node) % node_partitions) << ") " << std::endl << std::flush;
    }

    MPI_Barrier(mpi_comm);
    m_dont_need_graph();
    MPI_Barrier(mpi_comm);

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
          mpi_comm)) {

    if (m_mpi_rank == 0 && (loop_counter% 1000) == 0) {
      double cur_loop_time = MPI_Wtime();

      std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
        << " Node Turn: " << node_turn <<  " )Took: " << (cur_loop_time - last_loop_time) <<
        " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;

      last_loop_time = cur_loop_time;
    }
    loop_counter++;

    #if 0
    bool do_flush = true;
    #if FLUSH_ROUND_ROBIN
      const int processes_per_node = 24;
      static int check_id = m_mpi_rank  % processes_per_node;
      if (check_id-- == 0) {
        check_id = processes_per_node;
      } else {
        do_flush = false;
      }
    #elif FLUSH_BY_DIRTY_CHECK
      if (mpi_rank == 0) {
        do_flush = check_dirty_pages();
      }
      MPI_Bcast(&do_flush, 1, mpi_typeof(do_flush), 0, mpi_comm);
    #endif

    if (do_flush) {
      flush_advise_vector(m_owned_targets);
    }
    #endif

    // Generate Edges to Send
    std::vector<std::pair<uint64_t, uint64_t> > to_recv_edges_low;

    // Vector used to pass number of high edges
    std::vector<
      boost::container::map<uint64_t, uint64_t> > maps_to_send(m_mpi_size);
    int maps_to_send_element_count = 0;
    {
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_low;
      to_send_edges_low.reserve(edge_chunk_size);

      // for (size_t i=0; unsorted_itr != unsorted_itr_end && i < edge_chunk_size;
      //       ++unsorted_itr) {

      while (unsorted_itr != unsorted_itr_end &&
           to_send_edges_low.size()<edge_chunk_size) {
        // Get next edge
        const auto edge = *unsorted_itr;
        ++unsorted_itr;
        {
          const int owner = unsorted_itr->second % mpi_size;
          if ( (owner % processes_per_node) % node_partitions != node_turn) {
            continue;
          }
        }
        edge_counter++;


        if (global_hub_set.count(unsorted_itr->first) == 0) {
          to_send_edges_low.push_back(*unsorted_itr);
          //++i;
        } else if(global_hub_set.count(unsorted_itr->first)) {
          // This edge's source is a hub
          // 1) Increment the high edge count for the owner of the edge's dest
          tmp_high_count_per_rank[unsorted_itr->second % mpi_size]++;

          // 2) Increment the owner's count of edges for this hub.
          const int owner = unsorted_itr->second % mpi_size;
          if (owner == mpi_rank) {
            const uint64_t ver_id = unsorted_itr->first;

            const uint64_t new_source_id = m_map_delegate_locator[ver_id].local_id();
            assert(new_source_id < m_delegate_info.size()-1);
            m_delegate_info[new_source_id]++;
          } else {
            int c = maps_to_send.at(owner)[unsorted_itr->first]++;
            if (c == 0) {
              maps_to_send_element_count++;
            }
          }
        } else {
          assert(false);
        }
      }  // for

      // Exchange Edges/Recieve edges
      edge_source_partitioner paritioner(mpi_size);
      mpi_all_to_all_better(to_send_edges_low, to_recv_edges_low, paritioner,
          mpi_comm);

      // Send the hub edge count to the relevent nodes.
      send_high_info(mpi_comm, maps_to_send, maps_to_send_element_count);
    }

    std::sort(to_recv_edges_low.begin(), to_recv_edges_low.end());


#ifdef DEBUG
    // Sanity Check to make sure we recieve the correct edges
    for(size_t i=0; i<to_recv_edges_low.size(); ++i) {
      auto edge =  to_recv_edges_low[i];
      assert(int(edge.first % mpi_size) == mpi_rank);
      assert(global_hub_set.count(edge.first) == 0);
    }
#endif

    // Loop over recieved edges, appending them to the low CSR
    auto itr_end = to_recv_edges_low.end();
    for (auto itr = to_recv_edges_low.begin(); itr != itr_end; itr++) {

      auto edge = *itr;
      uint64_t new_vertex_id = local_source_id(m_mpi_size)(edge);
      assert(m_mpi_rank == int(edge.first % m_mpi_size));

      uint64_t temp_offset = (m_owned_info_tracker[new_vertex_id])++;
      uint64_t loc = temp_offset + m_owned_info[new_vertex_id].low_csr_idx;


      assert(loc <  m_owned_info[new_vertex_id+1].low_csr_idx);
      assert(!m_owned_targets[loc].is_valid());

      m_owned_targets[loc] = label_to_locator(edge.second);
    }  // for over recieved egdes
  }  // while global iterator range not empty
  }  // for node partition

  if (m_mpi_rank == 0) {
    double cur_loop_time = MPI_Wtime();
    std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
        << ")Took: " << (cur_loop_time - last_loop_time) <<
        " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;
  }

  edges_high_count = 0;
  for (size_t i = 0; i < m_delegate_info.size(); i++) {
    edges_high_count += m_delegate_info[i];
  }


#if DEBUG
  assert(m_delegate_info[m_delegate_info.size()-1] == 0);

  // Sync The high counts.
  std::vector<uint64_t> high_count_per_rank;
  mpi_all_reduce(tmp_high_count_per_rank, high_count_per_rank,
      std::plus<uint64_t>(), mpi_comm);

  uint64_t sanity_check_high_edge_count = high_count_per_rank[m_mpi_rank];
  assert(edges_high_count == sanity_check_high_edge_count);


#endif

  flush_advise_vector_dont_need(m_owned_targets);
  flush_advise_vector_dont_need(m_owned_info);
  flush_advise_vector_dont_need(m_owned_info_tracker);

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "partition_low_degree time = " << time_end - time_start
        << std::endl << std::flush;
    io_info_temp->log_diff();
  }
}  // partition_low_degree

/**
 * This function has each node exchange with one another their delegate vertex
 * incoming edge count.
 *
 * @param mpi_comm: the mpi communication group
 * @param maps_to_send: a vector of maps that map delegate_id to edge count.
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
send_high_info(MPI_Comm mpi_comm, std::vector< boost::container::map<
  uint64_t, uint64_t> >&maps_to_send, int maps_to_send_element_count) {

  int to_send_pos = 0;
  std::vector<uint64_t> to_send(
      maps_to_send_element_count*2, 0);
  std::vector<int> to_send_count(m_mpi_size, 0);

  assert(maps_to_send.size() == m_mpi_size);
  for (size_t i = 0; i < maps_to_send.size(); i++) {
    for (auto itr = maps_to_send[i].begin(); itr != maps_to_send[i].end(); itr++) {
      assert(to_send_pos < to_send.size());
      to_send[to_send_pos++] = itr->first;
      to_send[to_send_pos++] = itr->second;
    }
    to_send_count[i] = maps_to_send[i].size()*2;
  }

  std::vector<uint64_t> to_recv;
  std::vector<int> out_recvcnts;

  mpi_all_to_all(to_send, to_send_count,to_recv, out_recvcnts, mpi_comm);

  for (size_t i = 0; i < to_recv.size(); i++) {
    const uint64_t ver_id = to_recv[i++];
    const uint64_t delegate_dest_count = to_recv[i];

    const uint64_t new_source_id = m_map_delegate_locator[ver_id].local_id();
    assert(new_source_id < m_delegate_info.size()-1);
    m_delegate_info[new_source_id] += delegate_dest_count;
  }
}  // send_high_info


/**
 * This function iterates (3) over the edges and if an edge has a delegate
 * vertex as a source sends it to the apprioriate node based on the edge's
 * destination. If this node is giving edges to another node, then when this node has
 * enough edges for a delegate, extra edges are sent to a node that needs edges
 * for that delegate.
 *
 * @param mpi_comm: The mpi communication group
 * @param unsorted_itr: An iterator over edges
 * @param unsorted_itr_end: The end of the iterator over edges
 * @param global_hub_set: A set of all global hubs
 * @param transfer_info: A map of delegate_id to a deque of OverflowSendInfo
 * used to determine where overflowed edges go.
 */
template <typename SegmentManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegmentManager>::
partition_high_degree(MPI_Comm mpi_comm, InputIterator orgi_unsorted_itr,
    InputIterator unsorted_itr_end,
    boost::unordered_set<uint64_t>& global_hub_set,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info) {
  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

   IOInfo *io_info_temp = new IOInfo();
  if (mpi_rank == 0) {
    std::cout << "Starting partition_high_degree" << std::endl << std::flush;
  }
  // Initates the paritioner, which determines where overflowed edges go
  high_edge_partitioner paritioner(mpi_size, mpi_rank, &transfer_info);

  #ifdef DEBUG
    assert(m_delegate_degree.size() == m_delegate_info.size()-1);
    for (auto itr = m_delegate_degree.begin(); itr != m_delegate_degree.end();
        itr ++){
      assert((*itr) == 0);
    }

  #endif

  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;
  double last_loop_time = MPI_Wtime();


  // Scratch vector use for storing edges to send
  std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_high;
  to_send_edges_high.reserve(edge_chunk_size);

for (int node_turn = 0; node_turn < node_partitions; node_turn++) {

  if (m_mpi_rank == 0) {
    std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
      << " Node Turn: " << node_turn <<  " ): "<<
      " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;
  }

  MPI_Barrier(mpi_comm);
  m_dont_need_graph();
  MPI_Barrier(mpi_comm);

  InputIterator unsorted_itr = orgi_unsorted_itr;

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
        mpi_comm)) {


    if (m_mpi_rank == 0 && (loop_counter % 1000) == 0) {
      double cur_loop_time = MPI_Wtime();

      std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
        << " Node Turn: " << node_turn <<  " )Took: " << (cur_loop_time - last_loop_time) <<
        " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;

      last_loop_time = cur_loop_time;
    }
    loop_counter++;


#if 0
    bool do_flush = true;
    #ifdef FLUSH_ROUND_ROBIN
      const int processes_per_node = 24;
      static int check_id = m_mpi_rank  % processes_per_node;
      if (check_id-- == 0) {
        check_id = processes_per_node;
      } else {
        do_flush = false;
      }
    #elif defined FLUSH_BY_DIRTY_CHECK
      if (mpi_rank == 0) {
        do_flush = check_dirty_pages();
      }
      MPI_Bcast(&do_flush, 1, mpi_typeof(do_flush), 0, mpi_comm);
    #endif

    if (do_flush) {
      flush_advise_vector(m_delegate_targets);
    }
#endif


    while (unsorted_itr != unsorted_itr_end &&
           to_send_edges_high.size()<edge_chunk_size) {
      // Get next edge
      const auto edge = *unsorted_itr;
      ++unsorted_itr;

      {
          const int owner = unsorted_itr->second % mpi_size;
          if (owner % processes_per_node % node_partitions != node_turn) {
            continue;
          }
      }

      edge_counter++;

      if (global_hub_set.count(edge.first) == 1) {
        // If the edge's source is a hub node
        const uint64_t source_id = edge.first;
        uint64_t new_source_id = m_map_delegate_locator[source_id].local_id();
        assert(new_source_id >=0 && new_source_id < m_delegate_info.size()-1);

        // Send the edge if we don't own it or if we own it but have no room.
        to_send_edges_high.push_back(std::make_pair(new_source_id, edge.second));
      }  // end if is a hub
    }  // end while

    // Exchange edges we generated that we don't need with the other nodes and
    // recieve edges we may need
    // // Scratch vector use for storing recieved edges.
    std::vector< std::pair<uint64_t, uint64_t> > to_recv_edges_high;
    mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high, paritioner,
        mpi_comm);

    // Empty the vector
    {
      std::vector<std::pair<uint64_t, uint64_t> > temp;
      to_send_edges_high.swap(temp);
      to_send_edges_high.reserve(edge_chunk_size);
    }
    assert(to_send_edges_high.size() == 0);
    std::sort(to_recv_edges_high.begin(), to_recv_edges_high.end());

    for (size_t i=0; i<to_recv_edges_high.size(); ++i) {
      // Iterate over recieved edges, addiing them using similar logic from
      // above
      const auto edge = to_recv_edges_high[i];
      const uint64_t new_source_id = edge.first;
      assert(new_source_id >=0 && new_source_id < m_delegate_info.size()-1);

      uint64_t place_pos = m_delegate_degree[new_source_id];
      place_pos += m_delegate_info[new_source_id];

      if (place_pos == m_delegate_info[new_source_id+1]) {
        // We have no room for this node, so lets send it to a node that has
        // room.
        assert(transfer_info.size() > 0);
        assert(transfer_info.count(new_source_id) != 0);
        to_send_edges_high.push_back(edge);
      }
      else {
        assert(place_pos < m_delegate_info[new_source_id+1]);
        assert(place_pos < m_delegate_targets.size());

        uint64_t new_target_label = edge.second;
        m_delegate_targets[place_pos] = label_to_locator(new_target_label);
        assert(m_delegate_targets[place_pos].m_owner_dest < m_mpi_size);
        m_delegate_degree[new_source_id]++;

        if (owner_dest_id(m_mpi_size)(edge) != m_mpi_rank) {
          assert(transfer_info.size() == 0);
        }

      }  // else we have room
    }  // for edges recieved

  }  // end while get next edge
  } // For over node partitions
  if (m_mpi_rank == 0) {
    double cur_loop_time = MPI_Wtime();
    std::cout << "\t( Loops: " << loop_counter << "Edges: " << edge_counter
        << ")Took: " << (cur_loop_time - last_loop_time) <<
        " Dirty Pages: " << get_dirty_pages() << "kb." << std::endl << std::flush;
  }

  {//
  // Exchange edges we generated that we don't need with the other nodes and
    // recieve edges we may need
    // // Scratch vector use for storing recieved edges.
    std::vector< std::pair<uint64_t, uint64_t> > to_recv_edges_high;
    mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high, paritioner,
        mpi_comm);

    // Empty the vector
    {
      std::vector<std::pair<uint64_t, uint64_t> > temp;
      to_send_edges_high.swap(temp);
    }

    std::sort(to_recv_edges_high.begin(), to_recv_edges_high.end());
    for (size_t i=0; i<to_recv_edges_high.size(); ++i) {
      // Iterate over recieved edges, addiing them using similar logic from
      // above
      const auto edge = to_recv_edges_high[i];
      const uint64_t new_source_id = edge.first;
      assert(new_source_id >=0 && new_source_id < m_delegate_info.size()-1);

      uint64_t place_pos = m_delegate_degree[new_source_id];
      place_pos += m_delegate_info[new_source_id];

      if (place_pos == m_delegate_info[new_source_id+1]) {
        // We have no room for this node, so lets send it to a node that has
        // room. But this is the last round, so an error must have occurd
        assert(false);
      }
      else {
        assert(place_pos < m_delegate_info[new_source_id+1]);
        assert(place_pos < m_delegate_targets.size());

        uint64_t new_target_label = edge.second;
        m_delegate_targets[place_pos] = label_to_locator(new_target_label);
        assert(m_delegate_targets[place_pos].m_owner_dest < m_mpi_size);
        m_delegate_degree[new_source_id]++;

        if (owner_dest_id(m_mpi_size)(edge) != m_mpi_rank) {
          assert(transfer_info.size() == 0);
        }

      }  // else there is have room
    }  // for edges recieved
  }

#ifdef DEBUG
  int64_t recv_count2 = 0;
  for (size_t i = 0; i < m_delegate_degree.size(); i++) {
     recv_count2 += m_delegate_degree[i];
  }

  int64_t difference = recv_count2 - m_delegate_targets.size();

  for (size_t i = 0; i < m_delegate_info.size()-1; i++) {
    size_t pos = m_delegate_info[i];
    for (size_t j = pos; j < m_delegate_info[i+1]; j++) {
      assert(m_delegate_targets[j].m_owner_dest != 0xFFFFF);
      assert(m_delegate_targets[j].m_local_id != 0x7FFFFFFFFF);
    }
  }
#endif

  flush_advise_vector_dont_need(m_delegate_targets);
  flush_advise_vector_dont_need(m_delegate_info);
  flush_advise_vector_dont_need(m_delegate_degree);

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "partition_high_degree time = " << time_end - time_start
        << std::endl << std::flush;
    io_info_temp->log_diff();
  }
}  // partition_high_degre

/**
 * Builds a delegate_partitioned_graph with from and unsorted sequence of edges.
 *
 * @param sm       Pointer to segment manager
 * @param mpi_comm MPI communicator
 * @param Container input edges to partition
 * @param delegate_degree_threshold Threshold used to assign delegates
*/

template <typename SegmentManager>
template <typename Container>
delegate_partitioned_graph<SegmentManager>::
delegate_partitioned_graph(const SegmentAllocator<void>& seg_allocator,
                           MPI_Comm mpi_comm,
                           Container& edges, uint64_t max_vertex,
                           uint64_t delegate_degree_threshold,
                           std::function<void()> dont_need_graph)
    : m_global_edge_count(edges.size()),
      m_max_vertex(std::ceil(double(max_vertex) / double(m_mpi_size))),
      m_local_outgoing_count(seg_allocator),
      m_local_incoming_count(seg_allocator),
      m_owned_info(seg_allocator),
      m_owned_info_tracker(seg_allocator),
      m_owned_targets(seg_allocator),
      m_delegate_degree_threshold(delegate_degree_threshold),
      m_delegate_info(seg_allocator),
      m_delegate_degree(seg_allocator),
      m_delegate_label(seg_allocator),
      m_delegate_targets(seg_allocator),
      m_map_delegate_locator(100, boost::hash<uint64_t>(),
          std::equal_to<uint64_t>(), seg_allocator),
      m_controller_locators(seg_allocator)
       {

  MPI_Barrier(mpi_comm);
  m_mpi_comm = mpi_comm;
  m_dont_need_graph = dont_need_graph;
  double time_start = MPI_Wtime();
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank) );
  IOInfo *io_info_temp = new IOInfo();
  if (m_mpi_rank == 0) {
    std::cout << "Starting delegate_partitioned_graph: " << std::endl << std::flush;
  }
  assert(sizeof(vertex_locator) == 8);

  m_max_vertex = std::ceil(double(max_vertex) / double(m_mpi_size));

  m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0));
  // flush dont need
  flush_advise_vector_dont_need(m_owned_info);
  m_owned_info_tracker.resize(m_max_vertex+2, 0);
  flush_advise_vector_dont_need(m_owned_info_tracker);

  m_local_outgoing_count.resize(m_max_vertex+1, 0);
  flush_vector(m_local_outgoing_count);
  m_local_incoming_count.resize(m_max_vertex+1, 0);
  flush_vector(m_local_incoming_count);


  boost::unordered_set<uint64_t> global_hubs;

  if (m_mpi_rank == 0) {
    std::cout << "Initial Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;;
  }

  // Count Degree Information
  // For each owned vertex
  //  -count number of outgoing edges
  //  -count number of incoming edges
  // Generate global hubs information
  count_edge_degrees(mpi_comm, edges.begin(), edges.end(), global_hubs,
      delegate_degree_threshold);

  if (m_mpi_rank == 0) {
    std::cout << "After Count Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
  }


  // Using the above information construct the hub information, allocate space
  // for the low CSR
  if (m_mpi_rank == 0) {
    std::cout << "initialize_edge_storage " << std::endl << std::flush;
  }

  initialize_edge_storage(global_hubs, delegate_degree_threshold);
  flush_advise_vector_dont_need(m_local_outgoing_count);
  flush_advise_vector_dont_need(m_local_incoming_count);

  if (m_mpi_rank == 0) {
    std::cout << "After Initialize Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
  }

  MPI_Barrier(mpi_comm);
  if (m_mpi_rank == 0) {
    std::cout << "initialize_edge_storage complete" << std::endl << std::flush;
  }


  // Iterate (1) through the edges, sending all edges with a low degree source
  // vertex to the node that owns thats vertex
  // At the same time count the number of high edges and exchange that
  // information with the releveant node's owner.
  uint64_t edges_high_count;
  partition_low_degree_count_high(mpi_comm, edges.begin(), edges.end(),
      global_hubs, delegate_degree_threshold, edges_high_count);

  if (m_mpi_rank == 0) {
    std::cout << "After Low Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
  }
#if 0
  #define DELEGATE_PARTITION_GRAPH_IPP_PRINT_EDGE_COUNTS_TRANSFERS 1
  std::string temp = "[" + std::to_string(m_mpi_rank)+ "]\n";
  for (int i = 0; i < m_delegate_info.size(); i++) {
     temp += "[" + std::to_string(i) + ", "
            + std::to_string(m_delegate_info[i]) + "] ";
  }
  temp += "\n\tLow Edge Count:\t" + std::to_string(m_owned_targets.size());
  temp += "\tHigh Edge Count:\t" + std::to_string(edges_high_count);
  temp += "\tTotal:\t" + std::to_string(edges_high_count+m_owned_targets.size());
  temp += "\n";
#endif

// Calculate the overflow schedule, storing it into the transfer_info object
std::map< uint64_t, std::deque<OverflowSendInfo> > transfer_info;
calculate_overflow(mpi_comm, edges_high_count, m_owned_targets.size(),
    transfer_info);

if (m_mpi_rank == 0) {
  std::cout << "After Low Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
}

#ifdef DELEGATE_PARTITION_GRAPH_IPP_PRINT_EDGE_COUNTS_TRANSFERS
  for (int i = 0; i < m_delegate_info.size(); i++) {
    temp += "[" + std::to_string(i) + ", "
            + std::to_string(m_delegate_info[i]) + "] ";
  }
  temp += "\n\tLow Edge Count:\t" + std::to_string(m_owned_targets.size());
  temp += "\tHigh Edge Count:\t" + std::to_string(edges_high_count);
  temp += "\tTotal:\t" + std::to_string(edges_high_count+m_owned_targets.size());
  temp += "\n";

  // auto itr = transfer_info.begin();
  // auto itr_end = transfer_info.end();

  // for (; itr != itr_end; itr++) {
  //   auto list = (*itr).second;
  //   temp += "\t" + std::to_string((*itr).first) + ": ";
  //   for (size_t i = 0; i < list.size(); i++) {
  //     temp += "[" + std::to_string(list[i].to_send_id) + ", " + std::to_string(list[i].to_send_count) + "] ";
  //   }
  //   temp += "\n";
  // }

  for (int i = 0; i < m_mpi_size; i++) {
    if (i == m_mpi_rank) {
      std::cout << temp << std::flush;
    }
    for (int j = 0; j < 1000; j++) {
      if (i == m_mpi_rank) {
        std::cout << std::flush;
      }
      MPI_Barrier(mpi_comm);
    }
  }


#endif

  // Allocate and initilize the delegate edge CSR table and its index.
  initialize_delegate_target(edges_high_count); //flush/dont need the intenral vector

if(m_mpi_rank == 0) {
  std::cout << "After initilize delegate Dirty Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
}
  // Partition high degree, using overflow schedule
  partition_high_degree(mpi_comm, edges.begin(), edges.end(), global_hubs,
    transfer_info); //flush/dont need the intenral vector

if(m_mpi_rank == 0) {
  std::cout << "After partition high Pages:" <<  get_dirty_pages() << " kb" << std::endl << std::flush;
}
  MPI_Barrier(mpi_comm);
  double time_end = MPI_Wtime();
  if(m_mpi_rank == 0) {
    std::cout << "delegate_partitioned_graph time = " << time_end - time_start << std::endl << std::flush;
    io_info_temp->log_diff(true);
  }

  // We don't need this anymore, so free the resource
  free_edge_container<Container>(edges);


  uint64_t low_local_size      = m_owned_targets.size();
  uint64_t high_local_size     = m_delegate_targets.size();
  uint64_t total_local_size    = low_local_size + high_local_size;
  //
  // all-reduce hub degree
  {
    std::vector<uint64_t> my_hub_degrees(m_delegate_degree.begin(), m_delegate_degree.end());
    std::vector<uint64_t> tmp_hub_degrees;
    if(my_hub_degrees.size() > 0) {
      mpi_all_reduce(my_hub_degrees, tmp_hub_degrees, std::plus<uint64_t>(), mpi_comm);
      m_delegate_degree.clear();
      m_delegate_degree.insert(m_delegate_degree.end(),tmp_hub_degrees.begin(), tmp_hub_degrees.end());
    }
  }
  assert(m_delegate_degree.size() == m_delegate_label.size());

  //
  // Verify CSR integration properlly tagged owned delegates
  for (auto itr = m_map_delegate_locator.begin();
      itr != m_map_delegate_locator.end(); ++itr) {
    uint64_t label = itr->first;
    vertex_locator locator = itr->second;

    uint64_t local_id = label / uint64_t(m_mpi_size);
    if (label % uint64_t(m_mpi_size) == uint64_t(m_mpi_rank)) {
      assert(m_owned_info[local_id].is_delegate == 1);
      assert(m_owned_info[local_id].delegate_id == locator.local_id());
    }
  }


  //
  // Build controller lists
  for (size_t i=0; i < m_delegate_degree.size(); ++i) {
    if (int(i % m_mpi_size) == m_mpi_rank) {
      m_controller_locators.push_back(vertex_locator(true, i, m_mpi_rank));
    }
  }

  /*if(m_mpi_rank == 0) {
    for(size_t i=0; i<m_delegate_degree.size(); ++i) {
      std::cout << "Hub label = " << m_delegate_label[i] << ", degree = " << m_delegate_degree[i] << std::endl << std::flush;
    }
  }*/

  //
  //Print out debugging info
  uint64_t low_max_size       = mpi_all_reduce(low_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);
  uint64_t high_max_size      = mpi_all_reduce(high_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);
  //guint64_t overflow_max_size  = mpi_all_reduce(overflow_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);
  uint64_t total_max_size     = mpi_all_reduce(total_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);

  uint64_t low_sum_size       = mpi_all_reduce(low_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t high_sum_size      = mpi_all_reduce(high_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  //uint64_t overflow_sum_size  = mpi_all_reduce(overflow_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t total_sum_size     = mpi_all_reduce(total_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);

  uint64_t local_count_del_target = 0;
  for(uint64_t i=0; i<m_owned_targets.size(); ++i) {
    if(m_owned_targets[i].is_delegate()) ++local_count_del_target;
  }
  uint64_t total_count_del_target = mpi_all_reduce(local_count_del_target, std::plus<uint64_t>(), MPI_COMM_WORLD);

  if(m_mpi_rank == 0) {
    std::cout << "Max Vertex Id = " << max_vertex << std::endl << std::flush;
    std::cout << "Count of hub vertices = " << global_hubs.size() << std::endl << std::flush;
    std::cout << "Total percentage good hub edges = " << double(high_sum_size) / double(total_sum_size) * 100.0 << std::endl << std::flush;
    std::cout << "total count del target = " << total_count_del_target << std::endl << std::flush;
    std::cout << "Total percentage of localized edges = " << double(high_sum_size + total_count_del_target) / double(total_sum_size) * 100.0 << std::endl << std::flush;
    std::cout << "Global number of edges = " << total_sum_size << std::endl << std::flush;
    std::cout << "Number of small degree = " << low_sum_size << std::endl << std::flush;
    std::cout << "Number of hubs = " << high_sum_size << std::endl << std::flush;
    //std::cout << "Number of overfow = " << overflow_sum_size << std::endl << std::flush;
    std::cout << "oned imbalance = " << double(low_max_size) / double(low_sum_size/m_mpi_size) << std::endl << std::flush;
    std::cout << "hubs imbalance = " << double(high_max_size) / double(high_sum_size/m_mpi_size) << std::endl << std::flush;
    // if(overflow_sum_size > 0) {
    //   std::cout << "overflow imbalance = " << double(overflow_max_size) / double(overflow_sum_size/m_mpi_size) << std::endl << std::flush;
    // }
    std::cout << "TOTAL imbalance = " << double(total_max_size) / double(total_sum_size/m_mpi_size) << std::endl << std::flush;
  }
}

/**
 * @param  locator vertex_locator to convert
 * @return vertex label
 */
template <typename SegmentManager>
inline
uint64_t
delegate_partitioned_graph<SegmentManager>::
locator_to_label(delegate_partitioned_graph<SegmentManager>::vertex_locator
                  locator) const {
  uint64_t res;
  if(locator.is_delegate()) {
    res = m_delegate_label[locator.local_id()];
  } else {
    res = uint64_t(locator.local_id()) *
         uint64_t(m_mpi_size) +
         uint64_t(locator.owner());
  }

  return res;

}

/**
 * @param  label vertex label to convert
 * @return locator for the label
 */
template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_locator
delegate_partitioned_graph<SegmentManager>::
label_to_locator(uint64_t label) const {
  typename boost::unordered_map< uint64_t, vertex_locator,
              boost::hash<uint64_t>, std::equal_to<uint64_t>,
              SegmentAllocator< std::pair<uint64_t,vertex_locator> >
             >::const_iterator itr = m_map_delegate_locator.find(label);

  if(itr == m_map_delegate_locator.end()) {
    uint32_t owner    = label % uint64_t(m_mpi_size);
    uint64_t local_id = label / uint64_t(m_mpi_size);
    return vertex_locator(false, local_id, owner);
  }
  return itr->second;
}

/**
 * @details Gather operations performed when at least one process has
 *         found new local hubs
 * @param  local_hubs            set of local hubs
 * @param  global_hubs           set of global hubs to be updated
 * @param  found_new_hub_locally true, if new local hub has been found
 */
template <typename SegmentManager>
inline void
delegate_partitioned_graph<SegmentManager>::
sync_global_hub_set(const boost::unordered_set<uint64_t>& local_hubs,
                         boost::unordered_set<uint64_t>& global_hubs,
                         bool local_change, MPI_Comm mpi_comm) {
  uint32_t new_hubs = mpi_all_reduce(uint32_t(local_change),
                                     std::plus<uint32_t>(), mpi_comm);

  if(new_hubs > 0) {
    std::vector<uint64_t> vec_local_hubs(local_hubs.begin(), local_hubs.end());
    std::vector<uint64_t> vec_global_hubs;
    // global gather local hubs
    mpi_all_gather(vec_local_hubs, vec_global_hubs, mpi_comm);
    // Insert gathered global hubs to set
    global_hubs.insert(vec_global_hubs.begin(), vec_global_hubs.end());
  }
}

/**
 * @param  locator Vertex locator
 * @return Begin Edge Iterator
 */
template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator
delegate_partitioned_graph<SegmentManager>::
edges_begin(delegate_partitioned_graph<SegmentManager>::vertex_locator
             locator) const {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_info.size()-1);
    return edge_iterator(locator, m_delegate_info[locator.local_id()], this);
  }
  assert(locator.owner() == m_mpi_rank);
  assert(locator.local_id() < m_owned_info.size());
  return edge_iterator(locator,
                       m_owned_info[locator.local_id()].low_csr_idx,
                       this);
}

/**
 * @param  locator Vertex locator
 * @return End Edge Iterator
 */
template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator
delegate_partitioned_graph<SegmentManager>::
edges_end(delegate_partitioned_graph<SegmentManager>::vertex_locator
            locator) const {
  if(locator.is_delegate()) {
    assert(locator.local_id()+1 < m_delegate_info.size());
    return edge_iterator(locator, m_delegate_info[locator.local_id() + 1], this);
  }
  assert(locator.owner() == m_mpi_rank);
  assert(locator.local_id()+1 < m_owned_info.size());
  return edge_iterator(locator, m_owned_info[locator.local_id() + 1].low_csr_idx, this);
}

/**
 * @param  locator Vertex locator
 * @return Vertex degree
 */
template <typename SegmentManager>
inline
uint64_t
delegate_partitioned_graph<SegmentManager>::
degree(delegate_partitioned_graph<SegmentManager>::vertex_locator
        locator) const {
  uint64_t local_id = locator.local_id();
  if(locator.is_delegate()) {
    assert(local_id < m_delegate_degree.size());
    return m_delegate_degree[local_id];
  }
  assert(local_id + 1 < m_owned_info.size());
  return m_owned_info[local_id+1].low_csr_idx -
         m_owned_info[local_id].low_csr_idx;
}

/**
 * @param  locator Vertex locator
 * @return Vertex degree
 */
template <typename SegmentManager>
inline
uint64_t
delegate_partitioned_graph<SegmentManager>::
local_degree(delegate_partitioned_graph<SegmentManager>::vertex_locator
              locator) const {
  uint64_t local_id = locator.local_id();
  if(locator.is_delegate()) {
    assert(local_id + 1 < m_delegate_info.size());
    return m_delegate_info[local_id + 1] - m_delegate_info[local_id];
  }
  assert(local_id + 1 < m_owned_info.size());
  return m_owned_info[local_id+1].low_csr_idx -
         m_owned_info[local_id].low_csr_idx;
}


template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::
vertices_begin() const {
  return vertex_iterator(0,this);
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::
vertices_end() const {
  return vertex_iterator(m_owned_info.size()-1,this);
}

template <typename SegmentManager>
inline
bool
delegate_partitioned_graph<SegmentManager>::
is_label_delegate(uint64_t label) const {
  return m_map_delegate_locator.count(label) > 0;
}

template <typename SegmentManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegmentManager>::template vertex_data<
  T, SegManagerOther>*
delegate_partitioned_graph<SegmentManager>::
create_vertex_data(SegManagerOther* segment_manager_o,
    const char *obj_name) const {

  typedef typename delegate_partitioned_graph<SegmentManager>::template vertex_data<
  T, SegManagerOther> mytype;

  if (obj_name == nullptr) {
    return segment_manager_o->template construct<mytype>(bip::anonymous_instance)
        (m_owned_info.size(), m_delegate_info.size(), segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
        (m_owned_info.size(), m_delegate_info.size(), segment_manager_o);
  }
}

/**
 * @param   init initial value for each vertex
 */
template <typename SegmentManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegmentManager>::template vertex_data<
  T, SegManagerOther>*
delegate_partitioned_graph<SegmentManager>::
create_vertex_data(const T& init, SegManagerOther* segment_manager_o,
    const char *obj_name) const {

  typedef typename delegate_partitioned_graph<SegmentManager>::template vertex_data<
  T, SegManagerOther> mytype;

  if (obj_name == nullptr) {
    return segment_manager_o->template construct<mytype>(bip::anonymous_instance)
            (m_owned_info.size(), m_delegate_info.size(), init,
              segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
            (m_owned_info.size(), m_delegate_info.size(), init,
              segment_manager_o);
  }

}

template <typename SegmentManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegmentManager>::template edge_data<T, SegManagerOther>*
delegate_partitioned_graph<SegmentManager>::
create_edge_data(SegManagerOther* segment_manager_o,
    const char *obj_name) const {
  typedef typename delegate_partitioned_graph<SegmentManager>::template
                      edge_data<T, SegManagerOther> mytype;

  if (obj_name == nullptr) {
    return segment_manager_o->template construct<mytype>(bip::anonymous_instance)
            (m_owned_targets.size(), m_delegate_targets.size(),
              segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
            (m_owned_targets.size(), m_delegate_targets.size(),
              segment_manager_o);
  }
}

/**
 * @param   init initial value for each vertex
 */
template <typename SegmentManager>
template <typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::edge_data<T, SegManagerOther> *
delegate_partitioned_graph<SegmentManager>::
create_edge_data(const T& init, SegManagerOther * segment_manager_o,
    const char *obj_name) const {

  typedef delegate_partitioned_graph<SegmentManager>::
                      edge_data<T, SegManagerOther> mytype;

  if (obj_name == nullptr) {
    return segment_manager_o->template construct<mytype>(bip::anonymous_instance)
            (m_owned_targets.size(), m_delegate_targets.size(), init,
              segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
            (m_owned_targets.size(), m_delegate_targets.size(), init,
              segment_manager_o);
  }


}

///////////////////////////////////////////////////////////////////////////////
//                           Vertex Locator                                  //
///////////////////////////////////////////////////////////////////////////////
/**
 * @class  delegate_partitioned_graph::vertex_locator
 * @details Here are some very important details.
 */
/**
 *
 */
template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vertex_locator::
vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest) {
  m_is_bcast     = 0;
  m_is_intercept = 0;

  if (is_delegate) {
    m_is_delegate = true;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    assert(m_is_delegate == true
        && m_local_id    == local_id
        && m_owner_dest  == owner_dest);
  } else {
    m_is_delegate = false;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    assert(m_is_delegate == false
        && m_owner_dest  == owner_dest
        && m_local_id    == local_id);
  }
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::vertex_locator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::vertex_locator x) const {
  return m_is_delegate  == x.m_is_delegate
      && m_is_bcast     == x.m_is_bcast
      && m_is_intercept == x.m_is_intercept
      && m_owner_dest   == x.m_owner_dest
      && m_local_id     == x.m_local_id;
}



////////////////////////////////////////////////////////////////////////////////
//                               Edge Iterator                                //
////////////////////////////////////////////////////////////////////////////////
/**
 * \class delegate_partitioned_graph::edge_iterator
 * \details Put details here for class
 */
/**
 * @
 */
template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::edge_iterator::
edge_iterator(vertex_locator source,
              uint64_t edge_offset,
              const delegate_partitioned_graph* const pgraph)
  : m_source(source)
  , m_edge_offset(edge_offset)
  , m_ptr_graph(pgraph) { }

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator&
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator++() {
  ++m_edge_offset;
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::edge_iterator
delegate_partitioned_graph<SegmentManager>::edge_iterator::operator++(int) {
  edge_iterator to_return = *this;
  ++m_edge_offset;
  return to_return;
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::edge_iterator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x) const {
    assert(m_source      == x.m_source);
    assert(m_ptr_graph   == x.m_ptr_graph);
    return m_edge_offset == x.m_edge_offset;
}

template <typename SegmentManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& y) {
  return x.is_equal(y);

}

template <typename SegmentManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::edge_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_locator
delegate_partitioned_graph<SegmentManager>::edge_iterator::target() const {
  if(m_source.is_delegate()) {
    assert(m_edge_offset < m_ptr_graph->m_delegate_targets.size());
    assert(m_ptr_graph->m_delegate_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
    return m_ptr_graph->m_delegate_targets[m_edge_offset];
  }
  assert(m_edge_offset < m_ptr_graph->m_owned_targets.size());
  assert(m_ptr_graph->m_owned_targets[m_edge_offset].m_owner_dest <
          m_ptr_graph->m_mpi_size);
  return m_ptr_graph->m_owned_targets[m_edge_offset];
}

////////////////////////////////////////////////////////////////////////////////
//                             Vertex Iterator                                //
////////////////////////////////////////////////////////////////////////////////

template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
vertex_iterator(uint64_t index, const delegate_partitioned_graph<SegmentManager>*  pgraph)
  : m_ptr_graph(pgraph)
  , m_owned_vert_index(index) {
  update_locator();
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator&
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator++() {
  ++m_owned_vert_index;
  update_locator();
  return *this;
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::vertex_iterator::operator++(int) {
  vertex_iterator to_return = *this;
  ++m_owned_vert_index;
  update_locator();
  return to_return;
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
is_equal(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x) const {
  assert(m_ptr_graph        == x.m_ptr_graph);
  return m_owned_vert_index == x.m_owned_vert_index;
}

template <typename SegmentManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& y) {
  return x.is_equal(y);

}

template <typename SegmentManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegmentManager>::vertex_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegmentManager>
inline void
delegate_partitioned_graph<SegmentManager>::vertex_iterator::
update_locator() {
  for(; m_owned_vert_index < m_ptr_graph->m_owned_info.size()
        && m_ptr_graph->m_owned_info[m_owned_vert_index].is_delegate == true;
        ++ m_owned_vert_index);
  if(m_owned_vert_index < m_ptr_graph->m_owned_info.size()) {
    assert(m_ptr_graph->m_owned_info[m_owned_vert_index].is_delegate == false);
    uint32_t owner = m_ptr_graph->m_mpi_rank;
    m_locator = vertex_locator(false, m_owned_vert_index, owner);
  }
}

////////////////////////////////////////////////////////////////////////////////
//                                vert_info                                   //
////////////////////////////////////////////////////////////////////////////////

template <typename SegmentManager>
inline
delegate_partitioned_graph<SegmentManager>::vert_info::
vert_info(bool in_is_delegate, uint64_t in_delegate_id, uint64_t in_low_csr_idx)
  : is_delegate(in_is_delegate)
  , delegate_id(in_delegate_id)
  , low_csr_idx(in_low_csr_idx) {
  assert(is_delegate == in_is_delegate);
  assert(delegate_id == in_delegate_id);
  assert(low_csr_idx == in_low_csr_idx);
  assert(sizeof(vert_info) == 8);
}

////////////////////////////////////////////////////////////////////////////////
//                                vertex_data                                 //
////////////////////////////////////////////////////////////////////////////////
template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::vertex_data<T,SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, SegManagerOther* sm)
  : m_owned_vert_data(sm->template get_allocator<T>())
  , m_delegate_data(sm->template get_allocator<T>()) {
  m_owned_vert_data.resize(owned_data_size);
  m_delegate_data.resize(delegate_size);
  }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
  : m_owned_vert_data(owned_data_size, init, sm->template get_allocator<T>())
  , m_delegate_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
T&
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::
operator[](const vertex_locator& locator) {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
const T&
delegate_partitioned_graph<SegmentManager>::vertex_data<T, SegManagerOther>::operator[](const vertex_locator& locator) const {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

////////////////////////////////////////////////////////////////////////////////
//                                edge_data                                 //
////////////////////////////////////////////////////////////////////////////////
template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::edge_data<T,SegManagerOther>::
edge_data(uint64_t owned_size, uint64_t delegate_size, SegManagerOther* sm)
  : m_owned_edge_data(sm->template get_allocator<T>())
  , m_delegate_edge_data(sm->template get_allocator<T>()) {
  m_owned_edge_data.resize(owned_size);
  m_delegate_edge_data.resize(delegate_size);
  }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegmentManager>::edge_data<T, SegManagerOther>::
edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
  : m_owned_edge_data(owned_size, init, sm->template get_allocator<T>())
  , m_delegate_edge_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
T&
delegate_partitioned_graph<SegmentManager>::edge_data<T, SegManagerOther>::
operator[](const edge_iterator& itr) {
  if(itr.m_source.is_delegate()) {
    assert(itr.m_edge_offset < m_delegate_edge_data.size());
    return m_delegate_edge_data[itr.m_edge_offset];
  }
  assert(itr.m_edge_offset < m_owned_edge_data.size());
  return m_owned_edge_data[itr.m_edge_offset];
}

template <typename SegmentManager>
template<typename T, typename SegManagerOther>
const T&
delegate_partitioned_graph<SegmentManager>::edge_data<T, SegManagerOther>::operator[](const edge_iterator& itr) const {
  if(itr.m_source.is_delegate()) {
    assert(itr.m_edge_offset < m_delegate_edge_data.size());
    return m_delegate_edge_data[itr.m_edge_offset];
  }
  assert(itr.m_edge_offset < m_owned_edge_data.size());
  return m_owned_edge_data[itr.m_edge_offset];
}
} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
