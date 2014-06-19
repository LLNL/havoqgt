#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED

/**
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */





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

class edge_source_partitioner {
 public:
  explicit edge_source_partitioner(int p):m_mpi_size(p) { }
  int operator()(std::pair<uint64_t, uint64_t> i) const {
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
  T operator()(std::pair<T, T> i) const { return i.first % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class owner_dest_id {
 public:
  explicit owner_dest_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(std::pair<T, T> i) const { return i.second % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};

class get_owner_id {
 public:
  explicit get_owner_id(int p):m_mpi_size(p) {}
  template<typename T>
  T operator()(T i) const { return i % m_mpi_size; }
 private:
  uint64_t m_mpi_size;
};


template <typename SegementManager>
void
delegate_partitioned_graph<SegementManager>::
send_vertex_info(MPI_Comm mpi_comm, uint64_t& high_vertex_count,
  uint64_t delegate_degree_threshold, std::vector<
  boost::container::map< uint64_t, std::pair<uint64_t, uint64_t> >  >&
  maps_to_send) {
  // Send to Each Node in turn.
  for (int i = 0; i < m_mpi_size; ++i) {
    // start inner loop at i to avoid re-exchanging data
    for (int j = i; j < m_mpi_size; ++j) {
      int swap_id = -1;
      if (i == m_mpi_rank && j == m_mpi_rank) {
        continue;
      } else if (i == m_mpi_rank) {
        swap_id = j;
      } else if (j == m_mpi_rank) {
        swap_id = i;
      } else {
        continue;
      }

      assert(swap_id != m_mpi_rank);




      MPI_Status status;
      uint64_t send_count = maps_to_send[swap_id].size() * 3;
      uint64_t recv_count = send_count;

      // First determine the max buffer that we need
      // This way we only have to allocate it once
      CHK_MPI(MPI_Sendrecv_replace(&recv_count, 1,
          MPI_UNSIGNED_LONG, swap_id, 0, swap_id, 0, mpi_comm, &status) );

      std::vector<uint64_t> temp_buff(std::max(send_count, recv_count), 0);

      auto itr = maps_to_send[swap_id].begin();
      for (uint64_t pos = 0; itr != maps_to_send[swap_id].end(); itr++) {
        temp_buff[pos++] = (*itr).first;
        temp_buff[pos++] = (*itr).second.first;
        temp_buff[pos++] = (*itr).second.second;
      }

      assert(temp_buff.size() == std::max(send_count, recv_count));

      // std::cout << "[ " << m_mpi_rank << "]" << " Sending " << send_count <<
      //  " to " << swap_id << ". Recieving " << recv_count << " from " <<
      //  swap_id << ". " << std::endl;


      CHK_MPI(MPI_Sendrecv_replace(&(temp_buff[0]), temp_buff.size(),
            MPI_UNSIGNED_LONG, swap_id, 0, swap_id, 0, mpi_comm, &status));


      // Now update the recieved vector
      for (size_t j = 0; j < recv_count; ) {
        const uint64_t local_id = temp_buff[j++];
        const uint64_t source_count = temp_buff[j++];
        const uint64_t dest_count = temp_buff[j++];

        assert(local_id < m_local_degree_count.size());

        // If its not currently a high vertex but by adding this it becomes one
        // then increment high_vertex_count
        if (m_local_degree_count[local_id].second < delegate_degree_threshold
          && m_local_degree_count[local_id].second + dest_count >=
          delegate_degree_threshold) {
          high_vertex_count++;
        }
        m_local_degree_count[local_id].first += source_count;
        m_local_degree_count[local_id].second += dest_count;
      }  // for each recieved element.
    }  // for over nodes
  }  // for over nodes
}  // send_vertex_inf

template <typename SegementManager>
void
delegate_partitioned_graph<SegementManager>::
send_high_info(MPI_Comm mpi_comm, std::vector< boost::container::map<
  uint64_t, uint64_t> >&maps_to_send) {
  // Iterate over nodes
  for (int i = 0; i < m_mpi_size; ++i) {
    // start inner loop at i to avoid re-exchanging data
    for (int j = i; j < m_mpi_size; ++j) {
      int swap_id = -1;
      if (i == m_mpi_rank && j == m_mpi_rank) {
        continue;
      } else if (i == m_mpi_rank) {
        swap_id = j;
      } else if (j == m_mpi_rank) {
        swap_id = i;
      } else {
        continue;
      }

      assert(swap_id != m_mpi_rank);

      MPI_Status status;
      uint64_t send_count = maps_to_send[swap_id].size() * 2;
      uint64_t recv_count = send_count;

      // First determine the max buffer that we need
      // This way we only have to allocate it once
      CHK_MPI(MPI_Sendrecv_replace(&recv_count, 1,
          MPI_UNSIGNED_LONG, swap_id, 0, swap_id, 0, mpi_comm, &status) );
      std::vector<uint64_t> temp_buff(std::max(send_count, recv_count));

      auto itr = maps_to_send[swap_id].begin();
      for (uint64_t pos = 0; itr != maps_to_send[swap_id].end(); itr++) {
        temp_buff[pos++] = (*itr).first;
        temp_buff[pos++] = (*itr).second;
      }

      CHK_MPI(MPI_Sendrecv_replace(&(temp_buff[0]), temp_buff.size(),
            MPI_UNSIGNED_LONG, swap_id, 0, swap_id, 0, mpi_comm, &status));

      // Now update the received vector
      for (size_t j = 0; j < recv_count; ) {
        const uint64_t ver_id = temp_buff[j++];
        const uint64_t delegate_dest_count = temp_buff[j++];

        const uint64_t source_id = m_map_delegate_locator[ver_id].local_id();
        assert(source_id < m_delegate_info.size());
        m_delegate_info[source_id] += delegate_dest_count;
      }  // for each received element.
    }  // for over nodes
  }  // for over nodes

  auto itr = maps_to_send[m_mpi_rank].begin();
  for (uint64_t pos = 0; itr != maps_to_send[m_mpi_rank].end(); itr++) {
    const uint64_t ver_id = (*itr).first;
    const uint64_t delegate_dest_count = (*itr).second;

    const uint64_t new_source_id = m_map_delegate_locator[ver_id].local_id();
    assert(new_source_id < m_delegate_info.size());
    m_delegate_info[new_source_id] += delegate_dest_count;
  }
}  // send_high_info

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
count_degrees(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hubs,
                 uint64_t delegate_degree_threshold) {
  double time_start = MPI_Wtime();
  using boost::container::map;
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  if (mpi_rank == 0) {
    std::cout << "Starting:  count_degrees" << std::endl;
  }


  uint64_t high_vertex_count(0);

  // Loop until no processor is producing edges
  while (!detail::global_iterator_range_empty(unsorted_itr,
        unsorted_itr_end, mpi_comm)) {
    std::vector<
      boost::container::map< uint64_t, std::pair<uint64_t, uint64_t> >
    > maps_to_send(m_mpi_size);

    // Generate Enough information to send
    const size_t threshold = 1024;
    for (size_t i = 0; i < threshold && unsorted_itr != unsorted_itr_end; i++) {
      // First we get source info
      uint64_t local_id = local_source_id(m_mpi_size)(*unsorted_itr);
      uint64_t owner    = owner_source_id(m_mpi_size)(*unsorted_itr);

      if (owner == m_mpi_rank) {
        m_local_degree_count[local_id].first++;
      } else {
        maps_to_send[owner][local_id].first++;
      }

      // Then we get dest info
      local_id = local_dest_id(m_mpi_size)(*unsorted_itr);
      owner    = owner_dest_id(m_mpi_size)(*unsorted_itr);

      if (owner == m_mpi_rank) {
        m_local_degree_count[local_id].second++;
        if (m_local_degree_count[local_id].second == delegate_degree_threshold) {
          high_vertex_count++;
        }
      } else {
        maps_to_send[owner][local_id].second++;
      }
      unsorted_itr++;
    }  // for until threshold is reached

    send_vertex_info(mpi_comm, high_vertex_count, delegate_degree_threshold,
        maps_to_send);

  }  // while more edges

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "count_degree time = " << time_end-time_start << std::endl;
  }

  // Next
  // Figure out how many edges have a have a high degree source
  // Figure out how many edges have a low degree source
  std::vector<uint64_t> temp_local_hubs;
  temp_local_hubs.reserve(high_vertex_count);
  uint64_t sanity_high_vertex_count = 0;
  for (int i = 0; i < m_local_degree_count.size(); i++) {
    const uint64_t outgoing = m_local_degree_count[i].first;
    const uint64_t incoming = m_local_degree_count[i].second;

    if (incoming >= delegate_degree_threshold) {
      const uint64_t global_id = (i * m_mpi_size) + m_mpi_rank;
      assert(global_id != 0);
      temp_local_hubs.push_back(global_id);
      sanity_high_vertex_count++;
    }
  }

  assert(sanity_high_vertex_count == high_vertex_count);

  {
    //next sync the hub lists
    std::vector<uint64_t> vec_global_hubs;
    mpi_all_gather(temp_local_hubs, vec_global_hubs, mpi_comm);
    // Insert gathered global hubs to set
    global_hubs.insert(vec_global_hubs.begin(), vec_global_hubs.end());
  }


  // if(m_mpi_rank == 0) {
  //  std::vector<uint64_t> vec_hubs(global_hubs.begin(), global_hubs.end());
  //  std::sort(vec_hubs.begin(), vec_hubs.end());
  //  std::cout << "[" << m_mpi_rank << "]::: ";
  //   for (auto i = vec_hubs.begin(); i != vec_hubs.end(); i++) {
  //     std::cout << "[" << *i << "] ";
  //   }
  //   std::cout << std::endl;
  // }
}  // count_degres

template <typename SegementManager>
void
delegate_partitioned_graph<SegementManager>::
initialize_low_edge_storage(boost::unordered_set<uint64_t>& global_hubs,
  uint64_t delegate_degree_threshold) {
  //
  // Setup Low CSR
  //
  m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0));
  uint64_t edge_count = 0;

  // Initilize the m_owned_info
  for (uint64_t vert_id = 0; vert_id < m_owned_info.size(); vert_id++) {
    const uint64_t outgoing = m_local_degree_count[vert_id].first;
    const uint64_t incoming = m_local_degree_count[vert_id].second;

    m_owned_info[vert_id] = vert_info(false, 0, edge_count);

    if (incoming >= delegate_degree_threshold) {
      const uint64_t global_id = (vert_id * m_mpi_size) + m_mpi_rank;
      assert(global_id != 0);
      if (global_id < m_max_vertex) {
        // IF vert_id == size-1 then the above will be true
        // And this assert will hit incorrectly
        assert(global_hubs.count(global_id) != 0);
      }
    } else {

      edge_count += outgoing;
    }
  }
  m_owned_targets.resize(edge_count, vertex_locator());


  //
  // Setup and Compute Hub information
  //
  std::vector<uint64_t> vec_sorted_hubs(global_hubs.begin(), global_hubs.end());
  std::sort(vec_sorted_hubs.begin(), vec_sorted_hubs.end());

  m_delegate_degree.resize(vec_sorted_hubs.size(), 0);
  m_delegate_label.resize(vec_sorted_hubs.size());

  for(size_t i=0; i<vec_sorted_hubs.size(); ++i) {
    uint64_t t_local_id = vec_sorted_hubs[i] / uint64_t(m_mpi_size);
    uint32_t t_owner = uint32_t(vec_sorted_hubs[i] % uint32_t(m_mpi_size));
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
  }  // for items in vec_sorted_hubs
}

template <typename SegementManager>
void
delegate_partitioned_graph<SegementManager>::
allocate_high_edge_storage() {
  //
  // Setup High CSR
  //

  m_delegate_info.resize(m_map_delegate_locator.size()+1, 0);

}

template <typename SegementManager>
void
delegate_partitioned_graph<SegementManager>::
initialize_delegate_info(uint64_t edges_high_count) {
  //
  // Setup High CSR
  // Currently each position holds the number of edges.
  // But we need to turn that into index, so we replace it with a rolling total
  // then increment it by the num of edes
  // TODO(optimization): just add rolling total, then when filliung the targets
  // subtract one. We fill it in reverse, which means we dont need the temporary
  // vector to keep track of offsets.
  m_delegate_targets.resize(edges_high_count);
  uint64_t edge_count = 0;
  for (int i=0; i < m_delegate_info.size(); i++) {
    uint64_t num_edges = m_delegate_info[i];
    m_delegate_info[i] = edge_count;
    edge_count += num_edges;
    assert(edge_count <= edges_high_count);

  }
  assert(edges_high_count == m_delegate_targets.size());


}



template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
partition_low_degree_count_high(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 uint64_t delegate_degree_threshold,
                 uint64_t &edges_high_count) {

  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);

  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

  if (mpi_rank == 0) {
    std::cout << "Starting:  partition_low_degree" << std::endl;
  }

  // Temp Vector for storing offsets
  std::vector<uint64_t> m_owned_info_tracker(m_owned_info.size(), 0);

  // Used to store high_edge count
  std::vector<uint64_t> tmp_high_count_per_rank(mpi_size,0);

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
          mpi_comm)) {
    // Generate Edges to Send
    std::vector<std::pair<uint64_t, uint64_t> > to_recv_edges_low;

    // Vector used to pass number of high edges
    std::vector<
      boost::container::map<uint64_t, uint64_t> > maps_to_send(m_mpi_size);
    {
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_low;
      to_send_edges_low.reserve(16*1024);

      for (size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024;
            ++unsorted_itr) {
        if (global_hub_set.count(unsorted_itr->first) == 0) {
          to_send_edges_low.push_back(*unsorted_itr);
          ++i;
        }

        if(global_hub_set.count(unsorted_itr->first)) {
          // This edge's source is a hub
          // 1) Increment the high edge count for the owner of the edge's dest
          tmp_high_count_per_rank[unsorted_itr->second % mpi_size]++;

          // 2) Increment the owner's count of edges for this hub.
          maps_to_send[unsorted_itr->second % mpi_size][unsorted_itr->first]++;
        }

      }

      // Exchange Edges/Recieve edges
      mpi_all_to_all_better(to_send_edges_low, to_recv_edges_low,
            edge_source_partitioner(mpi_size), mpi_comm);

      // Send the hub edge count to the relevent nodes.
      send_high_info(mpi_comm, maps_to_send);
    }



    // Sanity Check to make sure we recieve the correct edges
    for(size_t i=0; i<to_recv_edges_low.size(); ++i) {
      auto edge =  to_recv_edges_low[i];
      assert(edge.first % mpi_size == mpi_rank);
      assert(global_hub_set.count(edge.first) == 0);
    }

    // Loop over recieved edges, appending them to the low CSR
    auto itr_end = to_recv_edges_low.end();
    for (auto itr = to_recv_edges_low.begin(); itr != itr_end; itr++) {

      auto edge = *itr;
      uint64_t new_vertex_id = local_source_id(m_mpi_size)(edge);
      assert(m_mpi_rank == edge.first % m_mpi_size);

      uint64_t temp_offset = (m_owned_info_tracker[new_vertex_id])++;
      uint64_t loc = temp_offset + m_owned_info[new_vertex_id].low_csr_idx;


      assert(loc <  m_owned_info[new_vertex_id+1].low_csr_idx);
      assert(!m_owned_targets[loc].is_valid());

      m_owned_targets[loc] = label_to_locator(edge.second);
    }  // for over recieved eges

  }  // while global iterator range not empty

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "partition_low_degree time = " << time_end - time_start
        << std::endl;
  }

  // Sync The high counts.
  std::vector<uint64_t> high_count_per_rank;
  mpi_all_reduce(tmp_high_count_per_rank, high_count_per_rank,
      std::plus<uint64_t>(), mpi_comm);

  edges_high_count = high_count_per_rank[m_mpi_rank];
  uint64_t sanity_check_high_edge_count = 0;
  for (int i = 0; i < m_delegate_info.size(); i++) {
    sanity_check_high_edge_count += m_delegate_info[i];
  }
  assert(edges_high_count == sanity_check_high_edge_count);
}  // partition_low_degree

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
partition_high_degree(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set) {
  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

  std::vector<uint64_t> m_delegate_info_offset(m_delegate_info.size(), 0);

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
        mpi_comm)) {

    std::vector< std::pair<uint64_t, uint64_t> > to_recv_edges_high;
    std::vector<  std::pair<int, std::pair<uint64_t, uint64_t> >  >
        to_recv_overflow;

    {  // send/recv block
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_high;
      to_send_edges_high.reserve(16*1024);

      for (size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024;
            ++unsorted_itr) {
        if (global_hub_set.count(unsorted_itr->first) == 1) {
          // IF it is a hub node
          ++i;
          to_send_edges_high.push_back(*unsorted_itr);
        }

      }
      mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high,
            edge_target_partitioner(mpi_size), mpi_comm);
    }  // end send/recv block

    for (size_t i=0; i<to_recv_edges_high.size(); ++i) {

      uint64_t source_id = to_recv_edges_high[i].first;
      uint64_t new_source_id = m_map_delegate_locator[source_id].local_id();
      m_delegate_degree[new_source_id]++;

      uint64_t place_pos = (m_delegate_info_offset[new_source_id])++;
      place_pos += m_delegate_info[new_source_id];
      assert(place_pos < m_delegate_targets.size());
      uint64_t new_target_label = to_recv_edges_high[i].second;

      m_delegate_targets[place_pos] = label_to_locator(new_target_label);
    }  // for edges recieved

  }  // end while get next edge
  double time_end = MPI_Wtime();

  if (mpi_rank == 0) {
    std::cout << "partition_high_degree time = " << time_end - time_start
        << std::endl;
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

template <typename SegementManager>
template <typename Container>
delegate_partitioned_graph<SegementManager>::
delegate_partitioned_graph(const SegmentAllocator<void>& seg_allocator,
                           MPI_Comm mpi_comm,
                           Container& edges, uint64_t max_vertex,
                           uint64_t delegate_degree_threshold
                           )
    : m_global_edge_count(edges.size()),
      m_local_degree_count(seg_allocator),
      m_owned_info(seg_allocator),
      m_owned_targets(seg_allocator),
      m_delegate_info(seg_allocator),
      m_delegate_degree(seg_allocator),
      m_delegate_label(seg_allocator),
      m_map_delegate_locator(100, boost::hash<uint64_t>(),
          std::equal_to<uint64_t>(), seg_allocator),
      m_controller_locators(seg_allocator),
      m_delegate_targets(seg_allocator),
      m_delegate_degree_threshold(delegate_degree_threshold) {

  MPI_Barrier(mpi_comm);
  double time_start = MPI_Wtime();
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank) );

  boost::unordered_set<uint64_t> global_hubs;

  assert(sizeof(vertex_locator) == 8);



  m_max_vertex = std::ceil(double(max_vertex) / double(m_mpi_size));

  m_local_degree_count.resize(m_max_vertex+1, std::make_pair(0,0));
  //
  // Count Degree Information
  // For each owned vertex
  //  -count number of outgoing edges
  //  -count number of incoming edges
  // Generate global hubs information

  count_degrees(mpi_comm, edges.begin(), edges.end(), global_hubs,
      delegate_degree_threshold);

  //
  // Using the above information construct the hub information, allocate space
  // for the low CSR
  initialize_low_edge_storage(global_hubs, delegate_degree_threshold);

  allocate_high_edge_storage();


  uint64_t edges_high_count;
  partition_low_degree_count_high(mpi_comm, edges.begin(), edges.end(),
      global_hubs, delegate_degree_threshold, edges_high_count);


  initialize_delegate_info(edges_high_count);

  //
  // Partition high degree, using overflow schedule
  partition_high_degree(mpi_comm, edges.begin(), edges.end(), global_hubs);



  MPI_Barrier(mpi_comm);
  double time_end = MPI_Wtime();
  if(m_mpi_rank == 0) {
    std::cout << "Partition time = " << time_end - time_start << std::endl;
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
    if (label % uint64_t(m_mpi_size) == m_mpi_rank) {
      assert(m_owned_info[local_id].is_delegate == 1);
      assert(m_owned_info[local_id].delegate_id == locator.local_id());
    }
  }


  //
  // Build controller lists
  for (size_t i=0; i < m_delegate_degree.size(); ++i) {
    if (i % m_mpi_size == m_mpi_rank) {
      m_controller_locators.push_back(vertex_locator(true, i, m_mpi_rank));
    }
  }

  /*if(m_mpi_rank == 0) {
    for(size_t i=0; i<m_delegate_degree.size(); ++i) {
      std::cout << "Hub label = " << m_delegate_label[i] << ", degree = " << m_delegate_degree[i] << std::endl;
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
    std::cout << "Max Vertex Id = " << max_vertex << std::endl;
    std::cout << "Count of hub vertices = " << global_hubs.size() << std::endl;
    std::cout << "Total percentage good hub edges = " << double(high_sum_size) / double(total_sum_size) * 100.0 << std::endl;
    std::cout << "total count del target = " << total_count_del_target << std::endl;
    std::cout << "Total percentage of localized edges = " << double(high_sum_size + total_count_del_target) / double(total_sum_size) * 100.0 << std::endl;
    std::cout << "Global number of edges = " << total_sum_size << std::endl;
    std::cout << "Number of small degree = " << low_sum_size << std::endl;
    std::cout << "Number of hubs = " << high_sum_size << std::endl;
    //std::cout << "Number of overfow = " << overflow_sum_size << std::endl;
    std::cout << "oned imbalance = " << double(low_max_size) / double(low_sum_size/m_mpi_size) << std::endl;
    std::cout << "hubs imbalance = " << double(high_max_size) / double(high_sum_size/m_mpi_size) << std::endl;
    // if(overflow_sum_size > 0) {
    //   std::cout << "overflow imbalance = " << double(overflow_max_size) / double(overflow_sum_size/m_mpi_size) << std::endl;
    // }
    std::cout << "TOTAL imbalance = " << double(total_max_size) / double(total_sum_size/m_mpi_size) << std::endl;
  }
}

/**
 * @param  locator vertex_locator to convert
 * @return vertex label
 */
template <typename SegementManager>
inline
uint64_t
delegate_partitioned_graph<SegementManager>::
locator_to_label(delegate_partitioned_graph<SegementManager>::vertex_locator
                  locator) const {
  if(locator.is_delegate()) {
    return m_delegate_label[locator.local_id()];
  }
  return uint64_t(locator.local_id()) *
         uint64_t(m_mpi_size) +
         uint64_t(locator.owner());
}

/**
 * @param  label vertex label to convert
 * @return locator for the label
 */
template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_locator
delegate_partitioned_graph<SegementManager>::
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
template <typename SegementManager>
inline void
delegate_partitioned_graph<SegementManager>::
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
template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::edge_iterator
delegate_partitioned_graph<SegementManager>::
edges_begin(delegate_partitioned_graph<SegementManager>::vertex_locator
             locator) const {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_info.size());
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
template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::edge_iterator
delegate_partitioned_graph<SegementManager>::
edges_end(delegate_partitioned_graph<SegementManager>::vertex_locator
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
template <typename SegementManager>
inline
uint64_t
delegate_partitioned_graph<SegementManager>::
degree(delegate_partitioned_graph<SegementManager>::vertex_locator
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
template <typename SegementManager>
inline
uint64_t
delegate_partitioned_graph<SegementManager>::
local_degree(delegate_partitioned_graph<SegementManager>::vertex_locator
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


template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_iterator
delegate_partitioned_graph<SegementManager>::
vertices_begin() const {
  return vertex_iterator(0,this);
}

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_iterator
delegate_partitioned_graph<SegementManager>::
vertices_end() const {
  return vertex_iterator(m_owned_info.size()-1,this);
}

template <typename SegementManager>
inline
bool
delegate_partitioned_graph<SegementManager>::
is_label_delegate(uint64_t label) const {
  return m_map_delegate_locator.count(label) > 0;
}

template <typename SegementManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegementManager>::template vertex_data<
  T, SegManagerOther>*
delegate_partitioned_graph<SegementManager>::
create_vertex_data(SegManagerOther* segment_manager_o,
    const char *obj_name) const {

  typedef typename delegate_partitioned_graph<SegementManager>::template vertex_data<
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
template <typename SegementManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegementManager>::template vertex_data<
  T, SegManagerOther>*
delegate_partitioned_graph<SegementManager>::
create_vertex_data(const T& init, SegManagerOther* segment_manager_o,
    const char *obj_name) const {

  typedef typename delegate_partitioned_graph<SegementManager>::template vertex_data<
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

template <typename SegementManager>
template <typename T, typename SegManagerOther>
typename delegate_partitioned_graph<SegementManager>::template edge_data<T, SegManagerOther>*
delegate_partitioned_graph<SegementManager>::
create_edge_data(SegManagerOther* segment_manager_o,
    const char *obj_name) const {
  typedef typename delegate_partitioned_graph<SegementManager>::template
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
template <typename SegementManager>
template <typename T, typename SegManagerOther>
delegate_partitioned_graph<SegementManager>::edge_data<T, SegManagerOther> *
delegate_partitioned_graph<SegementManager>::
create_edge_data(const T& init, SegManagerOther * segment_manager_o,
    const char *obj_name) const {

  typedef delegate_partitioned_graph<SegementManager>::
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
template <typename SegementManager>
inline
delegate_partitioned_graph<SegementManager>::vertex_locator::
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

template <typename SegementManager>
inline bool
delegate_partitioned_graph<SegementManager>::vertex_locator::
is_equal(const typename delegate_partitioned_graph<SegementManager>::vertex_locator x) const {
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
template <typename SegementManager>
inline
delegate_partitioned_graph<SegementManager>::edge_iterator::
edge_iterator(vertex_locator source,
              uint64_t edge_offset,
              const delegate_partitioned_graph* const pgraph)
  : m_source(source)
  , m_edge_offset(edge_offset)
  , m_ptr_graph(pgraph) { }

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::edge_iterator&
delegate_partitioned_graph<SegementManager>::edge_iterator::operator++() {
  ++m_edge_offset;
  return *this;
}

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::edge_iterator
delegate_partitioned_graph<SegementManager>::edge_iterator::operator++(int) {
  edge_iterator to_return = *this;
  ++m_edge_offset;
  return to_return;
}

template <typename SegementManager>
inline bool
delegate_partitioned_graph<SegementManager>::edge_iterator::
is_equal(const typename delegate_partitioned_graph<SegementManager>::edge_iterator& x) const {
    assert(m_source      == x.m_source);
    assert(m_ptr_graph   == x.m_ptr_graph);
    return m_edge_offset == x.m_edge_offset;
}

template <typename SegementManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegementManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegementManager>::edge_iterator& y) {
  return x.is_equal(y);

}

template <typename SegementManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegementManager>::edge_iterator& x,
           const typename delegate_partitioned_graph<SegementManager>::edge_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_locator
delegate_partitioned_graph<SegementManager>::edge_iterator::target() const {
  if(m_source.is_delegate()) {
    assert(m_edge_offset < m_ptr_graph->m_delegate_targets.size());
    return m_ptr_graph->m_delegate_targets[m_edge_offset];
  }
  assert(m_edge_offset < m_ptr_graph->m_owned_targets.size());
  return m_ptr_graph->m_owned_targets[m_edge_offset];
}

////////////////////////////////////////////////////////////////////////////////
//                             Vertex Iterator                                //
////////////////////////////////////////////////////////////////////////////////

template <typename SegementManager>
inline
delegate_partitioned_graph<SegementManager>::vertex_iterator::
vertex_iterator(uint64_t index, const delegate_partitioned_graph<SegementManager>*  pgraph)
  : m_ptr_graph(pgraph)
  , m_owned_vert_index(index) {
  update_locator();
}

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_iterator&
delegate_partitioned_graph<SegementManager>::vertex_iterator::operator++() {
  ++m_owned_vert_index;
  update_locator();
  return *this;
}

template <typename SegementManager>
inline
typename delegate_partitioned_graph<SegementManager>::vertex_iterator
delegate_partitioned_graph<SegementManager>::vertex_iterator::operator++(int) {
  vertex_iterator to_return = *this;
  ++m_owned_vert_index;
  update_locator();
  return to_return;
}

template <typename SegementManager>
inline bool
delegate_partitioned_graph<SegementManager>::vertex_iterator::
is_equal(const typename delegate_partitioned_graph<SegementManager>::vertex_iterator& x) const {
  assert(m_ptr_graph        == x.m_ptr_graph);
  return m_owned_vert_index == x.m_owned_vert_index;
}

template <typename SegementManager>
inline bool
operator==(const typename delegate_partitioned_graph<SegementManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegementManager>::vertex_iterator& y) {
  return x.is_equal(y);

}

template <typename SegementManager>
inline bool
operator!=(const typename delegate_partitioned_graph<SegementManager>::vertex_iterator& x,
           const typename delegate_partitioned_graph<SegementManager>::vertex_iterator& y) {
  return !(x.is_equal(y));
}

template <typename SegementManager>
inline void
delegate_partitioned_graph<SegementManager>::vertex_iterator::
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

template <typename SegementManager>
inline
delegate_partitioned_graph<SegementManager>::vert_info::
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
template <typename SegementManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegementManager>::vertex_data<T,SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, SegManagerOther* sm)
  : m_owned_vert_data(sm->template get_allocator<T>())
  , m_delegate_data(sm->template get_allocator<T>()) {
  m_owned_vert_data.resize(owned_data_size);
  m_delegate_data.resize(delegate_size);
  }

template <typename SegementManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegementManager>::vertex_data<T, SegManagerOther>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
  : m_owned_vert_data(owned_data_size, init, sm->template get_allocator<T>())
  , m_delegate_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename SegementManager>
template<typename T, typename SegManagerOther>
T&
delegate_partitioned_graph<SegementManager>::vertex_data<T, SegManagerOther>::
operator[](const vertex_locator& locator) {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

template <typename SegementManager>
template<typename T, typename SegManagerOther>
const T&
delegate_partitioned_graph<SegementManager>::vertex_data<T, SegManagerOther>::operator[](const vertex_locator& locator) const {
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
template <typename SegementManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegementManager>::edge_data<T,SegManagerOther>::
edge_data(uint64_t owned_size, uint64_t delegate_size, SegManagerOther* sm)
  : m_owned_edge_data(sm->template get_allocator<T>())
  , m_delegate_edge_data(sm->template get_allocator<T>()) {
  m_owned_edge_data.resize(owned_size);
  m_delegate_edge_data.resize(delegate_size);
  }

template <typename SegementManager>
template<typename T, typename SegManagerOther>
delegate_partitioned_graph<SegementManager>::edge_data<T, SegManagerOther>::
edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init, SegManagerOther* sm)
  : m_owned_edge_data(owned_size, init, sm->template get_allocator<T>())
  , m_delegate_edge_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename SegementManager>
template<typename T, typename SegManagerOther>
T&
delegate_partitioned_graph<SegementManager>::edge_data<T, SegManagerOther>::
operator[](const edge_iterator& itr) {
  if(itr.m_source.is_delegate()) {
    assert(itr.m_edge_offset < m_delegate_edge_data.size());
    return m_delegate_edge_data[itr.m_edge_offset];
  }
  assert(itr.m_edge_offset < m_owned_edge_data.size());
  return m_owned_edge_data[itr.m_edge_offset];
}

template <typename SegementManager>
template<typename T, typename SegManagerOther>
const T&
delegate_partitioned_graph<SegementManager>::edge_data<T, SegManagerOther>::operator[](const edge_iterator& itr) const {
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
