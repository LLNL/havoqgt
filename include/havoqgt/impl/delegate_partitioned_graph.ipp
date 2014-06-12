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
  source_partitioner(int p):m_mpi_size(p) { }
  int operator()(uint64_t i) const { return i % m_mpi_size; }
private:
  int m_mpi_size;
};

class edge_source_partitioner {
public:
  edge_source_partitioner(int p):m_mpi_size(p) { }
  int operator()(std::pair<uint64_t,uint64_t> i) const { return i.first % m_mpi_size; }
private:
  int m_mpi_size;
};

class edge_target_partitioner {
public:
  edge_target_partitioner(int p):m_mpi_size(p) { }
  int operator()(std::pair<uint64_t,uint64_t> i) const { return i.second % m_mpi_size; }
private:
  int m_mpi_size;
};

class dest_pair_partitioner {
public:
  template<typename T>
  int operator()(std::pair<int,T> i) const { return i.first; }
};

class local_source_id {
 public:
  local_source_id(int p):m_mpi_size(p) { }
  template<typename T>
  T operator()(T i) const { return i / m_mpi_size; }
  template<typename T>
  T operator()(std::pair<T,T> i) const { return i.first / m_mpi_size; }

 private:
  uint64_t m_mpi_size;
};

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
count_high_degree_transpose(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 std::vector<uint64_t>& high_count_per_rank)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );
  std::vector<uint64_t> tmp_high_count_per_rank(mpi_size,0);

  for(; unsorted_itr != unsorted_itr_end; ++unsorted_itr) {
    if(global_hub_set.count(unsorted_itr->first)) {
      tmp_high_count_per_rank[unsorted_itr->second % mpi_size]++;
    }
    //transpose version
    if(global_hub_set.count(unsorted_itr->second)) {
      tmp_high_count_per_rank[unsorted_itr->first % mpi_size]++;
    }
  }
  mpi_all_reduce(tmp_high_count_per_rank, high_count_per_rank, std::plus<uint64_t>(), mpi_comm);
}

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
count_low_degree(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 std::vector<uint64_t>& vertex_low_degree_count,
                 uint64_t delegate_degree_threshold) {
  double time_start = MPI_Wtime();
  using boost::container::map;
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );
  if (mpi_rank == 0) {
    std::cout << "Starting:  count_low_degree" << std::endl;
  }
  map<uint64_t,uint64_t> local_degree_count;



  uint64_t max_vertex_seen = 0;

  //
  // loop until all ranks have finished
  uint64_t gi=0;
  uint64_t gi_limit=4;
  while (!detail::global_iterator_range_empty(unsorted_itr,
          unsorted_itr_end, mpi_comm)) {
    // 1D exchange
    std::vector<uint64_t> to_recv_1d;
    std::vector<uint64_t> to_exchange_1d(0);
    ++gi;
    gi_limit *=2 ;
    for (uint64_t i=0;
          gi!=gi_limit && i<4*4096 && unsorted_itr != unsorted_itr_end;
          ++i, ++unsorted_itr, ++gi) {

      if (global_hub_set.count(unsorted_itr->first) == 0) {
        to_exchange_1d.push_back(unsorted_itr->first);
      }
      //transpose version
      if (global_hub_set.count(unsorted_itr->second) == 0) {
        to_exchange_1d.push_back(unsorted_itr->second);
      }

      max_vertex_seen = std::max(unsorted_itr->first, max_vertex_seen);
      max_vertex_seen = std::max(unsorted_itr->second, max_vertex_seen);
    }

    mpi_all_to_all_better(to_exchange_1d, to_recv_1d,
        source_partitioner(mpi_size), mpi_comm);
    std::vector<uint64_t>().swap(to_exchange_1d); //release memory

    //
    // Process new partitioned sources
    std::vector<uint64_t> new_hubs;
    for (size_t i=0; i<to_recv_1d.size(); ++i) {
      uint64_t source = to_recv_1d[i];
      assert(source % mpi_size == mpi_rank);
      if (global_hub_set.count(source) == 0) {
        ++local_degree_count[source];
        if (local_degree_count[source] == delegate_degree_threshold) {
          new_hubs.push_back(source);
        }

      }  // if in global hub set
    }  // for each vertex recieved

    //
    // Exchange newly found hubs
    std::vector<uint64_t> new_global_hubs;
    mpi_all_gather(new_hubs, new_global_hubs, mpi_comm);
    global_hub_set.insert(new_global_hubs.begin(), new_global_hubs.end());
  }

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "count_low_degree time = " << time_end-time_start << std::endl;
  }

  MPI_Allreduce(&max_vertex_seen, &m_max_vertex, 1, MPI_UNSIGNED_LONG,
        MPI_MAX, mpi_comm);

  m_max_vertex = local_source_id(m_mpi_size)(m_max_vertex);

  m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0));




  uint64_t cur_vertex_id = 0;
  m_owned_info[cur_vertex_id] = vert_info(false, 0, 0);


  auto degree_count_itr = local_degree_count.begin();
  auto degree_count_itr_end = local_degree_count.end();

  uint64_t edge_count = 0;
  for (; degree_count_itr != degree_count_itr_end; degree_count_itr++) {
    assert( (*degree_count_itr).first > cur_vertex_id);

    uint64_t raw_vertex_id = (*degree_count_itr).first;
    uint64_t new_vertex_id = local_source_id(m_mpi_size)(raw_vertex_id);
    uint64_t num_edges = (*degree_count_itr).second;

    if (num_edges >= delegate_degree_threshold) {
      num_edges = 0; // ie the edges are not stored in the low CSR
    }

    while (new_vertex_id > cur_vertex_id) {
      ++cur_vertex_id;
      m_owned_info[cur_vertex_id] = vert_info(false, 0, edge_count);
    }
    edge_count += num_edges;
  }
  cur_vertex_id++;
  for (; cur_vertex_id < m_owned_info.size(); cur_vertex_id++) {
    m_owned_info[cur_vertex_id] = vert_info(false, 0, edge_count);
  }

  // Allocate CSR Space based on the information gathered.
  m_owned_targets.resize(edge_count, vertex_locator());

}  // count_low_degree

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
partition_low_degree(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 std::deque<std::pair<uint64_t, uint64_t> >& edges_low)
{
  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);
  std::vector<uint64_t> m_owned_info_tracker(m_owned_info.size(), 0);

  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

  if (mpi_rank == 0) {
    std::cout << "Starting:  partition_low_degree" << std::endl;
  }

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
          mpi_comm)) {
    std::vector<std::pair<uint64_t, uint64_t> > to_recv_edges_low;
    {
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_low;
      to_send_edges_low.reserve(16*1024);
      for (size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024;
            ++unsorted_itr) {
        if (global_hub_set.count(unsorted_itr->first) == 0) {
          to_send_edges_low.push_back(*unsorted_itr);
          ++i;
        }
        //transpose version
        if(global_hub_set.count(unsorted_itr->second) == 0) {
          to_send_edges_low.push_back(std::make_pair(unsorted_itr->second,
                                                     unsorted_itr->first));
          ++i;
        }
      }
      mpi_all_to_all_better(to_send_edges_low, to_recv_edges_low,
            edge_source_partitioner(mpi_size), mpi_comm);
    }

    for(size_t i=0; i<to_recv_edges_low.size(); ++i) {
      assert(to_recv_edges_low[i].first % mpi_size == mpi_rank);
    }

    auto itr = to_recv_edges_low.begin();
    auto itr_end = to_recv_edges_low.end();
    while (itr != itr_end) {

      uint64_t raw_vertex_id = (*itr).first;
      uint64_t new_vertex_id = local_source_id(m_mpi_size)(raw_vertex_id);
      uint64_t dest = (*itr).second;

      itr++;

      uint64_t temp_offset = (m_owned_info_tracker[new_vertex_id])++;
      uint64_t place_pos = temp_offset + m_owned_info[new_vertex_id].low_csr_idx;

      assert(place_pos <  m_owned_info[new_vertex_id+1].low_csr_idx);
      assert(!m_owned_targets[place_pos].is_valid());

      m_owned_targets[place_pos] = label_to_locator(dest);
    }

    edges_low.insert(edges_low.end(), to_recv_edges_low.begin(),
          to_recv_edges_low.end());
  }  // while global iterator range not empty

  double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "partition_low_degree time = " << time_end - time_start
    		<< std::endl;
  }
}  // partition_low_degree

template <typename SegementManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegementManager>::
partition_high_degree(MPI_Comm mpi_comm,
                 InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 std::deque<std::pair<uint64_t, uint64_t> >& edges_high,
                 std::deque<std::pair<uint64_t, uint64_t> >& edges_high_overflow,
                 std::map<int, uint64_t>& overflow_schedule)
{
  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );

  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
  			mpi_comm)) {

    std::vector< std::pair<uint64_t, uint64_t> > to_recv_edges_high;
    std::vector<  std::pair<int, std::pair<uint64_t, uint64_t> >  >
    		to_recv_overflow;
    {
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_high;
      to_send_edges_high.reserve(16*1024);

      for (size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024;
      			++unsorted_itr) {
        if (global_hub_set.count(unsorted_itr->first) == 1) {
        	// IF it is a hub node
          ++i;
          to_send_edges_high.push_back(*unsorted_itr);
        }

        // transpose version
        if (global_hub_set.count(unsorted_itr->second) == 1) {
          ++i;
          to_send_edges_high.push_back(std::make_pair(unsorted_itr->second,
                                                   unsorted_itr->first));
        }
      }
      mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high,
      			edge_target_partitioner(mpi_size), mpi_comm);
    }

    for (size_t i=0; i<edges_high.size(); ++i) {
      assert(edges_high[i].second % mpi_size == mpi_rank);
    }

    {
      // Copy high edges to either edges_high or overflow
      std::vector<  std::pair< int, std::pair<uint64_t, uint64_t> >  >
      		to_send_overflow;

      for (size_t i=0; i<to_recv_edges_high.size(); ++i) {
        if (overflow_schedule.empty()) {
          edges_high.push_back(to_recv_edges_high[i]);
        } else {
          int dest = overflow_schedule.begin()->first;
          overflow_schedule[dest]--;

          if(overflow_schedule[dest] == 0) {
            overflow_schedule.erase(dest);
          }
          to_send_overflow.push_back(std::make_pair(dest,
          			to_recv_edges_high[i]));
        }
      }

      uint64_t global_overflow_size = mpi_all_reduce(
      			uint64_t(to_send_overflow.size()), std::plus<uint64_t>(), mpi_comm);

      if(global_overflow_size > 0) {
        mpi_all_to_all_better(to_send_overflow, to_recv_overflow,
        		dest_pair_partitioner(), mpi_comm);
      }
    }
    for(size_t i=0; i<to_recv_overflow.size(); ++i) {
      edges_high_overflow.push_back(to_recv_overflow[i].second);
    }
  }  // end while get next edge
  double time_end = MPI_Wtime();

  if (mpi_rank == 0) {
  	std::cout << "partition_high_degree time = " << time_end - time_start
  			<< std::endl;
  }
}



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
                           Container& edges,
                           uint64_t delegate_degree_threshold
                           )
    : m_owned_info(seg_allocator),
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

  std::deque< std::pair<uint64_t, uint64_t> > edges_low, edges_high, edges_high_overflow;
  boost::unordered_set<uint64_t> global_hubs;

  assert(sizeof(vertex_locator) == 8);

  std::vector<uint64_t> vertex_low_degree_count(1024);

  //
  // Count low degree & find delegates
  count_low_degree(mpi_comm, edges.begin(), edges.end(), global_hubs,
      vertex_low_degree_count, delegate_degree_threshold);

  //
  // Partition low degree vertices
  partition_low_degree(mpi_comm, edges.begin(), edges.end(), global_hubs, edges_low);


  std::vector<uint64_t> low_count_per_rank;

  assert( m_owned_targets.size() == edges_low.size() );
  mpi_all_gather(uint64_t(m_owned_targets.size()), low_count_per_rank, mpi_comm);

  //
  // Compute high degree vertices
  std::vector<uint64_t> high_count_per_rank;
  count_high_degree_transpose(mpi_comm, edges.begin(), edges.end(), global_hubs, high_count_per_rank);



  //
  // Compute Overflow schedule
  uint64_t global_edge_count = mpi_all_reduce(uint64_t(edges.size()*2),
        std::plus<uint64_t>(), mpi_comm);
  uint64_t target_edges_per_rank = global_edge_count / m_mpi_size;

  // std::string str_temp =
	 //  		"[" + std::to_string(m_mpi_rank)
	 //  		+ "] [H: " + std::to_string(high_count_per_rank[m_mpi_rank])
	 //  		+ "] [L: " + std::to_string(low_count_per_rank[m_mpi_rank])
	 //  		+ "] [T: " + std::to_string(high_count_per_rank[m_mpi_rank] +
	 //  																low_count_per_rank[m_mpi_rank])
	 //  		+ "] [E: " + std::to_string(high_count_per_rank[m_mpi_rank] +
	 //  																low_count_per_rank[m_mpi_rank] -
	 //  																target_edges_per_rank)
	 //  		+ "] [N: " + std::to_string(target_edges_per_rank -
	 //  																high_count_per_rank[m_mpi_rank] -
	 //  																low_count_per_rank[m_mpi_rank])
	 //  		+ "] \t Overflow values: ";

	// TODO: add code to calculate the number high edges I will recieve
	// The high_count_per_rank+ that will be the number of high_nodes
	// Which allows me to determine the number of edges I will have?

  std::map<int,uint64_t> overflow_schedule;
  uint64_t heavy_idx(0), light_idx(0);
  for (; heavy_idx < m_mpi_size && light_idx < m_mpi_size; ++heavy_idx) {
    uint64_t high_total_nodes = low_count_per_rank[heavy_idx] +
        high_count_per_rank[heavy_idx];

    while (high_total_nodes > target_edges_per_rank) {

      if (high_count_per_rank[heavy_idx] == 0) {
        break;
      }
      else if (low_count_per_rank[light_idx] + high_count_per_rank[light_idx]
                < target_edges_per_rank) {


        uint64_t num_high_nodes = high_count_per_rank[heavy_idx];
        uint64_t extra_nodes = high_total_nodes - target_edges_per_rank;
        // We can only send the extra high_nodes, so if we are over the
        // threshold but have no high nodes we can't send any.
        // Likewise if we have extra high nodes we don't want to send more than
        // the extra.
        uint64_t max_to_offload = std::min(num_high_nodes, extra_nodes);

        // At the same time we don't want to send more nodes then they need.
        uint64_t max_to_receive = target_edges_per_rank -
              high_count_per_rank[light_idx] - low_count_per_rank[light_idx];

        uint64_t to_move = std::min(max_to_offload, max_to_receive);

        high_count_per_rank[heavy_idx] -= to_move; // adjust the node count on
        high_count_per_rank[light_idx] += to_move; // each
        high_total_nodes -= to_move;  // keep this variable current

        if (heavy_idx == m_mpi_rank) { //if this is us, then we need to log what we are sending
          overflow_schedule[light_idx] += to_move;
        }
      } else { // light_idx does not need anymore nodes
        ++light_idx; // so check the next inline
        if (light_idx == m_mpi_size) {  // if there are no more
          break; // break the while loop, and the for loop's end condition is
                 // also met
        }
      } // else
    }  // while
  }  // for

 //  { // Debuging
	//   MPI_Barrier(mpi_comm); // TODO: remove with prin



	//   for (int i = 0; i < m_mpi_size; i++) {
	//     str_temp += std::to_string(overflow_schedule[i]) + ",";
	//   }

	//   std::cout  << str_temp << std::endl;
	// }

  //
  // Partition high degree, using overflow schedule
  partition_high_degree(mpi_comm, edges.begin(), edges.end(), global_hubs,
        edges_high, edges_high_overflow, overflow_schedule);


  std::cout  << str_temp << std::endl;

  MPI_Barrier(mpi_comm);


  double time_end = MPI_Wtime();
  if(m_mpi_rank == 0) {
    std::cout << "Partition time = " << time_end - time_start << std::endl;
  }

  free_edge_container<Container>(edges);


  //
  // Check correctness
  for(size_t i=0; i<edges_low.size(); ++i) {
    assert(edges_low[i].first % m_mpi_size == m_mpi_rank);
    assert(global_hubs.count(edges_low[i].first) == 0);
  }
  for(size_t i=0; i<edges_high.size(); ++i) {
    assert(edges_high[i].second % m_mpi_size == m_mpi_rank);
    assert(global_hubs.count(edges_high[i].first) > 0);
  }
  for(size_t i=0; i<edges_high_overflow.size(); ++i) {
    assert(global_hubs.count(edges_high_overflow[i].first) > 0);
  }

  uint64_t low_local_size      = edges_low.size();
  uint64_t high_local_size     = edges_high.size();
  uint64_t overflow_local_size = edges_high_overflow.size();
  uint64_t total_local_size    = low_local_size + high_local_size + overflow_local_size;

  //
  // Merge & sort edge lists
  edges_high.insert(edges_high.end(), edges_high_overflow.begin(), edges_high_overflow.end());
  {
    std::deque< std::pair<uint64_t,uint64_t> > empty(0);
    edges_high_overflow.swap(empty);
  }

  // This sort pred
  struct sort_pred {
    bool operator()(const std::pair<uint64_t,uint64_t> &left, const std::pair<uint64_t,uint64_t> &right) {
        if (left.first == right.first) {
          uint32_t o_left    = left.second % uint64_t(mpi_size);
          uint32_t o_right    = right.second % uint64_t(mpi_size);

          if (o_left == o_right) {
            uint64_t lid_left = left.second / uint64_t(mpi_size);
            uint64_t lid_right = right.second / uint64_t(mpi_size);

            return lid_left < lid_right;
          }
          return o_left < o_right;
        }
        return left.first < right.first;
    }
    uint64_t mpi_size;
  };

  sort_pred pred;
  pred.mpi_size = m_mpi_size;
  std::sort(edges_low.begin(), edges_low.end(), pred);
  std::sort(edges_high.begin(), edges_high.end(), pred);

  //We sort before to put it in the same order as the edges_low list
  //So we sort before we update owned_targets to delegegate versions
  for (uint64_t i=0; i<= m_max_vertex; ++i) {
    auto itr_start = m_owned_targets.begin() + m_owned_info[i  ].low_csr_idx;
    auto itr_end   = m_owned_targets.begin() + m_owned_info[i+1].low_csr_idx;

    assert(m_owned_info[i].low_csr_idx <= m_owned_info[i+1].low_csr_idx);
    assert( m_owned_info[i+1].low_csr_idx <= m_owned_targets.size());
    std::sort(itr_start, itr_end);
  }

  //
  // Compute Hub information
  std::vector<uint64_t> vec_sorted_hubs(global_hubs.begin(), global_hubs.end());
  m_delegate_degree.resize(vec_sorted_hubs.size(), 0);
  m_delegate_label.resize(vec_sorted_hubs.size());
  std::sort(vec_sorted_hubs.begin(), vec_sorted_hubs.end());

  for(size_t i=0; i<vec_sorted_hubs.size(); ++i) {
    uint64_t t_local_id = vec_sorted_hubs[i] / uint64_t(m_mpi_size);
    uint32_t t_owner = uint32_t(vec_sorted_hubs[i] % uint32_t(m_mpi_size));

    vertex_locator old_ver_loc(false, t_local_id, t_owner);
    vertex_locator new_ver_loc(true, i, t_owner);

    m_map_delegate_locator[vec_sorted_hubs[i]] = new_ver_loc;
    m_delegate_label[i] = vec_sorted_hubs[i];

    //Update the low CSR to indicate a global hub.
    auto itr_start = m_owned_targets.begin();
    auto itr_end   = m_owned_targets.end();
    itr_start = std::find(itr_start, itr_end, old_ver_loc);
    while (itr_start != itr_end) {
      *itr_start = new_ver_loc;
      itr_start = std::find(itr_start+1, itr_end, old_ver_loc);
    }

    //
    // Tag owned delegates
    //
    if (t_owner == m_mpi_rank) {
      m_owned_info[t_local_id].is_delegate = 1;
      m_owned_info[t_local_id].delegate_id = i;
    }
  }  // for items in vec_sorted_hubs


  MPI_Barrier(mpi_comm); // TODO(steven): delete this



    // Compares the generated Low CSR  with the one that would have been made
    // using the edge_low vector.
    // If they are equal then the Low CSR code should be correct
  {

    if(m_owned_targets.size() != edges_low.size() ) {
      std::cout << m_mpi_rank << ": m_owned_targets capacity is " <<
          m_owned_targets.size() << " but should be " << edges_low.size() <<
          "." << std::endl;
    }

    uint64_t max_local_id = local_source_id(m_mpi_size)(edges_low.back());
    uint64_t max_owned_delegate_id = edges_high.size() > 0 ?
        local_source_id(m_mpi_size)(edges_high.back()) : 0;

    // max_local_id = std::max(max_local_id, max_owned_delegate_id) +2 ;

    // if(m_owned_info.capacity() != max_local_id ) {
    //   std::cout << m_mpi_rank << ": m_owned_info capacity is " <<
    //       m_owned_info.capacity() << " but should be " << max_local_id <<
    //       "." << std::endl;
    // }

    uint64_t cur_source_id = 0;
    for (uint64_t i=0; i<edges_low.size(); ++i) {
      uint64_t new_source_id = local_source_id(m_mpi_size)(edges_low[i]);
      uint64_t new_target_label = edges_low[i].second;

      while (new_source_id > cur_source_id) {
        ++cur_source_id;
        if (m_owned_info[cur_source_id].is_delegate != 1) {
          // This check is necessary as some have been upgraded to delegate
          // already, which occurs after this point.
          assert(m_owned_info[cur_source_id] == vert_info(false,0,i));
        }
      }

      assert(m_owned_info[new_source_id].low_csr_idx <= i);
      assert(m_owned_info[new_source_id+1].low_csr_idx > i);
      assert(label_to_locator(new_target_label) == m_owned_targets[i]);
    }
  }  // end check low-degree CSR

  //
  // Build high-degree CSR
  {
    m_delegate_targets.resize(edges_high.size());
    m_delegate_info.resize(m_map_delegate_locator.size()+1, m_delegate_targets.size());
    uint64_t cur_source_id = 0;
    m_delegate_info[cur_source_id] = 0;
    for (uint64_t i=0; i<edges_high.size(); ++i) {
      uint64_t new_source_id = m_map_delegate_locator[edges_high[i].first].local_id();
      m_delegate_degree[new_source_id]++;

      uint64_t new_target_label = edges_high[i].second;

      while (new_source_id > cur_source_id) {
        ++cur_source_id;
        m_delegate_info[cur_source_id] = i;
      }
      m_delegate_targets[i] = label_to_locator(new_target_label);
    }
  }
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
  for(size_t i=0; i<m_delegate_degree.size(); ++i) {
    if(i%m_mpi_size == m_mpi_rank) {
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
  uint64_t overflow_max_size  = mpi_all_reduce(overflow_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);
  uint64_t total_max_size     = mpi_all_reduce(total_local_size, std::greater<uint64_t>(), MPI_COMM_WORLD);

  uint64_t low_sum_size       = mpi_all_reduce(low_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t high_sum_size      = mpi_all_reduce(high_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t overflow_sum_size  = mpi_all_reduce(overflow_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t total_sum_size     = mpi_all_reduce(total_local_size, std::plus<uint64_t>(), MPI_COMM_WORLD);

  uint64_t local_count_del_target = 0;
  for(uint64_t i=0; i<m_owned_targets.size(); ++i) {
    if(m_owned_targets[i].is_delegate()) ++local_count_del_target;
  }
  uint64_t total_count_del_target = mpi_all_reduce(local_count_del_target, std::plus<uint64_t>(), MPI_COMM_WORLD);

  if(m_mpi_rank == 0) {
    std::cout << "Count of hub vertices = " << global_hubs.size() << std::endl;
    std::cout << "Total percentage good hub edges = " << double(high_sum_size) / double(total_sum_size) * 100.0 << std::endl;
    std::cout << "total count del target = " << total_count_del_target << std::endl;
    std::cout << "Total percentage of localized edges = " << double(high_sum_size + total_count_del_target) / double(total_sum_size) * 100.0 << std::endl;
    std::cout << "Global number of edges = " << total_sum_size << std::endl;
    std::cout << "Number of small degree = " << low_sum_size << std::endl;
    std::cout << "Number of hubs = " << high_sum_size << std::endl;
    std::cout << "Number of overfow = " << overflow_sum_size << std::endl;
    std::cout << "oned imbalance = " << double(low_max_size) / double(low_sum_size/m_mpi_size) << std::endl;
    std::cout << "hubs imbalance = " << double(high_max_size) / double(high_sum_size/m_mpi_size) << std::endl;
    if(overflow_sum_size > 0) {
      std::cout << "overflow imbalance = " << double(overflow_max_size) / double(overflow_sum_size/m_mpi_size) << std::endl;
    }
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
