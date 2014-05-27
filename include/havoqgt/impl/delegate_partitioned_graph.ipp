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

template <typename InputIterator>
void 
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

template <typename InputIterator>
void 
count_low_degree(MPI_Comm mpi_comm, 
                 InputIterator unsorted_itr, 
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 uint64_t delegate_degree_threshold) 
{
  double time_start = MPI_Wtime();
  using boost::unordered_map;
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );
  if(mpi_rank == 0)
  std::cout << "Starting:  count_low_degree" << std::endl;
  unordered_map<uint64_t,uint64_t> local_degree_count;

  //
  // loop until all ranks have finished
  uint64_t gi=0;
  while(!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end, mpi_comm)) {
    // 1D exchange
    std::vector<uint64_t> to_recv_1d;
    std::vector<uint64_t> to_exchange_1d(0);
    ++gi;
    for(uint64_t i=0; gi!=8 & gi!=16 && gi!=32 && gi!=64 && gi!=128 && gi != 256 && gi != 512 
                      && gi != 1024 && gi != 4096 && i<4*4096 && unsorted_itr != unsorted_itr_end; 
                      ++i, ++unsorted_itr, ++gi) {
      if(global_hub_set.count(unsorted_itr->first) == 0) {
        to_exchange_1d.push_back(unsorted_itr->first);
      }
      //transpose version
      if(global_hub_set.count(unsorted_itr->second) == 0) {
        to_exchange_1d.push_back(unsorted_itr->second);
      }
    }
    mpi_all_to_all_better(to_exchange_1d, to_recv_1d, source_partitioner(mpi_size), mpi_comm);
    std::vector<uint64_t>().swap(to_exchange_1d); //release memory

    //
    // Process new partitioned sources
    std::vector<uint64_t> new_hubs;
    for(size_t i=0; i<to_recv_1d.size(); ++i) {
      uint64_t source = to_recv_1d[i];
      assert(source % mpi_size == mpi_rank);
      if(global_hub_set.count(source) == 0) {
        ++local_degree_count[source];
        if(local_degree_count[source] == delegate_degree_threshold) {
          new_hubs.push_back(source);
        }
      }
    }

    //
    // Exchange newly found hubs
    std::vector<uint64_t> new_global_hubs;
    mpi_all_gather(new_hubs, new_global_hubs, mpi_comm);
    global_hub_set.insert(new_global_hubs.begin(), new_global_hubs.end());
  }
  double time_end = MPI_Wtime();
  if(mpi_rank == 0)
  std::cout << "count_low_degree time = " << time_end - time_start << std::endl;
}

template <typename InputIterator>
void partition_low_degree(MPI_Comm mpi_comm, 
                 InputIterator unsorted_itr, 
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 std::deque<std::pair<uint64_t, uint64_t> >& edges_low) 
{
  double time_start = MPI_Wtime();
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) );
  if(mpi_rank == 0)
    std::cout << "Starting:  partition_low_degree" << std::endl;
  while(!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end, mpi_comm)) {
    std::vector<std::pair<uint64_t, uint64_t> > to_recv_edges_low;
    {
    std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_low;
    to_send_edges_low.reserve(16*1024);
    for(size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024; ++unsorted_itr) {
      if(global_hub_set.count(unsorted_itr->first) == 0) {
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
    mpi_all_to_all_better(to_send_edges_low, to_recv_edges_low, edge_source_partitioner(mpi_size), mpi_comm);
    }
    for(size_t i=0; i<to_recv_edges_low.size(); ++i) {
      assert(to_recv_edges_low[i].first % mpi_size == mpi_rank);
    }
    edges_low.insert(edges_low.end(), to_recv_edges_low.begin(), to_recv_edges_low .end());
  }
  double time_end = MPI_Wtime();
  if(mpi_rank == 0)
  std::cout << "partition_low_degree time = " << time_end - time_start << std::endl;
}

template <typename InputIterator>
void partition_high_degree(MPI_Comm mpi_comm, 
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

  while(!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end, mpi_comm)) {
    
    std::vector<std::pair<uint64_t, uint64_t> > to_recv_edges_high;
    std::vector<std::pair<int, std::pair<uint64_t, uint64_t> > > to_recv_overflow;
    {
      std::vector<std::pair<uint64_t, uint64_t> > to_send_edges_high;
      to_send_edges_high.reserve(16*1024);
      for(size_t i=0; unsorted_itr != unsorted_itr_end && i<16*1024; ++unsorted_itr) {
        if(global_hub_set.count(unsorted_itr->first) == 1) {
          ++i;
          to_send_edges_high.push_back(*unsorted_itr);
        }
        //transpose version
        if(global_hub_set.count(unsorted_itr->second) == 1) {
          ++i;
          to_send_edges_high.push_back(std::make_pair(unsorted_itr->second,
                                                   unsorted_itr->first));
        }
      }
      mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high, edge_target_partitioner(mpi_size), mpi_comm);
    }
    for(size_t i=0; i<edges_high.size(); ++i) {
      assert(edges_high[i].second % mpi_size == mpi_rank);
    }
    {
      // Copy high edges to either edges_high or overflow
      std::vector<std::pair<int, std::pair<uint64_t, uint64_t> > > to_send_overflow;
      for(size_t i=0; i<to_recv_edges_high.size(); ++i) {
        if(overflow_schedule.empty()) {
          edges_high.push_back(to_recv_edges_high[i]);
        } else {
          int dest = overflow_schedule.begin()->first;
          overflow_schedule[dest]--;
          if(overflow_schedule[dest] == 0) {
            overflow_schedule.erase(dest);
          }
          to_send_overflow.push_back(std::make_pair(dest, to_recv_edges_high[i]));
        }
      }
      uint64_t global_overflow_size = mpi_all_reduce(uint64_t(to_send_overflow.size()), std::plus<uint64_t>(), mpi_comm);
      if(global_overflow_size > 0) {
        mpi_all_to_all_better(to_send_overflow, to_recv_overflow, dest_pair_partitioner(), mpi_comm);
      }
    }
    for(size_t i=0; i<to_recv_overflow.size(); ++i) {
      edges_high_overflow.push_back(to_recv_overflow[i].second);
    }
  }
  double time_end = MPI_Wtime();
  if(mpi_rank == 0)
  std::cout << "partition_high_degree time = " << time_end - time_start << std::endl;
}



/**
 * Builds a delegate_partitioned_graph with from and unsorted sequence of edges.
 *
 * @param sm       Pointer to segment manager
 * @param mpi_comm MPI communicator
 * @param Container input edges to partition
 * @param delegate_degree_threshold Threshold used to assign delegates
 */
template <typename Arena>
template <typename Container>
delegate_partitioned_graph<Arena>::
delegate_partitioned_graph(Arena& arena, 
                           MPI_Comm mpi_comm, 
                           Container& edges, 
                           uint64_t delegate_degree_threshold) 
    : m_mpi_comm(mpi_comm),
      m_owned_info(arena.template get_allocator<vert_info>()),
      m_owned_targets(arena.template get_allocator<vertex_locator>()),
      m_delegate_info(arena.template get_allocator<uint64_t>()),
      m_delegate_degree(arena.template get_allocator<uint64_t>()),
      m_delegate_label(arena.template get_allocator<uint64_t>()),
      m_map_delegate_locator(100, boost::hash<uint64_t>(), std::equal_to<uint64_t>(), arena.template get_allocator<std::pair<uint64_t,vertex_locator> >()),
      m_delegate_targets(arena.template get_allocator<vertex_locator>()),
      m_delegate_degree_threshold(delegate_degree_threshold) {
  MPI_Barrier(m_mpi_comm);
  double time_start = MPI_Wtime();
  CHK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &m_mpi_size) );
  CHK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &m_mpi_rank) );
  
  std::deque< std::pair<uint64_t, uint64_t> > edges_low, edges_high, edges_high_overflow;
  boost::unordered_set<uint64_t> global_hubs;

  assert(sizeof(vertex_locator) == 8);
  
   

  //
  // Count low degree & find delegates
  count_low_degree(mpi_comm, edges.begin(), edges.end(), global_hubs, delegate_degree_threshold);

  //
  // Partition low degree vertices
  partition_low_degree(mpi_comm, edges.begin(), edges.end(), global_hubs, edges_low);
  std::vector<uint64_t> low_count_per_rank;
  mpi_all_gather(uint64_t(edges_low.size()), low_count_per_rank, m_mpi_comm);

  //
  // Compute high degree vertices
  std::vector<uint64_t> high_count_per_rank;
  count_high_degree_transpose(mpi_comm, edges.begin(), edges.end(), global_hubs, high_count_per_rank);
  //
  // Compute Overflow schedule
  uint64_t global_edge_count = mpi_all_reduce(uint64_t(edges.size()*2), std::plus<uint64_t>(), m_mpi_comm);
  uint64_t target_edges_per_rank = global_edge_count / m_mpi_size;
  std::map<int,uint64_t> overflow_schedule;
  uint64_t heavy_idx(0), light_idx(0);
  for(; heavy_idx < m_mpi_size && light_idx < m_mpi_size; ++heavy_idx) {
    while(low_count_per_rank[heavy_idx] + high_count_per_rank[heavy_idx] > target_edges_per_rank) {
      if(low_count_per_rank[light_idx] + high_count_per_rank[light_idx] < target_edges_per_rank) {
        if(high_count_per_rank[heavy_idx] == 0) break; //can't move more
        uint64_t max_to_offload = std::min(high_count_per_rank[heavy_idx], high_count_per_rank[heavy_idx]+low_count_per_rank[heavy_idx]-target_edges_per_rank);
        uint64_t max_to_receive = target_edges_per_rank - high_count_per_rank[light_idx] - low_count_per_rank[light_idx];
        uint64_t to_move = std::min(max_to_offload, max_to_receive);
        high_count_per_rank[heavy_idx]-=to_move;
        high_count_per_rank[light_idx]+=to_move;
        if(heavy_idx == m_mpi_rank) {
          overflow_schedule[light_idx]+=to_move;
        }
      } else {
        ++light_idx;
        if(light_idx == m_mpi_size) break;
      }
    }
  }

  //
  // Partition high degree, using overflow schedule
  partition_high_degree(mpi_comm, edges.begin(), edges.end(), global_hubs, edges_high, edges_high_overflow, overflow_schedule);

  MPI_Barrier(m_mpi_comm);
  double time_end = MPI_Wtime();
  if(m_mpi_rank == 0) {
    std::cout << "Partition time = " << time_end - time_start << std::endl;
  }
  {
    Container empty(0);
    edges.swap(empty);
  }
  

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
  std::sort(edges_low.begin(), edges_low.end());
  std::sort(edges_high.begin(), edges_high.end());

  //
  // Compute Hub information
  std::vector<uint64_t> vec_sorted_hubs(global_hubs.begin(), global_hubs.end());
  m_delegate_degree.resize(vec_sorted_hubs.size(), 0);
  m_delegate_label.resize(vec_sorted_hubs.size());
  std::sort(vec_sorted_hubs.begin(), vec_sorted_hubs.end());
  for(size_t i=0; i<vec_sorted_hubs.size(); ++i) {
    m_map_delegate_locator[vec_sorted_hubs[i]] = vertex_locator(true, i, 
                                            uint32_t(vec_sorted_hubs[i] % uint32_t(m_mpi_size)));
    m_delegate_label[i] = vec_sorted_hubs[i];
  }

  //
  // Build low-degree CSR
  {
    m_owned_targets.resize(edges_low.size());
    uint64_t max_local_id = edges_low.back().first / m_mpi_size;
    uint64_t max_owned_delegate_id = edges_high.size() > 0 ? edges_high.back().first / m_mpi_size : 0;
    max_local_id = std::max(max_local_id, max_owned_delegate_id); // not sure if this is correct??
    m_owned_info.resize(max_local_id+2, vert_info(false,0,edges_low.size()));
    uint64_t cur_source_id = 0;
    m_owned_info[cur_source_id] = vert_info(false,0,0);
    for(uint64_t i=0; i<edges_low.size(); ++i) {
      uint64_t new_source_id = edges_low[i].first / m_mpi_size;
      uint64_t new_target_label = edges_low[i].second;
      while(new_source_id > cur_source_id) {
        ++cur_source_id;
        m_owned_info[cur_source_id] = vert_info(false,0,i);
      }
      m_owned_targets[i] = label_to_locator(new_target_label);
    }
  }

  //
  // Build high-degree CSR
  {
    m_delegate_targets.resize(edges_high.size());
    m_delegate_info.resize(m_map_delegate_locator.size()+1, m_delegate_targets.size());
    uint64_t cur_source_id = 0;
    m_delegate_info[cur_source_id] = 0;
    for(uint64_t i=0; i<edges_high.size(); ++i) {
      uint64_t new_source_id = m_map_delegate_locator[edges_high[i].first].local_id();
      m_delegate_degree[new_source_id]++;
      uint64_t new_target_label = edges_high[i].second;
      while(new_source_id > cur_source_id) {
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
      mpi_all_reduce(my_hub_degrees, tmp_hub_degrees, std::plus<uint64_t>(), m_mpi_comm);
      m_delegate_degree.clear();
      m_delegate_degree.insert(m_delegate_degree.end(),tmp_hub_degrees.begin(), tmp_hub_degrees.end());
    }
  }
  assert(m_delegate_degree.size() == m_delegate_label.size());

  //
  // Tag owned delegates
  for(typename boost::unordered_map<uint64_t, vertex_locator, boost::hash<uint64_t>, 
          std::equal_to<uint64_t>, typename Arena::template allocator<std::pair<uint64_t,vertex_locator> >::type >::iterator itr = m_map_delegate_locator.begin();
      itr != m_map_delegate_locator.end(); ++itr) {
    uint64_t label = itr->first;
    vertex_locator locator = itr->second;
    //if(locator.owner() == m_mpi_rank) {
    if(label % uint64_t(m_mpi_size) == m_mpi_rank) {
      uint64_t local_id = label / uint64_t(m_mpi_size);
      m_owned_info[local_id].is_delegate = 1;
      m_owned_info[local_id].delegate_id = locator.local_id();
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
template <typename Arena>
inline
uint64_t 
delegate_partitioned_graph<Arena>::
locator_to_label(delegate_partitioned_graph<Arena>::vertex_locator locator) const {
  if(locator.is_delegate()) {
    return m_delegate_label[locator.local_id()];
  } 
  return uint64_t(locator.local_id()) * uint64_t(m_mpi_size) + uint64_t(locator.owner()); 
}

/**
 * @param  label vertex label to convert
 * @return locator for the label
 */
template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_locator 
delegate_partitioned_graph<Arena>::
label_to_locator(uint64_t label) const {
  typename boost::unordered_map<uint64_t, vertex_locator, boost::hash<uint64_t>, 
          std::equal_to<uint64_t>, typename Arena::template allocator<std::pair<uint64_t,vertex_locator> >::type >::const_iterator itr 
      = m_map_delegate_locator.find(label);
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
template <typename Arena>
inline void 
delegate_partitioned_graph<Arena>::
sync_global_hub_set(const boost::unordered_set<uint64_t>& local_hubs, 
                         boost::unordered_set<uint64_t>& global_hubs,
                         bool local_change) {
  uint32_t new_hubs = mpi_all_reduce(uint32_t(local_change), 
                                     std::plus<uint32_t>(), m_mpi_comm);

  if(new_hubs > 0) {
    std::vector<uint64_t> vec_local_hubs(local_hubs.begin(), local_hubs.end());
    std::vector<uint64_t> vec_global_hubs;
    // global gather local hubs
    mpi_all_gather(vec_local_hubs, vec_global_hubs, m_mpi_comm);
    // Insert gathered global hubs to set
    global_hubs.insert(vec_global_hubs.begin(), vec_global_hubs.end());
  }
}

/**
 * @param  locator Vertex locator
 * @return Begin Edge Iterator
 */
template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::edge_iterator 
delegate_partitioned_graph<Arena>::
edges_begin(delegate_partitioned_graph<Arena>::vertex_locator locator) const {
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
template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::edge_iterator 
delegate_partitioned_graph<Arena>::
edges_end(delegate_partitioned_graph<Arena>::vertex_locator locator) const {
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
template <typename Arena>
inline
uint64_t 
delegate_partitioned_graph<Arena>::
degree(delegate_partitioned_graph<Arena>::vertex_locator locator) const {
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
template <typename Arena>
inline
uint64_t 
delegate_partitioned_graph<Arena>::
local_degree(delegate_partitioned_graph<Arena>::vertex_locator locator) const {
  uint64_t local_id = locator.local_id();
  if(locator.is_delegate()) {
    assert(local_id + 1 < m_delegate_info.size());
    return m_delegate_info[local_id + 1] - m_delegate_info[local_id];
  }
  assert(local_id + 1 < m_owned_info.size());
  return m_owned_info[local_id+1].low_csr_idx - 
         m_owned_info[local_id].low_csr_idx;
}


template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_iterator 
delegate_partitioned_graph<Arena>::
vertices_begin() const { 
  return vertex_iterator(0,this); 
}

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_iterator 
delegate_partitioned_graph<Arena>::
vertices_end() const { 
  return vertex_iterator(m_owned_info.size()-1,this); 
}

template <typename Arena>
inline
bool 
delegate_partitioned_graph<Arena>::
is_label_delegate(uint64_t label) const { 
  return m_map_delegate_locator.count(label) > 0; 
}

template <typename Arena>
template <typename T, typename Arena2>
typename delegate_partitioned_graph<Arena>::template vertex_data<T, typename Arena2::segment_manager>* 
delegate_partitioned_graph<Arena>::
create_vertex_data(Arena2& arena) const {

  return arena.template construct<vertex_data<T, typename Arena2::segment_manager> >(bip::anonymous_instance)(m_owned_info.size(), m_delegate_info.size(), arena.get_segment_manager());
}

/**
 * @param   init initial value for each vertex
 */
template <typename Arena>
template <typename T, typename Arena2>
delegate_partitioned_graph<Arena>::vertex_data<T, typename Arena2::segment_manager>* 
delegate_partitioned_graph<Arena>::
create_vertex_data(const T& init, Arena2& arena) const {
  return arena.template construct<vertex_data<T, typename Arena2::segment_manager> >(bip::anonymous_instance)(m_owned_info.size(), m_delegate_info.size(), init,  arena.get_segment_manager());
}

template <typename Arena>
template <typename T, typename Arena2>
typename delegate_partitioned_graph<Arena>::template edge_data<T, typename Arena2::segment_manager>* 
delegate_partitioned_graph<Arena>::
create_edge_data(Arena2& arena) const {
  return arena.template construct<edge_data<T, typename Arena2::segment_manager> >(bip::anonymous_instance)(m_owned_targets.size(), m_delegate_targets.size(), arena.get_segment_manager());
}

/**
 * @param   init initial value for each vertex
 */
template <typename Arena>
template <typename T, typename Arena2>
delegate_partitioned_graph<Arena>::edge_data<T, typename Arena2::segment_manager>* 
delegate_partitioned_graph<Arena>::
create_edge_data(const T& init, Arena2& arena) const {
  return arena.template construct<edge_data<T, typename Arena2::segment_manager> >(bip::anonymous_instance)(m_owned_targets.size(), m_delegate_targets.size(), init,  arena.get_segment_manager());
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
template <typename Arena>
inline
delegate_partitioned_graph<Arena>::vertex_locator::
vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest) {
  m_is_bcast     = 0;
  m_is_intercept = 0;
  if(is_delegate) {
    m_is_delegate = true;
    m_local_id    = local_id;
    m_owner_dest  = owner_dest;
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

template <typename Arena>
inline bool
delegate_partitioned_graph<Arena>::vertex_locator::
is_equal(const typename delegate_partitioned_graph<Arena>::vertex_locator x) const {
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
template <typename Arena>
inline
delegate_partitioned_graph<Arena>::edge_iterator::
edge_iterator(vertex_locator source, 
              uint64_t edge_offset,
              const delegate_partitioned_graph* const pgraph)
  : m_source(source)
  , m_edge_offset(edge_offset)
  , m_ptr_graph(pgraph) { }

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::edge_iterator&
delegate_partitioned_graph<Arena>::edge_iterator::operator++() {
  ++m_edge_offset;
  return *this;
}

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::edge_iterator
delegate_partitioned_graph<Arena>::edge_iterator::operator++(int) {
  edge_iterator to_return = *this;
  ++m_edge_offset;
  return to_return;
}

template <typename Arena>
inline bool
delegate_partitioned_graph<Arena>::edge_iterator::
is_equal(const typename delegate_partitioned_graph<Arena>::edge_iterator& x) const {
    assert(m_source      == x.m_source);
    assert(m_ptr_graph   == x.m_ptr_graph);
    return m_edge_offset == x.m_edge_offset;
}

template <typename Arena>
inline bool
operator==(const typename delegate_partitioned_graph<Arena>::edge_iterator& x,
           const typename delegate_partitioned_graph<Arena>::edge_iterator& y) {
  return x.is_equal(y);

}

template <typename Arena>
inline bool
operator!=(const typename delegate_partitioned_graph<Arena>::edge_iterator& x,
           const typename delegate_partitioned_graph<Arena>::edge_iterator& y) {
  return !(x.is_equal(y));
}

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_locator
delegate_partitioned_graph<Arena>::edge_iterator::target() const {
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

template <typename Arena>
inline
delegate_partitioned_graph<Arena>::vertex_iterator::
vertex_iterator(uint64_t index, const delegate_partitioned_graph<Arena>*  pgraph)
  : m_ptr_graph(pgraph)
  , m_owned_vert_index(index) {
  update_locator();
}

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_iterator& 
delegate_partitioned_graph<Arena>::vertex_iterator::operator++() {
  ++m_owned_vert_index;
  update_locator();
  return *this;
}

template <typename Arena>
inline
typename delegate_partitioned_graph<Arena>::vertex_iterator 
delegate_partitioned_graph<Arena>::vertex_iterator::operator++(int) {
  vertex_iterator to_return = *this;
  ++m_owned_vert_index;
  update_locator();
  return to_return;
}

template <typename Arena>
inline bool
delegate_partitioned_graph<Arena>::vertex_iterator::
is_equal(const typename delegate_partitioned_graph<Arena>::vertex_iterator& x) const {
  assert(m_ptr_graph        == x.m_ptr_graph);
  return m_owned_vert_index == x.m_owned_vert_index;
}

template <typename Arena>
inline bool
operator==(const typename delegate_partitioned_graph<Arena>::vertex_iterator& x,
           const typename delegate_partitioned_graph<Arena>::vertex_iterator& y) {
  return x.is_equal(y);

}

template <typename Arena>
inline bool
operator!=(const typename delegate_partitioned_graph<Arena>::vertex_iterator& x,
           const typename delegate_partitioned_graph<Arena>::vertex_iterator& y) {
  return !(x.is_equal(y));
}

template <typename Arena>
inline void
delegate_partitioned_graph<Arena>::vertex_iterator::
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

template <typename Arena>
inline
delegate_partitioned_graph<Arena>::vert_info::
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
template <typename Arena>
template<typename T, typename SegmentManager>
delegate_partitioned_graph<Arena>::vertex_data<T,SegmentManager>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, SegmentManager* sm)
  : m_owned_vert_data(sm->template get_allocator<T>())
  , m_delegate_data(sm->template get_allocator<T>()) {
  m_owned_vert_data.resize(owned_data_size);
  m_delegate_data.resize(delegate_size);
  }

template <typename Arena>
template<typename T, typename SegmentManager>
delegate_partitioned_graph<Arena>::vertex_data<T, SegmentManager>::
vertex_data(uint64_t owned_data_size, uint64_t delegate_size, const T& init, SegmentManager* sm)
  : m_owned_vert_data(owned_data_size, init, sm->template get_allocator<T>())
  , m_delegate_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename Arena>
template<typename T, typename SegmentManager>
T& 
delegate_partitioned_graph<Arena>::vertex_data<T, SegmentManager>::
operator[](const vertex_locator& locator) {
  if(locator.is_delegate()) {
    assert(locator.local_id() < m_delegate_data.size());
    return m_delegate_data[locator.local_id()];
  }
  assert(locator.local_id() < m_owned_vert_data.size());
  return m_owned_vert_data[locator.local_id()];
}

template <typename Arena>
template<typename T, typename SegmentManager>
const T& 
delegate_partitioned_graph<Arena>::vertex_data<T, SegmentManager>::operator[](const vertex_locator& locator) const {
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
template <typename Arena>
template<typename T, typename SegmentManager>
delegate_partitioned_graph<Arena>::edge_data<T,SegmentManager>::
edge_data(uint64_t owned_size, uint64_t delegate_size, SegmentManager* sm)
  : m_owned_edge_data(sm->template get_allocator<T>())
  , m_delegate_edge_data(sm->template get_allocator<T>()) {
  m_owned_edge_data.resize(owned_size);
  m_delegate_edge_data.resize(delegate_size);
  }

template <typename Arena>
template<typename T, typename SegmentManager>
delegate_partitioned_graph<Arena>::edge_data<T, SegmentManager>::
edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init, SegmentManager* sm)
  : m_owned_edge_data(owned_size, init, sm->template get_allocator<T>())
  , m_delegate_edge_data(delegate_size, init, sm->template get_allocator<T>()) { }

template <typename Arena>
template<typename T, typename SegmentManager>
T& 
delegate_partitioned_graph<Arena>::edge_data<T, SegmentManager>::
operator[](const edge_iterator& itr) {
  if(itr.m_source.is_delegate()) {
    assert(itr.m_edge_offset < m_delegate_edge_data.size());
    return m_delegate_edge_data[itr.m_edge_offset];
  }
  assert(itr.m_edge_offset < m_owned_edge_data.size());
  return m_owned_edge_data[itr.m_edge_offset];
}

template <typename Arena>
template<typename T, typename SegmentManager>
const T& 
delegate_partitioned_graph<Arena>::edge_data<T, SegmentManager>::operator[](const edge_iterator& itr) const {
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
