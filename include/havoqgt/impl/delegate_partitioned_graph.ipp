
/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Re-written by Steven Feldman <feldman12@llnl.gov>.
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
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */

#include <type_traits>

#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED

namespace havoqgt {
namespace mpi {
/**
 * Builds a delegate_partitioned_graph with from and unsorted sequence of edges.
 *
 * @param sm       Pointer to segment manager
 * @paramm_mpi_comm MPI communicator
 * @param Container input edges to partition
 * @param delegate_degree_threshold Threshold used to assign delegates
*/

template <typename SegmentManager>
template <typename Container>
delegate_partitioned_graph<SegmentManager>::
delegate_partitioned_graph(const SegmentAllocator<void>& seg_allocator,
                           MPI_Comm mpi_comm,
                           Container& edges,
                           uint64_t max_vertex,
                           uint64_t delegate_degree_threshold,
                           uint64_t _node_partitions,
                           uint64_t _chunk_size,
                           ConstructionState stop_after
                           )
    : m_mpi_comm(mpi_comm),
      m_global_edge_count(edges.size()),
      m_local_outgoing_count(seg_allocator),
      m_local_incoming_count(seg_allocator),
      m_owned_info(seg_allocator),
      m_owned_info_tracker(seg_allocator),
      // m_owned_targets(seg_allocator),
      m_delegate_degree_threshold(delegate_degree_threshold),
      m_delegate_info(seg_allocator),
      m_delegate_degree(seg_allocator),
      m_delegate_label(seg_allocator),
      // m_delegate_targets(seg_allocator),
      //m_map_delegate_locator(100, boost::hash<uint64_t>(),
      //    std::equal_to<uint64_t>(), seg_allocator),
      m_map_delegate_locator(seg_allocator),
      m_controller_locators(seg_allocator),
      m_has_edge_data(false) 
{
  
  // this is just a dummy edge_data object 
  delegate_partitioned_graph<SegmentManager>::
    edge_data<uint8_t, SegmentAllocator<uint8_t>> _edge_data = 
      delegate_partitioned_graph<SegmentManager>::
      edge_data<uint8_t, SegmentAllocator<uint8_t>>(seg_allocator);
 
  CHK_MPI( MPI_Comm_size(m_mpi_comm, &m_mpi_size) );
  CHK_MPI( MPI_Comm_rank(m_mpi_comm, &m_mpi_rank) );

  processes_per_node = havoqgt_env()->node_local_comm().size();
  node_partitions = std::min(_node_partitions,uint64_t(processes_per_node));
  edge_chunk_size = _chunk_size;

  m_global_max_vertex = max_vertex;
  m_max_vertex = (std::ceil(double(max_vertex) / double(m_mpi_size)));

  LogStep logstep_main("Delegate Partitioning", m_mpi_comm, m_mpi_rank); 
 
////////////////////////////////////////////////////////////////////////////////
///// Meta data phase of graph construction
////////////////////////////////////////////////////////////////////////////////

  {
    LogStep logstep("Allocating 4 Arrays of length max local vertex.",
      m_mpi_comm, m_mpi_rank);
    m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0));
    m_owned_info_tracker.resize(m_max_vertex+2, 0);
    m_local_outgoing_count.resize(m_max_vertex+1, 0);
    m_local_incoming_count.resize(m_max_vertex+1, 0);
    MPI_Barrier(m_mpi_comm);
  }

  boost::unordered_set<uint64_t> global_hubs;

  {
    LogStep logstep("count_edge_degree", m_mpi_comm, m_mpi_rank);
    count_edge_degrees(edges.begin(), edges.end(), global_hubs,
      delegate_degree_threshold);
    if (m_mpi_rank == 0)
      std::cout << "\tNumber of Delegates: " << global_hubs.size() << std::endl;
    MPI_Barrier(m_mpi_comm);
  }

  {
      LogStep logstep("initialize_low_meta_data", m_mpi_comm, m_mpi_rank);
      initialize_low_meta_data(global_hubs);
        MPI_Barrier(m_mpi_comm);
  }

  {
      LogStep logstep("initialize_high_meta_data", m_mpi_comm, m_mpi_rank);
      initialize_high_meta_data(global_hubs);
        MPI_Barrier(m_mpi_comm);
  }

  {
    LogStep logstep("count_high_degree_edges", m_mpi_comm, m_mpi_rank);
    count_high_degree_edges(edges.begin(), edges.end(), global_hubs);
    MPI_Barrier(m_mpi_comm);
  }

  m_graph_state = MetaDataGenerated;
  if (m_graph_state != stop_after) {
    complete_construction(seg_allocator, mpi_comm, edges, _edge_data);
  }
   
}

template <typename SegmentManager>
template <typename Container, typename edge_data_type>
delegate_partitioned_graph<SegmentManager>::
delegate_partitioned_graph(const SegmentAllocator<void>& seg_allocator,
                           MPI_Comm mpi_comm,
                           Container& edges, 
                           uint64_t max_vertex,
                           uint64_t delegate_degree_threshold,
                           uint64_t _node_partitions,
                           uint64_t _chunk_size,  
                           edge_data_type& _edge_data,
                           bool _has_edge_data,
                           ConstructionState stop_after
                           )
    : m_mpi_comm(mpi_comm),
      m_global_edge_count(edges.size()),
      m_local_outgoing_count(seg_allocator),
      m_local_incoming_count(seg_allocator),
      m_owned_info(seg_allocator),
      m_owned_info_tracker(seg_allocator),
      // m_owned_targets(seg_allocator),
      m_delegate_degree_threshold(delegate_degree_threshold),
      m_delegate_info(seg_allocator),
      m_delegate_degree(seg_allocator),
      m_delegate_label(seg_allocator),
      // m_delegate_targets(seg_allocator),
      //m_map_delegate_locator(100, boost::hash<uint64_t>(),
      //    std::equal_to<uint64_t>(), seg_allocator),
      m_map_delegate_locator(seg_allocator),
      m_controller_locators(seg_allocator), 
      m_has_edge_data(_has_edge_data)
      {
 
    typedef typename Container::value_type edge_type;
    typedef typename edge_data_type::value_type edge_data_value_type;
  
    static_assert(std::is_same<edge_type, std::tuple<uint64_t, uint64_t, 
                               edge_data_value_type> >::value, 
                               "Container type mismatch error!");

  CHK_MPI( MPI_Comm_size(m_mpi_comm, &m_mpi_size) );
  CHK_MPI( MPI_Comm_rank(m_mpi_comm, &m_mpi_rank) );
  
  processes_per_node = havoqgt_env()->node_local_comm().size();
  node_partitions = std::min(_node_partitions,uint64_t(processes_per_node));
  edge_chunk_size = _chunk_size;

  m_global_max_vertex = max_vertex;
  m_max_vertex = (std::ceil(double(max_vertex) / double(m_mpi_size)));

  LogStep logstep_main("Delegate Partitioning", m_mpi_comm, m_mpi_rank);

////////////////////////////////////////////////////////////////////////////////
/// Meta data phase of graph construction
////////////////////////////////////////////////////////////////////////////////

  {
    LogStep logstep("Allocating 4 Arrays of length max local vertex.",
      m_mpi_comm, m_mpi_rank);
    m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0));
    m_owned_info_tracker.resize(m_max_vertex+2, 0);
    m_local_outgoing_count.resize(m_max_vertex+1, 0);
    m_local_incoming_count.resize(m_max_vertex+1, 0);
    MPI_Barrier(m_mpi_comm);
  }

  boost::unordered_set<uint64_t> global_hubs;
 
  {
    LogStep logstep("count_edge_degree", m_mpi_comm, m_mpi_rank);
    count_edge_degrees(edges.begin(), edges.end(), global_hubs,
      delegate_degree_threshold);
    if (m_mpi_rank == 0)
      std::cout << "\tNumber of Delegates: " << global_hubs.size() << std::endl;
    MPI_Barrier(m_mpi_comm);
  }

  {
      LogStep logstep("initialize_low_meta_data", m_mpi_comm, m_mpi_rank);
      initialize_low_meta_data(global_hubs);
        MPI_Barrier(m_mpi_comm);
  }

  {
      LogStep logstep("initialize_high_meta_data", m_mpi_comm, m_mpi_rank);
      initialize_high_meta_data(global_hubs);
        MPI_Barrier(m_mpi_comm);
  }

  {
    LogStep logstep("count_high_degree_edges", m_mpi_comm, m_mpi_rank);
    count_high_degree_edges(edges.begin(), edges.end(), global_hubs);
    MPI_Barrier(m_mpi_comm);
  }

  m_graph_state = MetaDataGenerated;
  if (m_graph_state != stop_after) {
    complete_construction(seg_allocator, mpi_comm, edges, _edge_data);
  }
 
}

////////////////////////////////////////////////////////////////////////////////
/// Edge storage allocation  phase of graph construction
////////////////////////////////////////////////////////////////////////////////
template <typename SegmentManager>
template <typename Container, typename edge_data_type>
void
delegate_partitioned_graph<SegmentManager>::
complete_construction(const SegmentAllocator<void>& seg_allocator,
  MPI_Comm mpi_comm, Container& edges, edge_data_type& _edge_data) {
  m_mpi_comm = mpi_comm;
  int temp_mpi_rank, temp_mpi_size;
  CHK_MPI( MPI_Comm_size(m_mpi_comm, &temp_mpi_size) );
  CHK_MPI( MPI_Comm_rank(m_mpi_comm, &temp_mpi_rank) );
  assert(temp_mpi_rank == m_mpi_rank);
  assert(temp_mpi_size == temp_mpi_size);

  std::map< uint64_t, std::deque<OverflowSendInfo> > transfer_info;
  switch (m_graph_state) {
    case MetaDataGenerated:
      {
        LogStep logstep("calculate_overflow", m_mpi_comm, m_mpi_rank);
        calculate_overflow(transfer_info);
            MPI_Barrier(m_mpi_comm);
      }

      {
        LogStep logstep("initialize_edge_storage", m_mpi_comm, m_mpi_rank);
        initialize_edge_storage(seg_allocator);      

        if (m_has_edge_data) {
          assert(_edge_data != NULL); 
          _edge_data.resize(*this);                
        }  

          MPI_Barrier(m_mpi_comm);
      }

      m_graph_state = EdgeStorageAllocated;
      

////////////////////////////////////////////////////////////////////////////////
/// Edge partitioning phase of graph construction
////////////////////////////////////////////////////////////////////////////////
  case EdgeStorageAllocated:
    {
      LogStep logstep("partition_low_degree", m_mpi_comm, m_mpi_rank);
      partition_low_degree(edges, _edge_data);
        MPI_Barrier(m_mpi_comm);
    }
    m_graph_state = LowEdgesPartitioned;

    {
      LogStep logstep("partition_high_degree", m_mpi_comm, m_mpi_rank);
      partition_high_degree(edges, transfer_info, _edge_data);
        MPI_Barrier(m_mpi_comm);
    }
    m_graph_state = HighEdgesPartitioned;

////////////////////////////////////////////////////////////////////////////////
/// Controller construction phase of graph construction
////////////////////////////////////////////////////////////////////////////////
  // all-reduce hub degree
  case HighEdgesPartitioned:
    {
      LogStep logstep("all-reduce hub degree", m_mpi_comm, m_mpi_rank);
      if(m_delegate_degree.size() > 0) {
        mpi_all_reduce_inplace(m_delegate_degree, std::plus<uint64_t>(), m_mpi_comm);
      }
        MPI_Barrier(m_mpi_comm);
    }
    assert(m_delegate_degree.size() == m_delegate_label.size());

    //
    // Build controller lists
    {
      LogStep logstep("Build controller lists", m_mpi_comm, m_mpi_rank);
      const int controllers = (m_delegate_degree.size() / m_mpi_size) + 1;
      m_controller_locators.reserve(controllers);
      for (size_t i=0; i < m_delegate_degree.size(); ++i) {
        if (int(i % m_mpi_size) == m_mpi_rank) {
          m_controller_locators.push_back(vertex_locator(true, i, m_mpi_rank));
        }
      }
        MPI_Barrier(m_mpi_comm);
    }

    //
    // Verify CSR integration properlly tagged owned delegates
    {
      LogStep logstep("Verify CSR integration", m_mpi_comm, m_mpi_rank);
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
        MPI_Barrier(m_mpi_comm);
    }
    m_graph_state = GraphReady;
  }  // switch
////////////////////////////////////////////////////////////////////////////////
/// End of graph construction
////////////////////////////////////////////////////////////////////////////////
};

/**
 * This function iterates (1) through the edges and calculates the following:
 *
 *   m_local_incoming_count: For each vertedx that is assigned to this node it
 *   determines the number of incoming edges
 *
 *   high_vertex_count: the number of high edges generated b
 *
 *  @paramm_mpi_comm: MPI communication
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
count_edge_degrees(InputIterator unsorted_itr, InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hubs,
                 uint64_t delegate_degree_threshold) {
  using boost::container::map;

  uint64_t high_vertex_count(0);

  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;

  double start_time, last_loop_time;
  start_time = last_loop_time = MPI_Wtime();


  // Loop until no processor is producing edges
  while (!detail::global_iterator_range_empty(unsorted_itr,
        unsorted_itr_end, m_mpi_comm)) {
    if (m_mpi_rank == 0 && (loop_counter% 1000) == 0) {
      double curr_time = MPI_Wtime();

      std::cout << "\t["
        << "Total Loops: " << loop_counter << ", "
        << "Total Edges: " << edge_counter
        <<  "] Time " << (curr_time - last_loop_time) << " second, "
        << " Total Time " << (curr_time - start_time) << " second, "
        << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
        << std::flush;


      last_loop_time = curr_time;
    }
    loop_counter++;

    std::vector<
      boost::container::map< int, std::pair<uint64_t, uint64_t> >
    > maps_to_send(m_mpi_size);
    int maps_to_send_element_count = 0;

    // Generate Enough information to send
    for (size_t i = 0; i < edge_chunk_size && unsorted_itr != unsorted_itr_end; i++) {
      edge_counter++;
      auto edge = std::make_pair(std::get<0>(*unsorted_itr), std::get<1>(*unsorted_itr));

      // Update this vertex's outgoing edge count (first member of the pair)
      uint64_t local_id = local_source_id(m_mpi_size)(edge);
      int owner    = owner_source_id(m_mpi_size)(edge);
      if (owner == m_mpi_rank) {
        m_local_outgoing_count[local_id]++;
        if (m_local_outgoing_count[local_id] == delegate_degree_threshold) {
          high_vertex_count++;
        }
      } else {
        int c = maps_to_send.at(owner)[local_id].first++;
        if (c == 0 && maps_to_send.at(owner)[local_id].second ==0) {
          maps_to_send_element_count++;
        }
      }

      // Update the vertex's incoming edge count (second member of the pair)
      local_id = local_dest_id(m_mpi_size)(edge);
      owner    = owner_dest_id(m_mpi_size)(edge);
      if (owner == m_mpi_rank) {
        m_local_incoming_count[local_id]++;
        // if (m_local_incoming_count[local_id] == delegate_degree_threshold) {
        //   high_vertex_count++;
        // }
      } else {
        int c = maps_to_send.at(owner)[local_id].second++;
        if (c == 0 && maps_to_send.at(owner)[local_id].first ==0) {
          maps_to_send_element_count++;
        }
      }

      unsorted_itr++;
    }  // for until threshold is reached


    //Send Vertex degree count information to other nodes.
    send_vertex_info(high_vertex_count, delegate_degree_threshold,
        maps_to_send, maps_to_send_element_count);

  }  // while more edges
  mpi_yield_barrier(m_mpi_comm);
  if (m_mpi_rank == 0 ) {
    double curr_time = MPI_Wtime();
    loop_counter++;
    std::cout << "\t["
      << "Total Loops: " << loop_counter << ", "
      << "Total Edges: " << edge_counter
      <<  "] Time " << (curr_time - last_loop_time) << " second, "
      << " Total Time " << (curr_time - start_time) << " second, "
      << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
      << std::flush;
  }
  mpi_yield_barrier(m_mpi_comm);


  // Now, the m_local_incoming_count contains the total incoming and outgoing
  // edges for each vertex owned by this node.
  // Using this information we identify the hubs.
  std::vector<uint64_t> temp_hubs;
  temp_hubs.reserve(high_vertex_count);
  for (size_t i = 0; i < m_local_incoming_count.size(); i++) {
    const uint64_t outgoing = m_local_outgoing_count[i];
    // const uint64_t incoming = m_local_incoming_count[i];

    if (outgoing >= delegate_degree_threshold) {
      const uint64_t global_id = (i * m_mpi_size) + m_mpi_rank;
      assert(global_id != 0);
      temp_hubs.push_back(global_id);
    } else {
      m_edges_low_count +=  m_local_outgoing_count[i];
    }
  }

  assert(temp_hubs.size() == high_vertex_count);

  // Gather the hub lists and add them to the map,
  std::vector<uint64_t> vec_global_hubs;
  mpi_yield_barrier(m_mpi_comm);
  mpi_all_gather(temp_hubs, vec_global_hubs, m_mpi_comm);
  // if (m_mpi_rank == 0) {
  //   std::cout << "[";
  //   for (size_t i = 0; i < vec_global_hubs.size(); i++) {
  //     std::cout << vec_global_hubs[i] << ",";
  //   }
  //   std::cout << "]" << std::endl;

  // }
  // Insert gathered global hubs to set
  global_hubs.insert(vec_global_hubs.begin(), vec_global_hubs.end());


}  // count_edge_degrees

/**
 * This function is used to send/recv information about vertexes during the
 * count_edge_degrees function.
 *
 * @paramm_mpi_comm: them_mpi_communication group
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
send_vertex_info(uint64_t& high_vertex_count, uint64_t delegate_degree_threshold,
  std::vector< boost::container::map< int, std::pair<uint64_t, uint64_t> >  >&
  maps_to_send, int maps_to_send_element_count) {


  int to_send_pos = 0;
  std::vector<uint64_t> to_send;//(maps_to_send_element_count*3, 0);
  std::vector<int> to_send_count(m_mpi_size, 0);

  assert(maps_to_send.size() == m_mpi_size);
  for (size_t i = 0; i < maps_to_send.size(); i++) {
    for (auto itr = maps_to_send[i].begin(); itr != maps_to_send[i].end(); itr++) {
      //assert(to_send_pos < to_send.size());
      std::pair<int, std::pair<uint64_t, uint64_t>> triple = (*itr);
      //to_send[to_send_pos++] = uint64_t(triple.first);
      //to_send[to_send_pos++] = triple.second.first;
      //to_send[to_send_pos++] = triple.second.second;
      
      to_send.push_back(uint64_t(triple.first));
      to_send.push_back(triple.second.first);
      to_send.push_back(triple.second.second);
    }
    to_send_count[i] = maps_to_send[i].size()*3;
  }

  std::vector<uint64_t> to_recv;
  std::vector<int> out_recvcnts;

  mpi_yield_barrier(m_mpi_comm);
  mpi_all_to_all(to_send, to_send_count,to_recv, out_recvcnts, m_mpi_comm);

  for (size_t k = 0; k < to_recv.size(); ) {
    const uint64_t local_id = to_recv[k++];
    const uint64_t source_count = to_recv[k++];
    const uint64_t dest_count = to_recv[k++];
    assert(local_id < m_local_incoming_count.size());

    // If its not currently a high vertex but by adding this it becomes one
    // then increment high_vertex_count
    // if (m_local_incoming_count[local_id] < delegate_degree_threshold
    //   && m_local_incoming_count[local_id] + dest_count >=
    //   delegate_degree_threshold) {

    //   high_vertex_count++;
    // }
    if (m_local_outgoing_count[local_id] < delegate_degree_threshold
	&& m_local_outgoing_count[local_id] + source_count >= // changed from dest_count to source_count
      delegate_degree_threshold) {

      high_vertex_count++;
    }

    m_local_outgoing_count[local_id] += source_count;
    m_local_incoming_count[local_id] += dest_count;
  }  // for each recieved element.

}  // send_vertex_info



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
initialize_low_meta_data(boost::unordered_set<uint64_t>& global_hubs) {


  // Allocate the index into the low edge csr.
  // +2: because the m_max_vertex is indexible and the last position must hold
  // the number of low edges.
  // m_owned_info.resize(m_max_vertex+2, vert_info(false, 0, 0)); moved to
  // constructor


  // Initilize the m_owned_info, by iterating through owned vertexes and
  //  if it is not a hub, then it incremenets the edge count by the number of
  //  outgoing edges.
  uint64_t edge_count = 0;
  for (uint64_t vert_id = 0; vert_id < m_owned_info.size(); vert_id++) {
    const uint64_t outgoing = m_local_outgoing_count[vert_id];
    const uint64_t incoming = m_local_incoming_count[vert_id];

    m_owned_info[vert_id] = vert_info(false, 0, edge_count);

    if (outgoing < m_delegate_degree_threshold) {
      edge_count += outgoing;
    } else {
      #ifdef DEBUG_DPG
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
}

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
initialize_high_meta_data(boost::unordered_set<uint64_t>& global_hubs) {
  //
  // Setup and Compute Hub information
  //
  for (int i = 0; i < processes_per_node; i++) {
    if (m_mpi_rank == 0) {
      std::cout << "\t Setup and Compute Hub information: " << i << "/"
        << processes_per_node << "." << std::endl << std::flush;
    }
    if (i == m_mpi_rank % processes_per_node) {
      std::vector<uint64_t> vec_sorted_hubs(global_hubs.begin(), global_hubs.end());
      std::sort(vec_sorted_hubs.begin(), vec_sorted_hubs.end());

      // Allocate space for the delegate csr index.
      // This is initlized during the paritioning of the low edges and then adjusted
      // in initialize_delegate_target

      #if 1
        m_delegate_info.resize(global_hubs.size()+1);
        m_delegate_degree.resize(vec_sorted_hubs.size());

      #else

      {
        size_t loop_limit = global_hubs.size()+1;
        m_delegate_info.reserve(loop_limit);
        for (size_t j = 0; j < loop_limit; ) {
          for (size_t k = 0; k < 10000000; k++) {
            m_delegate_info.emplace_back();
            j++;
            if (j >= loop_limit) {
              break;
            }
          }
        }
      }

      // Allocates and initilize the delegate (AKA hub) vertex infromation
      {
        size_t loop_limit = vec_sorted_hubs.size();
        m_delegate_degree.reserve(loop_limit);
        for (size_t j = 0; j < loop_limit; ) {
          for (size_t k = 0; k < 10000000; k++) {
            m_delegate_degree.emplace_back();
            j++;
            if (j >= loop_limit) {
              break;
            }
          }
        }
      }
      #endif

      // Loop over the hub vertexes:
      //  initilizing the delegate_degree tracking structures
      for(size_t i=0; i<vec_sorted_hubs.size(); ++i) {
        uint64_t t_local_id = vec_sorted_hubs[i] / uint64_t(m_mpi_size);
        int t_owner = uint32_t(vec_sorted_hubs[i] % uint32_t(m_mpi_size));
        vertex_locator new_ver_loc(true, i, t_owner);

        m_map_delegate_locator[vec_sorted_hubs[i]] = new_ver_loc;
        m_delegate_label.push_back(vec_sorted_hubs[i]);

        //
        // Tag owned delegates
        //
        if (t_owner == m_mpi_rank) {
          m_owned_info[t_local_id].is_delegate = 1;
          m_owned_info[t_local_id].delegate_id = i;
        }
      }  // for over vec_sorted_hubs

    }  // if i == mpi_rank % processes_per_node
    mpi_yield_barrier(m_mpi_comm);
  }
}

template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
initialize_edge_storage(const SegmentAllocator<void>& seg_allocator) {
  // Allocate the low edge csr to accommdate the number of edges
  // This will be filled by the partion_low_edge function
  for (int i = 0; i < processes_per_node; i++) {
    if (m_mpi_rank == 0) {
      std::cout << "\tResizing m_owned_targets and m_delegate_targets: "
        << i << "/" << processes_per_node << "." << std::endl << std::flush;
    }

    if (i == m_mpi_rank % processes_per_node) {
      // m_owned_targets.resize(m_edges_low_count);
      // m_delegate_targets.resize(m_edges_high_count);
      m_delegate_targets_size = m_edges_high_count;
      m_owned_targets_size = m_edges_low_count;

      SegmentManager *segment_manager = seg_allocator.get_segment_manager();

      {
        //void * temp = segment_manager->allocate(
        //   m_delegate_targets_size * sizeof(vertex_locator));
        //m_delegate_targets = reinterpret_cast<vertex_locator *>(temp);
        m_delegate_targets = (vertex_locator*) segment_manager->allocate(m_delegate_targets_size * sizeof(vertex_locator));
      }

      {
        //void * temp = segment_manager->allocate(
        //  m_owned_targets_size * sizeof(vertex_locator));
       // m_owned_targets = reinterpret_cast<vertex_locator *>(temp);
        m_owned_targets = (vertex_locator*) segment_manager->allocate(m_owned_targets_size * sizeof(vertex_locator));
      }

      // Currently, m_delegate_info holds the count of high degree edges
      // assigned to this node for each vertex.
      // Below converts it into an index into the m_delegate_targets array
      assert( m_delegate_info[ m_delegate_info.size()-1] == 0);
      uint64_t edge_count = 0;
      for (size_t i=0; i < m_delegate_info.size(); i++) {
        uint64_t num_edges = m_delegate_info[i];
        m_delegate_info[i] = edge_count;
        edge_count += num_edges;
        assert(edge_count <= m_edges_high_count);
      }

    }  // If this processes
    mpi_yield_barrier(m_mpi_comm);
  }  // for over processes

};

/**
 * This is a hack to determine the tuple index at compile time.
 * std::get<index>(tupel) wont compile if the value of "index" is not known at 
 * compile time. For std::pair or two element tuple, index is either 0 or 1. 
 */
constexpr std::size_t get_edge_data_tuple_index(
  bool is_edge_pair_type = false, bool is_edge_two_tuple_type = false) {
  return ((is_edge_pair_type == true) || (is_edge_two_tuple_type == true)) ? 0 : 2; 
  // c++11: must be a single return statemnet 
}

/**
 * This function iterates (2) through the edges and sends the low degree edges
 * to the nodes that own them.
 *
 * At the same time it tracks the number of outgoing edges for each delegate
 * vertex and exchanges that information with the other nodes.
 *
 */
template <typename SegmentManager>
template <typename Container, typename edge_data_type>
void
delegate_partitioned_graph<SegmentManager>::
partition_low_degree(Container& unsorted_edges, edge_data_type& _edge_data) {
  typedef typename Container::value_type edge_type;
  typedef std::pair<uint64_t, uint64_t> edge_pair_type; 
  constexpr bool is_edge_pair_type = std::is_same<edge_type, edge_pair_type>::value;
  typedef std::tuple<uint64_t, uint64_t> edge_two_tuple_type;
  constexpr bool is_edge_two_tuple_type = std::is_same<edge_type, edge_two_tuple_type>::value; 
  
  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;

  double start_time, last_loop_time, last_part_time;
  start_time = last_loop_time = last_part_time = MPI_Wtime();


  for (size_t node_turn = 0; node_turn < node_partitions; node_turn++) {

    auto unsorted_itr     = unsorted_edges.begin();
    auto unsorted_itr_end = unsorted_edges.end();
    if (m_mpi_rank == 0) {
      double curr_time = MPI_Wtime();

      std::cout << "\t***[LP]["
        << "Partition Number: " << node_turn << "/" << node_partitions << ", "
        << "Total Loops: " << loop_counter << ", "
        << "Total Edges: " << edge_counter
        <<  "] Time " << (curr_time - last_part_time) << " second, "
        << " Total Time " << (curr_time - start_time) << " second, "
        << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
        << std::flush;
      last_part_time = curr_time;
    }
    MPI_Barrier(m_mpi_comm);


    while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
            m_mpi_comm)) {
      if (m_mpi_rank == 0 && (loop_counter% 500) == 0) {
        double curr_time = MPI_Wtime();

        std::cout << "\t[LP]["
          << "Partition Number: " << node_turn << "/" << node_partitions << ", "
          << "Total Loops: " << loop_counter << ", "
          << "Total Edges: " << edge_counter
          <<  "] Time " << (curr_time - last_loop_time) << " second, "
          << " Total Time " << (curr_time - start_time) << " second, "
          << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
          << std::flush;


        last_loop_time = curr_time;
      }
      loop_counter++;

      // Generate Edges to Send
      std::vector<edge_type> to_recv_edges_low;    
      
      {
        std::vector<edge_type> to_send_edges_low; 
        to_send_edges_low.reserve(edge_chunk_size);

        for (size_t i = 0;
          unsorted_itr != unsorted_itr_end && i < edge_chunk_size;
          ++unsorted_itr, ++i) {
          // Get next edge
          const auto edge = std::make_pair(std::get<0>(*unsorted_itr), std::get<1>(*unsorted_itr));

          {
            const int owner = edge.first %m_mpi_size;
            if ( (owner % processes_per_node) % node_partitions != node_turn ) {
              continue;
            }
          }
          edge_counter++;

          if (m_map_delegate_locator.count(edge.first) == 0) {
            to_send_edges_low.push_back(*unsorted_itr);
            //++i;
          } else {
            continue;
          }
        }  // for

        // Exchange Edges/Recieve edges
        edge_source_partitioner<edge_type> paritioner(m_mpi_size);
        mpi_yield_barrier(m_mpi_comm);
        mpi_all_to_all_better(to_send_edges_low, to_recv_edges_low, paritioner,
            m_mpi_comm);
      }

      std::sort(to_recv_edges_low.begin(), to_recv_edges_low.end());

  #ifdef DEBUG_DPG
      // Sanity Check to make sure we recieve the correct edges
      for (size_t i = 0; i<to_recv_edges_low.size(); ++i) {
        auto edge =  std::make_pair(std::get<0>(to_recv_edges_low[i]), std::get<1>(to_recv_edges_low[i]));
        assert(int(edge.first) % m_mpi_size ==m_mpi_rank);
        assert(m_map_delegate_locator.count(edge.first) == 0);
      }
  #endif
      // Loop over recieved edges, appending them to the low CSR
      auto itr_end = to_recv_edges_low.end();
     
      for (auto itr = to_recv_edges_low.begin(); itr != itr_end; itr++) {

        auto edge = std::make_pair(std::get<0>(*itr), std::get<1>(*itr)); //std::get<0>(*itr);
         
        uint64_t new_vertex_id = local_source_id(m_mpi_size)(edge);
        assert(m_mpi_rank == int(edge.first % m_mpi_size));

        uint64_t temp_offset = (m_owned_info_tracker[new_vertex_id])++;
        uint64_t loc = temp_offset + m_owned_info[new_vertex_id].low_csr_idx;


        if (!(loc <  m_owned_info[new_vertex_id+1].low_csr_idx)) {
          std::cout << "Error rank: " << m_mpi_rank << " -- " <<
                   loc << " < " <<  m_owned_info[new_vertex_id+1].low_csr_idx
                   << ", new_vertex_id = " << new_vertex_id
                   << ", temp_offset = " << temp_offset
                   << ", edge = (" << edge.first << "," << edge.second << ")"
          << std::endl << std::flush;
          assert(false);
          exit(-1);
        }
        /// @todo this was tripping, is this old?   Could be left over from before when targets was vector based.
        //assert(!m_owned_targets[loc].is_valid());

        m_owned_targets[loc] = label_to_locator(edge.second);
        if(m_has_edge_data) {          
          _edge_data.m_owned_edge_data[loc] = std::get<get_edge_data_tuple_index(is_edge_pair_type, is_edge_two_tuple_type)>(*itr);         
        } 
      }  // for over recieved egdes
    }  // while global iterator range not empty
  }  // for node partition

  mpi_yield_barrier(m_mpi_comm);
  if (m_mpi_rank == 0) {
    double curr_time = MPI_Wtime();
    std::cout << "\t[LP]["
      << "Total Loops: " << loop_counter << ", "
      << "Total Edges: " << edge_counter
      <<  "] Time " << (curr_time - last_loop_time) << " second, "
      << " Total Time " << (curr_time - start_time) << " second, "
      << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
      << std::flush;
  }
  mpi_yield_barrier(m_mpi_comm);

}  // partition_low_degree



/**
 * This function iterates (through the edges and tracks the number of
 * outgoing edges for each delegate vertex and exchanges that information with
 * the other nodes.
 *
 */
template <typename SegmentManager>
template <typename InputIterator>
void
delegate_partitioned_graph<SegmentManager>::
count_high_degree_edges(InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set) {
  // Temp Vector for storing offsets
  // Used to store high_edge count
  #if DEBUG_DPG
    std::vector<uint64_t> tmp_high_count_per_rank(m_mpi_size, 0);
  #endif
  m_edges_high_count = 0;

  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;

  double start_time, last_loop_time;
  start_time = last_loop_time = MPI_Wtime();


  while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
          m_mpi_comm)) {
    if (m_mpi_rank == 0 && (loop_counter% 1000) == 0) {
      double curr_time = MPI_Wtime();

      std::cout << "\t[CH]["
        << "Total Loops: " << loop_counter << ", "
        << "Total Edges: " << edge_counter
        <<  "] Time " << (curr_time - last_loop_time) << " second, "
        << " Total Time " << (curr_time - start_time) << " second, "
        << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
        << std::flush;


      last_loop_time = curr_time;
    }
    loop_counter++;


    // Vector used to pass number of high edges
    std::vector<
      boost::container::map<uint64_t, uint64_t> > maps_to_send(m_mpi_size);
    int maps_to_send_element_count = 0;
    {
      for (size_t i=0; unsorted_itr != unsorted_itr_end && i < edge_chunk_size;
           ++unsorted_itr, ++i) {

        // Get next edge
        const auto edge = std::make_pair(std::get<0>(*unsorted_itr), std::get<1>(*unsorted_itr));
        edge_counter++;


        if (global_hub_set.count(edge.first) == 0) {
          continue;
        } else if(global_hub_set.count(edge.first)) {
          #if DEBUG_DPG
            // This edge's source is a hub
            // 1) Increment the high edge count for the owner of the edge's dest
            tmp_high_count_per_rank[edge.second %m_mpi_size]++;
          #endif

          // 2) Increment the owner's count of edges for this hub.
          const int owner = edge.second %m_mpi_size;
          if (owner == m_mpi_rank) {
            const uint64_t ver_id = edge.first;

            const uint64_t new_source_id = m_map_delegate_locator[ver_id].local_id();
            assert(new_source_id < m_delegate_info.size()-1);
            m_delegate_info[new_source_id]++;
            m_edges_high_count++;
          } else {
            int c = maps_to_send.at(owner)[edge.first]++;
            if (c == 0) {
              maps_to_send_element_count++;
              //i++;
            }
          }
        } else {
          assert(false);
        }
      }  // for
      // Send the hub edge count to the relevent nodes.
      send_high_info(maps_to_send, maps_to_send_element_count);
    }

  }  // while global iterator range not empty


  mpi_yield_barrier(m_mpi_comm);
  if (m_mpi_rank == 0) {
    double curr_time = MPI_Wtime();
    std::cout << "\t[CH]["
      << "Total Loops: " << loop_counter << ", "
      << "Total Edges: " << edge_counter
      <<  "] Time " << (curr_time - last_loop_time) << " second, "
      << " Total Time " << (curr_time - start_time) << " second, "
      << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
      << std::flush;
  }

  #if DEBUG_DPG
    std::cout << "Checking count_high_degree_edges:" << std::endl << std::flush;
    size_t l_edges_high_count = 0;
    for (size_t i = 0; i < m_delegate_info.size(); i++) {
      l_edges_high_count += m_delegate_info[i];
    }
    assert(l_edges_high_count == m_edges_high_count);

    // Sync The high counts.
    std::vector<uint64_t> high_count_per_rank;
    mpi_all_reduce(tmp_high_count_per_rank, high_count_per_rank,
        std::plus<uint64_t>(), m_mpi_comm);

    uint64_t sanity_check_high_edge_count = high_count_per_rank[m_mpi_rank];
    assert(m_edges_high_count == sanity_check_high_edge_count);
  #endif

}  // count_high_degree_edges


/**
 * This function has each node exchange with one another their delegate vertex
 * incoming edge count.
 *
 * @paramm_mpi_comm: the mpi communication group
 * @param maps_to_send: a vector of maps that map delegate_id to edge count.
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
send_high_info(std::vector< boost::container::map< uint64_t, uint64_t> >&
  maps_to_send, int maps_to_send_element_count) {

  int to_send_pos = 0;
  std::vector<uint64_t> to_send;//(maps_to_send_element_count*2, 0);
  std::vector<int> to_send_count(m_mpi_size, 0);

  assert(maps_to_send.size() == m_mpi_size);
  for (size_t i = 0; i < maps_to_send.size(); i++) {
    for (auto itr = maps_to_send[i].begin(); itr != maps_to_send[i].end(); itr++) {
      //assert(to_send_pos < to_send.size());
      //to_send[to_send_pos++] = itr->first;
      //to_send[to_send_pos++] = itr->second;
      to_send.push_back(itr->first);
      to_send.push_back(itr->second);
    }
    to_send_count[i] = maps_to_send[i].size()*2;
  }

  std::vector<uint64_t> to_recv;
  std::vector<int> out_recvcnts;

  mpi_yield_barrier(m_mpi_comm);
  mpi_all_to_all(to_send, to_send_count,to_recv, out_recvcnts, m_mpi_comm);

  for (size_t i = 0; i < to_recv.size(); i++) {
    const uint64_t ver_id = to_recv[i++];
    const uint64_t delegate_dest_count = to_recv[i];

    const uint64_t new_source_id = m_map_delegate_locator[ver_id].local_id();
    assert(new_source_id < m_delegate_info.size()-1);
    m_delegate_info[new_source_id] += delegate_dest_count;
    m_edges_high_count += delegate_dest_count;
  }
}  // send_high_info
   //

/**
 * This function initlizes the transfer_info object by determining which nodes
 * needs extra edges and which nodes have extra edges. This information is
 * exchanged so that nodes that need more edges will know how many extra edges
 * it will recieve.
 *
 * @paramm_mpi_comm: the mpi communication group
 * @param &owned_high_count: tracks the number of high edges owned by this node
 * it is updated when we give to or recieve edges from another node.
 * @param transfer_info: used to track to whome and how many edges are given
 * to another node
 */
template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
calculate_overflow(std::map< uint64_t, std::deque<OverflowSendInfo> >
    &transfer_info) {

  //
  // Get the number of high and low edges for each node.
  //
  mpi_yield_barrier(m_mpi_comm);
  std::vector<uint64_t> low_count_per_rank, high_count_per_rank;
  mpi_all_gather(m_edges_low_count, low_count_per_rank, m_mpi_comm);
  mpi_all_gather(m_edges_high_count, high_count_per_rank, m_mpi_comm);

  // Determine the total of edges accorss all nodes.
  const uint64_t owned_total_edges = m_edges_high_count + m_edges_low_count;
  uint64_t global_edge_count = mpi_all_reduce(owned_total_edges,
      std::plus<uint64_t>(), m_mpi_comm);

  // Determine the desired number of edges at each node.
  const uint64_t target_edges_per_rank = global_edge_count / m_mpi_size;

  uint64_t gave_edge_counter = 0;
  uint64_t recieve_edge_counter = 0;


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

          for (int i = 0; i < send_list.size();) {
            const uint64_t vert_id = send_list[i++];
            const uint64_t count = send_list[i++];
            #if 0
              std::cout << "[" << m_mpi_rank << "] giving " << vert_id << " +" << count
              << " to " << light_idx << ". "<< std::endl << std::flush;
            #endif

          }

          // Send the information
          int64_t send_len = send_list.size();
          CHK_MPI(MPI_Send(&send_len, 1, mpi_typeof(send_len), light_idx, 0, m_mpi_comm));
          CHK_MPI(MPI_Send(send_list.data(), send_len, mpi_typeof(to_move),
              light_idx, 0, m_mpi_comm));

          // Adjust the m_edges_high_count info.
          m_edges_high_count -= to_move;
          gave_edge_counter += to_move;
        } else if (light_idx == m_mpi_rank) {  // This node is reciving
          MPI_Status status;
          int64_t recv_length;

          // Recieve the information
          CHK_MPI(MPI_Recv(&recv_length, 1, mpi_typeof(recv_length), heavy_idx,
              0, m_mpi_comm, &status));
          std::vector<uint64_t> recv_list(recv_length);
          CHK_MPI(MPI_Recv(recv_list.data(), recv_length, mpi_typeof(to_move),
               heavy_idx, 0, m_mpi_comm, &status));


          // Update my delagate edge counts.
          uint64_t sanity_count = 0;
          for (int i = 0; i < recv_length;) {
            const uint64_t vert_id = recv_list[i++];
            const uint64_t count = recv_list[i++];
            m_delegate_info[vert_id] += count;
            sanity_count += count;

            #if 0
              std::cout << "[" << m_mpi_rank << "] recieved " << vert_id <<
               ", now has " << count << " edges for it." << std::endl << std::flush;
            #endif

          }

          // Adjust the m_edges_high_count info.
          m_edges_high_count += to_move;
          recieve_edge_counter += to_move;
          assert(sanity_count == to_move);
        } // else this node is not involved.
        mpi_yield_barrier(m_mpi_comm);
      } else {
        ++light_idx;
        if (light_idx == m_mpi_size) {
          break;
        }
      } // else
    }  // While
  }  // For

#ifdef DEBUG_DPG
  const uint64_t owned_total_edges2 = m_edges_high_count + m_edges_low_count;
  uint64_t sanity_global_edge_count = mpi_all_reduce(owned_total_edges2,
      std::plus<uint64_t>(), m_mpi_comm);
  assert(sanity_global_edge_count == global_edge_count);
  assert(m_edges_high_count == high_count_per_rank[m_mpi_rank]);
  uint64_t high_count_2 = 0;
  for (size_t i = 0; i < m_delegate_info.size()-1; i++) {
    high_count_2 += m_delegate_info[i];
  }

  assert(m_edges_high_count == high_count_2);

  std::vector<uint64_t> low_count_per_rank2, high_count_per_rank2;
  mpi_all_gather(m_edges_low_count, low_count_per_rank2, m_mpi_comm);
  mpi_all_gather(m_edges_high_count, high_count_per_rank2, m_mpi_comm);

  for (size_t i = 0; i < m_mpi_size; i++) {
    assert(low_count_per_rank2[i] == low_count_per_rank[i]);
    assert(high_count_per_rank2[i] == high_count_per_rank[i]);
  }

  #if 0
  for (int i = 0; i < m_mpi_size; i++) {
    if (i == m_mpi_rank) {
      std::cout << i << " has " << ((int64_t)owned_total_edges - (int64_t)target_edges_per_rank)
      << " extra edges. Givng:" << gave_edge_counter << " Recieving: "
      << recieve_edge_counter << "." << std::endl;
    }
    MPI_Barrier(m_mpi_comm);
  }
  #endif


#endif

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
 * This function iterates (3) over the edges and if an edge has a delegate
 * vertex as a source sends it to the apprioriate node based on the edge's
 * destination. If this node is giving edges to another node, then when this node has
 * enough edges for a delegate, extra edges are sent to a node that needs edges
 * for that delegate.
 *
 * @paramm_mpi_comm: The mpi communication group
 * @param unsorted_itr: An iterator over edges
 * @param unsorted_itr_end: The end of the iterator over edges
 * @param global_hub_set: A set of all global hubs
 * @param transfer_info: A map of delegate_id to a deque of OverflowSendInfo
 * used to determine where overflowed edges go.
 */
template <typename SegmentManager>
template <typename Container, typename edge_data_type>
void
delegate_partitioned_graph<SegmentManager>::
partition_high_degree(Container& unsorted_edges,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info, edge_data_type& _edge_data) {

  typedef typename Container::value_type edge_type;
  typedef std::pair<uint64_t, uint64_t> edge_pair_type;
  constexpr bool is_edge_pair_type = std::is_same<edge_type, edge_pair_type>::value;
  typedef std::tuple<uint64_t, uint64_t> edge_two_tuple_type;
  constexpr bool is_edge_two_tuple_type = std::is_same<edge_type, edge_two_tuple_type>::value;

  typedef typename edge_data_type::value_type edge_data_value_type;
  typedef std::tuple<uint64_t, uint64_t, edge_data_value_type> edge_type_comm; // used by the communication routines

  // Initates the paritioner, which determines where overflowed edges go
  high_edge_partitioner<edge_type_comm> paritioner(m_mpi_size, m_mpi_rank, &transfer_info);

  uint64_t loop_counter = 0;
  uint64_t edge_counter = 0;
  double start_time, last_loop_time, last_part_time;
  start_time = last_loop_time = last_part_time = MPI_Wtime();
  uint64_t gave_edge_counter = 0;

  // Scratch vector use for storing edges to send
  std::vector<edge_type_comm> to_send_edges_high;
  to_send_edges_high.reserve(edge_chunk_size);

  for (size_t node_turn = 0; node_turn < node_partitions; node_turn++) {

    if (m_mpi_rank == 0) {
      double curr_time = MPI_Wtime();

      std::cout << "\t***[HP]["
        << "Partition Number: " << node_turn << "/" << node_partitions << ", "
        << "Total Loops: " << loop_counter << ", "
        << "Total Edges: " << edge_counter
        <<  "] Time " << (curr_time - last_part_time) << " second, "
        << " Total Time " << (curr_time - start_time) << " second, "
        << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
        << std::flush;

      last_part_time = curr_time;
    }
    MPI_Barrier(m_mpi_comm);


    auto unsorted_itr     = unsorted_edges.begin();
    auto unsorted_itr_end = unsorted_edges.end();

    while (!detail::global_iterator_range_empty(unsorted_itr, unsorted_itr_end,
         m_mpi_comm)) {
      if (m_mpi_rank == 0 && (loop_counter % 500) == 0) {
        double curr_time = MPI_Wtime();

        std::cout << "\t[HP]["
          << "Partition Number: " << node_turn << "/" << node_partitions << ", "
          << "Total Loops: " << loop_counter << ", "
          << "Total Edges: " << edge_counter
          <<  "] Time " << (curr_time - last_loop_time) << " second, "
          << " Total Time " << (curr_time - start_time) << " second, "
          << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
          << std::flush;

        last_loop_time = curr_time;
      }
      loop_counter++;

      size_t i=0;
      while (unsorted_itr != unsorted_itr_end &&
             /*to_send_edges_high.size()*/i++ < edge_chunk_size) {
        // Get next edge
        const auto edge = std::make_pair(std::get<0>(*unsorted_itr), std::get<1>(*unsorted_itr));
        const auto edge_data = std::get<get_edge_data_tuple_index(is_edge_pair_type, is_edge_two_tuple_type)>(*unsorted_itr);
        ++unsorted_itr;

        {
            const int owner = edge.second %m_mpi_size;
            if (owner % processes_per_node % node_partitions != node_turn) {
              continue;
            }
        }

        edge_counter++;

        if (m_map_delegate_locator.count(edge.first) == 1) {
          // assert(global_hub_set.count(edge.first) > 0);
          // If the edge's source is a hub node
          const uint64_t source_id = edge.first;
          uint64_t new_source_id = m_map_delegate_locator[source_id].local_id();
          assert(new_source_id >=0 && new_source_id < m_delegate_info.size()-1);

          // Send the edge if we don't own it or if we own it but have no room.
          
          if (!m_has_edge_data) {            
              to_send_edges_high.push_back(std::forward_as_tuple(new_source_id, edge.second, (uint64_t)(0) ) );
          } else { 	
            to_send_edges_high.push_back(std::forward_as_tuple(new_source_id, edge.second, edge_data ) );    
          }
        }  // end if is a hub
        else {
          // assert(global_hub_set.count(edge.first) == 0);
        }

      }  // end while

      // Exchange edges we generated that we don't need with the other nodes and
      // recieve edges we may need
      // // Scratch vector use for storing recieved edges.
      std::vector<edge_type_comm> to_recv_edges_high;
      mpi_yield_barrier(m_mpi_comm);
      mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high, paritioner,
         m_mpi_comm);

      // Empty the vector
      {
        std::vector<edge_type_comm> temp; 
        to_send_edges_high.swap(temp);
      }
      to_send_edges_high.reserve(edge_chunk_size);

      assert(to_send_edges_high.size() == 0);
      std::sort(to_recv_edges_high.begin(), to_recv_edges_high.end());

      for (size_t i = 0; i < to_recv_edges_high.size(); ++i) {
        // Iterate over recieved edges, addiing them using similar logic from
        // above
        const auto edge = std::make_pair(std::get<0>(to_recv_edges_high[i]), std::get<1>(to_recv_edges_high[i]));
        const uint64_t new_source_id = edge.first;
        assert(new_source_id >=0 && new_source_id < m_delegate_info.size()-1);

        uint64_t place_pos = m_delegate_degree[new_source_id];
        place_pos += m_delegate_info[new_source_id];

        if (place_pos == m_delegate_info[new_source_id+1]) {
          // We have no room for this node, so lets send it to a node that has
          // room.
          assert(transfer_info.size() > 0);
          assert(transfer_info.count(new_source_id) != 0);
          to_send_edges_high.push_back(to_recv_edges_high[i]);
          gave_edge_counter++;
        }
        else {
          assert(place_pos < m_delegate_info[new_source_id+1]);
          assert(place_pos < m_delegate_targets_size);

          uint64_t new_target_label = edge.second;
          m_delegate_targets[place_pos] = label_to_locator(new_target_label);
          if(m_has_edge_data) {
	    _edge_data.m_delegate_edge_data[place_pos] = std::get<get_edge_data_tuple_index(is_edge_pair_type, is_edge_two_tuple_type)>(to_recv_edges_high[i]); 
          }
  
          assert(m_delegate_targets[place_pos].m_owner_dest < m_mpi_size);
          m_delegate_degree[new_source_id]++;

          if (owner_dest_id(m_mpi_size)(edge) != m_mpi_rank) {
            assert(transfer_info.size() == 0);
          }

        }  // else we have room
      }  // for edges recieved

    }  // end while get next edge
  } // For over node partitions
  mpi_yield_barrier(m_mpi_comm);
  if (m_mpi_rank == 0) {
    double curr_time = MPI_Wtime();
    std::cout << "\t[HP]["
      << "Total Loops: " << loop_counter << ", "
      << "Total Edges: " << edge_counter
      <<  "] Time " << (curr_time - last_loop_time) << " second, "
      << " Total Time " << (curr_time - start_time) << " second, "
      << "Dirty Pages: " << get_dirty_pages() << "kb." << std::endl
      << std::flush;
  }
  mpi_yield_barrier(m_mpi_comm);

  {//
  // Exchange edges we generated  with the other nodes and recieve edges we may need
    // // Scratch vector use for storing recieved edges.
    std::vector<edge_type_comm> to_recv_edges_high;
    mpi_yield_barrier(m_mpi_comm);
    mpi_all_to_all_better(to_send_edges_high, to_recv_edges_high, paritioner,
       m_mpi_comm);

    // Empty the vector
    {
      std::vector<edge_type_comm> temp;
      to_send_edges_high.swap(temp);
    }

    std::sort(to_recv_edges_high.begin(), to_recv_edges_high.end());
    for (size_t i=0; i<to_recv_edges_high.size(); ++i) {
      // Iterate over recieved edges, addiing them using similar logic from
      // above
      const auto edge = std::make_pair(std::get<0>(to_recv_edges_high[i]), std::get<1>(to_recv_edges_high[i])); //std::get<0>(to_recv_edges_high[i]);
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
        assert(place_pos < m_delegate_targets_size);

        uint64_t new_target_label = edge.second;
        m_delegate_targets[place_pos] = label_to_locator(new_target_label);
        if(m_has_edge_data) {
          _edge_data.m_delegate_edge_data[place_pos] = std::get<get_edge_data_tuple_index(is_edge_pair_type, is_edge_two_tuple_type)>(to_recv_edges_high[i]);
        }
 
        assert(m_delegate_targets[place_pos].m_owner_dest < m_mpi_size);
        m_delegate_degree[new_source_id]++;

        if (owner_dest_id(m_mpi_size)(edge) != m_mpi_rank) {
          assert(transfer_info.size() == 0);
        }

      }  // else there is have room
    }  // for edges recieved
  }

#ifdef DEBUG_DPG
  int64_t recv_count2 = 0;
  for (size_t i = 0; i < m_delegate_degree.size(); i++) {
     recv_count2 += m_delegate_degree[i];
  }

  int64_t difference = recv_count2 - m_delegate_targets_size;

  for (size_t i = 0; i < m_delegate_info.size()-1; i++) {
    size_t pos = m_delegate_info[i];
    for (size_t j = pos; j < m_delegate_info[i+1]; j++) {
      assert(m_delegate_targets[j].m_owner_dest != 0xFFFFF);
      assert(m_delegate_targets[j].m_local_id != 0x7FFFFFFFFF);
    }
  }
#endif

}  // partition_high_degre


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// MISC FUNCTIONS ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
    res = uint64_t(locator.local_id()) * uint64_t(m_mpi_size) +
         uint64_t(locator.owner());
  }
  return res;
}  // locator_to_label

/**
 * @param  label vertex label to convert
 * @return locator for the label
 */
template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_locator
delegate_partitioned_graph<SegmentManager>::
label_to_locator(uint64_t label) const {
  /*typename boost::unordered_map< uint64_t, vertex_locator,
              boost::hash<uint64_t>, std::equal_to<uint64_t>,
              SegmentAllocator< std::pair<uint64_t,vertex_locator> >
             >::const_iterator*/
 auto itr = m_map_delegate_locator.find(label);

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
  return edge_iterator(locator, m_owned_info[locator.local_id()].low_csr_idx,
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
  size_t last_index = m_owned_info.empty() ? 0 : m_owned_info.size()-1;
  while(last_index > 1 && (m_owned_info[last_index-1].low_csr_idx == m_owned_info[last_index].low_csr_idx)) {
    --last_index;
  }
  return vertex_iterator(last_index,this);
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::
delegate_vertices_begin() const {
  return vertex_iterator(0, this, true);
}

template <typename SegmentManager>
inline
typename delegate_partitioned_graph<SegmentManager>::vertex_iterator
delegate_partitioned_graph<SegmentManager>::
delegate_vertices_end() const {
  return vertex_iterator(m_delegate_info.size()-1, this, true);
}

template <typename SegmentManager>
inline
bool
delegate_partitioned_graph<SegmentManager>::
is_label_delegate(uint64_t label) const {
  return m_map_delegate_locator.count(label) > 0;
}

/**
 * @todo remove; make like vertex_data
 */ 
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
            (m_owned_targets_size, m_delegate_targets_size,
              segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
            (m_owned_targets_size, m_delegate_targets_size,
              segment_manager_o);
  }
}

/**
 * @param   init initial value for each vertex
 * @todo remove; make like vertex_data
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
            (m_owned_targets_size, m_delegate_targets_size, init,
              segment_manager_o);
  } else {
    return segment_manager_o->template construct<mytype>(obj_name)
            (m_owned_targets_size, m_delegate_targets_size, init,
              segment_manager_o);
  }
}

template <typename SegmentManager>
void
delegate_partitioned_graph<SegmentManager>::
print_graph_statistics() {
  /*if(m_mpi_rank == 0) {
    for(size_t i=0; i<m_delegate_degree.size(); ++i) {
      std::cout << "Hub label = " << m_delegate_label[i] << ", degree = " <<
        m_delegate_degree[i] << std::endl << std::flush;
    }
  }*/

  uint64_t low_local_size = m_owned_targets_size;
  uint64_t high_local_size = m_delegate_targets_size;
  uint64_t total_local_size = low_local_size + high_local_size;

  //
  // Print out debugging info
  if (m_mpi_rank == 0) {
    std::cout << "Reducing.." << std::endl;
  }
  uint64_t low_max_size = mpi_all_reduce(low_local_size,
    std::greater<uint64_t>(), m_mpi_comm);
  uint64_t high_max_size = mpi_all_reduce(high_local_size,
    std::greater<uint64_t>(), m_mpi_comm);
  uint64_t total_max_size = mpi_all_reduce(total_local_size,
    std::greater<uint64_t>(), m_mpi_comm);
  uint64_t low_sum_size = mpi_all_reduce(low_local_size,
    std::plus<uint64_t>(), m_mpi_comm);
  uint64_t high_sum_size = mpi_all_reduce(high_local_size,
    std::plus<uint64_t>(), m_mpi_comm);
  uint64_t total_sum_size = mpi_all_reduce(total_local_size,
    std::plus<uint64_t>(), m_mpi_comm);

  if (m_mpi_rank == 0) {
    std::cout << "counting..(" << m_owned_targets_size << ")" << std::endl;
  }
  uint64_t local_count_del_target = 0;
  for (uint64_t i = 0; i < m_owned_targets_size; ++i) {
    if (m_owned_targets[i].is_delegate())
      ++local_count_del_target;
    if (m_mpi_rank == 0 && i % 10000000 == 0) {
      std::cout <<  i << "..." << std::flush;
    }
  }

  if (m_mpi_rank == 0) {
    std::cout << "reducing(2).." << std::endl;
  }

  uint64_t total_count_del_target = mpi_all_reduce(local_count_del_target,
    std::plus<uint64_t>(), m_mpi_comm);

  if (m_mpi_rank == 0) {
    std::cout
      << "========================================================" << std::endl
      << "Graph Statistics" << std::endl
      << "========================================================" << std::endl
      << "\tMax Local Vertex Id = " << m_max_vertex << std::endl
      << "\tHub vertices = " << m_map_delegate_locator.size() << std::endl
      << "\tTotal percentage good hub edges = " <<
      double(high_sum_size) / double(total_sum_size) * 100.0 << std::endl
      << "\tTotal count del target = " << total_count_del_target << std::endl
      << "\tTotal percentage of localized edges = " <<
          double(high_sum_size + total_count_del_target) /
          double(total_sum_size) * 100.0 << std::endl
      << "\tGlobal number of edges = " << total_sum_size << std::endl
      << "\tNumber of small degree = " << low_sum_size << std::endl
      << "\tNumber of hubs = " << high_sum_size << std::endl
      << "\toned imbalance = "
      << double(low_max_size) / double(low_sum_size/m_mpi_size) << std::endl
      << "\thubs imbalance = "
      << double(high_max_size) / double(high_sum_size/m_mpi_size) << std::endl
      << "\tTOTAL imbalance = "
      << double(total_max_size) / double(total_sum_size/m_mpi_size) << std::endl
      << "\tNode Partition = " << node_partitions << std::endl
      << "\tEdge Chunk Size = " << edge_chunk_size << std::endl
      << "\tProcesses_per_node = " << processes_per_node << std::endl
      << "\tDebuging = " << IS_DEBUGING << std::endl
      << "========================================================" << std::endl
      << "End Graph Statistics" << std::endl
      << "========================================================" << std::endl
      << std::flush;
  }
}
}  // namespace mpi
}  // namespace havoqgt

#endif  // HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
