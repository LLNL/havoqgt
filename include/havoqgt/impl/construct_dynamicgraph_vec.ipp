#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED

/**
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */


#ifdef PROFILE_DETAIL
 #warning PROFILE_DETAIL is enabled.
#endif


namespace havoqgt {
namespace mpi {

  IOInfo::IOInfo()
  : read_total_mb_(0.0), written_total_mb_(0.0)
  {
    init();
  }

  void IOInfo::init() {
    read_total_mb_ = 0.0;
    written_total_mb_ = 0.0;
    get_status(read_init_mb_, written_init_mb_);
  }

  void IOInfo::get_status(int &r, int &w) {
    FILE *pipe;
    char str[250];
    pipe = popen("iostat -m | grep md0 2>&1", "r" );

    float temp;
    fscanf(pipe, "%s", str);
    fscanf(pipe, "%f %f %f", &temp, &temp, &temp);
    fscanf(pipe, "%d %d \n", &r, &w);
    pclose(pipe);
  };

  void IOInfo::log_diff(bool final = false) {
    int read, written;
    get_status(read, written);
    read    -= read_init_mb_;
    written -= written_init_mb_;

    read_total_mb_    += read;
    written_total_mb_ += written;

    std::cout << "MB Read:\t"     << read    << std::endl;
    std::cout << "MB Written:\t"  << written << std::endl;
    if (final) {
      std::cout << "Total MB Read:\t"    << read_total_mb_     << std::endl;
      std::cout << "Total MB Written:\t" << written_total_mb_  << std::endl;    
    }

  };



/**
 * Constructor : sets secment allocator
 *
 * @param seg_allocator       Reference to segment allocator
 */
template <typename SegementManager>
construct_dynamicgraph_vec<SegementManager>::
construct_dynamicgraph_vec(const SegmentAllocator<void>& seg_allocator, const int mode)
  : adjacency_matrix_vec_vec(seg_allocator)
  , adjacency_matrix_map_vec(seg_allocator)
  , data_structure_type(mode)
  , io_info() { }


/**
 * Add edges into vector-vector adjacency matrix using
 * boost:interprocess:vector with from and unsorted sequence of edges.
 *
 * @param edges               input edges to partition
 * @param max_vertex          Max vertex ID
*/
template <typename SegementManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph_vec<SegementManager>::
add_edges_adjacency_matrix_vector_vector(
  ManagedMappedFile& asdf,
  const SegmentAllocator<void>& seg_allocator,
  Container& edges)
{
  // XXX: make initializer
  uint64_vector_t init_vec(seg_allocator);
  if (adjacency_matrix_vec_vec.size() == 0)
      adjacency_matrix_vec_vec.resize(1, init_vec);
  // IOInfo *io_info = new IOInfo();

  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;

    while (adjacency_matrix_vec_vec.size() <= edge.first) {
      adjacency_matrix_vec_vec.resize(adjacency_matrix_vec_vec.size() * 2, init_vec);
    }
    
    adjacency_matrix_vec_vec[edge.first].push_back(edge.second);
  } 
  asdf.flush();
  double time_end = MPI_Wtime();

  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;  
  io_info.log_diff();

  //free_edge_container(edges);
}


template <typename SegmentManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph_vec<SegmentManager>::
add_edges_adjacency_matrix_map_vector(
  ManagedMappedFile& asdf, 
  const SegmentAllocator<void>& seg_allocator, 
  Container& edges)
{

  uint64_vector_t init_vec(seg_allocator);
  //IOInfo *io_info = new IOInfo();

  double time_start = MPI_Wtime();
  for (auto itr = edges.begin(); itr != edges.end(); itr++) {
    const auto edge = *itr;
    auto value = adjacency_matrix_map_vec.find(edge.first);
    if (value == adjacency_matrix_map_vec.end()) {
      uint64_vector_t vec(1, edge.second, seg_allocator);
      adjacency_matrix_map_vec.insert(map_value_t(edge.first, vec));
    } else {
      uint64_vector_t *adjacency_list_vec = &value->second;
      adjacency_list_vec->push_back(edge.second);
    }
  }  
  asdf.flush();
  double time_end = MPI_Wtime();  


  std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;
  io_info.log_diff();

}

template <typename SegementManager>
const int construct_dynamicgraph_vec<SegementManager>::USE_VEC_VEC_MATRIX = 1;
template <typename SegementManager>
const int construct_dynamicgraph_vec<SegementManager>::USE_MAP_VEC_MATRIX = 2;


} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
