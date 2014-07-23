#ifndef HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
#define HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED

/**
 * \file
 * Implementation of delegate_partitioned_graph and internal classes.
 */


#ifdef DEBUG
 #warning Debug MACRO is enabled.
#endif


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
                "\nTotal MB Written: " << (written - mb_written_) << std::endl;
    } else {
      std::cout << "\tMB Read: " << (read - mb_read_) <<
                "\n\tMB Written: " << (written - mb_written_) << std::endl;
    }
  };

 private:
  int mb_read_;
  int mb_written_;
};

/**
 * Construct a graph using boost:interprocess:vector with from and unsorted sequence of edges.
 *
 * @param edges               input edges to partition
 * @param max_vertex          Max vertex ID
*/
template <typename SegementManager>
template <typename ManagedMappedFile, typename Container>
void construct_dynamicgraph_vec<SegementManager>::
construct_adjacency_matrix_vector(ManagedMappedFile& asdf, Container& edges, uint64_t max_vertex)
{
  std::cout << "-- construct_adjacency_matrix_vector --" << std::endl;
  std::cout << "# of verticies = " << max_vertex << std::endl;
  std::cout << "# of edges = "     << edges.size() << std::endl;

  IOInfo *io_info_temp = new IOInfo();
  assert (adjacency_matrix.size() < max_vertex);
  double time_start = MPI_Wtime();
  for (auto it = edges.begin(); it != edges.end(); it++) {
    const auto edge = *it;  
    adjacency_matrix[edge.first].push_back(edge.second);
  }  
  asdf.flush();
  double time_end = MPI_Wtime();
  std::cout << "Execution time (sec.) = " << time_end - time_start << std::endl;
  
  io_info_temp->log_diff();

  free_edge_container(edges);

  std::cout << "-----------------------------------" << std::endl;
}

/**
 * Constructor : resize data space
 *
 * @param seg_allocator       Reference to segment allocator
 * @param max_vertex          Max vertex ID
*/

template <typename SegementManager>
construct_dynamicgraph_vec<SegementManager>::
construct_dynamicgraph_vec(const SegmentAllocator<void>& seg_allocator, uint64_t max_vertex )
  : adjacency_matrix(seg_allocator) 
{
  //IOInfo *io_info_temp = new IOInfo();

  //std::cout << "Resize adjacency_matrix: " << std::endl;
  uint64_vector_t ini_vec(seg_allocator);
  adjacency_matrix.resize(max_vertex+1, ini_vec);
  
  seg_allocator.get_segment_manager();

  //std::cout << "Max Vertex Id = " << max_vertex << std::endl;
  //io_info_temp->log_diff();
}


} // namespace mpi
} // namespace havoqgt


#endif //HAVOQGT_MPI_IMPL_DELEGATE_PARTITIONED_GRAPH_IPP_INCLUDED
