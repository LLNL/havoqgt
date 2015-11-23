#ifndef BFS_BENCH_HPP
#define BFS_BENCH_HPP

#include<vector>

#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/graphstore/csr/csr_graph.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#define BFS_USE_BITMAP 0


template<typename vertex_type>
struct trv_inf {
  const size_t num_vertices;
  const size_t num_edges;
  bool* is_visited;
  vertex_type* tree;
  size_t count_visited_vertices;
  size_t count_visited_edges;

  trv_inf(size_t max_vertex_id, size_t n_edges)
    : num_vertices(max_vertex_id + 1),
      num_edges(n_edges),
      is_visited(nullptr),
      tree(nullptr),
      count_visited_vertices(0),
      count_visited_edges(0)
  {
  }

  ~trv_inf()
  {
    delete[] is_visited;
    delete[] tree;
  }

  void init(const bool is_allocate_maps)
  {
    if (is_allocate_maps) {
      is_visited = new bool[num_vertices];
      std::cout << "Allocate is_visited:\t" <<  (double)num_vertices * sizeof(bool) / (1ULL<<30) << " GB" << std::endl;

      tree = new vertex_type[num_vertices];
      std::cout << "Allocate tree:\t" <<  (double)num_vertices * sizeof(vertex_type) / (1ULL<<30) << " GB" << std::endl;

      for (size_t i = 0 ; i < num_vertices; ++i) {
        is_visited[i] = false;
        tree[i] = std::numeric_limits<vertex_type>::max();
      }
    }
    count_visited_vertices = 0;
    count_visited_edges = 0;
  }


};


template <typename graphstore_type, typename vertex_type>
extern void run_bfs_sync(graphstore_type&,
                         trv_inf<vertex_type>&,
                         std::queue<vertex_type>&,
                         std::queue<vertex_type>&,
                         vertex_type&);

template <typename graphstore_type, typename vertex_type>
void run_bfs(graphstore_type& graph, size_t max_vertex_id, size_t num_edges, std::vector<vertex_type>& source_list)
{
  std::cout << "\n--- BFS ---" << std::endl;

  std::cout << "max_vertex_id:\t" << max_vertex_id << std::endl;
  std::cout << "num_edges:\t"     << num_edges     << std::endl;

  auto tic_init = graphstore::utility::duration_time();
  trv_inf<vertex_type> inf(max_vertex_id, num_edges);
  std::queue<vertex_type> frontier_queue;
  std::queue<vertex_type> next_queue;
  std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;


  for (int i = 0; i < source_list.size(); ++i) {
    std::cout << "BFS[" << i << "]: src=\t" << source_list[i] << std::endl;
    frontier_queue.push(source_list[i]);

    const auto tic = graphstore::utility::duration_time();
    run_bfs_sync<graphstore_type, vertex_type>(graph, inf, frontier_queue, next_queue, source_list[i]);
    double duration_sec = graphstore::utility::duration_time_sec(tic);

    std::cout << "-----------" << std::endl;
    std::cout << "#visited vertices:\t" << inf.count_visited_vertices << std::endl;
    std::cout << "#visited edges:\t"    << inf.count_visited_edges    << std::endl;
    std::cout << "BFS time (sec.):\t"  << duration_sec << std::endl;
    std::cout << "MegaTEPS:\t" << num_edges / duration_sec / (1ULL << 20) << std::endl;
  }
  std::cout << "BFS done." << std::endl;

}

template<typename vertex_type>
void generate_bfs_sources(const int num_sources,  const vertex_type max_vertex_id, std::vector<vertex_type>& source_list)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "generate sources using a seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<vertex_type> dis(0, max_vertex_id);
  for (size_t i = 0; i < num_sources; ++i) {
    source_list.push_back(dis(gen));
  }
}

template <typename vertex_id_type, typename graphstore_type, typename edgelist_type>
void construct_graph(segment_manager_type *const segment_manager,
                     graphstore_type& graph_store,
                     edgelist_type& edges,
                     const uint64_t chunk_size)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
  havoqgt::havoqgt_env()->world_comm().barrier();

  uint64_t loop_cnt = 0;
  bool global_is_finished = false;

  /// --- iterator and array for edgelist --- ///
  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  request_vector_type<vertex_id_type> update_request_vec;
  update_request_vec.reserve(chunk_size);

  const double time_start = MPI_Wtime();
  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n\n<< Loop no. " << loop_cnt << " >>" << std::endl;

    /// --- generate edges --- ///
    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);

    const unsigned char dummy = 0;
    for (auto request : update_request_vec) {
      graph_store.insert_edge(request.edge.first, request.edge.second, dummy);
    }
    if (mpi_rank == 0) sync_mmap();
    havoqgt::havoqgt_env()->world_comm().barrier();

    /// --- Has everyone finished ? --- ///
    const bool local_is_finished = (edges_itr == edges_itr_end);
    MPI_Allreduce(&local_is_finished, &global_is_finished, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    ++loop_cnt;
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// --- print results --- ///
  if (mpi_rank == 0) {
    const double time_end = MPI_Wtime();
    std::cout << "\n-- All edge updations done --" << std::endl;
    std::cout << "execution time (sec.) =\t" << time_end - time_start << std::endl;
    print_system_mem_usages();
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// --- print summary information --- ///
  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] Usage: segment size (GiB) =\t"<< get_segment_size(segment_manager) << std::endl;
    }
    havoqgt::havoqgt_env()->world_comm().barrier();
  }
}

#endif // BFS_BENCH_HPP
