#ifndef PAGE_RANK_BENCH_HPP
#define PAGE_RANK_BENCH_HPP

template <typename graphstore_type, typename vertex_type>
extern void run_page_rank_sync(graphstore_type&, const std::size_t, const double, const int);

template <typename graphstore_type, typename vertex_type>
void run_page_rank(graphstore_type& graph, const std::size_t max_vertex_id, const std::size_t num_edges, const double damping_factor, const int num_loop)
{
  std::cout << "\n--- Page Rank ---" << std::endl;

  std::cout << "max_vertex_id:\t" << max_vertex_id << std::endl;
  std::cout << "num_edges:\t"     << num_edges     << std::endl;

  const auto tic = graphstore::utility::duration_time();
  run_page_rank_sync<graphstore_type, vertex_type>(graph, max_vertex_id + 1, damping_factor, num_loop);
  double duration_sec = graphstore::utility::duration_time_sec(tic);

  std::cout << "-----------" << std::endl;
  std::cout << "Page Rank time (sec.):\t"  << duration_sec << std::endl;
  std::cout << "Page Rank done." << std::endl;

}

#endif // PAGE_RANK_BENCH_HPP
