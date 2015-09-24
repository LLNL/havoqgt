#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setfill, std::setw
#include <random>
#include <chrono>

#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/graphstore/csr/csr_graph.hpp>
#include <havoqgt/graphstore/graph_traversal/bfs.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

enum {
  kNumBFSLoop = 64
};

using index_type = uint64_t;
using vertex_type = uint64_t;
using graph_type = csr_graph_struct::csr_graph<havoqgt::parallel_edge_list_reader, index_type, vertex_type>;

template <typename gen_type, typename rnd_type>
void run_bfs(graph_type& graph, size_t max_vertex_id, size_t num_edges, gen_type& gen, rnd_type& dis)
{
  std::cout << "\n--- BFS ---" << std::endl;

  std::cout << "max_vertex_id:\t" << max_vertex_id << std::endl;
  std::cout << "num_edges:\t"     << num_edges     << std::endl;

#if USE_SRC_CANDIDATE_TBL
  bool* table = new bool[graph->num_vertices()];
  std::cout << "Allocated src_candidats_table: " << sizeof(bool) * graph->num_vertices() / (1ULL << 30) << " GB" << std::endl;
  find_startvertex_candidates(graph, table);
#endif

  for (int i = 0; i < kNumBFSLoop; ++i) {
    vertex_type src = dis(gen);
    std::cout << "BFS[" << i << "]: src=\t" << src << std::endl;

    graphstore::utility::print_time();
    bfs_sync<graph_type, vertex_type, false>(graph, src, max_vertex_id, num_edges);
    std::cout << "finish: ";
    graphstore::utility::print_time();

    std::cout << "\n" << std::endl;
  }
  std::cout << "BFS done." << std::endl;

}


int main(int argc, char* argv[])
{
  std::string fname_prefix;
  size_t num_files = 0;
  size_t max_vertex_id = 0;
  size_t num_edges = 0;
  bool is_load_grah_from_file;

  if (argc == 3 || argc == 6) {
    is_load_grah_from_file = static_cast<bool>(std::atoi(argv[1]));
    fname_prefix = argv[2];
    if (!is_load_grah_from_file) {
      num_files = std::atoll(argv[3]);
      max_vertex_id = std::atoll(argv[4]);
      num_edges = std::atoll(argv[5]);
    }
  } else {
    std::cout << "Invalid argments: " << std::endl;
    std::cout << "./a.out is_load_grah_from_file (True) fname_prefix" << std::endl;
    std::cout << "./a.out is_load_grah_from_file (False) fname_prefix num_files [max_vertex_id] [num_edges]" << std::endl;
    std::exit(1);
  }

  graph_type* graph;

  if (is_load_grah_from_file) {
    std::cout << "--- Constructing a graph from file ---" << std::endl;
    graph = new graph_type(fname_prefix);
    max_vertex_id = graph->num_vertices() - 1;
    num_edges = graph->num_edges();

  } else {
    std::vector<std::string> edgelis_fname_list;
    for (size_t i = 0; i < num_files; ++i) {
      std::stringstream fname;
      fname << fname_prefix << std::setfill('0') << std::setw(5) << i;
      std::cout << fname.str().c_str() << std::endl;
      edgelis_fname_list.push_back(fname.str().c_str());
    }

    std::cout << "\n--- Initializing edgelist ---" << std::endl;
    graphstore::utility::print_time();
    havoqgt::parallel_edge_list_reader* edge_list;
    edge_list = new havoqgt::parallel_edge_list_reader(edgelis_fname_list);
//    std::cout << "max_vertex_id:\t" << max_vertex_id << "" << std::endl;
//    std::cout << "#edges:\t" << num_edges << "" << std::endl;

    std::cout << "\n--- Constructing csr graph ---" << std::endl;
    graphstore::utility::print_time();
    graph = new graph_type(*edge_list, max_vertex_id, num_edges);

    delete edge_list;
  }

  /// ---------- Graph Traversal --------------- ///

  // obtain a seed from the system clock:
  //  std::random_device rd; /// not implemented in gccc ?
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << " use seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<vertex_type> dis(0, max_vertex_id);

  run_bfs(*graph, max_vertex_id, num_edges, gen, dis);

  delete graph;
}
