/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */
#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setfill, std::setw
#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include <functional>
#include <fcntl.h>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include <havoqgt/graphstore/graphstore_rhhda.hpp>
#include <havoqgt/graphstore/graph_traversal/bfs.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#include "dynamicgraphstore_bench.hpp"

#define VERBOSE 0

#define DEBUG_MODE 0
#if DEBUG_MODE
 std::ofstream ofs_edges;
#endif

 enum : size_t {
  midle_high_degree_threshold = 4
};

/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using vertex_meta_data_type = bool;
using edge_weight_type      = unsigned char;
using graphstore_type       = graphstore::graphstore_rhhda<
                                vertex_id_type,
                                vertex_meta_data_type,
                                edge_weight_type,
                                midle_high_degree_threshold>;

/// --- option variables --- ///
vertex_id_type max_vertex_id_ = 0;
size_t num_edges_ = 0;
std::string fname_segmentfile_;
std::vector<std::string> fname_edge_list_;
size_t segment_size_log2_ = 30;
std::vector<vertex_id_type> source_list_;


template <typename Edges>
void constract_graph(mapped_file_type& mapped_file,
                     segment_manager_type *const segment_manager,
                     graphstore_type& graph_store,
                     Edges& edges, const size_t chunk_size)
{
  std::cout << "-- Disp status of before generation --" << std::endl;
  print_usages(segment_manager);

  request_vector_type<vertex_id_type> update_request_vec = request_vector_type<vertex_id_type>();

  size_t count_inserted = 0;

  double construction_time = 0;
  size_t loop_cnt = 0;
  auto global_start = graphstore::utility::duration_time();
  for (auto edges_itr = edges.begin(), edges_itr_end = edges.end();
       edges_itr != edges_itr_end;
       ++edges_itr) {
    std::cout << "[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    update_request_vec.clear();
    generate_insertion_requests(edges_itr, edges_itr_end, chunk_size, update_request_vec, 0);

    unsigned char dummy_weight = 0;
    auto local_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      max_vertex_id_ = std::max(max_vertex_id_, edge.first);
      max_vertex_id_ = std::max(max_vertex_id_, edge.second);
      count_inserted += graph_store.insert_edge(edge.first, edge.second, dummy_weight);
    }
    std::cout << "shrink to fit low table" << std::endl;
    graph_store.shrink_to_fit_low_table();
    double t = graphstore::utility::duration_time_sec(local_start);
    construction_time += t;
    std::cout << "progress (sec.): " << t << std::endl;

    ++loop_cnt;
//    graph_store.print_all_elements_low();
//    graph_store.print_all_elements_mh();
  }
  std::cout << "sync mmap" << std::endl;
  flush_mmmap(mapped_file);
  std::cout << "sync di-mmap" << std::endl;
  sync_dimmap();
  const double whole_construction_time = graphstore::utility::duration_time_sec(global_start);

  num_edges_ = count_inserted;

  std::cout << "\n-- All edge updations done --" << std::endl;
  std::cout << "inserted edges : " << count_inserted << std::endl;
  std::cout << "construction time (insertion only) : " << construction_time << std::endl;
  std::cout << "whole construction time : " << whole_construction_time << std::endl;
  print_usages(segment_manager);

}


void run_bfs(graphstore_type& graph)
{
  std::cout << "\n--- BFS ---" << std::endl;

  std::cout << "max_vertex_id:\t" << max_vertex_id_ << std::endl;
  std::cout << "num_edges:\t"     << num_edges_     << std::endl;

  for (int i = 0; i < source_list_.size(); ++i) {
    std::cout << "BFS[" << i << "]: src=\t" << source_list_[i] << std::endl;

    graphstore::utility::print_time();
    bfs_sync<graphstore_type, vertex_id_type, 1>(graph, source_list_[i], max_vertex_id_, num_edges_);
    std::cout << "finish: ";
    graphstore::utility::print_time();
    std::cout << "\n" << std::endl;
  }
  std::cout << "BFS done." << std::endl;

}


void generate_source_list(const int num_sources)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "generate sources using a seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<vertex_id_type> dis(0, max_vertex_id_);
  for (size_t i = 0; i < num_sources; ++i) {
    source_list_.push_back(dis(gen));
  }
}


void parse_options(int argc, char **argv)
{

  char c;

  while ((c = getopt (argc, argv, "s:f:e:r:")) != -1) {
    switch (c) {
      case 's':
        segment_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'f':
        fname_segmentfile_ = optarg;
        break;

      case 'e':
      {
        std::string fname(optarg);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while (std::getline(fin, line)) {
          fname_edge_list_.push_back(line);
        }
        break;
      }

      case 'r':
      {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
            source_list_.push_back(boost::lexical_cast<size_t>(buf));
        break;
      }
    }
  }

}


int main(int argc, char** argv) {

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  parse_options(argc, argv);
  std::cout << "Midle-high degree threshold = " << midle_high_degree_threshold << std::endl;
  if (!fname_edge_list_.empty())
    std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
  for (auto itr : fname_edge_list_) {
    std::cout << "Load edge list from " << itr << std::endl;
  }


  /// --- create a segument file --- ///
  std::cout << "Delete segment file: " << fname_segmentfile_ << std::endl;
  boost::interprocess::file_mapping::remove(fname_segmentfile_.c_str());
  std::cout << "\n<<Construct segment>>" << std::endl;
  std::cout << "Create and map a segument file" << std::endl;
  uint64_t graph_capacity = std::pow(2, segment_size_log2_);
  mapped_file_type mapped_file = mapped_file_type(
                                   boost::interprocess::create_only,
                                   fname_segmentfile_.c_str(),
                                   graph_capacity);
  std::cout << "Call posix_fallocate\n";
  fallocate(fname_segmentfile_.c_str(), graph_capacity, mapped_file);

  /// --- Get a segument manager --- ///
  std::cout << "\n<Get a segment manager>" << std::endl;
  segment_manager_type* segment_manager = mapped_file.get_segment_manager();
  print_usages(segment_manager);


  /// --- Allocate graphstore_rhhda --- ///
  std::cout << "\n<Allocate graphstore_rhhda>" << std::endl;
  graphstore_type graph_store(segment_manager);

  /// --- Graph Construction --- ////
  std::cout << "\n<Construct graph>" << std::endl;
  havoqgt::havoqgt_init(&argc, &argv);
  {
    havoqgt::get_environment();

    havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_);
    constract_graph(
          mapped_file,
          segment_manager,
          graph_store,
          edgelist,
          static_cast<uint64_t>(std::pow(10, 6)));
  }


  /// ---------- Graph Traversal --------------- ///
  std::cout << "\n<Run BFS>" << std::endl;
  if (source_list_.empty())
    generate_source_list(4);
  run_bfs(graph_store);

  return 0;
}
