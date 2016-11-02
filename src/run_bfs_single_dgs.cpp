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
#include <havoqgt/distributed_db.hpp>

#include "dynamicgraphstore_bench.hpp" /// must include before the files below ??
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/baseline/baseline.hpp>
#include <havoqgt/graphstore/baseline/baseline_map.hpp>
#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>
#include "bfs_bench.hpp"


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = bool;

/// --- graph stores --- ///
using baseline_type = graphstore::graphstore_baseline<vertex_id_type,
                                                      vertex_property_type,
                                                      edge_property_type,
                                                      havoqgt::distributed_db::segment_manager_type>;

using baseline_map_type = graphstore::graphstore_baseline_map<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              havoqgt::distributed_db::segment_manager_type>;

 enum : size_t {
   middle_high_degree_threshold = 2 // must be more or equal than 1
 };
 using degawarerhh_type  = graphstore::degawarerhh<vertex_id_type,
                                                   vertex_property_type,
                                                   edge_property_type,
                                                   havoqgt::distributed_db::segment_manager_type,
                                                   middle_high_degree_threshold>;


template <typename graphstore_type, typename edgelist_type>
std::pair<vertex_id_type, size_t>
constract_graph(graphstore_type& graph_store,
                edgelist_type& edgelist,
                const size_t chunk_size)
{
  std::cout << "-- Disp a status of before construction --" << std::endl;
  std::cout << "segment size (GB); " << graphstore::utility::segment_size_gb(graph_store.get_segment_manager()) << std::endl;
  print_system_mem_usages();


  vertex_id_type max_vertex_id = 0;
  size_t num_edges = 0;

  size_t count_inserted = 0;
  size_t loop_cnt = 0;
  double construction_time = 0;

  auto edgelist_itr = edgelist.begin();
  auto edgelist_itr_end = edgelist.end();
  request_vector_type<vertex_id_type> update_request_vec = request_vector_type<vertex_id_type>();
  update_request_vec.reserve(chunk_size);

  auto global_start = graphstore::utility::duration_time();
  while (edgelist_itr != edgelist_itr_end) {
//    std::cout << "[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    generate_update_requests(edgelist_itr, edgelist_itr_end, update_request_vec, chunk_size);

    unsigned char dummy_weight = 0;
    auto local_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      max_vertex_id = std::max(max_vertex_id, edge.first);
      max_vertex_id = std::max(max_vertex_id, edge.second);
      count_inserted += graph_store.insert_edge(edge.first, edge.second, dummy_weight);
    }
    graphstore::utility::sync_files();
    double t = graphstore::utility::duration_time_sec(local_start);
    construction_time += t;
    std::cout << "progress (sec.): " << t << std::endl;

    ++loop_cnt;
  }
  const double whole_construction_time = graphstore::utility::duration_time_sec(global_start);
  assert(graph_store.num_edges() == count_inserted);

  std::cout << "\n-- All edge insertion done --" << std::endl;
  std::cout << "inserted edges : " << count_inserted << std::endl;
  std::cout << "construction time (insertion only) : " << construction_time << std::endl;
  std::cout << "whole construction time : " << whole_construction_time << std::endl;
  std::cout << "segment size (GB); " << graphstore::utility::segment_size_gb(graph_store.get_segment_manager()) << std::endl;
  print_system_mem_usages();

  return std::make_pair(max_vertex_id, num_edges);
}


template <typename graphstore_type, typename vertex_type>
  static void run_bfs_sync (graphstore_type& graphstore,
                            trv_inf<vertex_type>& inf,
                            std::queue<vertex_type>& frontier_queue,
                            std::queue<vertex_type>& next_queue,
                            vertex_type& root_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
#if BFS_USE_BITMAP
    inf.init(true);
    uint64_t* const visited = inf.visited;
    visited[bitmap_global_pos(root_vrtx)] |= 0x1 << bitmap_local_pos(root_vrtx);
#else
    inf.init(false);
    for (auto vert_itr = graphstore.vertices_begin(), end = graphstore.vertices_end();
         vert_itr != end;
         ++vert_itr) {
      vert_itr.property_data() = false;
    }
    graphstore.vertex_property_data(root_vrtx) = true;
#endif
    std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

    /// --- BFS main loop -------- ///
    size_t level = 0;
    while (true) {
      std::cout << "Lv. " << level << ", size of frontier queue =\t" << frontier_queue.size() << std::endl;

      /// loop for current frontier
      while (!frontier_queue.empty()) {
        vertex_type src = frontier_queue.front();
        frontier_queue.pop();
        ++inf.count_visited_vertices;

        size_t count_visited_edges = 0;
        /// push adjacent vertices to the next queue
        for (auto edge = graphstore.adjacent_edge_begin(src), end = graphstore.adjacent_edge_end(src);
             edge != end;
             ++edge) {
          const vertex_type& dst = edge.target_vertex();
#if BFS_USE_BITMAP
          const bool is_visited = visited[bitmap_global_pos(dst)] & (0x1 << bitmap_local_pos(dst));
#else
          const bool is_visited = graphstore.vertex_property_data(dst);
#endif
          if (!is_visited) {
            next_queue.push(dst);
#if BFS_USE_BITMAP
          visited[bitmap_global_pos(dst)] |= 0x1 << bitmap_local_pos(dst);
#else
          graphstore.vertex_property_data(dst) = true;
#endif
          }
          ++count_visited_edges;
        }
        inf.count_visited_edges += count_visited_edges;

      }  /// end of loop for a frontier

      if (next_queue.empty()) break; /// termination condition
      frontier_queue.swap(next_queue);
      ++level;
    } /// end of BFS loop
  }

/// Avoid linker errors with template function
template void run_bfs_sync<baseline_type, vertex_id_type>(baseline_type&,
                                                           trv_inf<vertex_id_type>&,
                                                           std::queue<vertex_id_type>&,
                                                           std::queue<vertex_id_type>&,
                                                           vertex_id_type&);

template void run_bfs_sync<degawarerhh_type, vertex_id_type>(degawarerhh_type&,
                                                         trv_inf<vertex_id_type>&,
                                                         std::queue<vertex_id_type>&,
                                                         std::queue<vertex_id_type>&,
                                                         vertex_id_type&);



/// --- option variables --- ///
std::string fname_segmentfile_ = "/dev/shm/graph";
std::vector<std::string> fname_edge_list_;
size_t segmentfile_init_size_gb_ = 30;
std::string graphstore_name_ = "DegAwareRHH";
std::vector<vertex_id_type> source_list_;

void parse_options(int argc, char **argv)
{

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  char c;

  while ((c = getopt (argc, argv, "S:o:E:r:g:")) != -1) {
    switch (c) {
      case 'S':
        segmentfile_init_size_gb_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'o':
        fname_segmentfile_ = optarg;
        break;

      case 'g':
        graphstore_name_ = optarg;
        break;

      case 'E':
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

  std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
  std::cout << "Initialize segment filse size (GB) = " << segmentfile_init_size_gb_ << std::endl;
  for (const auto itr : fname_edge_list_) {
    std::cout << "Load edge list from " << itr << std::endl;
  }

}

template<typename graphstore_type>
std::pair<vertex_id_type, size_t>
construct_graph(graphstore_type& graphstore,
                havoqgt::parallel_edge_list_reader& edgelist)
{
    return constract_graph(graphstore,
                           edgelist,
                           static_cast<size_t>(std::pow(10, 6)));
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
  havoqgt::get_environment();

  parse_options(argc, argv);


  /// --- init segment file --- ///
  uint64_t graph_capacity = segmentfile_init_size_gb_;
  havoqgt::distributed_db ddb(havoqgt::db_create(), fname_segmentfile_.c_str(), graph_capacity);
  print_system_mem_usages();


  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_);

  if (graphstore_name_ == "Baseline") {
    baseline_type graphstore(ddb.get_segment_manager());

    std::pair<vertex_id_type, size_t> ret = construct_graph(graphstore, edgelist);

    std::cout << "\n<Run BFS>" << std::endl;
    if (source_list_.empty())
      generate_bfs_sources(4, ret.first, source_list_);
    run_bfs(graphstore, ret.first, ret.second, source_list_);

  } else if (graphstore_name_ == "BaselineMap") {
    baseline_map_type graphstore(ddb.get_segment_manager());

    std::pair<vertex_id_type, size_t> ret = construct_graph(graphstore, edgelist);

    std::cout << "\n<Run BFS>" << std::endl;
    if (source_list_.empty())
      generate_bfs_sources(4, ret.first, source_list_);
    run_bfs(graphstore, ret.first, ret.second, source_list_);

  } else if (graphstore_name_ == "DegAwareRHH") {
    degawarerhh_type graphstore(ddb.get_segment_manager());

    std::pair<vertex_id_type, size_t> ret = construct_graph(graphstore, edgelist);

    std::cout << "\n<Run BFS>" << std::endl;
    if (source_list_.empty())
      generate_bfs_sources(4, ret.first, source_list_);
    run_bfs(graphstore, ret.first, ret.second, source_list_);

  } else {
    std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
    exit(1);
  }

  }
  havoqgt::havoqgt_finalize();

  return 0;
}
