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
#include "page_rank_bench.hpp"


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = graphstore::utility::packed_pair<double, double>;

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

  return std::make_pair(max_vertex_id, count_inserted);
}

template <typename graphstore_type>
void print_status(graphstore_type& graphstore)
{
  for (auto vrt_itr = graphstore.vertices_begin(), end = graphstore.vertices_end();
       vrt_itr != end;
       ++vrt_itr) {
    std::cout << vrt_itr.property_data().first << " ";
  }
  std::cout << std::endl;
}

template <typename graphstore_type, typename vertex_type>
static void run_page_rank_sync (graphstore_type& graphstore, const std::size_t num_vertices, const double damping_factor, const int num_loops)
{

  /// ---- init ---- ///
  auto tic_init = graphstore::utility::duration_time();
  for (auto vrt_itr = graphstore.vertices_begin(), end = graphstore.vertices_end();
       vrt_itr != end;
       ++vrt_itr) {
    vrt_itr.property_data() = vertex_property_type(1.0 / num_vertices, 0.0);
  }
  std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;
//  print_status(graphstore);

  for (int k = 0; k < num_loops; ++k) {
    auto tic_step = graphstore::utility::duration_time();
    for (auto vrt_itr = graphstore.vertices_begin(), vrt_end = graphstore.vertices_end();
         vrt_itr != vrt_end;
         ++vrt_itr) {
      const double pr = vrt_itr.property_data().first;
      const vertex_type& v = vrt_itr.source_vertex();
      const size_t degree = graphstore.degree(v);
      for (auto adj_itr = graphstore.adjacent_edge_begin(v), adj_end = graphstore.adjacent_edge_end(v);
           adj_itr != adj_end;
           ++adj_itr) {
        graphstore.vertex_property_data(adj_itr.target_vertex()).second += (pr / degree) * damping_factor;
      }
    }
    for (auto vrt_itr = graphstore.vertices_begin(), vrt_end = graphstore.vertices_end();
         vrt_itr != vrt_end;
         ++vrt_itr) {
      vrt_itr.property_data().first = vrt_itr.property_data().second + (1.0 - damping_factor) / num_vertices;
      vrt_itr.property_data().second = 0.0;
    }
    std::cout << "\n[" << k << "] : progress (sec.):\t"  << graphstore::utility::duration_time_sec(tic_step) << std::endl;
//    print_status(graphstore);
  }
}

/// Avoid linker errors with template function
template void run_page_rank_sync<baseline_type, vertex_id_type>(baseline_type&, const std::size_t, const double, const int);
template void run_page_rank_sync<baseline_map_type, vertex_id_type>(baseline_map_type&, const std::size_t, const double, const int);
template void run_page_rank_sync<degawarerhh_type, vertex_id_type>(degawarerhh_type&, const std::size_t, const double, const int);



/// --- option variables --- ///
std::string fname_segmentfile_ = "/dev/shm/graph";
std::vector<std::string> fname_edge_list_;
size_t segmentfile_init_size_gb_ = 30;
std::string graphstore_name_ = "DegAwareRHH";
int num_loop_ = 10;
double damping_factor_ = 1.0;

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

      case 'd':
        damping_factor_ = boost::lexical_cast<double>(optarg);
        break;

      case 'r':
        num_loop_ = boost::lexical_cast<int>(optarg);
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
    run_page_rank<baseline_type, vertex_id_type>(graphstore, ret.first, ret.second, damping_factor_, num_loop_);

  } else if (graphstore_name_ == "BaselineMap") {
    baseline_map_type graphstore(ddb.get_segment_manager());
    std::pair<vertex_id_type, size_t> ret = construct_graph(graphstore, edgelist);
    run_page_rank<baseline_map_type, vertex_id_type>(graphstore, ret.first, ret.second, damping_factor_, num_loop_);

  } else if (graphstore_name_ == "DegAwareRHH") {
    degawarerhh_type graphstore(ddb.get_segment_manager());
    std::pair<vertex_id_type, size_t> ret = construct_graph(graphstore, edgelist);
    run_page_rank<degawarerhh_type, vertex_id_type>(graphstore, ret.first, ret.second, damping_factor_, num_loop_);

  } else {
    std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
    exit(1);
  }

  }
  havoqgt::havoqgt_finalize();

  return 0;
}
