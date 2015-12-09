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

#include <havoqgt/mpi.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include "dynamicgraphstore_bench.hpp" /// must include before the files below ??
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/baseline.hpp>
#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>
#include <havoqgt/graphstore/dist_dynamic_graphstore.hpp>


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = unsigned char;
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                             vertex_property_type,
                                                             edge_property_type,
                                                             segment_manager_type>;

 enum : size_t {
   middle_high_degree_threshold = 2 // must be more or equal than 1
 };
 using degawarerhh_type  = graphstore::degawarerhh<vertex_id_type,
                                                       vertex_property_type,
                                                       edge_property_type,
                                                       segment_manager_type,
                                                       middle_high_degree_threshold>;


template <typename gstore_type>
using dist_gstore_type = dist_dynamic_graphstore<gstore_type>;

template <typename gstore_type>
using visitor_type = dg_visitor<dist_gstore_type<gstore_type>>;

template <typename gstore_type>
using dg_visitor_queue_type = havoqgt::mpi::visitor_queue<visitor_type<gstore_type>,
                                                          havoqgt::detail::visitor_priority_queue,
                                                          dist_gstore_type<gstore_type>>;


template <typename gstore_type, typename edgelist_type>
void constract_graph(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                graphstore::utility::interprocess_mmap_manager& mmap_manager,
                edgelist_type& edgelist,
                const size_t chunk_size)
{
  std::cout << "-- Disp status of before generation --" << std::endl;
  std::cout << "segment size (GB); " << mmap_manager.segment_size_gb() << std::endl;
  print_system_mem_usages();


  size_t loop_cnt = 0;
  double construction_time = 0;

  auto edgelist_itr = edgelist.begin();
  auto edgelist_itr_end = edgelist.end();
  request_vector_type<vertex_id_type> update_request_vec = request_vector_type<vertex_id_type>();
  update_request_vec.reserve(chunk_size);

  auto global_start = graphstore::utility::duration_time();
  while (edgelist_itr != edgelist_itr_end) {
    std::cout << "[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    generate_update_requests(edgelist_itr, edgelist_itr_end, update_request_vec, chunk_size);

    auto local_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      auto edge = request.edge;

      auto src = std::get<0>(edge);
      auto dst = std::get<1>(edge);

      typename visitor_type<gstore_type>::vertex_locator vl_src(src);
      typename visitor_type<gstore_type>::vertex_locator vl_dst(dst);

      visitor_type<gstore_type> vistor(vl_src, vl_dst, visitor_type<gstore_type>::ADD);

      dg_visitor_queue.queue_visitor(vistor);
    }
    double tp = graphstore::utility::duration_time_sec(local_start);
    std::cout << "start global construction (sec.): " << tp << std::endl;
    dg_visitor_queue.insert_edges();
    graphstore::utility::sync_files();
    double t = graphstore::utility::duration_time_sec(local_start);
    construction_time += t;
    std::cout << "progress (sec.): " << t << std::endl;

    ++loop_cnt;
  }
  const double whole_construction_time = graphstore::utility::duration_time_sec(global_start);

  std::cout << "\n-- All edge updations done --" << std::endl;
  std::cout << "construction time (insertion only) : " << construction_time << std::endl;
  std::cout << "whole construction time : " << whole_construction_time << std::endl;
  std::cout << "segment size (GB); " << mmap_manager.segment_size_gb() << std::endl;
  print_system_mem_usages();

}



/// --- option variables --- ///
uint64_t vertex_scale_               = 18;
uint64_t edge_factor_                = 16;
uint64_t segmentfile_init_size_log2_ = 30;
bool     is_delete_segmentfile_on_exit_ = false;
uint64_t chunk_size_log10_           = 6;
int edges_delete_ratio_              = 0;
std::string fname_segmentfile_;
std::string graphstore_name_;
std::vector<std::string> fname_edge_list_;


void parse_options(int argc, char **argv)
{

  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  char c;
  while ((c = getopt (argc, argv, "s:e:dc:o:S:g:E:")) != -1) {
    switch (c) {
      case 's':
        vertex_scale_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'e':
        edge_factor_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'o':
      {
        std::stringstream fname_local_segmentfile;
        fname_local_segmentfile << optarg << "_" << mpi_rank;
        fname_segmentfile_ = fname_local_segmentfile.str();
        break;
      }

      case 'S':
        segmentfile_init_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'g':
        graphstore_name_ = optarg;
        break;

      case 'c':
        chunk_size_log10_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'd':
        is_delete_segmentfile_on_exit_ = true;
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
          /// Note: parallel_edge_list_reader will assign files to multiple process
          fname_edge_list_.push_back(line);
        }
        break;
      }

    }
  }

  if (mpi_rank == 0) {
    std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
    std::cout << "Initialize segment filse size (log2) = " << segmentfile_init_size_log2_ << std::endl;
    std::cout << "Delete on Exit = " << is_delete_segmentfile_on_exit_ << std::endl;
    std::cout << "Chunk size (log10) = " << chunk_size_log10_ << std::endl;
    std::cout << "Edges Delete Ratio = " << edges_delete_ratio_ << std::endl;
    if (fname_edge_list_.empty()) {
      std::cout << "Building RMAT graph Scale: " << vertex_scale_ << std::endl;
      std::cout << "Building RMAT graph Edge factor: " << edge_factor_ << std::endl;
    } else {
      for (const auto itr : fname_edge_list_)
        std::cout << mpi_rank << " : " << "Load edge list from " << itr << std::endl;
    }

#if DEBUG_MODE
    {
    if (mpi_size > 1) {
        assert(false);
      }
      std::stringstream fname;
      fname << fname_segmentfile_ << ".debug_edges_raw";
      ofs_edges.open(fname.str());
    }
#endif
  }

}


template <typename gstore_type>
void run_benchmark_rmat(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                        graphstore::utility::interprocess_mmap_manager& mmap_manager,
                        size_t vertex_scale,
                        size_t num_edges,
                        size_t chunk_size)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  uint64_t num_edges_per_rank = num_edges / mpi_size;
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
    vertex_scale, num_edges_per_rank,
    0.57, 0.19, 0.19, 0.05, true, false);

  constract_graph(dg_visitor_queue, mmap_manager, rmat, chunk_size);
}


template <typename gstore_type>
void run_benchmark_edgefile(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                            graphstore::utility::interprocess_mmap_manager& mmap_manager,
                            std::vector<std::string>& fname_edge_list,
                            size_t chunk_size)
{
  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list);

  constract_graph(dg_visitor_queue, mmap_manager, edgelist, chunk_size);
}



int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
    havoqgt::get_environment();
    havoqgt::havoqgt_env()->world_comm().barrier();

    parse_options(argc, argv);
    uint64_t num_vertices = (1ULL << vertex_scale_);
    havoqgt::havoqgt_env()->world_comm().barrier();


    /// --- init segment file --- ///
    uint64_t graph_capacity = std::pow(2, segmentfile_init_size_log2_);
    graphstore::utility::interprocess_mmap_manager::delete_file(fname_segmentfile_);
    graphstore::utility::interprocess_mmap_manager mmap_manager(fname_segmentfile_, graph_capacity);
    print_system_mem_usages();


    /// --- allocate a graphstore and start a benchmark --- ///
    if (graphstore_name_ == "Baseline") {

      baseline_type graph_store(mmap_manager.get_segment_manager());
      dist_gstore_type<baseline_type> dist_graph(&graph_store);
      visitor_type<baseline_type>::set_graph_ref(&dist_graph);
      dg_visitor_queue_type<baseline_type> dg_visitor_queue(&dist_graph);


      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(dg_visitor_queue,
                           mmap_manager,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(dg_visitor_queue,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
    } else if (graphstore_name_ == "DegAwareRHH") {

      degawarerhh_type graph_store(mmap_manager.get_segment_manager());
      dist_gstore_type<degawarerhh_type> dist_graph(&graph_store);
      visitor_type<degawarerhh_type>::set_graph_ref(&dist_graph);
      dg_visitor_queue_type<degawarerhh_type> dg_visitor_queue(&dist_graph);


      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(dg_visitor_queue,
                           mmap_manager,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(dg_visitor_queue,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
    } else {
      std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
      exit(1);
    }

  }
  havoqgt::havoqgt_finalize();

  return 0;
}
