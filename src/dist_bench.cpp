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





 /// --- option variables --- ///
 uint64_t vertex_scale_               = 18;
 uint64_t edge_factor_                = 16;
 uint64_t segmentfile_init_size_log2_ = 30;
 bool     is_delete_segmentfile_on_exit_ = false;
 uint64_t chunk_size_log10_           = 6;
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
         fname_segmentfile_ = optarg;
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
     if (fname_edge_list_.empty()) {
       std::cout << "Building RMAT graph Scale: " << vertex_scale_ << std::endl;
       std::cout << "Building RMAT graph Edge factor: " << edge_factor_ << std::endl;
     } else {
       for (const auto itr : fname_edge_list_)
         std::cout << mpi_rank << " : " << "Load edge list from " << itr << std::endl;
     }
   }

 }


template <typename gstore_type, typename edgelist_type>
void constract_graph(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                     gstore_type& gstore,
                     graphstore::utility::interprocess_mmap_manager& mmap_manager,
                     edgelist_type& edges,
                     const size_t chunk_size)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
  havoqgt::havoqgt_env()->world_comm().barrier();

  if (mpi_rank == 0) std::cout << "-- statuses of before generation --" << std::endl;
  if (mpi_rank == 0) print_system_mem_usages();
  havoqgt::havoqgt_env()->world_comm().barrier();
  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] segment size (GB) =\t"<< mmap_manager.segment_size_gb() << std::endl;
    }
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// --- variables for analysys --- //
  uint64_t loop_cnt = 0;
  bool global_is_finished = false;

  /// --- iterator and array for edgelist --- ///
  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  request_vector_type<vertex_id_type> update_request_vec;
  update_request_vec.reserve(chunk_size);

#if DEBUG_MODE
  std::stringstream fname;
  fname << fname_segmentfile_ << ".debug_edges_raw_" << mpi_rank;
  ofs_edges.open(fname.str());
#endif


  const double whole_start = MPI_Wtime();
  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n\n<< Loop no. " << loop_cnt << " >>" << std::endl;

    /// --- generate edges --- ///
    /// \brief generate_update_requests
    if (mpi_rank == 0) std::cout << "-- generate requests --" << std::endl;
    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);
    havoqgt::havoqgt_env()->world_comm().barrier();

    if (mpi_rank == 0) std::cout << "\n-- process requests --" << std::endl;
    double time_start = MPI_Wtime();

    /// --- queue requests; local construction --- ///
    for (auto request : update_request_vec) {
      auto edge = request.edge;

      auto src = std::get<0>(edge);
      auto dst = std::get<1>(edge);

      typename visitor_type<gstore_type>::vertex_locator vl_src(src);
      typename visitor_type<gstore_type>::vertex_locator vl_dst(dst);

      visitor_type<gstore_type> vistor(vl_src, vl_dst, visitor_type<gstore_type>::ADD);

      dg_visitor_queue.queue_visitor(vistor);
#if DEBUG_MODE
        ofs_edges << src << " " << dst << " 0" << "\n";
#endif
    }
    const double time_local_end = MPI_Wtime();

    /// --- global construction --- ///
    dg_visitor_queue.insert_edges();

    /// --- sync --- ///
    if (mpi_rank == 0) graphstore::utility::sync_files();
    const double time_sync_end = MPI_Wtime();

    /// --- print a progress report --- ///
    const double time_end = MPI_Wtime();
    havoqgt::havoqgt_env()->world_comm().barrier();
    if (mpi_rank == 0) {
      std::cout << "\n-- results --" << std::endl;
      std::cout <<" exec_time (local, sync),\t segment_size(GB)" << std::endl;
    }
    for (int i = 0; i < mpi_size; ++i) {
      if (i == mpi_rank) {
        std::cout << "prg [" << mpi_rank << "] "
                  << (time_end - time_start) << " ( "
                  << (time_local_end - time_start) << " , "
                  << (time_sync_end  - time_local_end) << " ),\t"
                  << mmap_manager.segment_size_gb() << std::endl;
      }
      havoqgt::havoqgt_env()->world_comm().barrier();
    }
    if (mpi_rank == 0) print_system_mem_usages();

#if VERBOSE
    if (mpi_rank == 0) std::cout << "\n-- graph store's status --" << std::endl;
    for (int i = 0; i < mpi_size; ++i) {
      if (i == mpi_rank) {
        std::cout << "[" << mpi_rank << "]" << std::endl;
        gstore.print_status(0);
      }
      havoqgt::havoqgt_env()->world_comm().barrier();
    }
#endif

    ++loop_cnt;

    /// --- Has everyone finished ? --- ///
    const bool local_is_finished = (edges_itr == edges_itr_end);
    MPI_Allreduce(&local_is_finished, &global_is_finished, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
    havoqgt::havoqgt_env()->world_comm().barrier();
  }
  havoqgt::havoqgt_env()->world_comm().barrier();
  const double whole_end = MPI_Wtime();

  if (mpi_rank == 0) {
    std::cout << "\n-- All edge updations done --" << std::endl;
    std::cout << "whole construction time : " << (whole_end - whole_start) << std::endl;
    print_system_mem_usages();
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// --- print summary information --- ///
  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] Usage: segment size (GiB) =\t"<< mmap_manager.segment_size_gb() << std::endl;
    }
    havoqgt::havoqgt_env()->world_comm().barrier();
  }

}


template <typename gstore_type>
void run_benchmark_rmat(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                        gstore_type& gstore,
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

  constract_graph(dg_visitor_queue, gstore, mmap_manager, rmat, chunk_size);
}


template <typename gstore_type>
void run_benchmark_edgefile(dg_visitor_queue_type<gstore_type>& dg_visitor_queue,
                            gstore_type& gstore,
                            graphstore::utility::interprocess_mmap_manager& mmap_manager,
                            std::vector<std::string>& fname_edge_list,
                            size_t chunk_size)
{
  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list);

  constract_graph(dg_visitor_queue, gstore, mmap_manager, edgelist, chunk_size);
}


template <typename gstore_type>
void dump_all_edges(gstore_type& gstore)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  std::cout << "dumping all elements for debug" << std::endl;
  std::stringstream ofname;
  ofname << fname_segmentfile_ << ".debug_edges_graph_" << mpi_rank;
  std::ofstream of(ofname.str());
  gstore.fprint_all_elements(of);
  of.close();
  std::cout << "done" << std::endl;
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    havoqgt::get_environment();
    havoqgt::havoqgt_env()->world_comm().barrier();

    parse_options(argc, argv);
    havoqgt::havoqgt_env()->world_comm().barrier();
    uint64_t num_vertices = (1ULL << vertex_scale_);


    /// --- init segment file --- ///
    size_t graph_capacity = std::pow(2, segmentfile_init_size_log2_);
    std::stringstream fname_local_segmentfile;
    fname_local_segmentfile << fname_segmentfile_ << "_" << mpi_rank;
    graphstore::utility::interprocess_mmap_manager::delete_file(fname_local_segmentfile.str());
    graphstore::utility::interprocess_mmap_manager mmap_manager(fname_local_segmentfile.str(), graph_capacity);
    havoqgt::havoqgt_env()->world_comm().barrier();
    if (mpi_rank == 0) print_system_mem_usages();


    /// --- allocate a graphstore and start a benchmark --- ///
    if (graphstore_name_ == "Baseline") {

      baseline_type gstore(mmap_manager.get_segment_manager());
      dist_gstore_type<baseline_type> dist_graph(&gstore);
      visitor_type<baseline_type>::set_graph_ref(&dist_graph);
      dg_visitor_queue_type<baseline_type> dg_visitor_queue(&dist_graph);


      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(dg_visitor_queue,
                           gstore,
                           mmap_manager,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(dg_visitor_queue,
                               gstore,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
#if DEBUG_MODE
      dump_all_edges(gstore);
#endif
    } else if (graphstore_name_ == "DegAwareRHH") {

      degawarerhh_type gstore(mmap_manager.get_segment_manager());
      dist_gstore_type<degawarerhh_type> dist_graph(&gstore);
      visitor_type<degawarerhh_type>::set_graph_ref(&dist_graph);
      dg_visitor_queue_type<degawarerhh_type> dg_visitor_queue(&dist_graph);


      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(dg_visitor_queue,
                           gstore,
                           mmap_manager,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(dg_visitor_queue,
                               gstore,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
#if DEBUG_MODE
      dump_all_edges(gstore);
#endif
    } else {
      std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
      exit(1);
    }
  }
  havoqgt::havoqgt_finalize();

  return 0;
}
