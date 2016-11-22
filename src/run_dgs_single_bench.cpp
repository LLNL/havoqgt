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
#include <havoqgt/distributed_db.hpp>

#include "dynamicgraphstore_bench.hpp" /// must include before the files below ??
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/baseline/baseline.hpp>
#include <havoqgt/graphstore/baseline/baseline_map.hpp>
#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = unsigned char;
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                             vertex_property_type,
                                                             edge_property_type,
                                                             havoqgt::distributed_db::segment_manager_type>;

using baselinemap_type      = graphstore::graphstore_baseline_map<vertex_id_type,
                                                                  vertex_property_type,
                                                                  edge_property_type,
                                                                  havoqgt::distributed_db::segment_manager_type>;

enum : size_t {
 middle_high_degree_threshold = 2 // must be more than 1
};
using degawarerhh_type  = graphstore::degawarerhh<vertex_id_type,
                                                     vertex_property_type,
                                                     edge_property_type,
                                                     havoqgt::distributed_db::segment_manager_type,
                                                     middle_high_degree_threshold>;




/// --- option variables --- ///
uint64_t vertex_scale_               = 18;
uint64_t edge_factor_                = 16;
uint64_t segmentfile_size_gb_        = 1;
uint64_t chunk_size_log10_           = 6;
std::string fname_segmentfile_       = "/dev/shm/segment_file";
std::string graphstore_name_         = "DegAwareRHH";
std::vector<std::string> fname_edge_list_;
bool is_edgelist_with_delete_       = false;

void usage()  {
 if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
   std::cerr << "Usage: \n"
        << " -s <int>      - the logarithm base two of the number of vertices\n"
        << " -e <int>      - edge factor\n"
        << " -o <string>   - base filename to create segmentfiles\n"
        << " -S <int>      - the size of segmentfile per MPI rank in GB\n"
        << " -g <string>   - the name of graph store (DegAwareRHH, Baseline or BaselineMap)\n"
        << " -c <int>      - the logarithm with base ten of the chunk size (defaults is 6, i.g., chunk size is 10^6)\n"
        << " -E <string>   - the name of the file which has a list of edgelist file's makes\n"
        << " -d            - edgelist files have delete operations\n"
        << " -h            - print help and exit\n\n";
 }
}

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
 while ((c = getopt (argc, argv, "s:e:dc:o:S:g:E:zfh")) != -1) {
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
       segmentfile_size_gb_ = boost::lexical_cast<size_t>(optarg);
       break;

     case 'g':
       graphstore_name_ = optarg;
       break;

     case 'c':
       chunk_size_log10_ = boost::lexical_cast<size_t>(optarg);
       break;

     case 'd':
       is_edgelist_with_delete_ = true;
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

     case 'h':
       usage();
       break;
   }
 }

 if (mpi_rank == 0) {
   std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
   std::cout << "Segment filse size in GB = " << segmentfile_size_gb_ << std::endl;
   std::cout << "Chunk size (log10) = " << chunk_size_log10_ << std::endl;
   if (fname_edge_list_.empty()) {
     std::cout << "Building RMAT graph Scale: " << vertex_scale_ << std::endl;
     std::cout << "Building RMAT graph Edge factor (generate both direction): " << edge_factor_ << std::endl;
   } else {
     for (const auto itr : fname_edge_list_)
       std::cout << mpi_rank << " : " << "Load edge list from " << itr << std::endl;
   }
 }

}


template <typename gstore_type, typename edgelist_type>
void constract_graph(gstore_type& gstore,
                     edgelist_type& edges,
                     const size_t chunk_size)
{
  std::cout << "-- statuses of before generation --" << std::endl;
  print_system_mem_usages();
  std::cout << "segment size (GB) =\t"<< graphstore::utility::segment_size_gb(gstore.get_segment_manager()) << std::endl;

  /// --- variables for analysys --- //
  uint64_t loop_cnt = 0;

  /// --- iterator and array for edgelist --- ///
  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  request_vector_type<vertex_id_type> update_request_vec;
  update_request_vec.reserve(chunk_size);

#if DEBUG_MODE
  std::stringstream fname;
  fname << fname_segmentfile_ << ".debug_edges_raw";
  ofs_edges.open(fname.str());
#endif

  double total_graph_storing_time = 0;
  const auto whole_start = graphstore::utility::duration_time();
  while (edges_itr != edges_itr_end) {
    std::cout << "\n\n<< Loop no. " << loop_cnt << " >>" << std::endl;

    /// --- generate or load edge update requests --- ///
    std::cout << "-- generate or load edge update requests --" << std::endl;
    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);

    std::cout << "\n-- process requests --" << std::endl;
    const auto time_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      if (request.is_delete)
        gstore.erase_edge(request.edge.first, request.edge.second);
      else
        gstore.insert_edge(request.edge.first, request.edge.second, edge_property_type());
    }
    graphstore::utility::sync_files();
    const double time_update = graphstore::utility::duration_time_sec(time_start);

    /// --- print a progress report --- ///
    total_graph_storing_time += time_update;

    std::cout << "\n-- results --" << std::endl;
    std::cout <<" exec_time,\t segment_size(GB)" << std::endl;
    std::cout << "prg "
              << time_update << ",\t"
              << graphstore::utility::segment_size_gb(gstore.get_segment_manager()) << std::endl;
    print_system_mem_usages();

#if VERBOSE
    std::cout << "\n-- graph store's status --" << std::endl;
    std::cout << "num edges =\t" << gstore.num_edges() << std::endl;
    gstore.print_status(0);
#endif

    ++loop_cnt;
  }
  const double time_whole = graphstore::utility::duration_time_sec(whole_start);

  std::cout << "\n-- All edge updations done --" << std::endl;
  std::cout << "whole execution time (sec) : " << time_whole << std::endl;
  std::cout << "only graph construction time (sec) (not inluding edge generation/load time) : " << (total_graph_storing_time) << std::endl;
  print_system_mem_usages();

  /// --- print summary information --- ///
  std::cout << "Usage: segment size (GiB) =\t"<< graphstore::utility::segment_size_gb(gstore.get_segment_manager()) << std::endl;
  std::cout << "num edges =\t" << gstore.num_edges() << std::endl;
  gstore.print_status(0);

}


template <typename gstore_type>
void run_benchmark_rmat(gstore_type& gstore,
                        const size_t vertex_scale,
                        const size_t num_edges,
                        const size_t chunk_size)
{
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(0) * 3ULL,
    vertex_scale, num_edges,
    0.57, 0.19, 0.19, 0.05, true, true);

  constract_graph(gstore, rmat, chunk_size);
}


template <typename gstore_type>
void run_benchmark_edgefile(gstore_type& gstore,
                            std::vector<std::string>& fname_edge_list,
                            const size_t chunk_size,
                            const bool is_edgelist_with_delete)
{
  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list, is_edgelist_with_delete);

  constract_graph(gstore, edgelist, chunk_size);
}


template <typename gstore_type>
void dump_all_edges(gstore_type& gstore)
{
  std::cout << "dumping all elements for debug" << std::endl;
  std::stringstream ofname;
  ofname << fname_segmentfile_ << ".debug_edges_graph";
  std::ofstream of(ofname.str());
  gstore.fprint_all_elements(of);
  of.close();
  std::cout << "done" << std::endl;
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    assert(havoqgt::havoqgt_env()->world_comm().size() == 1);
    havoqgt::get_environment();

    parse_options(argc, argv);
    uint64_t num_vertices = (1ULL << vertex_scale_);


    /// --- init segment file --- ///
    havoqgt::distributed_db ddb(havoqgt::db_create(), fname_segmentfile_.c_str(), segmentfile_size_gb_);

    print_system_mem_usages();

    /// --- allocate a graphstore and start a benchmark --- ///
    if (graphstore_name_ == "Baseline") {

      baseline_type gstore(ddb.get_segment_manager());

      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(gstore,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(gstore,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_),
                               is_edgelist_with_delete_);
      }
#if DEBUG_MODE
      dump_all_edges(gstore);
#endif
    } else if (graphstore_name_ == "BaselineMap") {

      baselinemap_type gstore(ddb.get_segment_manager());

      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(gstore,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(gstore,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_),
                               is_edgelist_with_delete_);
      }
#if DEBUG_MODE
      dump_all_edges(gstore);
#endif
    } else if (graphstore_name_ == "DegAwareRHH") {

      degawarerhh_type gstore(ddb.get_segment_manager());

      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(gstore,
                           vertex_scale_,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(gstore,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_),
                               is_edgelist_with_delete_);
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
