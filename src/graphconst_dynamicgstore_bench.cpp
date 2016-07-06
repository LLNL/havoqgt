/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */


#include <iostream>
#include <fstream>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/offset_ptr.hpp>
#include <boost/interprocess/containers/set.hpp>
#include <havoqgt/distributed_db.hpp>

/// must include the files below with this order ??
#include "dynamicgraphstore_bench.hpp"
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/baseline/baseline.hpp>
#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = bool;
using baseline_type       = graphstore::graphstore_baseline<vertex_id_type,
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
        fname_segmentfile_ = optarg;
        break;

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
        is_delete_segmentfile_on_exit_ = true; /// TODO: this option is not working
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


template <typename vertex_id_type, typename graphstore_type, typename edgelist_type>
void apply_edge_update_requests(graphstore_type& graph_store,
                                graphstore::utility::interprocess_mmap_manager& mmap_manager,
                                edgelist_type& edges,
                                const uint64_t chunk_size)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
  havoqgt::havoqgt_env()->world_comm().barrier();

  if (mpi_rank == 0) std::cout << "-- statuses of before generation --" << std::endl;
  if (mpi_rank == 0) print_system_mem_usages();
  havoqgt::havoqgt_env()->world_comm().barrier();
  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] segment size (GiB) =\t"<< mmap_manager.segment_size_gb() << std::endl;
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

  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n\n<< Loop no. " << loop_cnt << " >>" << std::endl;

    /// --- generate edges --- ///
    /// \brief generate_update_requests
    if (mpi_rank == 0) std::cout << "-- generate requests --" << std::endl;
    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);
    havoqgt::havoqgt_env()->world_comm().barrier();
    if (mpi_rank == 0) std::cout << "\n-- process requests --" << std::endl;

    /// --- update edges --- ///
    size_t count_inserted = 0;
    size_t count_deleted = 0;
    size_t count_insert_req = 0;
    size_t count_delete_req = 0;

    const double time_start = MPI_Wtime();
    const unsigned char dummy = 0;
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      if (request.is_delete) {
        count_deleted += graph_store.erase_edge(edge.first, edge.second);
        ++count_delete_req;
      } else {
        count_inserted += graph_store.insert_edge(edge.first, edge.second, dummy);
        ++count_insert_req;
      }
#if DEBUG_MODE
        ofs_edges << edge.first << " " << edge.second << " " << request.is_delete << "\n";
#endif
    }
    /// --- sync --- ///
    const double time_start_sync = MPI_Wtime();
    graphstore::utility::sync_files();
    const double time_end_sync = MPI_Wtime();

    /// this is a temp implementation
    if (loop_cnt % 10 == 0) {
      graph_store.opt();
    }

    /// --- print results --- ///
    const double time_end = MPI_Wtime();
    havoqgt::havoqgt_env()->world_comm().barrier();
    if (mpi_rank == 0) {
      std::cout << "\n-- results --" << std::endl;
      std::cout << " inserted_edges,\t"
                <<" deleted_edge,\t"
                <<" exec_time,\t"
                <<" segment_size(GiB)" << std::endl;
    }
    for (int i = 0; i < mpi_size; ++i) {
      if (i == mpi_rank) {
        std::cout << "prg [" << mpi_rank << "] "
                  << count_insert_req        << " ( " << count_inserted                    << " ),\t"
                  << count_delete_req        << " ( " << count_deleted                     << " ),\t"
                  << (time_end - time_start) << " ( " << (time_end_sync - time_start_sync) << " ),\t"
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
        graph_store.print_status(0);
      }
      havoqgt::havoqgt_env()->world_comm().barrier();
    }
#endif

    ++loop_cnt;

    /// --- Has everyone finished ? --- ///
    const bool local_is_finished = (edges_itr == edges_itr_end);
    MPI_Allreduce(&local_is_finished, &global_is_finished, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  }
  havoqgt::havoqgt_env()->world_comm().barrier();


  if (mpi_rank == 0) {
    std::cout << "\n-- All edge updations done --" << std::endl;
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

#if DEBUG_MODE
    {
      std::cout << "dumping all elements for debug" << std::endl;
      std::stringstream ofname;
      ofname << graphstore_name_ << ".debug_edges_graph";
      std::ofstream of(ofname.str());
      graph_store.fprint_all_elements(of);
      std::cout << "done" << std::endl;
    }
#endif

}


template<typename graphstore_type>
void run_benchmark_rmat(graphstore_type& graph_store,
                        graphstore::utility::interprocess_mmap_manager& mmap_manager,
                        size_t num_edges,
                        size_t chunk_size)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  uint64_t num_edges_per_rank = num_edges / mpi_size;
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
    vertex_scale_, num_edges_per_rank,
    0.57, 0.19, 0.19, 0.05, true, false);
  apply_edge_update_requests<vertex_id_type>(graph_store,
                                             mmap_manager,
                                             rmat,
                                             chunk_size);
}

template<typename graphstore_type>
void run_benchmark_edgefile(graphstore_type& graph_store,
                            graphstore::utility::interprocess_mmap_manager& mmap_manager,
                            std::vector<std::string>& fname_edge_list,
                            size_t chunk_size)
{
  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list);
  apply_edge_update_requests<vertex_id_type>(graph_store,
                                             mmap_manager,
                                             edgelist,
                                             chunk_size);
}


int main(int argc, char** argv)
{

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
      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(graph_store,
                           mmap_manager,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(graph_store,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
    } else if (graphstore_name_ == "DegAwareRHH") {
      degawarerhh_type graph_store(mmap_manager.get_segment_manager());
      if (fname_edge_list_.empty()) {
        run_benchmark_rmat(graph_store,
                           mmap_manager,
                           num_vertices * edge_factor_,
                           std::pow(10, chunk_size_log10_));
      } else {
        run_benchmark_edgefile(graph_store,
                               mmap_manager,
                               fname_edge_list_,
                               std::pow(10, chunk_size_log10_));
      }
    } else {
      std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
      exit(1);
    }

    havoqgt::havoqgt_env()->world_comm().barrier();
  }
  havoqgt::havoqgt_finalize();

  return 0;
}
