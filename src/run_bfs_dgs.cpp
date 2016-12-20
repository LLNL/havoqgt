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
#include <havoqgt/graphstore/dist_dynamic_graphstore.hpp>

/// --- visitor class for BFS --- ///
template <typename graphstore_type>
class visitor_bfs {
 public:
  using vertex_locator = typename graphstore_type::vertex_locator;
  using self_type = visitor_bfs<graphstore_type>;

  // Default constructor.
  visitor_bfs() :
      vertex() {  }

  // Baseline constructor.
  explicit visitor_bfs(vertex_locator _vertex) :
      vertex(_vertex) {  }


  bool pre_visit() const {
    return true;
  }


  template<typename VisitorQueueHandle>
  bool visit(graphstore_type& graph, VisitorQueueHandle vis_queue) const {
    if(graph.vertex_property_data(vertex)) {
      return false; /// already visited
    }
    graph.vertex_property_data(vertex) = true;
//    std::cout << "vist " << vertex.id() << std::endl;
    for (auto edge = graph.adjacent_edge_begin(vertex), end = graph.adjacent_edge_end(vertex);
         edge != end;
         ++edge) {
      const vertex_locator dst = edge.target_vertex();
//       std::cout << "push " << dst.id() << std::endl;
      self_type new_visitor(dst);
      vis_queue->queue_visitor(new_visitor);
    }
    return false; // or true?
  }

  friend inline bool operator > (const self_type& v1,
                                 const self_type& v2) {
    // TODO
//    return v1.vertex > v2.vertex;
  }

  // Instance variables.
  vertex_locator vertex;
} __attribute__((packed));

/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = bool;
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


template <typename gstore_type>
using dist_gstore_type = graphstore::dist_dynamic_graphstore<gstore_type>;

template <typename gstore_type>
using const_visitor_type = dg_visitor<dist_gstore_type<gstore_type>>;

template <typename gstore_type>
using dgs_const_visitor_queue_type = havoqgt::mpi::visitor_queue<const_visitor_type<gstore_type>,
                                                          havoqgt::detail::visitor_priority_queue,
                                                          dist_gstore_type<gstore_type>>;


template <typename gstore_type>
using bfs_visitor_type = visitor_bfs<dist_gstore_type<gstore_type>>;

template <typename gstore_type>
using dgs_bfs_visitor_queue_type = havoqgt::mpi::visitor_queue<bfs_visitor_type<gstore_type>,
                                                          havoqgt::detail::visitor_priority_queue,
                                                          dist_gstore_type<gstore_type>>;


/// --- option variables --- ///
uint64_t segmentfile_size_gb_        = 1;
uint64_t chunk_size_log10_           = 6;
std::string fname_segmentfile_       = "/dev/shm/segment_file";
std::string graphstore_name_         = "DegAwareRHH";
bool is_no_sync_                     = false;
std::vector<std::string> fname_edge_list_;
bool is_edgelist_with_delete_       = false;

void usage()  {
 if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
   std::cerr << "Usage: \n"
        << " -o <string>   - base filename to create segmentfiles\n"
        << " -S <int>      - the size of segmentfile per MPI rank in GB\n"
        << " -g <string>   - the name of graph store (DegAwareRHH, Baseline or BaselineMap)\n"
        << " -c <int>      - the logarithm with base ten of the chunk size (defaults is 6, i.g., chunk size is 10^6)\n"
        << " -E <string>   - the name of the file which has a list of edgelist file's makes\n"
        << " -d            - edgelist files have delete operations\n"
        << " -n            - don't call sync(2) every chunk\n"
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
 while ((c = getopt (argc, argv, "dc:o:S:g:E:nzfh")) != -1) {
   switch (c) {
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

     case 'n':
       is_no_sync_ = true;
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
   for (const auto itr : fname_edge_list_)
     std::cout << mpi_rank << " : " << "Load edge list from " << itr << std::endl;
 }

}



template <typename gstore_type, typename edgelist_type>
void constract_graph(dgs_const_visitor_queue_type<gstore_type>& dg_visitor_queue,
                     gstore_type& gstore,
                     edgelist_type& edges,
                     const size_t chunk_size)
{
  const int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  const int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
  const int lc_rank  = havoqgt::havoqgt_env()->node_local_comm().rank();
  havoqgt::havoqgt_env()->world_comm().barrier();


  /// --- variables for analysys --- //
  uint64_t loop_cnt = 0;
  bool global_is_finished = false;


  /// --- iterator and array for edgelist --- ///
  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  request_vector_type<vertex_id_type> update_request_vec;
  update_request_vec.reserve(chunk_size);

  double total_graph_storing_time = 0;
  const double whole_start = MPI_Wtime();
  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n\n<< Loop no. " << loop_cnt << " >>" << std::endl;

    /// --- generate or load edge update requests --- ///
    generate_update_requests(edges_itr, edges_itr_end, update_request_vec, chunk_size);
    havoqgt::havoqgt_env()->world_comm().barrier();

    const double time_start = MPI_Wtime();
    dg_visitor_queue.dynamic_graphconst(&update_request_vec);

    /// --- sync --- ///
    if (lc_rank == 0 && !is_no_sync_) graphstore::utility::sync_files();

    /// --- print a progress report --- ///
    havoqgt::havoqgt_env()->world_comm().barrier();
    const double time_end = MPI_Wtime();
    total_graph_storing_time += time_end - time_start;

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
    std::cout << "whole execution time (sec) : " << (whole_end - whole_start) << std::endl;
    std::cout << "only graph construction time (sec) (not inluding edge generation/load time) : " << (total_graph_storing_time) << std::endl;
    print_system_mem_usages();
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// --- print summary information --- ///
  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] Usage: segment size (GiB) =\t"<< graphstore::utility::segment_size_gb(gstore.get_segment_manager()) << std::endl;
      std::cout << "num edges =\t" << gstore.num_edges() << std::endl;
      gstore.print_status(0);
    }
    havoqgt::havoqgt_env()->world_comm().barrier();
  }

}

template<typename graphstore_type>
void run_benchmark()
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();

  /// --- init segment file --- ///
  size_t graph_capacity_gb_per_rank = segmentfile_size_gb_;
  havoqgt::distributed_db ddb(havoqgt::db_create(), fname_segmentfile_.c_str(), graph_capacity_gb_per_rank);
  havoqgt::havoqgt_env()->world_comm().barrier();
  if (mpi_rank == 0) print_system_mem_usages();

  /// --- init graphstore --- ///
  graphstore_type gstore(ddb.get_segment_manager());
  dist_gstore_type<graphstore_type> dist_graph(&gstore);

  /// --- init edgelist --- ///
  havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_, is_edgelist_with_delete_);

  /// --- construct a graph --- ///
  const_visitor_type<graphstore_type>::set_graph_ref(&dist_graph);
  dgs_const_visitor_queue_type<graphstore_type> dg_visitor_queue(&dist_graph);
  constract_graph(dg_visitor_queue, gstore, edgelist, std::pow(10, chunk_size_log10_));
  havoqgt::havoqgt_env()->world_comm().barrier();

  /// ---- run BFS --- ///
  dgs_bfs_visitor_queue_type<graphstore_type> dgs_bfs_visitor_queue(&dist_graph);
  const double time_start = MPI_Wtime();
  dgs_bfs_visitor_queue.traversal_no_delegate(bfs_visitor_type<graphstore_type>(dist_graph.label_to_locator(0)));
  havoqgt::havoqgt_env()->world_comm().barrier();
  const double time_end = MPI_Wtime();
  if (mpi_rank == 0) {
    std::cout << "BFS time (sec): " << time_end - time_start << std::endl;
  }
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    havoqgt::get_environment();
    havoqgt::havoqgt_env()->world_comm().barrier();

    parse_options(argc, argv);
    havoqgt::havoqgt_env()->world_comm().barrier();

//    /// --- init segment file --- ///
//    size_t graph_capacity_gb_per_rank = segmentfile_size_gb_;
//    havoqgt::distributed_db ddb(havoqgt::db_create(), fname_segmentfile_.c_str(), graph_capacity_gb_per_rank);

//    havoqgt::havoqgt_env()->world_comm().barrier();
//    if (mpi_rank == 0) print_system_mem_usages();

//    /// --- init edgelist --- ///
//    havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_, is_edgelist_with_delete_);

    /// --- allocate a graphstore and start a benchmark --- ///
    if (graphstore_name_ == "Baseline") {

      run_benchmark<baseline_type>();
//      baseline_type gstore(ddb.get_segment_manager());
//      dist_gstore_type<baseline_type> dist_graph(&gstore);

//      const_visitor_type<baseline_type>::set_graph_ref(&dist_graph);
//      dgs_const_visitor_queue_type<baseline_type> dg_visitor_queue(&dist_graph);
//      constract_graph(dg_visitor_queue, gstore, edgelist, std::pow(10, chunk_size_log10_));

//      dgs_bfs_visitor_queue_type<baseline_type> dgs_bfs_visitor_queue(&dist_graph);
//      dgs_bfs_visitor_queue.traversal_no_delegate(bfs_visitor_type<baseline_type>(dist_graph.label_to_locator(0)));

    } else if (graphstore_name_ == "BaselineMap") {

      run_benchmark<baselinemap_type>();
//      baselinemap_type gstore(ddb.get_segment_manager());
//      dist_gstore_type<baselinemap_type> dist_graph(&gstore);

//      const_visitor_type<baselinemap_type>::set_graph_ref(&dist_graph);
//      dgs_const_visitor_queue_type<baselinemap_type> dg_visitor_queue(&dist_graph);
//      constract_graph(dg_visitor_queue, gstore, edgelist, std::pow(10, chunk_size_log10_));

//      dgs_bfs_visitor_queue_type<baselinemap_type> dgs_bfs_visitor_queue(&dist_graph);
//      dgs_bfs_visitor_queue.traversal_no_delegate(bfs_visitor_type<baselinemap_type>(dist_graph.label_to_locator(0)));

    } else if (graphstore_name_ == "DegAwareRHH") {

      run_benchmark<degawarerhh_type>();
//      degawarerhh_type gstore(ddb.get_segment_manager());
//      dist_gstore_type<degawarerhh_type> dist_graph(&gstore);

//      const_visitor_type<degawarerhh_type>::set_graph_ref(&dist_graph);
//      dgs_const_visitor_queue_type<degawarerhh_type> dg_visitor_queue(&dist_graph);
//      constract_graph(dg_visitor_queue, gstore, edgelist, std::pow(10, chunk_size_log10_));

//      dgs_bfs_visitor_queue_type<degawarerhh_type> dgs_bfs_visitor_queue(&dist_graph);
//      dgs_bfs_visitor_queue.traversal_no_delegate(bfs_visitor_type<degawarerhh_type>(dist_graph.label_to_locator(0)));

    } else {
      std::cout << "Wrong graphstore name : " << graphstore_name_ << std::endl;
      exit(1);
    }
  }
  havoqgt::havoqgt_finalize();

  return 0;
}
