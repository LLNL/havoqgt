
#include <iostream>
#include <fstream>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/sssp_dynamic.hpp>
#include <havoqgt/graphstore/dist_dynamic_graphstore.hpp>

using mapped_file_type     = boost::interprocess::managed_mapped_file;
using segment_manager_type = havoqgt::distributed_db::segment_manager_type;

using vertex_id_type        = uint64_t;
using edge_property_type    = int;
using vertex_property_type  = uint32_t;

/// --- Baseline model --- ///
#include <havoqgt/graphstore/baseline/baseline.hpp>
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              segment_manager_type>;

/// --- DegAwareRHH model --- ///
#include <havoqgt/graphstore/degawarerhh/degawarerhh.hpp>
enum : size_t {
  middle_high_degree_threshold = 2 // must be more or equal than 1
};
using degawarerhh_type      = graphstore::degawarerhh<vertex_id_type,
                                                      vertex_property_type,
                                                      edge_property_type,
                                                      segment_manager_type,
                                                      middle_high_degree_threshold>;

/// --- chooes graphstore type ---- ///
/// you don't need to change any code other than here to change the graphstore model, Baseline or DegAwareRHH.
using graphstore_type       = degawarerhh_type; // baseline_type or degawarerhh_type

/// --- wrapper class to easy fit into the visitor program --- ///
using dist_graphstore_type  = graphstore::dist_dynamic_graphstore<graphstore_type>;

void usage()  {
  if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
         << " -c            - Prints out count of vertices with a given cost.\n"
         << " -v            - Verifies correctness of result.\n"
         << " -s            - Source vertex id.\n"
         << " -g <string>   - output graph base filename (required)\n"
         << " -e <string>   - filename that has a list of edgelist files (required)\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_file,
                    std::vector<std::string>& edgelist_files,
                    bool* verify, bool* level_count, vertex_id_type* source_v) {
  if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_graph_filename_ = false;
  bool found_edgelist_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "cvg:e:hs:")) != -1) {
     switch (c) {
       case 'c':
         *level_count = true;
         break;
       case 'v':
         *verify = true;
         break;
       case 'h':
         prn_help = true;
         break;
       case 's':
         *source_v = atoll(optarg);
         break;
       case 'g':
         found_graph_filename_ = true;
         graph_file = optarg;
         break;
       case 'e':
       {
         found_edgelist_filename = true;
         std::string fname(optarg);
         std::ifstream fin(fname);
         std::string line;
         if (!fin.is_open()) {
           std::cerr << fname << std::endl;
           HAVOQGT_ERROR_MSG("Unable to open a file");
         }
         while (std::getline(fin, line)) {
           edgelist_files.push_back(line);
         }
         break;
       }
       default:
         std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
         prn_help = true;
         break;
     }
   }
   if (prn_help || !found_graph_filename_ || !found_edgelist_filename) {
     usage();
     exit(-1);
   }
}


int main(int argc, char** argv) {
  double graph_capacity_gb_per_rank = 4.0;
  int mpi_rank(0), mpi_size(0);
  bool level_count = false;
  bool verify = false;
  vertex_id_type source_v = 0;

  havoqgt::havoqgt_init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    havoqgt::get_environment().print();
  }
  MPI_Barrier(MPI_COMM_WORLD);


  /// --- parse argments ---- ///
  std::string graph_file;
  std::vector<std::string> edgelist_files;
  parse_cmd_line(argc, argv, graph_file, edgelist_files, &verify, &level_count, &source_v);
  MPI_Barrier(MPI_COMM_WORLD);


  /// --- init segment file --- ///
  std::stringstream fname_local;
  fname_local << graph_file << "_" << mpi_rank;
  havoqgt::distributed_db ddb(havoqgt::db_create(), fname_local.str().c_str(), graph_capacity_gb_per_rank);

  /// --- get a segment_manager --- ///
  segment_manager_type* segment_manager = ddb.get_segment_manager();

  // Get edgelist.
  havoqgt::parallel_edge_list_reader edgelist(edgelist_files);

  /// --- allocate a graphstore --- ///
  graphstore_type graphstore(segment_manager);
  dist_graphstore_type dist_graphstore(&graphstore);


  // Run algorithm.
  havoqgt::mpi::single_source_shortest_path_dynamic<dist_graphstore_type, havoqgt::parallel_edge_list_reader, vertex_property_type>(&dist_graphstore, &edgelist, dist_graphstore_type::vertex_locator(source_v), verify);

  if (mpi_rank == 0) {
    std::cout << "Collecting stats...\n";
  }

  // Collect stats.
  uint64_t visited_total = 0;

  std::map<uint64_t, uint64_t> costs_count;
  for (auto v = graphstore.vertices_begin(); v != graphstore.vertices_end(); v++) {
    auto v_cost = v.property_data();
    costs_count[v_cost]++;
  }

  auto costs_count_itr = costs_count.begin();
  uint64_t largest_cost = 0;
  uint64_t num_costs = 0;
  while(!havoqgt::mpi::detail::global_iterator_range_empty(costs_count_itr, costs_count.end(), MPI_COMM_WORLD)) {
    uint64_t local_next_cost = costs_count_itr->first;
    uint64_t global_next_cost = havoqgt::mpi::mpi_all_reduce(local_next_cost, std::less<uint64_t>(), MPI_COMM_WORLD);
    uint64_t local_count = (local_next_cost == global_next_cost) ?  (costs_count_itr++)->second : 0;
    uint64_t global_cost_count = havoqgt::mpi::mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
    if(mpi_rank == 0) {
      std::cout << "Cost " << global_next_cost << ", size = " << global_cost_count << std::endl;
    }
    largest_cost = std::max(global_cost_count, largest_cost);
    num_costs++;
    visited_total += global_cost_count;
  }

  if (mpi_rank == 0) {
    std::cout << "Number of cost groups = " << num_costs <<  std::endl
              << "Visited total = " << visited_total << std::endl
              ;  // << "GC Time = " << time_end - time_start << std::endl;
  }

  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}

