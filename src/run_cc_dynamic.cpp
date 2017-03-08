
#include <iostream>
#include <fstream>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/connected_components_dynamic.hpp>
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
         << " -c            - Prints out count of each component.\n"
         << " -v            - Verifies correctness of result.\n"
         << " -g <string>   - output graph base filename (required)\n"
         << " -e <string>   - filename that has a list of edgelist files (required)\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_file,
                    std::vector<std::string>& edgelist_files,
                    bool* verify, bool* co_count) {
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
  while ((c = getopt(argc, argv, "cvg:e:h")) != -1) {
     switch (c) {
       case 'c':
         *co_count = true;
         break;
       case 'v':
         *verify = true;
         break;
       case 'h':
         prn_help = true;
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
  bool co_count = false;
  bool verify = false;

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
  parse_cmd_line(argc, argv, graph_file, edgelist_files, &verify, &co_count);
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
  havoqgt::mpi::connected_components_dynamic<dist_graphstore_type, havoqgt::parallel_edge_list_reader, vertex_property_type>(&dist_graphstore, &edgelist, verify);

  if (mpi_rank == 0) {
    std::cout << "Collecting stats...\n";
  }

  // Collect stats.
  uint64_t visited_total = 0;

  std::map<uint64_t, uint64_t> cc_count;
  for (auto v = graphstore.vertices_begin(); v != graphstore.vertices_end(); v++) {
    auto v_component = v.property_data();
    cc_count[v_component]++;
  }

  auto cc_count_itr = cc_count.begin();
  uint64_t largest_cc = 0;
  uint64_t num_ccs = 0;
  while(!havoqgt::mpi::detail::global_iterator_range_empty(cc_count_itr, cc_count.end(), MPI_COMM_WORLD)) {
    uint64_t local_next_cc = cc_count_itr->first;
    uint64_t global_next_cc = havoqgt::mpi::mpi_all_reduce(local_next_cc, std::less<uint64_t>(), MPI_COMM_WORLD);
    uint64_t local_count = (local_next_cc == global_next_cc) ?  (cc_count_itr++)->second : 0;
    uint64_t global_cc_count = havoqgt::mpi::mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
    if(mpi_rank == 0 ) {
      std::cout << "CC " << global_next_cc << ", size = " << global_cc_count << std::endl;
    }
    largest_cc = std::max(global_cc_count, largest_cc);
    num_ccs++;
    visited_total += global_cc_count;
  }

  if (mpi_rank == 0 && visited_total > 1) {
    std::cout << "Number of Components = " << num_ccs << ", largest = " << largest_cc << std::endl
              << "Visited total = " << visited_total << std::endl;
  }

  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}

