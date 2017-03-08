
#include <iostream>
#include <fstream>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/bfs_dynamic.hpp>
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
         << " -c            - Prints out count of each level.\n"
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
  havoqgt::mpi::breadth_first_search_dynamic<dist_graphstore_type, havoqgt::parallel_edge_list_reader, vertex_property_type>(&dist_graphstore, &edgelist, dist_graphstore_type::vertex_locator(source_v), verify);

  if (mpi_rank == 0) {
    std::cout << "Collecting stats...\n";
  }

  // Collect stats.
  uint64_t visited_total = 0;
  bool finished = false;
  uint32_t num_levels = 0;
  const int AGGR_SIZE = 100;

  std::vector<uint32_t> level_counts = *(new std::vector<uint32_t>());

  // Agregates 100 levels per all reduce.
  for (uint32_t level = 0; finished != true; level += AGGR_SIZE) {
    std::vector<uint64_t> local_count = *(new std::vector<uint64_t>());
    local_count.resize(AGGR_SIZE, 0);
    std::vector<uint64_t> global_count = *(new std::vector<uint64_t>());
    global_count.resize(AGGR_SIZE, 0);
    level_counts.resize(level_counts.size() + AGGR_SIZE, 0);

    for (auto v = graphstore.vertices_begin(); v != graphstore.vertices_end(); v++) {
      auto v_level = v.property_data();
      if (v_level >= level && v_level < (level + AGGR_SIZE)) {
        local_count[v_level % AGGR_SIZE]++;
      }
    }
    havoqgt::mpi::mpi_all_reduce(local_count, global_count,
                                 std::plus<uint64_t>(), MPI_COMM_WORLD);

    for (uint32_t i = 0; i < AGGR_SIZE; i++) {
      // note: level 0 should have 0 count.
      if (global_count[i] != 0 || (level == 0 && i == 0)) {
        level_counts[level + i] = global_count[i];

        visited_total += global_count[i];

        // Print out per-level count if desired.
        if (level_count && mpi_rank == 0) {
          if (level + i == 0) {
            std::cout << "Sanity check (should be zero): ";
          } else {
            std::cout << "Level " << level + i << ": ";
          }
          std::cout << global_count[i] << std::endl;
        }
      } else {
        num_levels = level + i - 1;  // One offset since 0 is not a level.
        finished = true;
        break;
      }
    }
  }

  if (mpi_rank == 0) {
    std::cout << "Number of levels = " << num_levels <<  std::endl
              << "Visited total = " << visited_total << std::endl
              ;  // << "GC Time = " << time_end - time_start << std::endl;
  }

  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}

