
#include <iostream>
#include <fstream>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include "../include/havoqgt/graph_colour_dynamic.hpp"


using mapped_file_type     = boost::interprocess::managed_mapped_file;
using segment_manager_type = boost::interprocess::managed_mapped_file::segment_manager;

using vertex_id_type        = uint64_t;
using edge_property_type    = int;
using vertex_property_type  = uint32_t;

#include <havoqgt/graphstore/baseline.hpp>
using graphstore_type       = graphstore::graphstore_baseline<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              segment_manager_type>;



void fallocate(const char* const fname, size_t size, mapped_file_type& asdf)
{
#ifdef __linux__
    // std::cout << "Call fallocate()" << std::endl;
    int fd  = open(fname, O_RDWR);
    assert(fd != -1);
    /// posix_fallocate dosen't work on XFS ?
    /// (dosen't actually expand the file size ?)
    int ret = fallocate(fd, 0, 0, size);
    assert(ret == 0);
    close(fd);
    asdf.flush();
#else
#warning fallocate() is not supported
#endif
}


void usage()  {
  if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
         << " -g <string>   - output graph base filename (required)\n"
         << " -e <string>   - filename that has a list of edgelist files (required)\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_file, std::vector<std::string>& edgelist_files) {
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
  while ((c = getopt(argc, argv, "g:e:h")) != -1) {
     switch (c) {
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
  const int SCALE = 32;
  int mpi_rank(0), mpi_size(0);

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
  parse_cmd_line(argc, argv, graph_file, edgelist_files);
  MPI_Barrier(MPI_COMM_WORLD);


  /// --- create a segument file --- ///
  std::stringstream fname_local;
  fname_local << graph_file << "_" << mpi_rank;
  boost::interprocess::file_mapping::remove(fname_local.str().c_str());
  uint64_t graphfile_init_size = std::pow(2, SCALE) / mpi_size;
  mapped_file_type mapped_file = mapped_file_type(boost::interprocess::create_only, fname_local.str().c_str(), graphfile_init_size);

  /// --- Call fallocate --- ///
  fallocate(fname_local.str().c_str(), graphfile_init_size, mapped_file);

  /// --- get a segment_manager --- ///
  segment_manager_type* segment_manager = mapped_file.get_segment_manager();

  // Get edgelist.
  havoqgt::parallel_edge_list_reader edgelist(edgelist_files);

  /// --- allocate a graphstore --- ///
  graphstore_type graphstore(segment_manager, &edgelist);

  // Run algorithm.
  havoqgt::mpi::graph_colour_dynamic<graphstore_type, vertex_property_type>(&graphstore);

  if (mpi_rank == 0) {
    std::cout << "Collecting stats...\n";
  }

  // Collect stats.
  uint64_t visited_total = 0;
  bool finished = false;
  uint32_t num_colours = 0;
  const int AGGR_SIZE = 100;

  std::vector<uint32_t> colour_counts = *(new std::vector<uint32_t>());

  // Agregates 100 colours per all reduce.
  for (uint32_t colour = 0; finished != true; colour += AGGR_SIZE) {
    std::vector<uint64_t> local_count = *(new std::vector<uint64_t>());
    local_count.resize(AGGR_SIZE, 0);
    std::vector<uint64_t> global_count = *(new std::vector<uint64_t>());
    global_count.resize(AGGR_SIZE, 0);
    colour_counts.resize(colour_counts.size() + AGGR_SIZE, 0);

    for (auto v = graphstore.vertices_begin(); v != graphstore.vertices_end(); v++) {
      auto v_col = v.property_data();
      if (v_col >= colour && v_col < (colour + AGGR_SIZE)) {
        local_count[v_col % AGGR_SIZE]++;
      }
    }
    havoqgt::mpi::mpi_all_reduce(local_count, global_count,
                                 std::plus<uint64_t>(), MPI_COMM_WORLD);

    for (uint32_t i = 0; i < AGGR_SIZE; i++) {
      // Note: colour 0 (colour = 0, i = 0) will have 0 count.
      if (global_count[i] != 0 || (colour == 0 && i == 0)) {
        colour_counts[colour + i] = global_count[i];
        visited_total += global_count[i];

        // Print out per-colour count if desired.
        if (true && mpi_rank == 0) {
          std::cout << "Colour " << colour + i << ": " << global_count[i]
                    << std::endl;
        }
      } else {
        num_colours = colour + i - 1;  // One offset since 0 is not a colour.
        finished = true;
        break;
      }
    }
  }

  if (mpi_rank == 0 && visited_total > 1) {
    std::cout << "Number of Colours = " << num_colours <<  std::endl
              << "Visited total = " << visited_total << std::endl
              ;  // << "GC Time = " << time_end - time_start << std::endl;
  }

  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}

