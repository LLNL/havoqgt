
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
    std::cout << "Call fallocate()" << std::endl;
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
  uint64_t graphfile_init_size = std::pow(2, 30) / mpi_size;
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

/*
  /// ------- delete edges and update vertices' property data ------- ///
  {
    for (int i = 0; i < 10; ++i) {
      vertex_id_type src = i % 2;
      vertex_id_type dst = i;

      vertex_property_type v_prop = graphstore.vertex_property_data(src);  // get a vertex property data or
      graphstore.vertex_property_data(src) = v_prop;                       // update a vertex property data.
                                                                           // Note that vertex_property_data() return a reference

      size_t erased_edges = graphstore.erase_edge(src, dst);               // erase edges
    }
  }


  /// ------- iterator over an adjacencylist ------- ///
  {
    vertex_id_type src_vrtx = 0;
    for (auto adj_edges = graphstore.adjacencylist(src_vrtx), end = graphstore.adjacencylist_end(src_vrtx);
         adj_edges != end;
         ++adj_edges) {
      std::cout << "destination vertex: " << adj_edges->first << ", weight: " << adj_edges->second << std::endl;
    }
  }
*/

  std::cout << "END" << std::endl;
  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}

