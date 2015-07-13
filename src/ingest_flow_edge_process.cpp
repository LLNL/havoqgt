#define MASTER_DO(task) if(mpi_rank == 0) { task }
#define MASTER_MSG(msg) if(mpi_rank == 0) { std::cout << msg << std::endl; }

#include <havoqgt/environment.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/ingest_flow_edge_list.hpp>
#include <havoqgt/read_edge_metadata.hpp>
#include <havoqgt/parallel_flow_edge_list_reader.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <unistd.h>

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef graph_type::edge_data<flow, segment_manager_t> edge_data_type;
void usage() {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -o <string> -d <int> [file ...]\n"
	      << "-o <string>    - output graph base filename (required)\n"
	      << "-m <string>    - output metadata base filename (required)\n"
	      << "-b <string>    - backup graph base filename\n"
	      << "-d <int>       - delegate threshold ( Default is 1048576 )\n"
	      << "-h             - print help and exit\n"
	      << "-p <int>       - number of Low & High partition passes (Default is 1)\n"
	      << "-f <float>     - Gigabytes reserved per rank (Default is 0.25)\n"
	      << "-c <int>       - Edge partitioning chunk size (Default is 8192)\n"
	      << "-u <bool>      - Treat edgelist as undirected (Default is 0)\n"
	      << "[file ...]     - list of edge list files to ingest\n\n";
  }
}

void parse_cmd_line(int argc, char **argv, std::string& output_filename, std::string& metadata_output_filename/* ++ */
		    ,std::string& backup_filename, uint64_t& delegate_threshold, std::vector<std::string>& input_filenames
		    ,double& gbyte_per_rank, uint64_t& partition_passes, uint64_t& chunk_size
		    ,bool& undirected) {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for(int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }
  
  bool found_output_filename = false;
  bool found_metadata_output_filename = false;
  delegate_threshold = 1048576;
  input_filenames.clear();
  gbyte_per_rank = 0.25; // change this and @todo track what this value is
  partition_passes = 1;
  chunk_size = 8*1024;
  undirected = false;

  char c;
  bool prn_help = false;
  
  while((c = getopt(argc, argv, "o:m:d:p:f:c:b:u:h ")) != -1) {
    switch (c) {
    case 'h':
      prn_help = true;
      break;
    case 'd':
      delegate_threshold = atoll(optarg);
      break;
    case 'o':
      found_output_filename = true;
      output_filename = optarg;
      break;
    case 'm':
      found_metadata_output_filename = true;
      metadata_output_filename = optarg;
      break;
    case 'b':
      backup_filename = optarg;
      break;
    case 'p':
      partition_passes = atoll(optarg);
      break;
    case 'f':
      gbyte_per_rank = atof(optarg);
      break;
    case 'c':
      chunk_size = atoll(optarg);
      break;
    case 'u':
      undirected = atoi(optarg);
      break;
    default:
      std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
      prn_help = true;
      break;
    }
  }

  if(prn_help || !found_output_filename || !found_metadata_output_filename ) {
    usage();
    exit(-1);
  }

  for(int index = optind; index < argc; ++index) {
    input_filenames.push_back( argv[index] );
  }
}

template <typename Graph, typename EdgeData, typename MetaData>
int edge_metadata_visitor<Graph, EdgeData, MetaData>::count = 0;

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  havoqgt_init(&argc, &argv);
  {
    std::string output_filename;
    std::string backup_filename;
    std::string metadata_output_filename;
    {
      int mpi_rank = havoqgt_env()->world_comm().rank();
      int mpi_size = havoqgt_env()->world_comm().size();
      havoqgt::get_environment();

      MASTER_DO(
	std::cout << "MPI Initialized with " << mpi_size << "ranks. " << std::endl;
	havoqgt::get_environment().print();
	);
      havoqgt_env()->world_comm().barrier();

      // identify what parameters are required here?
      std::vector<std::string>   input_filenames;
      uint64_t                   delegate_threshold;
      uint64_t                   partition_passes;
      double                     gbyte_per_rank;
      uint64_t                   chunk_size;
      bool                       undirected;

      parse_cmd_line(argc, argv, output_filename, metadata_output_filename, backup_filename, delegate_threshold,
		     input_filenames, gbyte_per_rank, partition_passes, chunk_size, undirected);
      MASTER_DO(
	std::cout << "Ingesting graph from " << input_filenames.size() << " files." << std::endl;
	);

      havoqgt::distributed_db ddb(havoqgt::db_create(), output_filename.c_str(), gbyte_per_rank);
      
      segment_manager_t* segment_manager = ddb.get_segment_manager();
      bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);
      
      havoqgt::parallel_flow_edge_list_reader pelr(input_filenames, undirected);

      MASTER_MSG("Generating new graph.");
      
      graph_type *graph = segment_manager->construct<graph_type>
	("graph_obj")
	(alloc_inst, havoqgt_env()->world_comm().comm(), pelr, pelr.max_vertex_id(), delegate_threshold,
	 partition_passes, chunk_size);
      
      havoqgt_env()->world_comm().barrier();
      
      MASTER_MSG("Graph Topology Ready. Integrating Meta Data.");
      
      havoqgt::distributed_db metadata_ddb(havoqgt::db_create(), metadata_output_filename.c_str(), gbyte_per_rank);
      segment_manager_t* metadata_segment_manager = metadata_ddb.get_segment_manager();

      //typedef graph_type::edge_data<flow, segment_manager_t> edge_data_type;
      edge_data_type *edge_metadata = metadata_segment_manager->construct<edge_data_type>("edge_data")(graph->num_owned_targets(), graph->num_delegate_targets(), metadata_segment_manager);

      havoqgt::ingest_flow_edge_list ifel(input_filenames);
      havoqgt_env()->world_comm().barrier();
      graph->print_graph_statistics();
      havoqgt_env()->world_comm().barrier();
      int count = 0;
      hmpi::read_edge_metadata(graph, edge_metadata, ifel.begin(), ifel.end(), count);
      havoqgt_env()->world_comm().barrier();
      MASTER_MSG("Reading of metadata completed.");

      int total_count = mpi_all_reduce(count, std::plus<int>(), havoqgt_env()->world_comm().comm());
      MASTER_DO(
		std::cout << "Total Number of vertices : " << total_count << std::endl;
		);
    }
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
