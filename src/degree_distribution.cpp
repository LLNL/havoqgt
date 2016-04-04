#include <iostream>
#include <sstream>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>

#include <vector>
#include <unordered_map>

using namespace havoqgt;

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

class histogram_tracker {
public:
  histogram_tracker() : max_deg(0) { }

  void operator()(uint64_t key) {
    if(histogram.find(key) == histogram.end()) histogram.insert( std::make_pair( key, 0 ) );
    histogram[key]++;
    max_deg = std::max( max_deg, key);
  }

  uint64_t get_max_deg() {
    return max_deg;
  }

  std::unordered_map<uint64_t, uint64_t> histogram;
  uint64_t max_deg;
};

int main(int argc, char** argv) {
  /*  if( argc != 5 ) {
    std::cout << "Usage: run_backup <Edgelist> <Metadata> <Backup EdgeList> <Backup Metadata>" << std::endl;
    exit(1);
    }*/
  std::string input_graph_filename = std::string( argv[1]);
  havoqgt_init(&argc, &argv);
  {
  havoqgt::get_environment();
  int mpi_rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
  histogram_tracker ht;
  havoqgt::distributed_db graph_ddb(havoqgt::db_open(), input_graph_filename.c_str());
  //  havoqgt::distributed_db metadata_ddb(havoqgt::db_open(), input_metadata_filename.c_str());
  
  graph_type *graph = graph_ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
  for( auto itr = graph->vertices_begin(); itr != graph->vertices_end(); ++itr) {
    graph_type::vertex_locator v = *itr;
    uint64_t edge_count = static_cast<uint64_t>( graph->edges_end(v) - graph->edges_begin(v) );
    ht(edge_count);
  }
  havoqgt_env()->world_comm().barrier();
  uint64_t global_max = mpi_all_reduce( ht.get_max_deg(), std::greater<uint64_t>(), MPI_COMM_WORLD);
  if( mpi_rank == 0) { std::cout << "Max: " << global_max << std::endl; }
  havoqgt_env()->world_comm().barrier();
 
  //get max vertex id
  for(uint64_t v = 1; v <= global_max; ) {
    std::vector<uint64_t> vec, out;
    uint64_t before = v;
    for(uint64_t k = 0; k < 1000000UL && v <= global_max; v++, k++) {
      if( ht.histogram.find(v) == ht.histogram.end())
	vec.push_back(0);
      else vec.push_back( ht.histogram[v] );
    }
    out.resize( vec.size() );
    mpi_all_reduce( vec, out, std::plus<uint64_t>(),MPI_COMM_WORLD );
    if( mpi_rank == 0 ) {
      uint64_t i = 0;
      for(; before < v; ++i, ++before) {
	if( out[i] != 0 )
	  std::cout << before << ";" << out[i] << std::endl;
      }
    }
    havoqgt_env()->world_comm().barrier();
  }
  havoqgt_env()->world_comm().barrier();
  }
  havoqgt::havoqgt_finalize();
  return 0;
}
