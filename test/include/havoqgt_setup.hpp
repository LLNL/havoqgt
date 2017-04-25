#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef typename graph_type::edge_iterator eitr_type;
typedef typename graph_type::vertex_iterator vitr_type;
typedef typename graph_type::vertex_locator vloc_type;

typedef double edge_data_type;
typedef std::tuple<std::pair<uint64_t, uint64_t>, edge_data_type> edge_type; 

