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


/// --- typenames --- ///
using vertex_id_type        = uint64_t;
using edge_property_type    = unsigned char;
using vertex_property_type  = unsigned char;
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                             vertex_property_type,
                                                             edge_property_type,
                                                             havoqgt::distributed_db::segment_manager_type>;

using baselinemap_type      = graphstore::graphstore_baseline_map<vertex_id_type,
                                                                  vertex_property_type,
                                                                  edge_property_type,
                                                                  havoqgt::distributed_db::segment_manager_type>;

enum : size_t { middle_high_degree_threshold = 8 }; // must be more than 1
using degawarerhh_type  = graphstore::degawarerhh<vertex_id_type,
                                                     vertex_property_type,
                                                     edge_property_type,
                                                     havoqgt::distributed_db::segment_manager_type,
                                                     middle_high_degree_threshold>;


/// --- Choose graph store type --- ///
/// baseline_type, baselinemap_type or degawarerhh_type
using gstore_type = degawarerhh_type;

using dist_gstore_type = graphstore::dist_dynamic_graphstore<gstore_type>;
using visitor_type = dg_visitor<dist_gstore_type>;
using dg_visitor_queue_type = havoqgt::mpi::visitor_queue<visitor_type,
                                                          havoqgt::detail::visitor_priority_queue,
                                                          dist_gstore_type>;

void usage()  {
 if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
   std::cerr << "Usage: \n"
        << " -o <string>   - base filename to create segmentfiles\n"
        << " -h            - print help and exit\n\n";
 }
}

void parse_options(int argc, char **argv, std::string& segment_file_prefix)
{
 int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
 int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

 if (mpi_rank == 0) {
   std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
 }

 char c;
 while ((c = getopt (argc, argv, "o:h")) != -1) {
   switch (c) {
     case 'o':
     {
       segment_file_prefix = optarg;
       break;
     }
     case 'h':
       usage();
       break;
   }
 }
}


void generate_edgelist(request_vector_type<vertex_id_type>& edgelist)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  uint64_t num_edges_per_rank = (1ULL << 17) * 16ULL / mpi_size;
  havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                    17, num_edges_per_rank,
                                    0.57, 0.19, 0.19, 0.05, true, true);
  edgelist.clear();
  for (const auto& edge : rmat) {
    edgelist.push_back(make_update_request(edge));
  }
}

using mpi_buf_vec_type = std::vector<request_vector_type<vertex_id_type>>;
void global_sort(request_vector_type<vertex_id_type>& edgelist)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  mpi_buf_vec_type send_buf_vec(mpi_size);
  mpi_buf_vec_type recev_buf_vec(mpi_size);

  for (const auto& edge : edgelist) {
    const int target_rank = std::get<0>(edge.edge) % mpi_rank;
    send_buf_vec[target_rank].push_back(edge);
  }
  edgelist.clear();
  havoqgt::havoqgt_env()->world_comm().barrier();

  havoqgt::mpi::mpi_all_to_all(send_buf_vec, recev_buf_vec, havoqgt::havoqgt_env()->world_comm().comm());
  havoqgt::havoqgt_env()->world_comm().barrier();

  for (auto& vec_per_rank : recev_buf_vec) {
    for (auto& edge : vec_per_rank) {
      edgelist.push_back(edge);
    }
  }
  sort_requests(edgelist);
}

void dump_all_edges(dist_gstore_type& dist_graph, request_vector_type<vertex_id_type>& buffer)
{
  buffer.clear();
  for (auto vrtx = dist_graph.vertices_begin(), vrtx_end = dist_graph.vertices_end();
       vrtx != vrtx_end;
       ++vrtx) {
    for (auto edge = dist_graph.adjacent_edge_begin(vrtx.source_vertex()), edge_end = dist_graph.adjacent_edge_end(vrtx.source_vertex());
         edge != edge_end;
         ++edge) {
      buffer.push_back(make_update_request(std::make_pair(vrtx.source_vertex().id(), edge.target_vertex().id())));
    }
  }
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    havoqgt::get_environment();
    havoqgt::havoqgt_env()->world_comm().barrier();

    /// --- init segment file --- ///
    std::string segment_file_prefix = "/dev/shm/gstore";
    parse_options(argc, argv, segment_file_prefix);
    size_t graph_capacity_gb_per_rank = 2;
    havoqgt::distributed_db ddb(havoqgt::db_create(), segment_file_prefix.c_str(), graph_capacity_gb_per_rank);
    havoqgt::havoqgt_env()->world_comm().barrier();

    gstore_type gstore(ddb.get_segment_manager());
    dist_gstore_type dist_graph(&gstore);
    visitor_type::set_graph_ref(&dist_graph);
    dg_visitor_queue_type dg_visitor_queue(&dist_graph);
    havoqgt::havoqgt_env()->world_comm().barrier();

    request_vector_type<vertex_id_type> edgelist;
    generate_edgelist(edgelist);
    havoqgt::havoqgt_env()->world_comm().barrier();

    /// --- start dynamic graph construction --- ///
    dg_visitor_queue.dynamic_graphconst(&edgelist);

    request_vector_type<vertex_id_type> stored_edgelist;
    dump_all_edges(dist_graph, stored_edgelist);
    sort_requests(stored_edgelist);
    global_sort(edgelist);

    assert(stored_edgelist.size() == edgelist.size());
    for (size_t i = 0; i < stored_edgelist.size(); ++i) {
      assert(stored_edgelist[i].edge == edgelist[i].edge);
    }
  }
}
