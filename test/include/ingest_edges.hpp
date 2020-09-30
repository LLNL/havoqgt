// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_TEST_INCLUDE_INGEST_EDGES_HPP
#define HAVOQGT_TEST_INCLUDE_INGEST_EDGES_HPP

#include <string>
#include <vector>

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <util.hpp>

namespace havoqgt::test {

// Graph type
typedef havoqgt::delegate_partitioned_graph<
    havoqgt::distributed_db::allocator<>>
    graph_type;

// Unweighted Edge data type
typedef uint8_t dummy_edge_data_type;
typedef havoqgt::distributed_db::allocator<dummy_edge_data_type>
    unweighted_edge_data_allocator_type;
typedef graph_type::edge_data<dummy_edge_data_type,
                              unweighted_edge_data_allocator_type>
    unweighted_edge_data_type;

/// \brief Ingest unweighted edge lists.
/// Edge lists must be unweighted.
void ingest_unweighted_edges(
    const std::vector<std::string> &edge_list_file_names,
    const std::string &             data_store_path,
    const uint64_t                  delegate_degree_threshold = 1024 * 10) {

  havoqgt::parallel_edge_list_reader<dummy_edge_data_type> pelr(
      edge_list_file_names, true);

  bool has_edge_data = pelr.has_edge_data();
  if (has_edge_data) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // Create a new graph
  havoqgt::distributed_db::remove(data_store_path);
  havoqgt::distributed_db   ddb(havoqgt::db_create(), data_store_path);

  unweighted_edge_data_type dummy_edge_data(ddb.get_allocator());

  // Allocates and constructs a graph
  graph_type *              graph =
      ddb.get_manager()->construct<graph_type>(metall::unique_instance)(
          ddb.get_allocator(), MPI_COMM_WORLD, pelr, pelr.max_vertex_id(),
          delegate_degree_threshold, 1, 8192, dummy_edge_data);

  if (!graph) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}
}  // namespace havoqgt::test

#endif  // HAVOQGT_TEST_INCLUDE_INGEST_EDGES_HPP
