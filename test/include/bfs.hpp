// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_TEST_INCLUDE_BFS_HPP
#define HAVOQGT_TEST_INCLUDE_BFS_HPP

#include <cmath>
#include <iostream>
#include <utility>

#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>
#include <ingest_edges.hpp>

namespace havoqgt::test {

// BFS data types
typedef graph_type::vertex_data<uint16_t, std::allocator<uint16_t>>
    bfs_level_data_type;
typedef graph_type::vertex_data<graph_type::vertex_locator,
                                std::allocator<graph_type::vertex_locator>>
    bfs_parent_data_type;

static constexpr const char *k_test_name      = "test_bfs";
static constexpr uint16_t k_unvisited_level = std::numeric_limits<uint16_t>::max();


/// \brief Run BFS from vertex 0
void run_bfs(const uint64_t max_vertex_id,
             const std::function<void(const graph_type &, const uint64_t,
                                      const bfs_level_data_type &,
                                      const bfs_parent_data_type &)>
                 &validate_bfs_result) {

  // Restore graph
  havoqgt::distributed_db ddb(havoqgt::db_open_read_only(),
                              gen_test_dir_path(k_test_name));
  auto graph = ddb.get_manager()->find<graph_type>(metall::unique_instance).first;

  // Allocates and init BFS data
  bfs_level_data_type bfs_level_data(*graph);
  bfs_level_data.reset(k_unvisited_level);
  bfs_parent_data_type bfs_parent_data(*graph);

  havoqgt::comm_world().barrier();

  // Run BFS
  auto source = graph->label_to_locator(0);
  havoqgt::breadth_first_search(graph, bfs_level_data, bfs_parent_data, source);

  // Validate BFS results
  validate_bfs_result(*graph, max_vertex_id, bfs_level_data, bfs_parent_data);
}

}  // namespace havoqgt::test

#endif  // HAVOQGT_TEST_INCLUDE_BFS_HPP
