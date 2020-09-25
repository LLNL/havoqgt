// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_TEST_INCLUDE_CC_HPP
#define HAVOQGT_TEST_INCLUDE_CC_HPP

#include <iostream>
#include <string>

#include <havoqgt/connected_components.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/mpi.hpp>
#include <ingest_edges.hpp>

namespace havoqgt::test {

typedef graph_type::vertex_data<graph_type::vertex_locator,
                                std::allocator<graph_type::vertex_locator>>
    cc_data_type;

static constexpr const char *k_test_name = "test_cc";

void run_cc(const uint64_t max_vertex_id,
            const std::function<void(const graph_type &, const uint64_t,
                                     const cc_data_type &)>
                &validate_cc_result) {
  // Restore graph
  havoqgt::distributed_db ddb(havoqgt::db_open_read_only(),
                              gen_test_dir_path(k_test_name));
  auto                    graph =
      ddb.get_manager()->find<graph_type>(metall::unique_instance).first;

  // Allocate and init CC data
  cc_data_type cc_data(*graph);
  havoqgt::comm_world().barrier();

  // Run CC
  connected_components(graph, cc_data);
  havoqgt::comm_world().barrier();

  // Validate CC
  validate_cc_result(*graph, max_vertex_id, cc_data);
}

}  // namespace havoqgt::test

#endif  // HAVOQGT_TEST_INCLUDE_CC_HPP
