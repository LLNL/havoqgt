/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License. http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by the United States Government or Lawrence Livermore National Security, LLC.
 * The views and opinions of authors expressed herein do not necessarily state
 * or reflect those of the United States Government or Lawrence Livermore
 * National Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

#ifndef HAVOQGT_TEST_INCLUDE_BFS_HPP
#define HAVOQGT_TEST_INCLUDE_BFS_HPP

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <havoqgt/breadth_first_search.hpp>
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

// Edge data type
typedef uint8_t edge_data_type;
typedef havoqgt::distributed_db::allocator<edge_data_type>
    edge_data_allocator_type;

// BFS data types
typedef graph_type::vertex_data<uint16_t, std::allocator<uint16_t>>
    bfs_level_data_type;
typedef graph_type::vertex_data<graph_type::vertex_locator,
                                std::allocator<graph_type::vertex_locator>>
    bfs_parent_data_type;

static constexpr const char *k_test_name      = "test_bfs";
static constexpr const char *k_graph_key_name = "graph_obj";
static constexpr uint16_t k_unvisited_level = std::numeric_limits<uint16_t>::max();

/// \brief Ingest edge lists
void ingest_edges(const std::vector<std::string> &edge_list_file_names,
                  const uint64_t delegate_degree_threshold = 1024 * 10) {
  havoqgt::parallel_edge_list_reader<edge_data_type> pelr(edge_list_file_names,
                                                          true);
  bool has_edge_data = pelr.has_edge_data();
  if (has_edge_data) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  havoqgt::distributed_db ddb(havoqgt::db_create(),
                              gen_test_dir_path(k_test_name));
  graph_type::edge_data<edge_data_type, edge_data_allocator_type>
              dummy_edge_data(ddb.get_allocator());
  graph_type *graph =
      ddb.get_manager()->construct<graph_type>(k_graph_key_name)(
          ddb.get_allocator(), MPI_COMM_WORLD, pelr, pelr.max_vertex_id(),
          delegate_degree_threshold, 1, 8192, dummy_edge_data);
  if (!graph) {
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
}

/// \brief Run BFS from vertex 0
void run_bfs(const uint64_t max_vertex_id,
             const std::function<void(const graph_type &, const uint64_t,
                                      const bfs_level_data_type &,
                                      const bfs_parent_data_type &)>
                 &validate_bfs_result) {
  havoqgt::distributed_db ddb(havoqgt::db_open_read_only(),
                              gen_test_dir_path(k_test_name));
  auto graph = ddb.get_manager()->find<graph_type>(k_graph_key_name).first;
  bfs_level_data_type bfs_level_data(*graph);
  bfs_level_data.reset(k_unvisited_level);
  bfs_parent_data_type bfs_parent_data(*graph);
  havoqgt::comm_world().barrier();

  auto source = graph->label_to_locator(0);
  havoqgt::breadth_first_search(graph, bfs_level_data, bfs_parent_data, source);

  validate_bfs_result(*graph, max_vertex_id, bfs_level_data, bfs_parent_data);
}

}  // namespace havoqgt::test

#endif  // HAVOQGT_TEST_INCLUDE_BFS_HPP
