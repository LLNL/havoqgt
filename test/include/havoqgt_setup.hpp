// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>

#include <havoqgt/parallel_edge_list_reader.hpp>

using namespace havoqgt;

typedef havoqgt::distributed_db::manager_type manager_t;
typedef havoqgt::delegate_partitioned_graph<havoqgt::distributed_db::allocator<>> graph_type;
typedef typename graph_type::edge_iterator eitr_type;
typedef typename graph_type::vertex_iterator vitr_type;
typedef typename graph_type::vertex_locator vloc_type;

typedef double edge_data_type;
typedef std::tuple<std::pair<uint64_t, uint64_t>, edge_data_type> edge_type; 

