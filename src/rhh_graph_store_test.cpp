// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <metall/metall.hpp>

#include <havoqgt/rmat_edge_generator.hpp>
#include <rhh/graph_store.hpp>
#include <rhh/simple_stl_graph_store.hpp>

using rhh_graph_type =
    rhh::graph_store<uint64_t, double, bool,
                     metall::manager::allocator_type<std::byte>>;

using stl_graph_type =
    rhh::simple_stl_graph_store<uint64_t, double, bool,
                                metall::manager::allocator_type<std::byte>>;

template <typename graph_type, typename edge_type>
void run_test(const std::string& path, std::vector<edge_type> edges,
              const std::function<std::vector<edge_type>(const graph_type&)>&
                  dump_edges) {
  std::shuffle(edges.begin(), edges.end(), std::random_device());

  metall::manager manager(metall::create_only, path.c_str());
  graph_type      graph(manager.get_allocator<std::byte>());

  // Insert edges
  std::cout << "Insert Edges" << std::endl;
  for (const auto& edge : edges) {
    graph.insert_edge(edge.first, edge.second);
  }

  // Dump edges
  std::cout << "Dump Edges" << std::endl;
  auto dumped_edges = dump_edges(graph);

  std::cout << "Sort Generated Edges" << std::endl;
  std::sort(edges.begin(), edges.end());

  std::cout << "Sort Dumped Edges" << std::endl;
  std::sort(dumped_edges.begin(), dumped_edges.end());

  std::cout << "Compare Edges" << std::endl;
  if (edges != dumped_edges) {
    std::cerr << "Two edge lists do not match" << std::endl;
    std::abort();
  }
}

int main() {
  std::size_t  scale     = 20;
  std::size_t  num_edges = (1ULL << scale) * 16;
  unsigned int rnd_seed  = std::random_device()();
  std::string  path("/tmp/havoqgt_graph_store_test");

  havoqgt::rmat_edge_generator edge_generator(rnd_seed, scale, num_edges, 0.57,
                                              0.19, 0.19, 0.05, true, true);
  using edge_type = havoqgt::rmat_edge_generator::value_type;
  std::vector<edge_type> edges;
  for (auto itr = edge_generator.begin(), end = edge_generator.end();
       itr != end; ++itr) {
    edges.emplace_back(*itr);
  }

  std::cout << "\nTest RHH ver" << std::endl;
  auto dump_rhh_graph =
      [](const rhh_graph_type& graph) -> std::vector<edge_type> {
    std::vector<edge_type> dumped_edges;
    for (auto vit = graph.vertices_begin(), vend = graph.vertices_end();
         vit != vend; ++vit) {
      std::size_t degree = 0;
      for (auto eitr = graph.edges_begin(vit), eend = graph.edges_end(vit);
           eitr != eend; ++eitr) {
        dumped_edges.emplace_back(*vit, *eitr);
        ++degree;
      }

      if (degree != graph.degree(*vit)) {
        std::cerr << "At vertex " << *vit
                  << ", degree values are not the same: " << degree << ", "
                  << graph.degree(*vit) << std::endl;
        std::abort();
      }
    }
    return dumped_edges;
  };
  run_test<rhh_graph_type, edge_type>(path, edges, dump_rhh_graph);

  std::cout << "\nTest STL ver" << std::endl;
  auto dump_stl_graph =
      [](const stl_graph_type& graph) -> std::vector<edge_type> {
    std::vector<edge_type> dumped_edges;
    for (auto vit = graph.vertices_begin(), vend = graph.vertices_end();
         vit != vend; ++vit) {
      std::size_t degree = 0;
      for (auto eitr = graph.edges_begin(vit), eend = graph.edges_end(vit);
           eitr != eend; ++eitr) {
        dumped_edges.emplace_back(vit->first, eitr->first);
        ++degree;
      }

      if (degree != graph.degree(vit->first)) {
        std::cerr << "At vertex " << vit->first
                  << ", degree values are not the same: " << degree << ", "
                  << graph.degree(vit->first) << std::endl;
        std::abort();
      }
    }
    return dumped_edges;
  };
  run_test<stl_graph_type, edge_type>(path, edges, dump_stl_graph);

  std::cout << "Success!!" << std::endl;

  return 0;
}