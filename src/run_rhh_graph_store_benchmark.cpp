// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT


/*
Benchmark results:

------
Flash, /dev/shm
------
STL ver
#of inserted edges      134217856
Total time (s)  176.268
#of inserted edges/s    761442

RHH ver
#of inserted edges      134217856
Total time (s)  94.673
#of inserted edges/s    1.4177e+06
*/

#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <havoqgt/rmat_edge_generator.hpp>
#include <metall/metall.hpp>
#include <rhh/graph_store.hpp>
#include <rhh/simple_stl_graph_store.hpp>

using rhh_graph_type =
    rhh::graph_store<uint64_t, double, bool,
                     metall::manager::allocator_type<std::byte>>;

using stl_graph_type =
    rhh::simple_stl_graph_store<uint64_t, double, bool,
                                metall::manager::allocator_type<std::byte>>;

inline std::chrono::high_resolution_clock::time_point elapsed_time_sec() {
  return std::chrono::high_resolution_clock::now();
}

inline double elapsed_time_sec(
    const std::chrono::high_resolution_clock::time_point& tic) {
  auto duration_time = std::chrono::high_resolution_clock::now() - tic;
  return static_cast<double>(
      std::chrono::duration_cast<std::chrono::microseconds>(duration_time)
          .count() /
      1e6);
}

template <typename edge_generator_t, typename graph_type>
void run_bench_core(const std::size_t chunk_size, const unsigned int rnd_seed,
               edge_generator_t* edge_generator, graph_type* graph) {
  using edge_type = typename edge_generator_t::value_type;
  std::vector<edge_type> edge_buf;
  std::size_t            num_inserted = 0;
  double                 total_time   = 0;

  auto edge_itr = edge_generator->begin();
  auto edge_end = edge_generator->end();

  std::mt19937_64 rnd_generator(rnd_seed);

  while (true) {
    // Generate edges
    {
      edge_buf.clear();
      while (edge_itr != edge_end) {
        edge_buf.push_back(*edge_itr);
        if (edge_buf.size() == chunk_size) break;
        ++edge_itr;
      }
      if (edge_buf.empty()) break;
      std::shuffle(edge_buf.begin(), edge_buf.end(), rnd_generator);
    }

    // Insert edges
    {
      const auto start_time = elapsed_time_sec();
      for (const auto& edge : edge_buf) {
        graph->insert_edge(edge.first, edge.second);
      }
      const auto elapsed_time = elapsed_time_sec(start_time);
      total_time += elapsed_time;
      num_inserted += edge_buf.size();
    }
  }
  std::cout << "#of inserted edges\t" << num_inserted << std::endl;
  std::cout << "Total time (s)\t" << total_time << std::endl;
  std::cout << "#of inserted edges/s\t" << num_inserted / total_time
            << std::endl;
}

void dump_edges(const rhh_graph_type& graph) {
  for (auto vit = graph.vertices_begin(), vend = graph.vertices_end();
       vit != vend; ++vit) {
    for (auto eitr = graph.edges_begin(vit), eend = graph.edges_end(vit);
         eitr != eend; ++eitr) {
      std::cout << *vit << "\t" << *eitr << std::endl;
    }
  }
}

void dump_edges(const stl_graph_type& graph) {
  for (auto vit = graph.vertices_begin(), vend = graph.vertices_end();
       vit != vend; ++vit) {
    for (auto eitr = graph.edges_begin(vit), eend = graph.edges_end(vit);
         eitr != eend; ++eitr) {
      std::cout << vit->first << "\t" << eitr->first << std::endl;
    }
  }
}

template <typename graph_type>
void run_bench(const std::size_t scale, const std::size_t num_edges,
               const std::size_t chunk_size, const unsigned int rnd_seed,
               const std::string& path) {
  metall::manager manager(metall::create_only, path.c_str());

  graph_type graph(manager.get_allocator<std::byte>());

  havoqgt::rmat_edge_generator edge_generator(rnd_seed, scale, num_edges, 0.57,
                                              0.19, 0.19, 0.05, true, true);

  run_bench_core(chunk_size, rnd_seed, &edge_generator, &graph);
  // dump_edges(graph);
}

int main() {
  std::size_t  scale      = 22;
  std::size_t  num_edges  = (1ULL << scale) * 16;
  std::size_t  chunk_size = 1ULL << 20ULL;
  unsigned int rnd_seed   = std::random_device()();

  std::cout << "\nSTL ver" << std::endl;
  run_bench<stl_graph_type>(scale, num_edges, chunk_size, rnd_seed, "/tmp/test");

  std::cout << "\nRHH ver" << std::endl;
  run_bench<rhh_graph_type>(scale, num_edges, chunk_size, rnd_seed, "/tmp/test");

  return 0;
}