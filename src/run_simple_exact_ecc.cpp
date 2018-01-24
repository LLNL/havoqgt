//
// Created by Iwabuchi, Keita on 1/22/18.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>


using id_type = uint64_t;
using level_type = uint16_t;
constexpr level_type k_unvisited_level = std::numeric_limits<level_type>::max();

using edge_list_type = std::vector<id_type>;
struct ecc_data_type {
  level_type lower{std::numeric_limits<level_type>::min()};
  level_type upper{std::numeric_limits<level_type>::max()};
};
struct vertex_data_type
{
  level_type level{k_unvisited_level};
  ecc_data_type ecc_data;
  edge_list_type edge_list;
};
using graph_type = std::unordered_map<id_type, vertex_data_type>;

size_t ingest_edge_list(const std::vector<std::string> &edge_list, graph_type &graph)
{
  size_t num_edges = 0;

  for (auto &edge_list_name : edge_list) {
    std::ifstream ifs(edge_list_name);
    assert(ifs.is_open());
    id_type source;
    id_type target;
    while (ifs >> source >> target) {
      graph[source].edge_list.emplace_back(target);
      ++num_edges;
    }
  }
  return num_edges;
}

void reset_bfs_data(graph_type &graph)
{
  for (auto &vertex : graph) {
    vertex.second.level = k_unvisited_level;
  }
}

level_type bfs(const id_type root, graph_type &graph)
{
  std::unordered_set<id_type> queue;
  queue.insert(root);
  graph.at(root).level = 0;

  std::cout << "Level\t" << "# vertices" << std::endl;
  level_type level = 0;
  while (!queue.empty()) {
    std::cout << level << " " << queue.size() << std::endl;
    std::unordered_set<id_type> next;
    for (const id_type source : queue) {
      for (const id_type target : graph.at(source).edge_list) {
        if (graph.at(target).level == k_unvisited_level) {
          graph.at(target).level = level + 1;
          next.insert(target);
        }
      } // Loop over a vertex
    } // Loop over queue
    queue = std::move(next);
    ++level;
  }

  return level - 1;
}

size_t bound_ecc(const id_type source, const level_type exact_ecc, graph_type &graph)
{
  graph.at(source).ecc_data.lower = graph.at(source).ecc_data.upper = exact_ecc;

  size_t num_bounded = 0;
  for (auto &vertex : graph) {
    vertex_data_type& vertex_data = vertex.second;
    ecc_data_type& ecc_data = vertex_data.ecc_data;
    if (vertex_data.level == k_unvisited_level) continue;
    if (ecc_data.lower == ecc_data.upper) continue;

    ecc_data.lower = std::max(ecc_data.lower, std::max(vertex_data.level, static_cast<uint16_t>(exact_ecc - vertex_data.level)));
    ecc_data.upper = std::min(ecc_data.upper, static_cast<uint16_t>(exact_ecc + vertex_data.level));

    if (ecc_data.lower == ecc_data.upper) ++num_bounded;
  }

  return num_bounded;
}

int main (int argc, char *argv[])
{
  std::vector<id_type> source_list;
  {
    std::stringstream sstrm(argv[1]);
    std::string buf;
    while (std::getline(sstrm, buf, ':')) {
      source_list.push_back(std::stoull(buf));
    }
  }

  std::vector<std::string> edge_list;
  for (int i = 2; i < argc; ++i) {
    edge_list.emplace_back(std::string(argv[i]));
  }

  graph_type graph;

  const size_t num_edges = ingest_edge_list(edge_list, graph);
  std::cout << "# vertices: " << graph.size() << std::endl;
  std::cout << "# edges: " << num_edges << std::endl;

  for (const id_type source : source_list) {
    std::cout << "\nSource: " << source << std::endl;
    reset_bfs_data(graph);
    const level_type exact_ecc = bfs(source, graph);
    const size_t num_bounded = bound_ecc(source, exact_ecc, graph);
    std::cout << "num_bounded: " << num_bounded << std::endl;
  }

  return 0;
}