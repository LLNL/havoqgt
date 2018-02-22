//
// Created by Iwabuchi, Keita on 2/22/18.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <algorithm>

using edge_list_t = std::vector<uint64_t>;
struct vertex_data_t
{
  uint16_t ecc;
  edge_list_t edge_list;
};
using graph_t = std::unordered_map<uint64_t, vertex_data_t>;

void read_exact_ecc_result(const std::string &vid_file, const std::string &ecc_file,
                           graph_t &graph)
{
  std::cout << "In " << __FUNCTION__ << std::endl;

  std::ifstream fin_vid(vid_file);
  std::ifstream fin_ecc(ecc_file);

  assert(fin_vid.is_open() && fin_ecc.is_open());

  while (true) {

    std::string buf1;
    std::string buf2;

    std::getline(fin_vid, buf1);
    std::getline(fin_ecc, buf2);

    if (buf1[0] == '\0' || buf2[0] == '\0') {
      assert(buf1[0] == '\0' && buf2[0] == '\0');
      break;
    }

    const uint64_t vid = std::stoull(buf1);
    const uint16_t ecc = std::stoi(buf2);


    assert(graph.count(vid) == 0);
    graph[vid].ecc = ecc;

    if (fin_vid.eof() || fin_ecc.eof()) {
      assert(fin_vid.eof() && fin_ecc.eof());
      break;
    }
  }
}

void read_edge_list(const std::vector<std::string> &edge_file_list, graph_t &graph)
{
  std::cout << "In " << __FUNCTION__ << std::endl;

  for (auto &edge_file : edge_file_list) {
    std::ifstream edge_list(edge_file);
    assert(edge_list.is_open());

    uint64_t src;
    uint64_t dst;
    while (edge_list >> src >> dst) {
      if (graph.count(src) == 0) continue; // skip unvisited vertices
      graph[src].edge_list.emplace_back(dst);
      assert(graph.count(dst) == 1);
      graph[dst].edge_list.emplace_back(src);
    }
  }
}

void remove_duplicated_edges(graph_t &graph)
{
  std::cout << "In " << __FUNCTION__ << std::endl;

  for (auto &vertex : graph) {
    auto& edge_list = vertex.second.edge_list;
    std::sort(edge_list.begin(), edge_list.end());
    auto last = std::unique(edge_list.begin(), edge_list.end());
    edge_list.erase(last, edge_list.end());
  }
}


void count_off_2_neighbor(graph_t &graph)
{
  std::cout << "In " << __FUNCTION__ << std::endl;

  size_t count = 0;
  for (auto &vertex : graph) {
    auto& edge_list = vertex.second.edge_list;
    size_t count_ninus_1 = 0;
    size_t count_plus_1 = 0;
    for (auto &trg : edge_list) {
      const int16_t diff = vertex.second.ecc - graph[trg].ecc;
      if (diff == 1) ++count_plus_1;
      if (diff == -1) ++count_ninus_1;
    }
    count = count_plus_1 * count_ninus_1;
  }

  std::cout << "count " << count << std::endl;
}

int main(int argc, char *argv[])
{
  std::string vid_file(argv[1]);
  std::string ecc_file(argv[2]);
  std::vector<std::string> edge_file_list;
  for (int i = 3; i < argc; ++i) {
    edge_file_list.emplace_back(std::string(argv[i]));
  }

  graph_t graph;

  read_exact_ecc_result(vid_file, ecc_file, graph);
  read_edge_list(edge_file_list, graph);
  remove_duplicated_edges(graph);

  count_off_2_neighbor(graph);

  return 0;
}