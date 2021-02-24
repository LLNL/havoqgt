// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <fstream>
#include <iostream>
#include <map>
#include <stack>
#include <string>
#include <vector>

using graph_t = std::map<uint64_t, std::vector<uint64_t>>;

void ingest_edges(const std::vector<std::string>& edgelists,
                  graph_t* const                  graph) {
  for (const auto& ef : edgelists) {
    std::ifstream ifs(ef);
    if (!ifs.is_open()) {
      std::cerr << "Failed to open " << ef << std::endl;
      std::abort();
    }

    uint64_t src;
    uint64_t dst;
    while (ifs >> src >> dst) {
      (*graph)[src].push_back(dst);
      (*graph)[dst].push_back(src);
    }
  }
}

int main() {
  graph_t graph;

  std::vector<std::string> edgelists;
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-0");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-1");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-2");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-3");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-4");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-5");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-6");
  edgelists.emplace_back("../test/datasets/edge_list_rmat_s17-7");
  ingest_edges(edgelists, &graph);

  std::cout << "Start CC" << std::endl;

  std::map<uint64_t, uint64_t> cc_data;
  for (const auto& root_data : graph) {  // Assumes that graph uses std::map
    const uint64_t       root = root_data.first;
    if (cc_data.count(root) > 0) continue; // Already visited
    std::stack<uint64_t> st;
    st.push(root);
    cc_data[root] = root;
    while (!st.empty()) {
      const auto src = st.top();
      st.pop();
      for (const auto& dst : graph[src]) {
        if (cc_data.count(dst) == 0) { // Not visited yet
          st.push(dst);
          cc_data[dst] = root;
        }
      }
    }
  }

  std::cout << "Write CC data" << std::endl;
  std::ofstream ofs("cc_data");
  const uint64_t max_vid = (--(cc_data.end()))->first;
  for (uint64_t vid = 0; vid <= max_vid; ++vid) {
    if (cc_data.count(vid) == 0)
      ofs << vid << " " << vid << "\n";
    else
      ofs << vid << " " << cc_data.at(vid) << "\n";
  }
  ofs.close();

  std::cout << "Finished" << std::endl;

  return 0;
}