/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <unordered_map>
#include <limits>
#include <utility>
#include <algorithm>

constexpr size_t max_vid = 1ULL << 20;

using vid_t = uint32_t;
using graph_t = std::vector<std::vector<vid_t>>;
using bitmap_t = std::bitset<max_vid + 1>;

void bfs(const graph_t &graph, bitmap_t &visited) {
  vid_t root = 0;
  while (root < graph.size()) {
    if (!graph[root].empty()) break;
    ++root;
  }


  std::cout << "\nBFS Root: " << root << std::endl;
  bitmap_t frontier;
  bitmap_t next;
  frontier.set(root);
  visited.set(root);
  size_t level = 0;

  std::cout << "Level \t Frontiers" << std::endl;
  while (frontier.any()) {
    std::cout << level << " \t " << frontier.count() << std::endl;
    next.reset();
    for (uint64_t src = 0; src < frontier.size(); ++src) {
      if (!frontier.test(src)) continue;
      for (const auto dst : graph[src]) {
        if (!visited.test(dst)) next.set(dst);
      } // edge list
    } // frontier
    visited |= next;
    std::swap(frontier, next);
    next.reset();
    ++level;
  } // BFS loop
}

std::pair<size_t, size_t> cal_range(const size_t first, const size_t last, const size_t myid, const size_t nthreads) {
  size_t len = last - first + 1;
  size_t chunk = len / nthreads;
  size_t r = len % nthreads;

  size_t start;
  size_t end;

  if (myid < r) {
    start = first + (chunk + 1) * myid;
    end = start + chunk;
  } else {
    start = first + (chunk + 1) * r + chunk * (myid - r);
    end = start + chunk - 1;
  }

  return std::make_pair(start, end);
}

inline uint32_t hash32(uint32_t a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

#include <ctime>
void print_time() {
  std::time_t result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));
}


int main(int argc, char **argv) {
  std::string out_path(argv[1]);
  std::vector<std::string> edge_list_file;
  for (int i = 2; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  graph_t graph(max_vid + 1);

  size_t count_edges = 0;
  uint64_t actual_max_vid = 0;

  print_time();
  for (const auto &f : edge_list_file) {
    std::ifstream ifs(f);
    std::cout << "Open " << f << std::endl;
    if (!ifs.is_open()) std::abort();

    uint64_t src;
    uint64_t dst;
    while (ifs >> src >> dst) {
      if (src > max_vid || dst > max_vid) {
        std::abort();
      }
      graph[src].push_back(dst);
      graph[dst].push_back(src);
      actual_max_vid = std::max(actual_max_vid, src);
      actual_max_vid = std::max(actual_max_vid, dst);
      ++count_edges;
    }
  }
  std::cout << "Edge loading done" << std::endl;
  std::cout << "Actual max id: " << actual_max_vid << std::endl;
  std::cout << "Edges: " << count_edges << std::endl;
  print_time();

  for (auto &edge_list : graph) {
    std::sort(edge_list.begin(), edge_list.end());
    auto end = std::unique(edge_list.begin(), edge_list.end());
    edge_list.resize(std::distance(edge_list.begin(), end));
  }
  std::cout << "\nRemoving duplicated edges done" << std::endl;
  print_time();

  bitmap_t visited;
  bfs(graph, visited);
  print_time();

  std::cout << "\nDump edge" << std::endl;
  for (size_t i = 0; i < edge_list_file.size(); ++i) {
    auto const pos = edge_list_file[i].find_last_of('/');
    const auto fname = edge_list_file[i].substr(pos); // will return something "/xxx" if the input is "/aa/bb/xxx"
    std::ofstream ofs(out_path + fname);

    vid_t start;
    vid_t end;
    std::tie(start, end) = cal_range(0, actual_max_vid, i, edge_list_file.size());
    std::cout << "[ " << start << " ~ " << end << " ] : " << out_path + fname << std::endl;
    for (vid_t src = start; src <= end; ++src) {
      if (!visited.test(src)) continue;
      for (const auto &dst : graph[src]) {
        if (src == dst) continue;
        if (src != dst && hash32(src) == hash32(dst)) std::abort();
        if (hash32(src) < hash32(dst)) ofs << src << " " << dst << "\n";
        // if (src < dst) ofs << src << " " << dst << "\n";
      }
    }
    ofs.close();
  }

  std::cout << "Done" << std::endl;
  print_time();

  return 0;
}