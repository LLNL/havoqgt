//
// Created by Iwabuchi, Keita on 4/7/18.
//

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

// g++ -std=c++11 -O3 -fopenmp ../src/run_2core.cpp -o run_2core

#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <unordered_map>
#include <limits>
#include <utility>
#include <algorithm>
#include <cassert>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr size_t max_vid = 1 << 20; //3563602788

using vid_t = uint32_t;
// using graph_t = std::vector<std::vector<vid_t>>;
using bitmap_t = std::bitset<max_vid + 1>;

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

template <size_t max_id, size_t num_bank = 256>
class graph_bank {
 public:
  using edge_list_t = std::vector<vid_t>;
  using graph_t = std::vector<edge_list_t>;
  using table_t = std::vector<graph_t>;
  static constexpr size_t k_each_size = ((max_id + 1) + num_bank - 1) / num_bank;

  graph_bank()
      : m_offset(num_bank),
        m_table(num_bank),
        m_locks() {
    for (size_t i = 0; i < num_bank; ++i) {
      m_offset[i] = i * k_each_size;
      m_table[i].resize(k_each_size);
      omp_init_lock(&m_locks[i]);
    }
  }

  edge_list_t &operator[](const uint64_t i) {
    return const_cast<edge_list_t &>(const_cast<const graph_bank &>(*this)[i]);
  }

  const edge_list_t &operator[](const uint64_t i) const {
    const uint64_t idx = index(i);
    return m_table[idx][i - m_offset[idx]];
  }

  void add(const vid_t src, const vid_t dst) {
    {
      const uint64_t idx = index(src);
      while (!omp_test_lock(&m_locks[idx]));
      m_table[idx][src - m_offset[idx]].push_back(dst);
      omp_unset_lock(&m_locks[idx]);
    }
  }

  static uint64_t index(const uint64_t i) {
    assert(i / k_each_size < num_bank);
    return i / k_each_size;
  }

  size_t size() const {
    return (max_id + 1);
  }

 private:
  std::vector<uint64_t> m_offset;
  table_t m_table;
  std::array<omp_lock_t, num_bank> m_locks;
};

size_t compute_2_core(const graph_bank<max_vid, 256> &graph,
                      const size_t actual_max_vid,
                      std::vector<uint64_t> &left_degree,
                      std::vector<uint16_t> &max_depth,
                      std::vector<uint64_t> &num_kill)
{
  size_t total_num_killed = 0;
  uint16_t depth = 1;
  std::cout << "depth num_killed" << std::endl;
  while (true) {
    bool cont = false;
    size_t num_killed = 0;
    for (vid_t vid = 0; vid < actual_max_vid + 1; ++vid) {
      if (left_degree[vid] == 1) {
        --left_degree[vid];
        ++num_killed;
        for (vid_t trg : graph[vid]) {
          if (left_degree[trg] > 0) {
            --left_degree[trg];
            num_kill[trg] += num_kill[vid] + 1;
            max_depth[trg] = depth;
          }
        }
        cont = true;
      }
    }
    if (!cont) break;
    std::cout << depth << " " << num_killed << std::endl;
    ++depth;
    total_num_killed += num_killed;
  }

  return total_num_killed;
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

void find_parent(const graph_bank<max_vid, 256> &graph,
                 const vid_t actual_max_vid,
                 const std::vector<uint64_t> &left_degree,
                 const std::vector<uint16_t> &max_depth,
                 const std::vector<uint64_t> &num_kill,
                 std::vector<vid_t> &parent)
{
  bitmap_t frontier;
  bitmap_t next;

  for (vid_t vid = 0; vid < actual_max_vid + 1; ++vid) {
    if (left_degree[vid] > 0 && num_kill[vid] > 0) {
      frontier.set(vid);
      parent[vid] = vid;
    }
  }

  while (frontier.any()) {
    next.reset();
    for (uint64_t src = 0; src < frontier.size(); ++src) {
      if (!frontier.test(src)) continue;
      for (const auto dst : graph[src]) {
        if (left_degree[dst] == 0 && max_depth[dst] < max_depth[src]) {
          next.set(dst);
          parent[dst] = parent[src];
        }
      } // edge list
    } // frontier
    std::swap(frontier, next);
    next.reset();
  } // BFS loop
}

int main(int argc, char **argv) {
  std::string out_path(argv[1]);
  std::vector<std::string> edge_list_file;
  for (int i = 2; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  size_t count_edges = 0;
  uint64_t actual_max_vid = 0;
  print_time();

  std::vector<uint64_t> num_edges(max_vid + 1, 0);
  {
#pragma omp parallel for
    for (size_t i = 0; i < edge_list_file.size(); ++i) {
      const auto &f = edge_list_file[i];
      std::ifstream ifs(f);
      std::cout << "Open " << f << std::endl;
      if (!ifs.is_open()) std::abort();

      uint64_t src;
      uint64_t dst;
      while (ifs >> src >> dst) {
        if (src > max_vid || dst > max_vid) {
          std::abort();
        }
#pragma omp atomic
        ++num_edges[src];
#pragma omp atomic
        ++num_edges[dst];
      }
    }
  }
  std::cout << "Read edge: " << std::accumulate(num_edges.cbegin(), num_edges.cend(), 0ULL) << std::endl;
  print_time();

  std::cout << "Alloc graph" << std::endl;
  graph_bank<max_vid, 256> graph;

  std::cout << "Extend edge list" << std::endl;
#pragma omp parallel for
  for (uint64_t i = 0; i < graph.size(); ++i) {
    graph[i].reserve(num_edges[i]);
  }
  print_time();

#pragma omp parallel for reduction(+:count_edges), reduction(max:actual_max_vid)
  for (size_t i = 0; i < edge_list_file.size(); ++i) {
    const auto &f = edge_list_file[i];
    std::ifstream ifs(f);
    std::cout << "Open " << f << std::endl;
    if (!ifs.is_open()) std::abort();

    uint64_t src;
    uint64_t dst;
    while (ifs >> src >> dst) {
      if (src > max_vid || dst > max_vid) {
        std::abort();
      }

      graph.add(src, dst);
      graph.add(dst, src);

      actual_max_vid = std::max(actual_max_vid, src);
      actual_max_vid = std::max(actual_max_vid, dst);
      ++count_edges;
    }
  }
  assert(count_edges * 2 == std::accumulate(num_edges.begin(), num_edges.end(), 0ULL));
  std::cout << "Edge loading done" << std::endl;
  std::cout << "Actual max id: " << actual_max_vid << std::endl;
  std::cout << "Edges: " << count_edges << std::endl;
  num_edges.resize(0);
  print_time();

  std::vector<uint64_t> degree(actual_max_vid + 1, 0);
  std::vector<uint16_t> max_depth(actual_max_vid + 1, 0);
  std::vector<uint64_t> num_kill(actual_max_vid + 1, 0);
  for (vid_t vid = 0; vid < actual_max_vid + 1; ++vid) {
    degree[vid] = graph[vid].size();
  }
  const size_t num_killed = compute_2_core(graph, actual_max_vid, degree, max_depth, num_kill);
  std::cout << "\n2 core num killed: " << num_killed << std::endl;
  print_time();

  std::vector<vid_t> parent(actual_max_vid + 1, std::numeric_limits<vid_t>::max());
  find_parent(graph, actual_max_vid, degree, max_depth, num_kill, parent);
  print_time();

  std::cout << "\nDump info" << std::endl;
  std::ofstream ofs(out_path);
  for (vid_t vid = 0; vid < actual_max_vid + 1; ++vid) {
    ofs << vid << " " << (degree[vid] > 0) << " " << max_depth[vid] << " " << num_kill[vid] << " " << parent[vid] << "\n";
  }

  std::cout << "Done" << std::endl;
  print_time();

  return 0;
}