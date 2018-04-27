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

// g++ -std=c++11 -O3 -fopenmp ../src/extract_2core_graph.cpp -o extract_2core_graph

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
#include <strstream>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr size_t max_vid = 1ULL << 20; //3563602788

int main(int argc, char **argv) {
  std::string table_path(argv[1]);
  std::string out_path(argv[2]);
  std::vector<std::string> edge_list_file;
  for (int i = 3; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  std::vector<bool> is_alive_table(max_vid + 1, false);
  std::cout << "Allocated table" << std::endl;
  {
    std::ifstream ifs(table_path);
    uint64_t vid;
    bool is_alive;
    int height;
    size_t num_child;
    uint64_t parent_id;

    while (ifs >> vid >> is_alive >> height >> num_child >> parent_id) {
      is_alive_table[vid] = is_alive;
    }
  }

  std::cout << "Constructed table" << std::endl;
  size_t count_dumped_edges = 0;
#pragma omp parallel for reduction(+:count_dumped_edges)
  for (size_t i = 0; i < edge_list_file.size(); ++i) {
    std::ifstream ifs(edge_list_file[i]);
    assert(ifs.is_open());

    auto const pos = edge_list_file[i].find_last_of('/');
    const auto fname = edge_list_file[i].substr(pos); // will return something "/xxx" if the input is "/aa/bb/xxx"
    std::ofstream ofs(out_path + fname);
    assert(ofs.is_open());
    std::cout << "Create " << out_path + fname << std::endl;

    uint64_t src;
    uint64_t dst;
    while (ifs >> src >> dst) {
      if (src > max_vid || dst > max_vid) {
        std::abort();
      }
      if (is_alive_table[src]) {
        ofs << src << " " << dst << "\n";
        ++count_dumped_edges;
      }
    }
    ofs.close();
  }

  std::cout << "count_dumped_edges: " << count_dumped_edges << std::endl;
  std::cout << "Done" << std::endl;

  return 0;
}