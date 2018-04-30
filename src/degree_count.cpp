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

// g++ -std=c++11 -O3 -fopenmp ../src/degree_count.cpp -o degree_count

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <strstream>

constexpr size_t max_vid = 3563602788;

int main(int argc, char **argv) {

  std::string out_path(argv[1]);
  std::vector<std::string> edge_list_file;
  for (int i = 2; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  std::vector<uint64_t> degree_table(max_vid + 1, 0);
  {
#pragma omp parallel for
    for (size_t i = 0; i < edge_list_file.size(); ++i) {
      const auto &f = edge_list_file[i];
      std::ifstream ifs(f);
      std::cout << "Open " << f << std::endl;
      if (!ifs.is_open()) std::abort();

      uint64_t src;
      uint64_t dst;
      while(ifs >> src >> dst) {
#pragma omp atomic
        ++degree_table[src];
#pragma omp atomic
        ++degree_table[dst];
      }
    }
  }
  std::cout << "Count done" << std::endl;

  std::ofstream ofs(out_path);
  for (uint64_t i = 0; i < degree_table.size(); ++i) {
    if (degree_table[i] == 0) continue;
    ofs << i << "\t" << degree_table[i] << "\n";
  }
  ofs.close();

  std::cout << "Done" << std::endl;

  return 0;
}