/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. 
 * Written by Roger Pearce <rpearce@llnl.gov>. 
 * LLNL-CODE-644630. 
 * All rights reserved.
 * 
 * This file is part of HavoqGT, Version 1.0. 
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


#include <external_memory_arena.hpp>
#include <adjacency_list_graph.hpp>
#include <breadth_first_search.hpp>
#include <iostream>

using namespace havoqgt::omp;

int main(int argc, char** argv) {

  if(argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <filename> <scale>" <<std::endl;
    return 1;
  }

  havoqgt::omp::init_environment();

  const char* fname = argv[1];
  int SCALE = atoi(argv[2]);
  std::cout << "Graph: " << fname << ", scale = " << SCALE << std::endl;

  external_memory_arena em(fname);

  typedef adjacency_list_graph<external_memory_arena> graph_type;
  graph_type graph(em, "rmat_adj_list");

  for(uint64_t i=0, raw_source=0; i<10; ++i) {
    //std::cout << "BFS from source " << i << std::endl;
    graph_type::vertex_data<uint8_t> bfs_level(&graph, 255);
    graph_type::vertex_descriptor source;
    do {
      source.local_id = raw_source++;
      source.thread_id = 0;
    } while(graph.num_edges(source) == 0);
    double time_start = omp_get_wtime();
    breadth_first_search(graph, source, bfs_level);
    double time_end = omp_get_wtime();

    uint64_t level_count = 0;
    uint64_t level = 0;
    do {
      level_count = bfs_level.count_equal(level);
      std::cout << "BFS Level " << level << " has " << level_count << " vertices." << std::endl;
      ++level;
    } while(level_count > 0);
    if(level > 1) {
      std::cout << "TEPS = " << (double(uint64_t(1) << uint64_t(SCALE)) * 16) / (time_end - time_start) << std::endl;
    }
  }

};
