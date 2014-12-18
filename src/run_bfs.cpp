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

#include <havoqgt/page_rank.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/upper_triangle_edge_generator.hpp>
#include <havoqgt/single_source_shortest_path.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <havoqgt/distributed_db.hpp>
#include <assert.h>

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>

#include <boost/interprocess/managed_heap_memory.hpp>

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

int main(int argc, char** argv) {
  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
    havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  std::string graph_input;

  if (argc < 2) {
    std::cerr << "usage: <graph input file name>"
      << " (argc:" << argc << " )." << std::endl;
    exit(-1);
  } else {
    int pos = 1;
    graph_input = argv[pos++];
  }

//  graph_input += "_" + std::to_string(mpi_rank);
//  if (mpi_rank == 0) {
//    std::cout << "[0]Graph input file = " << graph_input << std::endl;
//  }


  MPI_Barrier(MPI_COMM_WORLD);

  havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
  // graph_mapped_t graph_mapped_file(bip::open_read_only, graph_input.c_str());

  // boost::interprocess::mapped_region::advice_types rand_advice;
  // rand_advice = boost::interprocess::mapped_region::advice_types::advice_random;
  // bool assert_res = graph_mapped_file.advise(rand_advice);
  // assert(assert_res);


  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }
  //graph->print_graph_statistics();
  MPI_Barrier(MPI_COMM_WORLD);

  //
  // Calculate max degree
  /*uint64_t max_degree = 0;//graph->max_vertex_id();
  for (graph_type::vertex_iterator vitr = graph->vertices_begin();
        vitr != graph->vertices_end(); ++vitr) {
    max_degree = std::max(max_degree, graph->degree(*vitr));
  }
  for (graph_type::controller_iterator citr = graph->controller_begin();
        citr != graph->controller_end(); ++citr) {
    max_degree = std::max(max_degree, graph->degree(*citr));
  }
  uint64_t global_max_degree = havoqgt::mpi::mpi_all_reduce(max_degree,
        std::greater<uint64_t>(), MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "Max Local Degree = " << global_max_degree << std::endl;
  }
*/
  MPI_Barrier(MPI_COMM_WORLD);
 
  //   NEED TO FIX LOCAL DATA!!
  // BFS Experiments
  {
    #if 0
    typedef bip::managed_mapped_file bfs_mapped_t;
    uint64_t file_size = (173249084000.0/24.0);

    std::string bfs_filename =  graph_input + std::to_string(mpi_rank) + "_bfs";
    bfs_mapped_t bfs_mapped_data(bip::create_only, bfs_filename.c_str(),
        file_size);
    // assert_res = bfs_mapped_data.advise(rand_advice);
    // assert(assert_res);

    #else

    typedef bip::managed_heap_memory bfs_mapped_t;
    //uint64_t filesize = (21474836480/24.0);
    //uint64_t filesize = (4973120026ULL);
    uint64_t filesize = 256*1024*1024;
    bfs_mapped_t bfs_mapped_data(filesize);
    #endif


    typedef bfs_mapped_t::segment_manager bfs_segment_manager_t;



    graph_type::vertex_data<uint8_t, bfs_segment_manager_t >* bfs_level_data;
    graph_type::vertex_data<graph_type::vertex_locator, bfs_segment_manager_t >* bfs_parent_data;

    bfs_level_data = graph->create_vertex_data<uint8_t, bfs_segment_manager_t>(
          bfs_mapped_data.get_segment_manager(), "bfs_level_data");
    bfs_parent_data = graph->create_vertex_data<graph_type::vertex_locator, bfs_segment_manager_t>(
          bfs_mapped_data.get_segment_manager(), "bfs_parent_data");
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "BFS data allocated." << std::endl;
  }

    //  Run BFS experiments
    double time(0);
    int count(0);
    uint64_t isource = 0;
    for (int nrbfs = 0; nrbfs < 16; ++nrbfs) {
      graph_type::vertex_locator source = graph->label_to_locator(isource);
      uint64_t global_degree(0);
      do {
        uint64_t local_degree = 0;
        source = graph->label_to_locator(isource++);
        if (source.is_delegate()) {
          break;
        }
        if (uint32_t(mpi_rank) == source.owner()) {
          local_degree = graph->degree(source);
        }
        global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(),
            MPI_COMM_WORLD);
      } while (global_degree == 0);
      if (uint32_t(mpi_rank) == source.owner()) {
        std::cout << "Starting vertex = " << isource << std::endl;
        std::cout << "delegate? = " << source.is_delegate() << std::endl;
        std::cout << "local_id = " << source.local_id() << std::endl;
        std::cout << "degree = " << graph->degree(source) << std::endl;
      }

      bfs_level_data->reset(128);

      MPI_Barrier(MPI_COMM_WORLD);
      double time_start = MPI_Wtime();
      hmpi::breadth_first_search(graph, *bfs_level_data, *bfs_parent_data,
          source);
      MPI_Barrier(MPI_COMM_WORLD);
      double time_end = MPI_Wtime();

      uint64_t visited_total(0);
      for (uint64_t level = 0; level < 15; ++level) {
        uint64_t local_count(0);
        graph_type::vertex_iterator vitr;
        for (vitr = graph->vertices_begin();
             vitr != graph->vertices_end();
             ++vitr) {
          if ((*bfs_level_data)[*vitr] == level) {
            ++local_count;
          }
        }

        // Count the controllers!
        graph_type::controller_iterator citr;
        for (citr = graph->controller_begin();
             citr != graph->controller_end();
             ++citr) {
          if ((*bfs_level_data)[*citr] == level) {
            ++local_count;
          }
        }

        uint64_t global_count = mpi_all_reduce(local_count,
          std::plus<uint64_t>(), MPI_COMM_WORLD);
        visited_total += global_count;
        if (mpi_rank == 0 && global_count > 0) {
          std::cout << "Level " << level << ": " << global_count << std::endl;
        }
        if (global_count == 0) {
          break;
        }
      }  // end for level

      if (mpi_rank == 0) {
        if (visited_total > 1000) {
          std::cout
            << "Visited total = " << visited_total
            << ", percentage visited = "
            << double(visited_total) / double(graph->max_global_vertex_id()) * 100
            << "%" << std::endl
            << "BFS Time = " << time_end - time_start << std::endl;
          time += time_end - time_start;
          ++count;
        }
      }
    }
    if (mpi_rank == 0) {
      std::cout << "Count BFS = " << count << std::endl;
      std::cout << "AVERAGE BFS = " << time / double(count) << std::endl;
    }
  }  // End BFS Test
  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  if (mpi_rank == 0) {
    std::cout << "FIN." << std::endl;
    //print_system_info(false);
    // print_dmesg();
  }
  return 0;
}
