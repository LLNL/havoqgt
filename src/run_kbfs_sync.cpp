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

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/fixed_size_unordered_map.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <havoqgt/distributed_db.hpp>
#include <assert.h>

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>
#include <unordered_map>

#include <boost/interprocess/managed_heap_memory.hpp>

using namespace havoqgt;

#ifdef NUM_SOURCES
constexpr uint32_t k_num_sources = NUM_SOURCES;
#else
constexpr uint32_t k_num_sources = 8;
#endif
using level_t = uint16_t;

using graph_t = havoqgt::delegate_partitioned_graph<distributed_db::allocator<>>;
using kbfs_t = havoqgt::k_breadth_first_search<distributed_db::allocator<>, level_t, k_num_sources>;
// using kbfs_vertex_data_t = typename havoqgt::kbfs_t::vertex_data_t;


void usage()  {
  if(comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -i <string>   - input graph base filename (required)\n"
              << " -b <string>   - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
              << " -s <int>      - Source vertex of BFS (Default is 0)\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename, std::string& backup_filename, std::vector<uint64_t>& source_vertex_list) {
  if(comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;
  source_vertex_list.clear();

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:b:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's':
      {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_vertex_list.push_back(atoll(buf.c_str()));
        break;
      }
      case 'i':
        found_input_filename = true;
        input_filename = optarg;
        break;
      case 'b':
        backup_filename = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_input_filename || (source_vertex_list.size() != k_num_sources)) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {


  int mpi_rank(0), mpi_size(0);


  havoqgt::init(&argc, &argv);
  {
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      std::cout << "k_num_sources " << k_num_sources << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    std::string graph_input;
    std::string backup_filename;
    std::vector<uint64_t> source_vertex_list;

    parse_cmd_line(argc, argv, graph_input, backup_filename, source_vertex_list);

    if (source_vertex_list.empty()) {
      source_vertex_list.push_back(0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    distributed_db ddb(db_open_read_only(), graph_input);

    auto graph = ddb.get_manager()->find<graph_t>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    //graph->print_graph_statistics();
    MPI_Barrier(MPI_COMM_WORLD);


    // BFS Experiments
    {
      kbfs_t kbfs(*graph);
      MPI_Barrier(MPI_COMM_WORLD);

      //  Run BFS experiments
      double time(0);
      int count(0);
      std::vector<graph_t::vertex_locator> source_list;
      for (auto& source_vertex : source_vertex_list) {
        uint64_t isource = source_vertex;
        graph_t::vertex_locator source = graph->label_to_locator(isource);
        uint64_t global_degree(0);
        do {
          uint64_t local_degree = 0;
          source = graph->label_to_locator(isource);
          if (source.is_delegate()) {
            break;
          }
          if (uint32_t(mpi_rank) == source.owner()) {
            local_degree = graph->degree(source);
          }
          global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(),
                                         MPI_COMM_WORLD);
          if (global_degree == 0) ++isource;
        } while (global_degree == 0);
        if (uint32_t(mpi_rank) == source.owner()) {
          if (isource != source_vertex) {
            std::cout << "\nVertex " << source_vertex << " has a degree of 0.   New source vertex = " << isource
                      << std::endl;
          } else {
            std::cout << "\nStarting vertex = " << isource << std::endl;
          }
          std::cout << "delegate? = " << source.is_delegate() << std::endl;
          std::cout << "local_id = " << source.local_id() << std::endl;
          std::cout << "degree = " << graph->degree(source) << std::endl;
        }
        source_vertex = isource;
        source_list.push_back(source);
      }

      MPI_Barrier(MPI_COMM_WORLD);
      double time_start = MPI_Wtime();
      kbfs.run(source_list);
      MPI_Barrier(MPI_COMM_WORLD);
      double time_end = MPI_Wtime();


      if (mpi_rank == 0) {
        std::cout << "kBFS_time\t" << time_end - time_start << std::endl;
      }

    }  // End BFS
  }  // END Main MPI

  return 0;
}

