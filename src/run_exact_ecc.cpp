//
// Created by Iwabuchi, Keita on 12/14/17.
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

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <chrono>
#include <random>
#include <assert.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/fixed_size_unordered_map.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/exact_eccentricity.hpp>

// -------------------------------------------------------------------------------------------------------------- //
// Types
// -------------------------------------------------------------------------------------------------------------- //
using namespace havoqgt;

using graph_t = delegate_partitioned_graph<distributed_db::allocator<>>;

#ifdef NUM_SOURCES
constexpr uint32_t k_num_sources = NUM_SOURCES;
#else
constexpr uint32_t k_num_sources = 8;
#endif

using level_t = uint16_t;
using exact_eccentricity_t = exact_eccentricity<distributed_db::allocator<>, level_t, k_num_sources>;

// -------------------------------------------------------------------------------------------------------------- //
// Parse command line
// -------------------------------------------------------------------------------------------------------------- //
void usage()
{
  if (comm_world().rank() == 0) {
    std::cout << "Usage: env [USE_TAKE=1, USE_DS_FIX=1, USE_DS_ADP=1, or USE_DS_FIX_MS=1] ./run_exact_ecc -i <string>\n"
              << " -i <string>    - input graph base filename (required)\n"
              << " -b <string>    - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
              << " -a <string>    - colon separated source selection algorithms IDs for USE_DS_AD\n"
              << " -h             - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char **argv, std::string &input_filename,
                    std::string &backup_filename,
                    std::string &ecc_output_filename,
                    std::set<int>& use_algorithm,
                    std::string& two_core_info_filename)
{
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:b:e:a:c:h")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;

      case 'i':
        found_input_filename = true;
        input_filename = optarg;
        break;

      case 'b':
        backup_filename = optarg;
        break;

      case 'e':
        ecc_output_filename = optarg;
        break;

      case 'a': {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          use_algorithm.insert(std::stoi(buf.c_str()));
        break;
      }


      case 'c':
        two_core_info_filename = optarg;
        break;

      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_input_filename) {
    usage();
    exit(-1);
  }
}


// -------------------------------------------------------------------------------------------------------------- //
// select_non_zero_degree_source
// -------------------------------------------------------------------------------------------------------------- //
std::vector<typename graph_t::vertex_locator>
select_non_zero_degree_source(const graph_t *const graph, std::vector<uint64_t> &source_id_list)
{
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector<typename graph_t::vertex_locator> vertex_locator_list;
  for (auto id_itr = source_id_list.begin(), end = source_id_list.end(); id_itr != end; ++id_itr) {
    uint64_t source_id = *id_itr;
    typename graph_t::vertex_locator source = graph->label_to_locator(source_id);
    uint64_t global_degree(0);
    while (true) {
      uint64_t local_degree = 0;
      source = graph->label_to_locator(source_id);
      if (source.is_delegate()) {
        break;
      }
      if (uint32_t(mpi_rank) == source.owner()) {
        local_degree = graph->degree(source);
      }
      global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
      const bool already_exist = (std::find(source_id_list.begin(), id_itr, source_id) != id_itr);
      if (global_degree == 0 || already_exist) ++source_id;
      else break;
    }
    if (uint32_t(mpi_rank) == source.owner()) {
      if (source_id != *id_itr) {
        std::cout << "\nVertex " << *id_itr << " has a degree of 0.   New source vertex = " << source_id << std::endl;
      } else {
        std::cout << "\nStarting vertex = " << source_id << std::endl;
      }
      std::cout << "delegate? = " << source.is_delegate() << std::endl;
      std::cout << "local_id = " << source.local_id() << std::endl;
      std::cout << "degree = " << graph->degree(source) << std::endl;
    }
    vertex_locator_list.emplace_back(std::move(source));
    *id_itr = source_id;
  }

  return vertex_locator_list;
}

// -------------------------------------------------------------------------------------------------------------- //
// Main
// -------------------------------------------------------------------------------------------------------------- //
int main(int argc, char **argv)
{

  int mpi_rank(0), mpi_size(0);


  havoqgt::init(&argc, &argv);
  {
    // -------------------------------------------------------------------------------------------------------------- //
    //                                            Parse options & Prepare graph
    // -------------------------------------------------------------------------------------------------------------- //
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      std::cout << "k_num_sources: " << k_num_sources << std::endl;
      //print_system_info(false);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_filename;
    std::string ecc_output_filename;
    std::set<int> use_algorithm;
    std::string two_core_info_filename;

    parse_cmd_line(argc, argv, graph_input, backup_filename, ecc_output_filename, use_algorithm, two_core_info_filename);

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

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Touch graph
    // -------------------------------------------------------------------------------------------------------------- //
    {
      if (mpi_rank == 0) std::cout << "Touching graph" << std::endl;
      size_t sum = 0;
      for (auto itr = graph->vertices_begin(), end = graph->vertices_end();
          itr != end; ++itr) {
        sum += graph->degree(*itr);
      }
      for (auto itr = graph->controller_begin(), end = graph->controller_end();
           itr != end; ++itr) {
        sum += graph->degree(*itr);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0) std::cout << "Touching graph done" << std::endl;
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Compute exact ecc
    // -------------------------------------------------------------------------------------------------------------- //
    exact_eccentricity_t eeec(*graph, use_algorithm);
    eeec.load_2core_info(two_core_info_filename);

    const double total_start_time = MPI_Wtime();
    eeec.run();
    const double total_end_time = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Total execution time: " << total_end_time - total_start_time << std::endl;

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Compute diameter
    // -------------------------------------------------------------------------------------------------------------- //
    const uint16_t diameter = eeec.find_max_ecc();
    if (mpi_rank == 0) std::cout << "Diameter: " << diameter << std::endl;

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Dump exact ecc
    // -------------------------------------------------------------------------------------------------------------- //
    if (!ecc_output_filename.empty()) {
      eeec.dump_ecc(ecc_output_filename);
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpi_rank == 0) std::cout << "Dumped EEC in " << ecc_output_filename << std::endl;
    }
  }  // END Main MPI

  return 0;
}
