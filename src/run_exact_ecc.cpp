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

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/fixed_size_unordered_map.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/k_breadth_first_search_sync_level_per_source.hpp>
#include <havoqgt/exact_eccentricity.hpp>

// -------------------------------------------------------------------------------------------------------------- //
// Types
// -------------------------------------------------------------------------------------------------------------- //
using namespace havoqgt;

using segment_manager_t = havoqgt::distributed_db::segment_manager_type;
using graph_t = havoqgt::delegate_partitioned_graph<segment_manager_t>;

#ifdef NUM_SOURCES
constexpr int k_num_sources = NUM_SOURCES;
#else
constexpr int k_num_sources = 8;
#endif

using kbfs_t = havoqgt::kbfs_type<graph_t, k_num_sources>;
using eecc_t = havoqgt::eecc_type<graph_t, k_num_sources>;

// -------------------------------------------------------------------------------------------------------------- //
// Parse command line
// -------------------------------------------------------------------------------------------------------------- //
void usage()
{
  if (havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -i <string>    - input graph base filename (required)\n"
              << " -b <string>    - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
              << " -s <int>:<int> - Colon separated k source vertices (required)\n"
              << " -h             - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char **argv, std::string &input_filename,
                    std::string &backup_filename,
                    std::vector<uint64_t> &source_id_list)
{
  if (havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:b:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's': {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_id_list.push_back(std::stoull(buf.c_str()));
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
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_input_filename || source_id_list.size() < k_num_sources) {
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
// find_max_lower_upper
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
std::pair<level_t, level_t> find_max_lower_upper(iterator_t vitr, iterator_t end,
                                                 eecc_t::vertex_data &ecc_vertex_data)
{
  level_t max_lower = std::numeric_limits<level_t>::min();
  level_t max_upper = std::numeric_limits<level_t>::min();
  for (; vitr != end; ++vitr) {
    if (ecc_vertex_data.upper[*vitr] == std::numeric_limits<level_t>::max()) continue; // skip unvisited vertices
    max_lower = std::max(max_lower, ecc_vertex_data.lower[*vitr]);
    max_upper = std::max(max_upper, ecc_vertex_data.upper[*vitr]);
  }

  return std::make_pair(max_lower, max_upper);
}

// -------------------------------------------------------------------------------------------------------------- //
// find_max_ecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
level_t find_max_ecc(iterator_t vitr, iterator_t end, typename eecc_t::vertex_data &ecc_vertex_data)
{
  level_t max_ecc = std::numeric_limits<level_t>::min();
  for (; vitr != end; ++vitr) {
    max_ecc = std::max(max_ecc, ecc_vertex_data.lower[*vitr]);
  }
  return max_ecc;
}

// -------------------------------------------------------------------------------------------------------------- //
// Main
// -------------------------------------------------------------------------------------------------------------- //
int main(int argc, char **argv)
{

  int mpi_rank(0), mpi_size(0);


  havoqgt::havoqgt_init(&argc, &argv);
  {
    // -------------------------------------------------------------------------------------------------------------- //
    //                                            Parse options & Prepare graph
    // -------------------------------------------------------------------------------------------------------------- //
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    havoqgt::get_environment();

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      havoqgt::get_environment().print();
      std::cout << "k_num_sources " << k_num_sources << std::endl;
      //print_system_info(false);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_filename;
    std::vector<uint64_t> source_id_list;

    parse_cmd_line(argc, argv, graph_input, backup_filename, source_id_list);
    MPI_Barrier(MPI_COMM_WORLD);

    if (backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_t *graph = ddb.get_segment_manager()->find<graph_t>("graph_obj").first;
    assert(graph != nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }

    std::vector<typename graph_t::vertex_locator> source_locator_list = select_non_zero_degree_source(graph, source_id_list);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Selected initial vertices" << std::endl;
    }

    // -------------------------------------------------------------------------------------------------------------- //
    //                                        Compute exact ecc and diameter
    //-------------------------------------------------------------------------------------------------------------- //
    const double total_start_time = MPI_Wtime();
    kbfs_t::vertex_data kbfs_vertex_data(*graph);
    eecc_t::vertex_data ecc_vertex_data(*graph);
    if (mpi_rank == 0) std::cout << "Allocated vertex data: " << std::endl;

    // ------------------------------ Init ecc vertex data ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      ecc_vertex_data.init();
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) std::cout << "Init ecc vertex data: " << time_end - time_start << std::endl;
    }

    // ------------------------------ Compute exact ECC ------------------------------ //
    size_t count_iteration(0);
    while (true) {
      if (mpi_rank == 0)
        std::cout << "\n==================== " << count_iteration << " ====================" << std::endl;

      // ------------------------------ Select sources ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        if (count_iteration == 0) {
          // Do nothing
        } else {
          source_locator_list = select_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data,
                                                                      eecc_source_select_mode_tag::far_and_hdeg());
        }
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Sources: " << time_end - time_start << std::endl;
          for (auto locator : source_locator_list)
            std::cout << graph->locator_to_label(locator) << " ";
          std::cout << std::endl;
        }
      }

      // ------------------------------ Reset data for KBFS ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        kbfs_vertex_data.init();
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Reset data for KBFS: " << time_end - time_start << std::endl;
        }
      }

      // ------------------------------ KBFS ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        k_breadth_first_search_level_per_source<graph_t, k_num_sources>(graph, kbfs_vertex_data, source_locator_list);
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << "BFS Time: " << time_end - time_start << std::endl;
      }

      // ------------------------------ Compute exact ecc ------------------------------ //
      {
        const size_t num_remains = compute_eecc<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data, source_locator_list);
        if (num_remains == 0) break; // Terminal condition
      }

      {
        plun_single_degree_vertices<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data);
      }

      {
        std::vector<size_t> histgram = compute_distance_histgram<graph_t, k_num_sources>(graph, kbfs_vertex_data, ecc_vertex_data, 5);
        if (mpi_rank == 0) {
          std::cout << "distance score (upper - lower)" << std::endl;
          for (size_t i = 0; i < histgram.size(); ++i) {
            std::cout << "DS: " << i << "\t" << histgram[i] << std::endl;
          }
        }
      }

      // ------------------------------ Diameter calculation termination test ------------------------------ //
      {
        const auto ret1 = find_max_lower_upper(graph->vertices_begin(), graph->vertices_end(), ecc_vertex_data);
        const auto ret2 = find_max_lower_upper(graph->controller_begin(), graph->controller_end(), ecc_vertex_data);

        level_t max_lower = std::max(ret1.first, ret2.first);
        level_t max_upper = std::max(ret1.second, ret2.second);

        max_lower = mpi_all_reduce(max_lower, std::greater<level_t>(), MPI_COMM_WORLD);
        max_upper = mpi_all_reduce(max_upper, std::greater<level_t>(), MPI_COMM_WORLD);

        if (mpi_rank == 0) std::cout << "terminate test " << max_lower << " : " << max_upper << std::endl;
      }
      ++count_iteration;
    } // End exact ecc loop
    const double total_end_time = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Total execution time: " << total_end_time - total_start_time << std::endl;

    // ------------------------------ Compute diameter ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      const level_t max_ecc = std::max(find_max_ecc(graph->vertices_begin(), graph->vertices_end(), ecc_vertex_data),
                                       find_max_ecc(graph->controller_begin(), graph->controller_end(), ecc_vertex_data));
      level_t diameter;
      CHK_MPI(MPI_Reduce(&max_ecc, &diameter, 1, mpi_typeof(max_ecc), MPI_MAX, 0, MPI_COMM_WORLD));
      const double time_end = MPI_Wtime();

      if (mpi_rank == 0) {
        std::cout << "\nCompute diameter: " << time_end - time_start << std::endl;
        std::cout << "Diameter: " << diameter << std::endl;
      }
    } // End compute exact ecc & diameter

  }  // END Main MPI
  havoqgt_finalize();

  return 0;
}
