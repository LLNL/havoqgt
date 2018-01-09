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

#include <havoqgt/environment.hpp>
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
#include <array>
#include <chrono>
#include <random>

#include <boost/interprocess/managed_heap_memory.hpp>

// -------------------------------------------------------------------------------------------------------------- //
// Types
// -------------------------------------------------------------------------------------------------------------- //
using namespace havoqgt;

using segment_manager_t = havoqgt::distributed_db::segment_manager_type ;

using graph_t                       = havoqgt::delegate_partitioned_graph<segment_manager_t>;
using k_level_t                     = std::array<level_t, k_num_sources>;

// -------------------------------------------------------------------------------------------------------------- //
// Data structure
// -------------------------------------------------------------------------------------------------------------- //
struct kbfs_vertex_data_t
{
  using vertex_data_k_level_t         = typename graph_t::vertex_data<k_level_t, std::allocator<k_level_t>>;
  using vertex_data_k_visit_bitmap_t  = typename graph_t::vertex_data<visit_bitmap_t, std::allocator<visit_bitmap_t>>;
  using vertex_data_visit_flag_t      = typename graph_t::vertex_data<bool, std::allocator<bool>>;

  kbfs_vertex_data_t(const graph_t& graph)
    : level(graph),
      visit_bitmap(graph),
      next_visit_bitmap(graph),
      visit_flag(graph),
      next_visit_flag(graph) { }

  void init()
  {
    k_level_t zero_initialied;
    zero_initialied.fill(0);
    level.reset(zero_initialied);

    visit_bitmap.reset(visit_bitmap_t());
    next_visit_bitmap.reset(visit_bitmap_t());

    visit_flag.reset(false);
    next_visit_flag.reset(false);
  }

  vertex_data_k_level_t level;
  vertex_data_k_visit_bitmap_t visit_bitmap;
  vertex_data_k_visit_bitmap_t next_visit_bitmap;
  vertex_data_visit_flag_t visit_flag; // whether visited at the previous level
  vertex_data_visit_flag_t next_visit_flag;
};


struct ecc_vertex_data_t
{
  using vertex_data_ecc_t = typename graph_t::vertex_data<level_t, std::allocator<level_t>>;

  ecc_vertex_data_t(const graph_t &graph)
    : lower(graph),
      upper(graph) { }

  void init()
  {
    lower.reset(std::numeric_limits<level_t>::min());
    upper.reset(std::numeric_limits<level_t>::max());
  }

  vertex_data_ecc_t lower;
  vertex_data_ecc_t upper;
};


// -------------------------------------------------------------------------------------------------------------- //
// Parse command line
// -------------------------------------------------------------------------------------------------------------- //
void usage()  {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -i <string>   - input graph base filename (required)\n"
              << " -b <string>   - backup graph base filename.  If set, \"input\" graph will be deleted if it exists\n"
              << " -s <int>      - Source vertex of BFS (Default is 0)\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& input_filename,
                    std::string& backup_filename,
                    std::vector<uint64_t>& source_id_list) {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
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
      case 's':
      {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_id_list.push_back(std::stoull(buf.c_str()));
        assert(source_id_list.size() == k_num_sources);
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
  if (prn_help || !found_input_filename) {
    usage();
    exit(-1);
  }
}


// -------------------------------------------------------------------------------------------------------------- //
// select_source
// -------------------------------------------------------------------------------------------------------------- //
void select_source(const graph_t *const graph, std::vector<uint64_t> &source_id_list)
{

  std::vector<uint64_t> global_source_id_list;
  havoqgt::mpi_all_gather(source_id_list, global_source_id_list, MPI_COMM_WORLD);

  std::mt19937 generator(123);
  std::uniform_int_distribution<size_t> dist(0, global_source_id_list.size());
  std::random_shuffle(global_source_id_list.begin(), global_source_id_list.end(), [&](size_t i) {
    return dist(generator);
  });
  if (global_source_id_list.size() >= k_num_sources) global_source_id_list.resize(k_num_sources); // discard candidates

  source_id_list = std::move(global_source_id_list);
}

// -------------------------------------------------------------------------------------------------------------- //
// compute_ecc_k_source
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
void compute_ecc_k_source(iterator_t vitr, iterator_t end,
                          const kbfs_vertex_data_t &kbfs_vertex_data,
                          k_level_t &k_source_ecc)
{
  for (; vitr != end; ++vitr) {
    for (size_t i = 0; i < k_num_sources; ++i) {
      k_source_ecc[i] = std::max(kbfs_vertex_data.level[*vitr][i], k_source_ecc[i]);
    }
  }
}


// -------------------------------------------------------------------------------------------------------------- //
// bound_ecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
std::pair<size_t, size_t> bound_ecc(const graph_t *const graph,
                                    iterator_t vitr,
                                    iterator_t end,
                                    const kbfs_vertex_data_t &kbfs_vertex_data,
                                    const k_level_t &k_source_ecc,
                                    ecc_vertex_data_t &ecc_vertex_data,
                                    std::vector<uint64_t>& source_candidate_id_list)
{
  size_t num_completed_vertices(0);
  size_t num_non_completed_vertices(0);

  for (; vitr != end; ++vitr) {
    if (!(kbfs_vertex_data.visit_bitmap[*vitr].has_bit())) continue;

    auto& current_lower = ecc_vertex_data.lower[*vitr];
    auto& current_upper = ecc_vertex_data.upper[*vitr];
    if (current_lower == current_upper) continue; // Exact ecc has been already found

    bool completed = false;
    for (size_t k = 0; k < k_num_sources; ++k) {
      auto& level = kbfs_vertex_data.level[*vitr][k];

      current_lower = std::max(current_lower, std::max(level, static_cast<uint16_t>(k_source_ecc[k] - level)));
      current_upper = std::min(current_upper, static_cast<uint16_t>(k_source_ecc[k] + level));

      if (current_lower == current_upper) {
        completed = true;
        break;
      }
    }
    num_completed_vertices += completed;
    num_non_completed_vertices += !completed;

    // While updating ecc boundary, correct the candidate of next BFS sources
    if (!completed && source_candidate_id_list.size() < k_num_sources) {
      source_candidate_id_list.push_back(graph->locator_to_label(*vitr));
    }
  }

  return std::make_pair(num_completed_vertices, num_non_completed_vertices);
}

// -------------------------------------------------------------------------------------------------------------- //
// find_max_ecc
// -------------------------------------------------------------------------------------------------------------- //
template <typename iterator_t>
level_t find_max_ecc(iterator_t vitr, iterator_t end, ecc_vertex_data_t &ecc_vertex_data)
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
int main(int argc, char** argv) {

  int mpi_rank(0), mpi_size(0);


  havoqgt::havoqgt_init(&argc, &argv);
  {
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

    // ---------------------------------------- Parsing options ---------------------------------------- //
    std::string graph_input;
    std::string backup_filename;
    std::vector<uint64_t> source_id_list;

    parse_cmd_line(argc, argv, graph_input, backup_filename, source_id_list);

    MPI_Barrier(MPI_COMM_WORLD);

    // ---------------------------------------- Prepare graph ---------------------------------------- //
    if(backup_filename.size() > 0) {
      distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }

    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

    graph_t *graph = ddb.get_segment_manager()->find<graph_t>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    //graph->print_graph_statistics();
    MPI_Barrier(MPI_COMM_WORLD);

    const double total_start_time = MPI_Wtime();
    // ---------------------------------------- Compute diameter ---------------------------------------- //
    kbfs_vertex_data_t kbfs_vertex_data(*graph);
    ecc_vertex_data_t ecc_vertex_data(*graph);
    if (mpi_rank == 0) std::cout << "Allocated vertex data: " << std::endl;

    // ------------------------------ Init vertex data ------------------------------ //
    {
      const double time_start = MPI_Wtime();
      ecc_vertex_data.init();
      MPI_Barrier(MPI_COMM_WORLD);
      const double time_end = MPI_Wtime();
      if (mpi_rank == 0) std::cout << "Init vertex data: " << time_end - time_start << std::endl;
    }

    size_t current_iteration_number(0);
    while (true) {
      if (mpi_rank == 0) std::cout << "\n==================== " << current_iteration_number << " ====================" << std::endl;

      // ------------------------------ Init BFS ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        kbfs_vertex_data.init();
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << "BFS initialization: " << time_end - time_start << std::endl;
      }

      // ------------------------------ Select sources ------------------------------ //
      std::vector<typename graph_t::vertex_locator> source_locator_list;
      {
        const double time_start = MPI_Wtime();
        if (current_iteration_number > 0) {
          select_source(graph, source_id_list);
        }
        for (auto& source_id : source_id_list) {
          source_locator_list.push_back(graph->label_to_locator(source_id));
        }
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << "Select sources: " << time_end - time_start << std::endl;
      }

      // ------------------------------ KBFS ------------------------------ //
      {
        const double time_start = MPI_Wtime();
        const kbfs_result_t result = havoqgt::k_breadth_first_search_level_per_source(graph,
                                                                                      kbfs_vertex_data.level,
                                                                                      kbfs_vertex_data.visit_bitmap,
                                                                                      kbfs_vertex_data.next_visit_bitmap,
                                                                                      kbfs_vertex_data.visit_flag,
                                                                                      kbfs_vertex_data.next_visit_flag,
                                                                                      source_locator_list);
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();

        if (mpi_rank == 0) std::cout << "BFS Time: " << time_end - time_start << std::endl;
      }

      // ------------------------------ Find ECC for k sources ------------------------------ //
      k_level_t k_source_ecc;
      {
        const double time_start = MPI_Wtime();
        k_source_ecc.fill(0);

        compute_ecc_k_source(graph->vertices_begin(), graph->vertices_end(), kbfs_vertex_data, k_source_ecc);
        compute_ecc_k_source(graph->controller_begin(), graph->controller_end(), kbfs_vertex_data, k_source_ecc);
        mpi_all_reduce_inplace(k_source_ecc, std::greater<level_t>(), MPI_COMM_WORLD);
        for (size_t k = 0; k < source_id_list.size(); ++k) {
          if (source_locator_list[k].owner() == mpi_rank) {
            ecc_vertex_data.lower[source_locator_list[k]] = ecc_vertex_data.upper[source_locator_list[k]] = k_source_ecc[k];
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) std::cout << "Find ECC for k sources: " << time_end - time_start << std::endl;
        for (size_t k = 0; k < source_id_list.size(); ++k)
          std::cout << k_source_ecc[k] << " ";
        std::cout << std::endl;
      }

      // ------------------------------ Bound ECC ------------------------------ //
      {
        size_t num_completed_vertices(0);
        size_t num_non_completed_vertices(0);

        const double time_start = MPI_Wtime();
        source_id_list.clear();

        {
          auto ret = bound_ecc(graph,
                               graph->vertices_begin(),
                               graph->vertices_end(),
                               kbfs_vertex_data,
                               k_source_ecc,
                               ecc_vertex_data,
                               source_id_list);
          num_completed_vertices += ret.first;
          num_non_completed_vertices += ret.second;
        }

        {
          auto ret = bound_ecc(graph,
                               graph->controller_begin(),
                               graph->controller_end(),
                               kbfs_vertex_data,
                               k_source_ecc,
                               ecc_vertex_data,
                               source_id_list);
          num_completed_vertices += ret.first;
          num_non_completed_vertices += ret.second;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        size_t global_num_completed_vertices(0);
        CHK_MPI(MPI_Reduce(&num_completed_vertices,
                           &global_num_completed_vertices,
                           1,
                           mpi_typeof(num_completed_vertices),
                           MPI_SUM,
                           0,
                           MPI_COMM_WORLD));

        size_t global_num_non_completed_vertices(0);
        CHK_MPI(MPI_Reduce(&num_non_completed_vertices,
                           &global_num_non_completed_vertices,
                           1,
                           mpi_typeof(num_non_completed_vertices),
                           MPI_SUM,
                           0,
                           MPI_COMM_WORLD));

        const double time_end = MPI_Wtime();
        if (mpi_rank == 0) {
          std::cout << "Update ECC bounds: " << time_end - time_start << std::endl;
          std::cout << "# completed vertices: " << global_num_completed_vertices << std::endl;
          std::cout << "# non-completed vertices: " << global_num_non_completed_vertices << std::endl;
        }
      }

      // ---------------------------------------- Check termination ---------------------------------------- //
      {
        const char wk = static_cast<char>(source_id_list.size() > 0);
        const char global = mpi_all_reduce(wk, std::logical_or<char>(), MPI_COMM_WORLD);
        if (!global) break;
        MPI_Barrier(MPI_COMM_WORLD); // Just in case
      }
      ++current_iteration_number;
    } // End Exact ECC

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
    } // End Compute diameter

    const double total_end_time = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "Total execution time: " << total_end_time - total_start_time << std::endl;
  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}
