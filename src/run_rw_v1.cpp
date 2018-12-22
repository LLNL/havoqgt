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
#include <string>
#include <fstream>
#include <random>
#include <tuple>
#include <cassert>
#include <chrono>
#include <sstream>

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/higher_order_random_walk.hpp>

static constexpr int k_num_histories = 3; // 0:square; 1:triangle; 2:previous vertex

void parse_cmd_line(int argc, char **argv,
                    std::string &input_filename,
                    bool &personalized,
                    uint64_t &seed_vertex,
                    uint64_t &num_walkers,
                    uint64_t &die_rate,
                    uint64_t &max_walk_length,
                    bool &warp_to_seed,
                    std::array<int, k_num_histories> &closing_rates,
                    std::string &score_dump_file_prefix,
                    uint64_t &num_top) {
  if (havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  // Initial values
  input_filename.clear();
  bool found_input_filename = false;
  personalized = false;
  num_walkers = 1024;
  die_rate = 10;
  max_walk_length = 128;
  warp_to_seed = false;
  closing_rates.fill(100);
  closing_rates[closing_rates.size() - 1] = 0; // Does not return to the previous vertex
  score_dump_file_prefix.clear();
  num_top = 10;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:ps:w:d:l:rc:o:t:h")) != -1) {
    switch (c) {
      case 'h':prn_help = true;
        break;
      case 'i':found_input_filename = true;
        input_filename = optarg;
        break;
      case 'p': personalized = true;
        break;
      case 's': seed_vertex = std::stoll(optarg);
        break;
      case 'w':num_walkers = std::stoll(optarg);
        break;
      case 'd':die_rate = std::stoll(optarg);
        break;
      case 'l':max_walk_length = std::stoll(optarg);
        break;
      case 'r': warp_to_seed = true;
        break;
      case 'c': {
        std::stringstream ss(optarg);
        std::string buf;
        int count = 0;
        while (std::getline(ss, buf, ':')) {
          closing_rates[count] = std::stoi(buf);
          ++count;
        }
        break;
      }
      case 'o':score_dump_file_prefix = optarg;
        break;
      case 't': num_top = std::stoll(optarg);
        break;
      default:std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_input_filename) {
    // usage();
    exit(-1);
  }

  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  if (mpi_rank == 0) {
    std::cout << "input_filename: " << input_filename
              << "\npersonalized: " << personalized
              << "\nseed_vertex: " << seed_vertex
              << "\nnum_walkers: " << num_walkers
              << "\ndie_rate: " << die_rate
              << "\nmax_walk_length: " << max_walk_length
              << "\nwarp_to_seed: " << warp_to_seed
              << "\nscore_dump_file_prefix: " << score_dump_file_prefix
              << "\nnum_top: " << num_top
              << "\nclosing_rates: ";
    for (const auto item : closing_rates) std::cout << item << " ";
    std::cout << std::endl;
  }
}

template <typename graph_type>
uint64_t select_seed(const graph_type *const graph, const uint64_t candidate_vertex_id) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  uint64_t seed_id = candidate_vertex_id;
  uint64_t global_degree(0);

  while (true) {
    typename graph_type::vertex_locator seed = graph->label_to_locator(seed_id);
    if (seed.is_delegate()) {
      std::cerr << "Delegate vertex found" << std::endl;
      std::abort();
    }

    uint64_t local_degree = 0;
    if (uint32_t(mpi_rank) == seed.owner()) {
      local_degree = graph->degree(seed);
    }
    global_degree = havoqgt::mpi_all_reduce(local_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
    if (global_degree > 0) break;

    ++seed_id;
  }

  typename graph_type::vertex_locator seed = graph->label_to_locator(seed_id);
  if (uint32_t(mpi_rank) == seed.owner()) {
    if (seed_id == candidate_vertex_id) {
      std::cout << "Starting vertex = " << seed_id << std::endl;
    } else {
      std::cout << "Vertex " << candidate_vertex_id << " has a degree of 0." << std::endl;
      std::cout << "New seed vertex = " << seed_id << std::endl;
    }
    std::cout << "delegate? = " << seed.is_delegate() << std::endl;
    std::cout << "local_id = " << seed.local_id() << std::endl;
    std::cout << "degree = " << graph->degree(seed) << std::endl;
  }

  return seed_id;
}

void dump_score(const std::vector<std::pair<double, uint64_t>> &score,
                const uint64_t num_top,
                const std::string &dump_file_name) {
  std::ofstream ofs(dump_file_name);
  if (!ofs.is_open()) {
    std::cerr << "Can not open: " << dump_file_name << std::endl;
    std::abort();
  }
  for (int i = 0; i < std::min((uint64_t)num_top, (uint64_t)score.size()); ++i) {
    ofs << score[i].second << " : " << score[i].first << std::endl;
  }
  ofs.close();
}

template <typename graph_type>
void dump_all_score(const graph_type *const graph, const std::string &file_name,
                    const std::function<double(typename graph_type::vertex_locator)> &score_function) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::ofstream ofs(file_name + "_" + std::to_string(mpi_rank));
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    ofs << graph->locator_to_label(*vitr) << " " << std::fixed << score_function(*vitr) << "\n";
  }
  ofs.close();
}

int main(int argc, char **argv) {
  using graph_type = havoqgt::delegate_partitioned_graph<havoqgt::distributed_db::segment_manager_type>;
  using vertex_locator = typename graph_type::vertex_locator;

  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    havoqgt::get_environment();

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      havoqgt::get_environment().print();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    bool personalized;
    uint64_t seed_vertex_id;
    uint64_t num_walkers;
    uint64_t die_rate;
    uint64_t max_walk_length;
    bool warp_to_seed;
    std::array<int, k_num_histories> closing_rates;
    std::string score_dump_file_prefix;
    uint64_t num_top;

    parse_cmd_line(argc,
                   argv,
                   graph_input,
                   personalized,
                   seed_vertex_id,
                   num_walkers,
                   die_rate,
                   max_walk_length,
                   warp_to_seed,
                   closing_rates,
                   score_dump_file_prefix,
                   num_top);

    MPI_Barrier(MPI_COMM_WORLD);
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::size_t num_wk_per_process = (personalized) ? num_walkers : num_walkers / mpi_size;
    vertex_locator seed_vertex = graph->label_to_locator(seed_vertex_id);
    ho_rw::algo_data_type<graph_type, k_num_histories> algorithm_data(graph,
                                                                      personalized,
                                                                      seed_vertex,
                                                                      num_wk_per_process,
                                                                      die_rate,
                                                                      max_walk_length,
                                                                      warp_to_seed,
                                                                      closing_rates);
    run_rw(graph, algorithm_data);
    const auto top_scores = compute_top_scores(graph, algorithm_data, num_top);

    if (!score_dump_file_prefix.empty()) {
      if (mpi_rank == 0) std::cout << "\nDumping the top dead scores" << std::endl;
      if (mpi_rank == 0) dump_score(top_scores.first, num_top, score_dump_file_prefix + std::string("_dead_score"));
      MPI_Barrier(MPI_COMM_WORLD);

      if (mpi_rank == 0) std::cout << "\nDumping the top visit scores" << std::endl;
      if (mpi_rank == 0) dump_score(top_scores.second, num_top, score_dump_file_prefix + std::string("_visit_score"));
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) std::cout << "Finished rw v1" << std::endl;
  }

  return 0;
}
