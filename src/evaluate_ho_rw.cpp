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
#include <unordered_map>
#include <random>
#include <algorithm>
#include <iterator>

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/higher_order_random_walk.hpp>

static constexpr int k_num_histories = 3; // 0:square; 1:triangle; 2:previous vertex

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
    uint64_t die_rate{15};

    MPI_Barrier(MPI_COMM_WORLD);
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::unordered_map<uint64_t, uint64_t> true_community_table;
    std::unordered_map<uint64_t, std::size_t> community_size_table;
    std::vector<uint64_t> seed_vertex_id_list;
    {
      std::ifstream ifs("true");
      if (!ifs.is_open()) {
        std::abort();
      }
      uint64_t vid;
      uint64_t comm_id;
      while (ifs >> vid >> comm_id) {
        true_community_table[vid] = comm_id;
        seed_vertex_id_list.emplace_back(vid);
      }

      for (auto elem : true_community_table) {
        ++community_size_table[elem.second];
      }

      std::mt19937 g(123); // Note all processors must use the same seed
      std::shuffle(seed_vertex_id_list.begin(), seed_vertex_id_list.end(), g);
    }

    std::ofstream ofs;
    if (havoqgt::mpi_comm_rank() == 0) {
      ofs.open("evaluation_result.log");
      ofs << "seed\tcs\tct\twk\tdr\trv\trd" << std::endl;
    }

    for (auto seed_vertex_id : seed_vertex_id_list) {
      vertex_locator seed_vertex = graph->label_to_locator(seed_vertex_id);
      const uint64_t community_id = true_community_table[seed_vertex_id];
      const std::size_t community_size = community_size_table[community_id];

      for (int closing_rate = 0; closing_rate <= 100; closing_rate += 20) {
        std::array<int, k_num_histories> closing_rates{closing_rate, closing_rate, 0};
        for (std::size_t num_walkers = 10000; num_walkers <= 10000000ULL; num_walkers *= 10ULL) {
          ho_rw::algo_data_type<graph_type, k_num_histories> algorithm_data(graph,
                                                                            true, // personalized
                                                                            seed_vertex,
                                                                            num_walkers,
                                                                            die_rate, // dead rate
                                                                            std::numeric_limits<uint64_t>::max(),
                                                                            false, // warp to seed
                                                                            closing_rates);
          run_rw(graph, algorithm_data);
          const auto top_scores = compute_top_scores(graph, algorithm_data, community_size);

          std::vector<uint64_t> top_visit_score_vertices;
          for (auto elem : top_scores.first) top_visit_score_vertices.emplace_back(elem.second);
          assert(top_visit_score_vertices.size() == community_size);
          const auto recall_visit = ho_rw::compute_recall(top_visit_score_vertices,
                                                          community_id,
                                                          community_size_table[community_id],
                                                          true_community_table);

          std::vector<uint64_t> top_dead_score_vertices;
          for (auto elem : top_scores.second) top_dead_score_vertices.emplace_back(elem.second);
          assert(top_dead_score_vertices.size() == community_size);
          const auto recall_dead = ho_rw::compute_recall(top_dead_score_vertices,
                                                         community_id,
                                                         community_size_table[community_id],
                                                         true_community_table);

          if (havoqgt::mpi_comm_rank() == 0) {
            std::stringstream ss;
            ss << seed_vertex_id << "\t" << closing_rates[0] << "\t" << closing_rates[1] << "\t" << num_walkers << "\t"
               << die_rate
               << "\t" << recall_visit << "\t" << recall_dead;

            std::cout << ss.str() << "\n" << std::endl;

            ofs << ss.str() << std::endl;
            ofs.flush();
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) std::cout << "Finished evaluation" << std::endl;
  }

  return 0;
}
