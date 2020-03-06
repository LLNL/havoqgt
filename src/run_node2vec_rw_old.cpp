/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <iomanip>
#include <exception>
#include <tuple>

#include <boost/container/vector.hpp>
#include <boost/container/scoped_allocator.hpp>
#include <metall_utility/metall_mpi_adaptor.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/mpi.hpp>
#include <havoqgt/distributed_db.hpp>
#include <node2vec_rw/k_bfs.hpp>
#include <node2vec_rw/global_shuffle.hpp>
#include <node2vec_rw/hash.hpp>
#include <node2vec_rw/random.hpp>

#include <node2vec/node2vec_rw.hpp>

namespace bc = boost::container;

using graph_allocator_type = typename havoqgt::distributed_db::segment_manager_type::template allocator<void>::type;
using graph_type = havoqgt::delegate_partitioned_graph<graph_allocator_type>;
using vertex_type = typename graph_type::vertex_locator;
using level_type = uint16_t;
using edge_weight_type = double;
using edge_weight_data_type = typename graph_type::template edge_data<edge_weight_type,
                                                                      std::allocator<edge_weight_type>>;

enum walk_data {
  walk_id = 0,
  step_no = 1,
  visited_vertex = 2
};
using walk_data_type = std::tuple<uint64_t, uint16_t, vertex_type>;
using walk_history_type = bc::vector<walk_data_type, metall::manager::allocator_type<walk_data_type>>;

using train_data_list_type = bc::vector<std::pair<uint64_t, uint64_t>,
                                        metall::manager::allocator_type<std::pair<uint64_t, uint64_t>>>;

void parse_cmd_line(int argc, char **argv,
                    std::string *graph_filename,
                    std::string *backup_filename,
                    std::string *storage_dir_path,
                    std::vector<uint64_t> *source_id_list,
                    node2vec_rw::option *rw_option,
                    int *skip_gram_context_size,
                    std::string *out_file_name,
                    bool *generate_train_data) {
  int c;
  while ((c = getopt(argc, argv, "v:g:b:d:l:r:p:q:s:k:o:t")) != -1) {
    switch (c) {
      case 'v': {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_id_list->push_back(std::stoull(buf.c_str()));
        break;
      }

      case 'g':*graph_filename = optarg;
        break;

      case 'b':*backup_filename = optarg;
        break;

      case 'd':*storage_dir_path = optarg;
        break;

      case 'l':rw_option->length = std::stoull(optarg);
        break;

      case 'r': rw_option->num_walkers_per_vertex = std::stoll(optarg);
        break;

      case 'p':rw_option->p = std::stod(optarg);
        break;

      case 'q':rw_option->q = std::stod(optarg);
        break;

      case 's':rw_option->seed = std::stoul(optarg);
        break;

      case 'k':*skip_gram_context_size = std::stoul(optarg);
        break;

      case 'o':*out_file_name = optarg;
        break;

      case 'w':*generate_train_data = true;
        break;
    }
  }

  if (source_id_list->empty()) {
    std::cerr << "k-BFS source ID list cannot be empty" << std::endl;
    std::abort();
  }

  if (graph_filename->empty() && backup_filename->empty()) {
    std::cerr << "Graph file name and backup file name are empty" << std::endl;
    std::abort();
  }
}

template <typename k_bfs_level_diff_table_type, typename allocator_type>
auto run_walker(graph_type &graph,
                const edge_weight_data_type &edge_weight,
                const k_bfs_level_diff_table_type &k_bfs_level_diff_table,
                const node2vec_rw::option &rw_option,
                const allocator_type &allocator = allocator_type()) {
  auto walk_store = std::make_unique<walk_history_type>(allocator);

  node2vec_rw::run_node2vec_rw_visitor<graph_allocator_type, edge_weight_data_type, walk_history_type>
      (graph,
       edge_weight,
       k_bfs_level_diff_table,
       rw_option.p,
       rw_option.q,
       rw_option.length,
       rw_option.num_walkers_per_vertex,
       walk_store.get());

  return walk_store;
}

void dump_local_walk_list(const graph_type &graph,
                          const int walk_length,
                          const std::string &out_file_name,
                          const walk_history_type &walk_list) {
  int mpi_rank = -1;
  int mpi_size = -1;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::string fout_name(out_file_name + "-" + std::to_string(mpi_rank) + "-" + std::to_string(mpi_size));
  std::ofstream ofs(fout_name);

  if (mpi_rank == 0) {
    std::cout << "Writing to file..." << std::endl;
  }

  const std::size_t num_walks = walk_list.size() / walk_length;
  for (std::size_t w = 0; w < num_walks; ++w) {
    for (int i = 0; i < walk_length; ++i) {
      const auto v = graph.locator_to_label(std::get<walk_data::visited_vertex>(walk_list[w * walk_length + i]));
      ofs << v << " ";
    }
    ofs << "\n";
  }
  ofs.close();
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Writing to file done" << std::endl;
  }
}

template <typename allocator_type>
auto shuffle_walk_data(const walk_history_type &walk_dump, const allocator_type &allocator = allocator_type()) {
  int mpi_rank = -1;
  int mpi_size = -1;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  walk_history_type shuffled_walk_list(allocator);

  auto walk_data_partitioner = [mpi_size](const walk_data_type &history_data) -> int {
    return node2vec_rw::hash<std::tuple_element<walk_data::walk_id, walk_data_type>::type>{}
        (std::get<walk_data::walk_id>(history_data)) % mpi_size;
  };

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Global history shuffle start" << std::endl;
  }
  const double global_history_shuffle_time_start = MPI_Wtime();
  node2vec_rw::mpi_large_all_to_all(walk_dump.begin(), walk_dump.end(),
                                  walk_data_partitioner,
                                  MPI_COMM_WORLD,
                                  std::back_inserter(shuffled_walk_list));
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Global history shuffle took:\t" << MPI_Wtime() - global_history_shuffle_time_start << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Local history shuffle start" << std::endl;
  }
  const double local_history_shuffle_time_start = MPI_Wtime();
  std::sort(shuffled_walk_list.begin(), shuffled_walk_list.end(),
            [](const walk_data_type &lhd, const walk_data_type &rhd) {
              if (std::get<walk_data::walk_id>(lhd) != std::get<walk_data::walk_id>(rhd))
                return std::get<walk_data::walk_id>(lhd) < std::get<walk_data::walk_id>(rhd);
              return std::get<walk_data::step_no>(lhd) < std::get<walk_data::step_no>(rhd);
            });
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Local history shuffle took:\t" << MPI_Wtime() - local_history_shuffle_time_start << std::endl;
  }

  return shuffled_walk_list;
}

template <typename allocator_type>
void shuffle_train_data(const graph_type &graph,
                        const int context_size,
                        const int walk_length,
                        const std::size_t num_words_per_rank,
                        const std::string &out_file_name,
                        const walk_history_type &walk_list,
                        const allocator_type &allocator = allocator_type()) {
  int mpi_rank = -1;
  int mpi_size = -1;
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  train_data_list_type train_data_send_buf(allocator);

  const std::size_t num_walks = walk_list.size() / walk_length;
  assert(walk_list.size() % walk_length == 0);
  node2vec_rw::rand_1024 rand(std::random_device{}());
  for (std::size_t w = 0; w < num_walks; ++w) {

    for (int center_word_index = 0; center_word_index < walk_length; ++center_word_index) {
      std::uniform_int_distribution<int> dist(0, context_size - 1);
      const int offset = dist(rand);
      for (int c = offset; c < context_size * 2 + 1 - offset; ++c) {
        const int context_word_pos = center_word_index - context_size + c;

        // Chosen context position was too small or large to fit into 'sentence' array
        if (context_word_pos < 0) continue;
        if (context_word_pos >= walk_length) continue;

        const auto &context_word = std::get<walk_data::visited_vertex>(walk_list[w * walk_length + context_word_pos]);
        const auto &center_word = std::get<walk_data::visited_vertex>(walk_list[w * walk_length + center_word_index]);
        if (context_word == center_word) continue;

        train_data_send_buf.emplace_back(std::make_pair(graph.locator_to_label(context_word),
                                                        graph.locator_to_label(center_word)));

      }
    }
  }
  std::cout << "Unshuffled local train data size at " << mpi_rank << " \t " << train_data_send_buf.size() << std::endl;

  train_data_list_type train_data_recv_buf(allocator);

  auto train_data_partitioner = [num_words_per_rank](const train_data_list_type::value_type &train_data) -> int {
    return train_data.first / num_words_per_rank;
  };

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Train data global shuffle start" << std::endl;
  }
  const double global_train_data_shuffle_time_start = MPI_Wtime();
  node2vec_rw::mpi_large_all_to_all(train_data_send_buf.begin(), train_data_send_buf.end(),
                                  train_data_partitioner,
                                  MPI_COMM_WORLD,
                                  std::back_inserter(train_data_recv_buf));
  if (mpi_rank == 0) {
    std::cout << "Train data global shuffle took:\t" << MPI_Wtime() - global_train_data_shuffle_time_start << std::endl;
  }
  train_data_send_buf.clear();
  train_data_send_buf.shrink_to_fit();

  if (mpi_rank == 0) {
    std::cout << "Writing to file..." << std::endl;
  }
  std::string
      fout_name(out_file_name + "-" + std::to_string(mpi_rank) + "-" + std::to_string(mpi_size));
  std::ofstream ofs(fout_name);
  for (const auto &elem : train_data_recv_buf) {
    ofs << elem.first << "\t" << elem.second << "\n";
  }
  ofs.close();
  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Writing to file done" << std::endl;
  }
}

int main(int argc, char **argv) {

  havoqgt::init(&argc, &argv);
  {
    std::cout << std::fixed;
    std::cout << std::setprecision(1);

    int mpi_rank = -1;
    int mpi_size = -1;
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_filename;
    std::vector<uint64_t> source_id_list;
    node2vec_rw::option rw_option;
    std::string storage_dir_path = "/tmp/rw_datastore";
    int skip_gram_context_size = 10;
    std::string out_file_name;
    bool generate_train_data = false;
    parse_cmd_line(argc,
                   argv,
                   &graph_input,
                   &backup_filename,
                   &storage_dir_path,
                   &source_id_list,
                   &rw_option,
                   &skip_gram_context_size,
                   &out_file_name,
                   &generate_train_data);
    MPI_Barrier(MPI_COMM_WORLD);

    if (!backup_filename.empty()) {
      havoqgt::distributed_db::transfer(backup_filename.c_str(), graph_input.c_str());
    }
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    node2vec_rw::k_bfs<graph_allocator_type, level_type> k_bfs(*graph);
    std::vector<graph_type::vertex_locator> source_candidate_list;
    for (const auto vid : source_id_list) {
      source_candidate_list.emplace_back(graph->label_to_locator(vid));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const double k_bfs_time_start = MPI_Wtime();
    auto k_bfs_level_table = k_bfs.run(source_candidate_list);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "\nk-BFS took:\t" << MPI_Wtime() - k_bfs_time_start << std::endl;
    }

    const double k_bfs_level_diff_compute_time_start = MPI_Wtime();
    auto k_bfs_level_diff_table = node2vec_rw::construct_k_bfs_level_diff_table(*graph, k_bfs_level_table);
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "\nk-BFS level diff computing took:\t"
                << MPI_Wtime() - k_bfs_level_diff_compute_time_start << std::endl;
    }

    edge_weight_data_type edge_weight(*graph);
    edge_weight.reset(1.0f);

    metall::utility::metall_mpi_adaptor metall_adaptor(metall::create_only, storage_dir_path);
    auto &metall_manager = metall_adaptor.get_local_manager();

    const double rw_time_start = MPI_Wtime();
    auto walk_history_store = run_walker(*graph,
                                         edge_weight,
                                         k_bfs_level_diff_table,
                                         rw_option,
                                         metall_manager.get_allocator());
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "\nRandom walk took:\t" << MPI_Wtime() - rw_time_start << std::endl;
    }
    std::cout << "Local walk data dump size at " << mpi_rank << " \t " << walk_history_store->size() << std::endl;

    if (!out_file_name.empty()) {
      auto walk_list = shuffle_walk_data(*walk_history_store, metall_manager.get_allocator());

      if (generate_train_data) {
        walk_history_store.reset(nullptr);
        std::cout << "Shuffled walk data dump size at " << mpi_rank << " \t " << walk_list.size() << std::endl;

        const std::size_t num_words_per_rank = (graph->max_global_vertex_id() + mpi_size - 1) / mpi_size;
        shuffle_train_data(*graph,
                           skip_gram_context_size,
                           rw_option.length,
                           num_words_per_rank,
                           out_file_name,
                           walk_list,
                           metall_manager.get_allocator());
      } else {
        dump_local_walk_list(*graph, rw_option.length, out_file_name, walk_list);
      }
    }
  }
  // END Main MPI

  return 0;
}
