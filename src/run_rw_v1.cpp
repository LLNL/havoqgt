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

#define WRITE_LOG 0
#if WRITE_LOG
std::ofstream *ofs_log;
#define LOG_OUT (*ofs_log)
//#define LOG_OUT std::cerr
#endif

static constexpr int k_num_history = 3; // 0:square; 1:triangle; 2:previous vertex
using rw_traversal_history_type = std::array<uint64_t, k_num_history>;
static constexpr int num_top = 20;

template <typename graph_type, typename history_type, int num_history>
class rw_algorithm_v1 {
 private:
  using vertex_locator = typename graph_type::vertex_locator;
  using count_table_type = typename graph_type::template vertex_data<uint64_t, std::allocator<uint64_t>>;
  using history_select_prob_table_type = std::array<int, num_history>;

 public:
  explicit rw_algorithm_v1(graph_type *const graph,
                           const bool personalized,
                           const vertex_locator &seed,
                           const std::size_t num_local_walkers,
                           const uint64_t die_rate,
                           const uint64_t max_walk_length,
                           const bool warp_to_seed,
                           const history_select_prob_table_type &history_select_prob_table)
      : m_graph(graph),
        m_personalized(personalized),
        m_seed_locator(seed),
        m_num_local_walkers(num_local_walkers),
        m_die_rate(die_rate),
        m_max_walk_length(max_walk_length),
        m_num_launched_walkers(0),
        m_rnd_gen((uint32_t)std::chrono::system_clock::now().time_since_epoch().count()),
        m_warp_to_seed(warp_to_seed),
        m_num_dead_table(*graph),
        m_num_visits_table(*graph),
        m_history_select_prob_table(history_select_prob_table) {
    m_num_dead_table.reset(0);
    m_num_visits_table.reset(0);
  }

  bool rw_launcher(const vertex_locator vertex) {
    if (m_graph->degree(vertex) == 0) return false;

    // Personalized model
    if (m_personalized) {
      if (m_seed_locator == vertex && m_num_launched_walkers < m_num_local_walkers) {
        ++m_num_launched_walkers;
        return true;
      }
    } else {
      std::uniform_int_distribution<uint64_t> dis(0, m_graph->num_local_vertices() - 1);
      if (dis(m_rnd_gen) < m_num_local_walkers) {
        ++m_num_launched_walkers;
        return true;
      }
    }

    return false;
  }

  history_type new_history() {
    history_type history;
    for (int i = 0; i < num_history; ++i) history[i] = m_graph->locator_to_label(warp_machine());
    return history;
  }

  bool russian_roulette() {
    std::uniform_int_distribution<int> dis(0, 99);
    return (dis(m_rnd_gen) < m_die_rate);
  }

  vertex_locator warp_machine() {
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->max_global_vertex_id() - 1);
    return m_graph->label_to_locator(dis(m_rnd_gen));
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex) {
    assert(m_graph->degree(vertex) > 0);
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->degree(vertex) - 1);
    const uint64_t offset = dis(m_rnd_gen);
    auto edge = m_graph->edges_begin(vertex) + offset;
    return edge.target();
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex, const history_type &history) {
    for (int i = 0; i < num_history; ++i) {
      for (auto eitr = m_graph->edges_begin(vertex), end = m_graph->edges_end(vertex); eitr != end; ++eitr) {
        if (m_graph->locator_to_label(eitr.target()) != history[i]) continue;

        std::uniform_int_distribution<int> dis(0, 99);
        if (dis(m_rnd_gen) < m_history_select_prob_table[i]) return eitr.target();
      }
    }
    return neighbor_roulette(vertex);
  }

  uint64_t max_walk_length() const {
    return m_max_walk_length;
  }

  void increment_num_dead(const vertex_locator vertex) {
    ++m_num_dead_table[vertex];
  }

  uint64_t num_dead(const vertex_locator vertex) const {
    return m_num_dead_table[vertex];
  }

  void increment_num_visits(const vertex_locator vertex) {
    ++m_num_visits_table[vertex];
  }

  uint64_t num_visits(const vertex_locator vertex) const {
    return m_num_visits_table[vertex];
  }

  uint64_t num_launched_walkers() const {
    return m_num_launched_walkers;
  }

  bool personalized() const {
    return m_personalized;
  }

  vertex_locator seed() const {
    return m_seed_locator;
  }

  graph_type *graph() {
    return m_graph;
  }

 private:
  graph_type *m_graph;
  bool m_personalized;
  vertex_locator m_seed_locator;
  uint64_t m_num_local_walkers;
  uint64_t m_die_rate;
  uint64_t m_max_walk_length;
  uint64_t m_num_launched_walkers;
  std::mt19937 m_rnd_gen;
  bool m_warp_to_seed;
  count_table_type m_num_dead_table;
  count_table_type m_num_visits_table;
  history_select_prob_table_type m_history_select_prob_table;
};

template <typename graph_type, typename history_type, int num_history>
class rw_visitor {
 public:
  using vertex_locator = typename graph_type::vertex_locator;

  rw_visitor()
      : vertex(),
        walk_length(0),
        history() {}

  explicit rw_visitor(vertex_locator _vertex)
      : vertex(_vertex), walk_length(0), id(-1), history() {}

  rw_visitor(vertex_locator _vertex, int _walk_length, int _id, const history_type &_history)
      : vertex(_vertex), walk_length(_walk_length), id(_id), history(_history) {}

  template <typename visitor_queue_handle, typename alg_data_type>
  bool init_visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    bool launched = false;
    if (alg_data.personalized()) {
      while (alg_data.rw_launcher(vertex)) {
#if WRITE_LOG
        LOG_OUT << "Personalized start from " << g.locator_to_label(vertex) << std::endl;
#endif
        vis_queue->queue_visitor(rw_visitor(vertex,
                                            1,
                                            g.locator_to_label(vertex),
                                            alg_data.new_history()));
        launched = true;
      }
    } else {
      if (alg_data.rw_launcher(vertex)) {
#if WRITE_LOG
        LOG_OUT << "Start from " << g.locator_to_label(vertex) << std::endl;
#endif
        vis_queue->queue_visitor(rw_visitor(vertex,
                                            1,
                                            g.locator_to_label(vertex),
                                            alg_data.new_history()));
        launched = true;
      }
    }
    return launched;
  }

  template <typename alg_data_type>
  bool pre_visit(alg_data_type &alg_data) const {
    if (walk_length == alg_data.max_walk_length()) {
#if WRITE_LOG
      LOG_OUT << "Has reached visit limit " << serialize(*(alg_data.graph())) << "\n";
#endif
      return false;
    }

    return true;
  }

  template <typename visitor_queue_handle, typename alg_data_type>
  bool visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    alg_data.increment_num_visits(vertex);

    if (g.degree(vertex) == 0) {
      auto target = alg_data.warp_machine();
      rw_visitor new_visitor = rw_visitor(target, walk_length + 1, id, append_history(g.locator_to_label(vertex)));
#if WRITE_LOG
      LOG_OUT << "Dead end (degree 0) at " << serialize(g) << ". Restart from " << new_visitor.serialize(g)
              << "\n";
#endif
      vis_queue->queue_visitor(new_visitor);
      return true;
    }

    if (g.degree(vertex) == 1) {
      const auto vitr = g.vertices_begin();
      const auto previous_vertex = history[std::min(walk_length - 1, num_history - 1)];
      if (g.locator_to_label(*vitr) == previous_vertex) {
        auto target = alg_data.warp_machine();
        rw_visitor new_visitor = rw_visitor(target, walk_length + 1, id, append_history(g.locator_to_label(vertex)));
#if WRITE_LOG
        LOG_OUT << "Dead end (degree 1) at " << serialize(g)
                << ". Restart from " << new_visitor.serialize(g)
                << "\n";
#endif
        vis_queue->queue_visitor(new_visitor);
        return true;
      }
    }

    if (alg_data.russian_roulette()) { // Die at here and increment the die score
      alg_data.increment_num_dead(vertex);
#if WRITE_LOG
      LOG_OUT << "Increment die score at " << serialize(g) << "\n";
#endif
      return false;
    }

    auto target = alg_data.neighbor_roulette(vertex, history);
    rw_visitor new_visitor = rw_visitor(target, walk_length + 1, id, append_history(g.locator_to_label(vertex)));
#if WRITE_LOG
    LOG_OUT << "Visit " << serialize(g) << " -> " << new_visitor.serialize(g) << "\n";
#endif

    vis_queue->queue_visitor(new_visitor);
    return true;
  }

  friend inline
  bool operator>(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex != v2.vertex) && !(v1.vertex < v2.vertex);
  }

  friend inline bool operator<(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex < v2.vertex);
  }

  history_type append_history(const uint64_t new_vertex) const {
    history_type new_history;

    for (int i = 0; i < num_history - 1; ++i) new_history[i] = history[i + 1];
    new_history[num_history - 1] = new_vertex;

    return new_history;
  }

  std::string serialize(const graph_type &g) const {
    std::stringstream ss;
    ss << "id " << id << ", lenght " << walk_length << ", vertex " << g.locator_to_label(vertex) << ", history";
    for (int i = 0; i < num_history; ++i) {
      ss << " [" << i << "] " << history[i];
    }
    return ss.str();
  }

  vertex_locator vertex;
  int walk_length;
  int id;
  history_type history;
};

void parse_cmd_line(int argc, char **argv,
                    std::string &input_filename,
                    bool &personalized,
                    uint64_t &seed_vertex,
                    uint64_t &num_walkers,
                    uint64_t &die_rate,
                    uint64_t &max_walk_length,
                    bool &warp_to_seed,
                    std::array<int, k_num_history> &closing_rates,
                    std::string &score_dump_file_prefix) {
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

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:ps:w:d:l:rc:o:h")) != -1) {
    switch (c) {
      case 'h':prn_help = true;
        break;
      case 'i':found_input_filename = true;
        input_filename = optarg;
        break;
      case 'p': personalized = true;
        break;
      case 's': seed_vertex = atoll(optarg);
        break;
      case 'w':num_walkers = atoll(optarg);
        break;
      case 'd':die_rate = atoll(optarg);
        break;
      case 'l':max_walk_length = atoll(optarg);
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

template <typename graph_type, typename score_type>
void compute_global_top_scores(const graph_type *const graph,
                               const std::function<score_type(typename graph_type::vertex_locator)> &score_function) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::map<score_type, uint64_t> top_vertices;
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    const auto score = score_function(*vitr);
    // std::cout << graph->locator_to_label(*vitr) << " >> " << score << std::endl;
    if (top_vertices.size() < num_top && 0 < score) {
      top_vertices.emplace(score, graph->locator_to_label(*vitr));
    } else if (top_vertices.begin()->first < score) {
      top_vertices.emplace(score, graph->locator_to_label(*vitr));
      top_vertices.erase(top_vertices.begin());
    }
  }

  std::vector<score_type> local_top_score;
  std::vector<uint64_t> local_top_id;
  for (const auto elem : top_vertices) {
    local_top_score.emplace_back(elem.first);
    local_top_id.emplace_back(elem.second);
  }

  std::vector<score_type> global_top_score;
  std::vector<uint64_t> global_top_id;
  havoqgt::mpi_all_gather(local_top_score, global_top_score, MPI_COMM_WORLD);
  havoqgt::mpi_all_gather(local_top_id, global_top_id, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::vector<std::pair<score_type, uint64_t>> global_top;
    for (uint64_t i = 0; i < global_top_score.size(); ++i) {
      global_top.emplace_back(std::make_pair(global_top_score.at(i), global_top_id.at(i)));
    }
    std::partial_sort(global_top.begin(),
                      global_top.begin() + std::min((uint64_t)num_top, (uint64_t)global_top.size()),
                      global_top.end(),
                      [](const std::pair<score_type, uint64_t> &lh, std::pair<score_type, uint64_t> &rh) -> bool {
                        return (lh.first > rh.first);
                      });
    std::cout << "Vertex ID : Score" << std::endl;
    for (int i = 0; i < std::min((uint64_t)num_top, (uint64_t)global_top.size()); ++i) {
      std::cout << global_top[i].second << " : " << global_top[i].first << std::endl;
    }
  }
}

template <typename graph_type, typename score_type>
void dump_score(const graph_type *const graph, const std::string file_name,
                const std::function<score_type(typename graph_type::vertex_locator)> &score_function) {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::ofstream ofs(file_name + "_" + std::to_string(mpi_rank));
  for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
    ofs << graph->locator_to_label(*vitr) << " " << score_function(*vitr) << "\n";
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
    std::array<int, k_num_history> closing_rates;
    std::string score_dump_file_prefix;

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
                   score_dump_file_prefix);

    MPI_Barrier(MPI_COMM_WORLD);
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    vertex_locator seed_vertex = graph->label_to_locator(seed_vertex_id);
    const std::size_t num_wk_per_process = (personalized) ? num_walkers : num_walkers / mpi_size;
    rw_algorithm_v1<graph_type, rw_traversal_history_type, k_num_history> algorithm_data(graph,
                                                                                         personalized,
                                                                                         seed_vertex,
                                                                                         num_wk_per_process,
                                                                                         die_rate,
                                                                                         max_walk_length,
                                                                                         warp_to_seed,
                                                                                         closing_rates);
#if WRITE_LOG
    ofs_log = new std::ofstream(std::string("./log_") + std::to_string(mpi_rank));
#endif
    auto vq = havoqgt::create_visitor_queue<rw_visitor<graph_type, rw_traversal_history_type, k_num_history>,
                                            havoqgt::detail::visitor_priority_queue>(graph, algorithm_data);

    MPI_Barrier(MPI_COMM_WORLD);
    const double start_time = MPI_Wtime();
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD);
    const double end_time = MPI_Wtime();
#if WRITE_LOG
    ofs_log->close();
    delete ofs_log;
#endif

//    for (int i = 0 ; i < mpi_size; ++i)  {
//      if (i == mpi_rank)
//    for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
//      for (auto e = graph->edges_begin(*vitr), en = graph->edges_end(*vitr); e != en; ++e) {
//        assert(*vitr == e.source());
//        std::cout << graph->locator_to_label(*vitr) << " " << graph->locator_to_label(e.target()) << std::endl;
//      }
//    }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }

    const uint64_t global_num_walkers = havoqgt::mpi_all_reduce(algorithm_data.num_launched_walkers(),
                                                                std::plus<uint64_t>(),
                                                                MPI_COMM_WORLD);

    uint64_t local_total_dead = 0;
    for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
      local_total_dead += algorithm_data.num_dead(*vitr);
    }
    const uint64_t global_total_dead = havoqgt::mpi_all_reduce(local_total_dead, std::plus<uint64_t>(), MPI_COMM_WORLD);

    uint64_t local_total_visits = 0;
    for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
      local_total_visits += algorithm_data.num_visits(*vitr);
    }
    const uint64_t
        global_total_visits = havoqgt::mpi_all_reduce(local_total_visits, std::plus<uint64_t>(), MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      std::cout << "Execution time: " << end_time - start_time << std::endl;
      std::cout << "#launched wk: " << global_num_walkers << std::endl;
      std::cout << "sum of dead: " << global_total_dead << std::endl;
      std::cout << "sum of visits: " << global_total_visits << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const auto dead_score_func = [&algorithm_data, &global_total_dead](vertex_locator vertex) -> double {
      return static_cast<double>(algorithm_data.num_dead(vertex)) / global_total_dead;
    };

    const auto visit_score_func = [&algorithm_data, &global_total_visits](vertex_locator vertex) -> double {
      return static_cast<double>(algorithm_data.num_visits(vertex)) / global_total_visits;
    };

    if (mpi_rank == 0) std::cout << "\nDead score" << std::endl;
    compute_global_top_scores<graph_type, double>(graph, dead_score_func);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) std::cout << "\nVisit score" << std::endl;
    compute_global_top_scores<graph_type, double>(graph, visit_score_func);
    MPI_Barrier(MPI_COMM_WORLD);

    if (!score_dump_file_prefix.empty()) {
      if (mpi_rank == 0) std::cout << "\nDumping dead score" << std::endl;
      dump_score<graph_type, double>(graph,
                                     score_dump_file_prefix + std::string("_dead_score_") + std::to_string(global_num_walkers),
                                     dead_score_func);
      MPI_Barrier(MPI_COMM_WORLD);

      if (mpi_rank == 0) std::cout << "\nDumping visit score" << std::endl;
      dump_score<graph_type, double>(graph,
                                     score_dump_file_prefix + std::string("_visit_score_") + std::to_string(global_num_walkers),
                                     visit_score_func);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  return 0;
}
