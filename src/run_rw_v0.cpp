//
// Created by Iwabuchi, Keita on 10/29/18.
//

#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <tuple>
#include <cassert>

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

template <typename graph_type>
class rw_algorithm_v0 {
 private:
  using vertex_locator = typename graph_type::vertex_locator;

 public:
  explicit rw_algorithm_v0(graph_type *const graph,
                           const vertex_locator seed,
                           const vertex_locator goal,
                           const std::size_t max_num_walkers,
                           const uint64_t die_rate,
                           const uint64_t num_max_visits)
      : m_graph(graph),
        m_seed_locator(seed),
        m_goal_locator(goal),
        m_max_num_walkers(max_num_walkers),
        m_die_rate(die_rate),
        m_num_max_visits(num_max_visits),
        m_count_walkers(0),
        m_gen(123),
        m_count_visited_goal(0) {}

  bool rw_launcher(const vertex_locator vertex) {
    if (vertex == m_seed_locator && m_count_walkers < m_max_num_walkers) {
      ++m_count_walkers;
      return true;
    }
    return false;
  }

  bool arrived_goal(const vertex_locator vertex) {
    return m_goal_locator == vertex;
  }

  bool russian_roulette() {
    std::uniform_int_distribution<int> dis(0, 99);
    return (dis(m_gen) < m_die_rate);
  }

  vertex_locator warp_machine() {
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->max_global_vertex_id() - 1);
    return m_graph->label_to_locator(dis(m_gen));
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex) {
    assert(m_graph->degree(vertex) > 0);
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->degree(vertex) - 1);
    const uint64_t offset = dis(m_gen);
    auto edge = m_graph->edges_begin(vertex) + offset;
    return edge.target();
  }

  uint64_t max_visits() const {
    return m_num_max_visits;
  }

  void increment_arrived_goal() {
    ++m_count_visited_goal;
  }

  uint64_t num_arrived_goal() const {
    return m_count_visited_goal;
  }

 private:
  graph_type *m_graph;
  vertex_locator m_seed_locator;
  vertex_locator m_goal_locator;
  uint64_t m_max_num_walkers;
  uint64_t m_die_rate;
  uint64_t m_num_max_visits;
  uint64_t m_count_walkers;
  std::mt19937 m_gen;
  uint64_t m_count_visited_goal;
};

template <typename graph_type>
class rw_visitor {
 public:
  using vertex_locator = typename graph_type::vertex_locator;
  rw_visitor()
      : vertex(),
        num_visited(0) {}

  explicit rw_visitor(vertex_locator _vertex)
      : vertex(_vertex), num_visited(0) {}

  rw_visitor(vertex_locator _vertex, int _num_visited)
      : vertex(_vertex), num_visited(_num_visited) {}

  template <typename visitor_queue_handle, typename alg_data_type>
  bool init_visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    bool launched = false;
    while (alg_data.rw_launcher(vertex)) {
      vis_queue->queue_visitor(rw_visitor(alg_data.neighbor_roulette(vertex), 0));
      launched = true;
    }
    return launched;
  }

  template <typename alg_data_type>
  bool pre_visit(alg_data_type &alg_data) const {
    if (num_visited == alg_data.max_visits()) {
      // std::cout << "Reached visit limit" << std::endl;
      return false;
    }
    if (alg_data.arrived_goal(vertex)) {
      // std::cout << "Arrived goal" << std::endl;
      alg_data.increment_arrived_goal();
      return false;
    }
    assert(num_visited < alg_data.max_visits());

    return true;
  }

  template <typename visitor_queue_handle, typename alg_data_type>
  bool visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    vertex_locator target;
    if (alg_data.russian_roulette()) { // Die at here and warp to somewhere
      target = alg_data.warp_machine();
      // std::cout << "Warp" << std::endl;
    } else {
      target = alg_data.neighbor_roulette(vertex);
    }
    std::cout << g.locator_to_label(vertex) << " -> " << g.locator_to_label(target) << std::endl;

    rw_visitor new_visitor(target, num_visited + 1);
    vis_queue->queue_visitor(new_visitor);
    return true;
  }

  friend inline bool operator>(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex != v2.vertex) && !(v1.vertex < v2.vertex);
  }

  friend inline bool operator<(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex < v2.vertex);
  }

  vertex_locator vertex;
  int num_visited;
};

void parse_cmd_line(int argc, char **argv,
                    std::string &input_filename,
                    uint64_t &seed_vertex,
                    uint64_t &goal_vertex,
                    uint64_t &max_num_walkers,
                    uint64_t &die_rate,
                    uint64_t &num_max_visits) {
  if (havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_input_filename = false;
  seed_vertex = 0;
  goal_vertex = 100;
  max_num_walkers = 1024;
  die_rate = 10;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:g:w:d:v:h ")) != -1) {
    switch (c) {
      case 'h':prn_help = true;
        break;
      case 'i':found_input_filename = true;
        input_filename = optarg;
        break;
      case 's':seed_vertex = atoll(optarg);
        break;
      case 'g':goal_vertex = atoll(optarg);
        break;
      case 'w':max_num_walkers = atoll(optarg);
        break;
      case 'd':die_rate = atoll(optarg);
        break;
      case 'v':num_max_visits = atoll(optarg);
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
              << "\nseed_vertex: " << seed_vertex
              << "\ngoal_vertex: " << goal_vertex
              << "\nmax_num_walkers: " << max_num_walkers
              << "\ndie_rate: " << die_rate
              << "\nnum_max_visits: " << num_max_visits << std::endl;
  }
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
    uint64_t seed_vertex_id;
    uint64_t goal_vertex_id;
    uint64_t max_num_walkers;
    uint64_t die_rate;
    uint64_t num_max_visits;
    parse_cmd_line(argc, argv, graph_input, seed_vertex_id, goal_vertex_id, max_num_walkers, die_rate, num_max_visits);
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());

    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    vertex_locator seed_vertex(graph->label_to_locator(seed_vertex_id));
    if (graph->degree(seed_vertex) == 0) {
      std::cout << "Select another source vertex" << std::endl;
      for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
        if (graph->degree(*vitr) > 0) {
          seed_vertex = *vitr;
          break;
        }
      }
    }
    std::cout << "Source: " << graph->locator_to_label(seed_vertex)
              << ", Degree: " << graph->degree(seed_vertex) << std::endl;

    vertex_locator goal_vertex(graph->label_to_locator(goal_vertex_id));
    rw_algorithm_v0<graph_type> algorithm_data(graph, seed_vertex, goal_vertex, max_num_walkers, die_rate, num_max_visits);
    auto vq = havoqgt::create_visitor_queue<rw_visitor<graph_type>, havoqgt::detail::visitor_priority_queue>(graph,
                                                                                                             algorithm_data);
    const double start_time = MPI_Wtime();
    vq.init_visitor_traversal();
    const double end_time = MPI_Wtime();
    std::cout << "Execution time: " << end_time - start_time << std::endl;
    std::cout << "#reached wk: " << algorithm_data.num_arrived_goal() << std::endl;
  }
}