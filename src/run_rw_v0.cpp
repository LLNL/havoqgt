//
// Created by Iwabuchi, Keita on 10/29/18.
//

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

#define WRITE_LOG 1
#if WRITE_LOG
std::ofstream *ofs_log;
#endif

template <typename graph_type>
class rw_algorithm_v0 {
 private:
  using vertex_locator = typename graph_type::vertex_locator;
  using count_table_type = typename graph_type::template vertex_data<uint64_t, std::allocator<uint64_t>>;

 public:
  explicit rw_algorithm_v0(graph_type *const graph,
                           const vertex_locator seed,
                           const vertex_locator goal,
                           const std::size_t num_walkers,
                           const uint64_t die_rate,
                           const uint64_t max_walk_length,
                           const bool warp_to_seed)
      : m_graph(graph),
        m_seed_locator(seed),
        m_goal_locator(goal),
        m_num_walkers(num_walkers),
        m_die_rate(die_rate),
        m_max_walk_length(max_walk_length),
        m_count_walkers(0),
        m_gen(std::chrono::system_clock::now().time_since_epoch().count()),
        m_count_visited_goal(0),
        m_warp_to_seed(warp_to_seed),
        m_num_visited_wk_table(*graph) {
    m_num_visited_wk_table.reset(0);
  }

  bool rw_launcher(const vertex_locator vertex) {
    if (vertex == m_seed_locator && m_count_walkers < m_num_walkers) {
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
    if (m_warp_to_seed) {
      return m_seed_locator;
    } else {
      std::uniform_int_distribution<uint64_t> dis(0, m_graph->max_global_vertex_id() - 1);
      return m_graph->label_to_locator(dis(m_gen));
    }
  }

  vertex_locator neighbor_roulette(const vertex_locator vertex) {
    assert(m_graph->degree(vertex) > 0);
    std::uniform_int_distribution<uint64_t> dis(0, m_graph->degree(vertex) - 1);
    const uint64_t offset = dis(m_gen);
    auto edge = m_graph->edges_begin(vertex) + offset;
    return edge.target();
  }

  uint64_t max_walk_length() const {
    return m_max_walk_length;
  }

  void increment_arrived_goal() {
    ++m_count_visited_goal;
  }

  uint64_t num_arrived_goal() const {
    return m_count_visited_goal;
  }

  void increment_visited_walkers(const vertex_locator vertex) {
    ++m_num_visited_wk_table[vertex];
  }

  uint64_t num_visited_walkers(const vertex_locator vertex) const {
    return m_num_visited_wk_table[vertex];
  }

  graph_type *graph() {
    return m_graph;
  }

 private:
  graph_type *m_graph;
  vertex_locator m_seed_locator;
  vertex_locator m_goal_locator;
  uint64_t m_num_walkers;
  uint64_t m_die_rate;
  uint64_t m_max_walk_length;
  uint64_t m_count_walkers;
  std::mt19937 m_gen;
  uint64_t m_count_visited_goal;
  bool m_warp_to_seed;
  count_table_type m_num_visited_wk_table;
};

template <typename graph_type>
class rw_visitor {
 public:
  using vertex_locator = typename graph_type::vertex_locator;
  rw_visitor()
      : vertex(),
        walk_length(0) {}

  explicit rw_visitor(vertex_locator _vertex)
      : vertex(_vertex), walk_length(0), id(-1) {}

  rw_visitor(vertex_locator _vertex, int _walk_length, int _id)
      : vertex(_vertex), walk_length(_walk_length), id(_id) {}

  template <typename visitor_queue_handle, typename alg_data_type>
  bool init_visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    int count = 0;
    while (alg_data.rw_launcher(vertex)) {
      const auto target = alg_data.neighbor_roulette(vertex);
#if WRITE_LOG
      *ofs_log << "Start from " << g.locator_to_label(target) << std::endl;
#endif
      vis_queue->queue_visitor(rw_visitor(target, 0, count));
      ++count;
    }
    return count > 0;
  }

  template <typename alg_data_type>
  bool pre_visit(alg_data_type &alg_data) const {
    if (walk_length == alg_data.max_walk_length()) {
#if WRITE_LOG
      *ofs_log << "Has reached visit limit " << serialize(*(alg_data.graph())) << "\n";
#endif
      return false;
    }
    if (alg_data.arrived_goal(vertex)) {
      alg_data.increment_arrived_goal();
#if WRITE_LOG
      *ofs_log << "Has arrived goal " << serialize(*(alg_data.graph())) << "\n";
#endif
      return false;
    }
    assert(walk_length < alg_data.max_walk_length());

    return true;
  }

  template <typename visitor_queue_handle, typename alg_data_type>
  bool visit(graph_type &g, visitor_queue_handle vis_queue, alg_data_type &alg_data) const {
    rw_visitor new_visitor;

    if (g.degree(vertex) <= 1) {
      auto target = alg_data.warp_machine();
      new_visitor = rw_visitor(target, walk_length, id);
#if WRITE_LOG
      *ofs_log << "Dead end at " << serialize(g) << ". Restart from " << new_visitor.serialize(g)
               << "\n";
#endif

    } else {
      if (alg_data.russian_roulette()) { // Die at here and warp to somewhere
        auto target = alg_data.warp_machine();
        new_visitor = rw_visitor(target, walk_length + 1, id);
#if WRITE_LOG
        *ofs_log << "Warp " << serialize(g) << " -> " << new_visitor.serialize(g) << "\n";
#endif
      } else {
        auto target = alg_data.neighbor_roulette(vertex);
        new_visitor = rw_visitor(target, walk_length + 1, id);
#if WRITE_LOG
        *ofs_log << "Visit " << serialize(g) << " -> " << new_visitor.serialize(g) << "\n";
#endif
      }
      alg_data.increment_visited_walkers(vertex);
    }

    vis_queue->queue_visitor(new_visitor);
    return true;
  }

  friend inline bool operator>(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex != v2.vertex) && !(v1.vertex < v2.vertex);
  }

  friend inline bool operator<(const rw_visitor &v1, const rw_visitor &v2) {
    return (v1.vertex < v2.vertex);
  }

  std::string serialize(const graph_type& g) const {
    std::stringstream ss;
    ss << id  << " : " << walk_length << " : " << g.locator_to_label(vertex);
    return ss.str();
  }

  vertex_locator vertex;
  int walk_length;
  int id;
};

void parse_cmd_line(int argc, char **argv,
                    std::string &input_filename,
                    uint64_t &seed_vertex,
                    uint64_t &goal_vertex,
                    uint64_t &num_walkers,
                    uint64_t &die_rate,
                    uint64_t &max_walk_length,
                    bool &warp_to_seed) {
  if (havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  // Initial values
  bool found_input_filename = false;
  seed_vertex = 0;
  goal_vertex = 100;
  num_walkers = 1024;
  die_rate = 10;
  max_walk_length = 128;
  warp_to_seed = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "i:s:g:w:d:l:rh ")) != -1) {
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
      case 'w':num_walkers = atoll(optarg);
        break;
      case 'd':die_rate = atoll(optarg);
        break;
      case 'l':max_walk_length = atoll(optarg);
        break;
      case 'r': warp_to_seed = true;
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
              << "\nnum_walkers: " << num_walkers
              << "\ndie_rate: " << die_rate
              << "\nmax_walk_length: " << max_walk_length
              << "\nwarp_to_seed: " << warp_to_seed << std::endl;
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
    uint64_t num_walkers;
    uint64_t die_rate;
    uint64_t max_walk_length;
    bool warp_to_seed;
    parse_cmd_line(argc,
                   argv,
                   graph_input,
                   seed_vertex_id,
                   goal_vertex_id,
                   num_walkers,
                   die_rate,
                   max_walk_length,
                   warp_to_seed);

    MPI_Barrier(MPI_COMM_WORLD);
    havoqgt::distributed_db ddb(havoqgt::db_open(), graph_input.c_str());
    graph_type *graph = ddb.get_segment_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Graph Loaded Ready." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    seed_vertex_id = select_seed(graph, seed_vertex_id);

    vertex_locator seed_vertex(graph->label_to_locator(seed_vertex_id));
    vertex_locator goal_vertex(graph->label_to_locator(goal_vertex_id));
    rw_algorithm_v0<graph_type>
        algorithm_data(graph, seed_vertex, goal_vertex, num_walkers, die_rate, max_walk_length, warp_to_seed);
#if WRITE_LOG
    ofs_log = new std::ofstream(std::string("./log_") + std::to_string(mpi_rank));
#endif
    auto vq = havoqgt::create_visitor_queue<rw_visitor<graph_type>, havoqgt::detail::visitor_priority_queue>(graph,
                                                                                                             algorithm_data);

    MPI_Barrier(MPI_COMM_WORLD);
    const double start_time = MPI_Wtime();
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD);
    const double end_time = MPI_Wtime();
#if WRITE_LOG
    ofs_log->close();
    delete ofs_log;
#endif

    const uint64_t global_num_arrived_goal = havoqgt::mpi_all_reduce(algorithm_data.num_arrived_goal(),
                                                                     std::plus<uint64_t>(),
                                                                     MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Execution time: " << end_time - start_time << std::endl;
      std::cout << "#reached wk: " << global_num_arrived_goal << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::map<uint64_t, uint64_t> top_vertices;
    for (auto vitr = graph->vertices_begin(), end = graph->vertices_end(); vitr != end; ++vitr) {
      const uint64_t num_visited = algorithm_data.num_visited_walkers(*vitr);
      if (top_vertices.size() < 10) {
        top_vertices.emplace(num_visited, graph->locator_to_label(*vitr));
      } else if (top_vertices.begin()->first < num_visited) {
        top_vertices.emplace(num_visited, graph->locator_to_label(*vitr));
        top_vertices.erase(top_vertices.begin());
      }
    }

    if (mpi_rank == 0) std::cout << "#visited : Vertex ID" << std::endl;
    for (int i = 0; i < mpi_size; ++i) {
      if (mpi_rank == i) {
        std::cout << "--------------------" << std::endl;
        for (const auto elem : top_vertices) {
          std::cout << elem.first << " : " << elem.second << std::endl;
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  return 0;
}
