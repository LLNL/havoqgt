/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */


#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/map.hpp>
#include <boost/range/algorithm.hpp>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/offset_ptr.hpp>
#include <boost/interprocess/containers/set.hpp>

#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include <havoqgt/graphstore/graphstore_rhhda.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#include "dynamicgraphstore_bench.hpp"
#include "bfs_bench.hpp"

#define VERBOSE 0

#define DEBUG_MODE 0
#if DEBUG_MODE
std::ofstream ofs_edges;
#endif

/// --- typenames --- ///
using vertex_type        = uint64_t;

/// adjacency list
using edge_property_type    = unsigned char;
using edge_vec_element_type = std::pair<vertex_type, edge_property_type>;
using vec_allocator_type    = boost::interprocess::allocator<edge_vec_element_type, segment_manager_type>;
using edge_vec_type         = boost::interprocess::vector<edge_vec_element_type, vec_allocator_type>;

/// vertex table
using vertex_property_type  = bool;
using map_value_type        = std::pair<vertex_property_type, edge_vec_type>;
using map_element_type      = std::pair<const vertex_type, map_value_type>;
using map_allocator_type    = boost::interprocess::allocator<map_element_type, segment_manager_type>;
using graphstore_type       = boost::unordered_map<vertex_type,
                                                   map_value_type,
                                                   boost::hash<vertex_type>,
                                                   std::equal_to<vertex_type>,
                                                   map_allocator_type>;

/// --- global variables --- ///
vertex_type max_vertex_id_ = 0;
size_t num_edges_ = 0;

template <typename edges_type, typename allocator_type>
void constract_graph(mapped_file_type& mapped_file,
                     segment_manager_type *const segment_manager,
                     graphstore_type& graphstore,
                     allocator_type& allocator,
                     edges_type& edges, const size_t chunk_size)
{
  std::cout << "-- Disp status of before generation --" << std::endl;
  print_usages(segment_manager);

  request_vector_type<vertex_type> update_request_vec = request_vector_type<vertex_type>();

  uint64_t count_inserted = 0;
  uint64_t count_duplicated = 0;

  double construction_time = 0;
  size_t loop_cnt = 0;

  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();

    auto global_start = graphstore::utility::duration_time();
  while (edges_itr != edges_itr_end) {
    std::cout << "[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    update_request_vec.clear();
    generate_insertion_requests(edges_itr, edges_itr_end, chunk_size, update_request_vec, 0);

    auto local_start = graphstore::utility::duration_time();
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      max_vertex_id_ = std::max(max_vertex_id_, edge.first);
      max_vertex_id_ = std::max(max_vertex_id_, edge.second);

      edge_property_type dummy_weight = 0;

      /// ------- insertion core ------- ///
      auto value = graphstore.find(edge.first);
      if (value == graphstore.end()) { // new vertex
        const vertex_property_type vp = false;
        map_value_type map_value(vp,
                                 edge_vec_type(1, edge_vec_element_type(edge.second, dummy_weight), allocator));
        graphstore.insert(map_element_type(edge.first, map_value));
        ++count_inserted;
      } else {
        edge_vec_type& edge_vec = value->second.second;
        edge_vec_element_type trg_edge(edge.second, dummy_weight);
        if (boost::find<edge_vec_type>(edge_vec, trg_edge) != edge_vec.end() ) {
          ++count_duplicated;
          continue;
        }
        edge_vec.push_back(trg_edge);
        ++count_inserted;
      }
    } /// end of chunk
    double t = graphstore::utility::duration_time_sec(local_start);
    construction_time += t;
    std::cout << "progress (sec.): " << t << std::endl;

    ++loop_cnt;
  }
  // flush_mmmap(mapped_file);
  sync_dimmap();
  const double whole_construction_time = graphstore::utility::duration_time_sec(global_start);

  num_edges_ = count_inserted;

  std::cout << "\n-- All edge updations done --" << std::endl;
  std::cout << "inserted edges : " << count_inserted << std::endl;
  std::cout << "depulicated edges : " << count_duplicated << std::endl;
  std::cout << "construction time (insertion only) : " << construction_time << std::endl;
  std::cout << "whole construction time : " << whole_construction_time << std::endl;
  print_usages(segment_manager);

}


///
/// \brief The bfs_core<graphstore_type, vertex_type, _Tp3> class
/// bfs core for baseline model
template <typename graphstore_type, typename vertex_type>
void run_bfs_sync (graphstore_type& graph,
                   trv_inf<vertex_type>& inf,
                   std::queue<vertex_type>& frontier_queue,
                   std::queue<vertex_type>& next_queue,
                   vertex_type& start_vrtx)
  {

    /// ---- init inf ---- ///
    auto tic_init = graphstore::utility::duration_time();
    inf.init(BFS_USE_BITMAP);
#if BFS_USE_BITMAP
    inf.is_visited[start_vrtx] = true;
#else
    for (auto itr : graph) {
      itr.second.first = false;
    }
    graph.find(start_vrtx)->second.first = true;
#endif
    std::cout << "Init time (sec.):\t"  << graphstore::utility::duration_time_sec(tic_init) << std::endl;

    /// --- BFS main loop -------- ///
    size_t level = 0;
    while (true) {
      std::cout << "Lv. " << level << ", size of frontier queue =\t" << frontier_queue.size() << std::endl;

      /// loop for current frontier
      while (!frontier_queue.empty()) {
        vertex_type src = frontier_queue.front();
        frontier_queue.pop();
        ++inf.count_visited_vertices;

        /// push adjacent vertices to the next queue
        auto adjlist_vec = graph.find(src)->second.second;
        for (const auto edge : adjlist_vec) {
          const vertex_type dst = edge.first;
#if BFS_USE_BITMAP
          bool& is_visited = inf.is_visited[dst];
#else
          bool& is_visited = graph.find(dst)->second.first;
#endif
          if (!is_visited) {
            next_queue.push(dst);
            /// inf.tree[dst] = src;
            is_visited = true;
          }
          ++(inf.count_visited_edges);
        }
      }  /// end of loop for a frontier

      if (next_queue.empty()) break; /// termination condition
      frontier_queue.swap(next_queue);
      ++level;
    } /// end of BFS loop
  }

/// Avoid linker errors with template function
template void run_bfs_sync<graphstore_type, vertex_type>(graphstore_type&,
                                                         trv_inf<vertex_type>&,
                                                         std::queue<vertex_type>&,
                                                         std::queue<vertex_type>&,
                                                         vertex_type&);


/// --- option variables --- ///
std::string fname_segmentfile_;
std::vector<std::string> fname_edge_list_;
size_t segment_size_log2_ = 30;
std::vector<vertex_type> source_list_;

void parse_options(int argc, char **argv)
{

  char c;

  while ((c = getopt (argc, argv, "s:f:e:r:")) != -1) {
    switch (c) {
      case 's':
        segment_size_log2_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'f':
        fname_segmentfile_ = optarg;
        break;

      case 'e':
      {
        std::string fname(optarg);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while (std::getline(fin, line)) {
          fname_edge_list_.push_back(line);
        }
        break;
      }

      case 'r':
      {
        std::string buf;
        std::stringstream sstrm(optarg);
        while (std::getline(sstrm, buf, ':'))
          source_list_.push_back(boost::lexical_cast<size_t>(buf));
        break;
      }
    }
  }

}


int main(int argc, char** argv)
{

  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  parse_options(argc, argv);
  if (!fname_edge_list_.empty())
    std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
  for (auto itr : fname_edge_list_) {
    std::cout << "Load edge list from " << itr << std::endl;
  }


  /// --- create a segument file --- ///
  std::cout << "Delete segment file: " << fname_segmentfile_ << std::endl;
  boost::interprocess::file_mapping::remove(fname_segmentfile_.c_str());
  std::cout << "\n<<Construct segment>>" << std::endl;
  std::cout << "Create and map a segument file" << std::endl;
  uint64_t graph_capacity = std::pow(2, segment_size_log2_);
  mapped_file_type mapped_file = mapped_file_type(
                                   boost::interprocess::create_only,
                                   fname_segmentfile_.c_str(),
                                   graph_capacity);
  std::cout << "Call posix_fallocate\n";
  fallocate(fname_segmentfile_.c_str(), graph_capacity, mapped_file);

  /// --- Get a segument manager --- ///
  std::cout << "\n<Get a segment manager>" << std::endl;
  segment_manager_type* segment_manager = mapped_file.get_segment_manager();
  print_usages(segment_manager);


  /// --- Allocate graphstore_rhh_matrix --- ///
  std::cout << "\n<Allocate graphstore_rhh_matrix>" << std::endl;
  boost::interprocess::allocator<void, segment_manager_type> boost_seg_allocator(segment_manager);
  graphstore_type graphstore(boost_seg_allocator);


  /// --- Graph Construction --- ////
  std::cout << "\n<Construct graph>" << std::endl;
  havoqgt::havoqgt_init(&argc, &argv);
  {
    havoqgt::get_environment();

    havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_);
    constract_graph(
          mapped_file,
          segment_manager,
          graphstore,
          boost_seg_allocator,
          edgelist,
          static_cast<uint64_t>(std::pow(10, 6)));
  }


  /// ---------- Graph Traversal --------------- ///
  std::cout << "\n<Run BFS>" << std::endl;
  if (source_list_.empty())
    generate_source_list(4, max_vertex_id_, source_list_);
  run_bfs(graphstore, max_vertex_id_, num_edges_, source_list_);

  return 0;
}
