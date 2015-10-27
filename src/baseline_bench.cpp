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

#include "dynamicgraphstore_bench.hpp"

#define VERBOSE 0

#define DEBUG_MODE 0
#if DEBUG_MODE
std::ofstream ofs_edges;
#endif

/// --- typenames --- ///
using vertex_id_type         = uint64_t;

/// adjacency list
using edge_property_type    = unsigned char;
using edge_vec_element_type = std::pair<vertex_id_type, edge_property_type>;
using vec_allocator_type    = boost::interprocess::allocator<edge_vec_element_type, segment_manager_type>;
using edge_vec_type         = boost::interprocess::vector<edge_vec_element_type, vec_allocator_type>;

/// vertex table
using vertex_property_type  = bool;
using map_value_type        = std::pair<vertex_property_type, edge_vec_type>;
using map_element_type      = std::pair<const vertex_id_type, map_value_type>;
using map_allocator_type    = boost::interprocess::allocator<map_element_type, segment_manager_type>;
using graphstore_type       = boost::unordered_map<vertex_id_type,
                                                   map_value_type,
                                                   boost::hash<vertex_id_type>,
                                                   std::equal_to<vertex_id_type>,
                                                   map_allocator_type>;


template <typename Edges, typename allocator_type>
void apply_edges_update_requests(mapped_file_type& mapped_file,
                                 segment_manager_type *const segment_manager,
                                 graphstore_type& graphstore,
                                 allocator_type& allocator,
                                 Edges& edges,
                                 const size_t chunk_size,
                                 int edges_delete_ratio)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
  havoqgt::havoqgt_env()->world_comm().barrier();

  if (mpi_rank == 0) std::cout << "-- Disp status of before generation --" << std::endl;
  if (mpi_rank == 0) print_usages(segment_manager);
  havoqgt::havoqgt_env()->world_comm().barrier();

  request_vector_type<vertex_id_type> update_request_vec = request_vector_type<vertex_id_type>();

  size_t count_inserted = 0;
  size_t count_delete = 0;
  uint64_t loop_cnt = 0;

  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  bool global_is_finished = false;

  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    update_request_vec.clear();
    generate_insertion_requests(edges_itr, edges_itr_end, chunk_size, update_request_vec, edges_delete_ratio);
    havoqgt::havoqgt_env()->world_comm().barrier();

    uint64_t count_inserted = 0;
    uint64_t count_duplicated = 0;
    unsigned char dummy_weight = 0;

    const double time_start = MPI_Wtime();
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      if (request.is_delete) {
        assert(false);
      } else {

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
      }
    }

    // flush_mmmap(mapped_file);
    sync_dimmap();
    havoqgt::havoqgt_env()->world_comm().barrier();

    const double time_end = MPI_Wtime();

    if (mpi_rank == 0) std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;
    if (mpi_rank == 0) std::cout << "ins =\t" << count_inserted << std::endl;
    if (mpi_rank == 0) std::cout << "dup =\t" << count_duplicated << std::endl;
    if (mpi_rank == 0) print_usages(segment_manager);
    havoqgt::havoqgt_env()->world_comm().barrier();

#if VERBOSE
    for (int i = 0; i < mpi_size; ++i) {
      if (i == mpi_rank) {
        std::cout << "[" << mpi_rank << "]" << std::endl;
        graph_store.print_status();
      }
      havoqgt::havoqgt_env()->world_comm().barrier();
    }
#endif

    ++loop_cnt;

    const bool local_is_finished = (edges_itr == edges_itr_end);
    MPI_Allreduce(&local_is_finished, &global_is_finished, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  if (mpi_rank == 0) {
    std::cout << "\n-- All edge updations done --" << std::endl;
    print_usages(segment_manager);
  }
  havoqgt::havoqgt_env()->world_comm().barrier();

  for (int i = 0; i < mpi_size; ++i) {
    if (i == mpi_rank) {
      std::cout << "[" << mpi_rank << "] inserted edges : " << count_inserted << std::endl;
      std::cout << "[" << mpi_rank << "] deleted edges : " << count_delete << std::endl;
    }
  }

}


/// --- option variables --- ///
uint64_t vertex_scale_               = 18;
uint64_t edge_factor_                = 16;
uint64_t segmentfile_init_size_log2_ = 30;
bool     is_delete_segmentfile_on_exit_ = false;
uint64_t chunk_size_log10_           = 6;
int edges_delete_ratio_              = 0;
std::string fname_segmentfile_;
std::vector<std::string> fname_edge_list_;

void parse_options(int argc, char **argv)
{

  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
    havoqgt::get_environment().print();
  }

<<<<<<< HEAD
    uint64_t num_vertices = 1;
    uint64_t vert_scale;
    uint64_t edge_factor;
    uint32_t delete_file;
    uint64_t chunk_size_log10;
    uint64_t segmentfile_init_size_log2;
    uint64_t edges_delete_ratio = 0;
    std::string fname_output;
    std::vector<std::string> fname_edge_list;

    if (argc < 8) {
      std::cerr << "usage: <Scale> <Edge factor> <segmentfile name>"
                << " <segmentfile_init_size_log2 (log2)> <delete file on exit>"
                << " <chunk_size_log10> <edges_delete_ratio>"
                << " <edgelist file>"
                << " (argc:" << argc << " )." << std::endl;
      exit(-1);
    } else {
      int pos = 1;
      vert_scale            = boost::lexical_cast<uint64_t>(argv[pos++]);
      edge_factor           = boost::lexical_cast<uint64_t>(argv[pos++]);
      fname_output          = argv[pos++];
      segmentfile_init_size_log2 = boost::lexical_cast<uint64_t>(argv[pos++]);
      delete_file           = boost::lexical_cast<uint32_t>(argv[pos++]);
      chunk_size_log10      = boost::lexical_cast<uint64_t>(argv[pos++]);
      edges_delete_ratio    = boost::lexical_cast<uint64_t>(argv[pos++]);
      if (pos < argc) {
        std::string fname(argv[pos++]);
=======
  char c;
  while ((c = getopt (argc, argv, "s:e:dc:o:r:i:")) != -1) {
    switch (c) {
      case 's':
        vertex_scale_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'e':
        edge_factor_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'o':
        fname_segmentfile_ = optarg;
        break;

      case 'd':
        is_delete_segmentfile_on_exit_ = true;
        break;

      case 'c':
        chunk_size_log10_ = boost::lexical_cast<size_t>(optarg);
        break;

      case 'r':
        edges_delete_ratio_ = boost::lexical_cast<int>(optarg);
        break;

      case 'i':
      {
        std::string fname(optarg);
>>>>>>> origin/develop_keita
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        size_t num_files = 0;
        while (std::getline(fin, line)) {
          ++num_files;
        }
        const size_t num_local_files = (num_files + mpi_size - 1) / mpi_size;
        fin.close();
        fin.open(fname);
        size_t count = 0;
        while (std::getline(fin, line)) {
<<<<<<< HEAD
          if (num_local_files * mpi_rank <= count && count < num_local_files * (mpi_rank + 1)) {
            fname_edge_list.push_back(line);
=======
         if (num_local_files * mpi_rank <= count && count < num_local_files * (mpi_rank + 1)) {
            fname_edge_list_.push_back(line);
>>>>>>> origin/develop_keita
          }
          ++count;
        }
        break;
      }

    }
  }
}

int main(int argc, char** argv)
{

  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
    havoqgt::get_environment();
    havoqgt::havoqgt_env()->world_comm().size();

    parse_options(argc, argv);
    uint64_t num_vertices = (1ULL << vertex_scale_);

    if (mpi_rank == 0) {
<<<<<<< HEAD
      std::cout << "Segment file name = " << fname_output << std::endl;
      std::cout << "Initialize segment filse size (log2) = " << segmentfile_init_size_log2 << std::endl;
      std::cout << "Delete on Exit = " << delete_file << std::endl;
      std::cout << "Chunk size (log10) = " << chunk_size_log10 << std::endl;
      std::cout << "Edges Delete Ratio = " << edges_delete_ratio << std::endl;
      //      std::cout << "Midle-high degree threshold = " << midle_high_degree_threshold << std::endl;

      if (fname_edge_list.empty()) {
        std::cout << "Building RMAT graph Scale: " << vert_scale << std::endl;
        std::cout << "Building RMAT graph Edge factor: " << edge_factor << std::endl;
=======
      std::cout << "Segment file name = " << fname_segmentfile_ << std::endl;
      std::cout << "Initialize segment filse size (log2) = " << segmentfile_init_size_log2_ << std::endl;
      std::cout << "Delete on Exit = " << is_delete_segmentfile_on_exit_ << std::endl;
      std::cout << "Chunk size (log10) = " << chunk_size_log10_ << std::endl;
      std::cout << "Edges Delete Ratio = " << chunk_size_log10_ << std::endl;
      if (fname_edge_list_.empty()) {
        std::cout << "Building RMAT graph Scale: " << vertex_scale_ << std::endl;
        std::cout << "Building RMAT graph Edge factor: " << edge_factor_ << std::endl;
>>>>>>> origin/develop_keita
      } else {
        for (const auto itr : fname_edge_list_)
          std::cout << mpi_rank << " : " << "Load edge list from " << itr << std::endl;
      }

#if DEBUG_MODE
      {
        if (mpi_size > 1) {
          assert(false);
        }
        std::stringstream fname;
        fname << fname_segmentfile_ << ".debug_edges_raw";
        ofs_edges.open(fname.str());
      }
#endif
    }
    havoqgt::havoqgt_env()->world_comm().barrier();


    /// --- create a segument file --- ///
    std::stringstream fname;
    fname << fname_segmentfile_ << "_" << mpi_rank;
    boost::interprocess::file_mapping::remove(fname.str().c_str());
    if (mpi_rank == 0) {
      std::cout << "\n<<Construct segment>>" << std::endl;
      std::cout << "Create and map a segument file" << std::endl;
    }
    uint64_t segmentfile_init_size = std::pow(2, segmentfile_init_size_log2_) / mpi_size;
    mapped_file_type mapped_file = mapped_file_type(boost::interprocess::create_only, fname.str().c_str(), segmentfile_init_size);

#if 0
    boost::interprocess::mapped_region::advice_types advise = boost::interprocess::mapped_region::advice_types::advice_random;
    assert(mapped_file.advise(advise));
    std::cout << "Call adise_randam" << std::endl;
#endif
    if (mpi_rank == 0) {
      std::cout << "Call posix_fallocate\n";
    }
    fallocate(fname.str().c_str(), segmentfile_init_size, mapped_file);
    havoqgt::havoqgt_env()->world_comm().barrier();


    /// --- create a segument --- ///
    segment_manager_type* segment_manager = mapped_file.get_segment_manager();
    if (mpi_rank == 0) print_usages(segment_manager);
    havoqgt::havoqgt_env()->world_comm().barrier();


    std::cout << "Allocate graphstore baseline" << std::endl;
    boost::interprocess::allocator<void, segment_manager_type> allocator(segment_manager);
    graphstore_type graph_store(allocator);
    havoqgt::havoqgt_env()->world_comm().barrier();

    if (mpi_rank == 0) {
      std::cout << "\n<<Update edges>>" << std::endl;
    }
    if (fname_edge_list_.empty()) { /// generate edges
      uint64_t num_edges_per_rank = num_vertices * edge_factor_ / mpi_size;
      havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
<<<<<<< HEAD
                                        vert_scale, num_edges_per_rank,
                                        0.57, 0.19, 0.19, 0.05, true, false);
      apply_edges_update_requests(adjacency_matrix_map_vec,
                                  rmat,
=======
        vertex_scale_, num_edges_per_rank,
        0.57, 0.19, 0.19, 0.05, true, false);
      apply_edges_update_requests(mapped_file,
>>>>>>> origin/develop_keita
                                  segment_manager,
                                  graph_store,
                                  allocator,
                                  rmat,
                                  static_cast<uint64_t>(std::pow(10, chunk_size_log10_)),
                                  edges_delete_ratio_);
    } else {
      const double time_start = MPI_Wtime();
      havoqgt::parallel_edge_list_reader edgelist(fname_edge_list_);
      havoqgt::havoqgt_env()->world_comm().barrier();
      if (mpi_rank == 0) {
        std::cout << "TIME: Initializing a edge list reader (sec.) =\t" << MPI_Wtime() - time_start << std::endl;
      }
      apply_edges_update_requests(mapped_file,
                                  segment_manager,
                                  graph_store,
                                  allocator,
                                  edgelist,
                                  static_cast<uint64_t>(std::pow(10, chunk_size_log10_)),
                                  edges_delete_ratio_);
    }


#if DEBUG_MODE
    {
      std::cout << "dumping all elements for debug" << std::endl;
      std::stringstream ofname;
      ofname << fname_output << ".debug_edges_graph";
      std::ofstream  of(ofname.str());
      std::cout << "done" << std::endl;
    }
#endif

    if (is_delete_segmentfile_on_exit_) {
      boost::interprocess::file_mapping::remove(fname.str().c_str());
    }

    havoqgt::havoqgt_env()->world_comm().barrier();

  }
}
