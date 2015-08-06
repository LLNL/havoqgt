/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#include "dynamicgraphstore_bench.hpp"

#define VERBOSE 0

#define DEBUG_MODE 0
#if DEBUG_MODE
std::ofstream ofs_edges;
#endif

enum : size_t {
  midle_high_degree_threshold = 3
};

using vertex_meta_data_type = unsigned char;
using edge_weight_type = unsigned char;
using graphstore_type  = graphstore::graphstore_rhhda<vertex_id_type, vertex_meta_data_type, edge_weight_type, static_cast<size_t>(midle_high_degree_threshold)>;

template <typename Edges>
void apply_edges_update_requests(graphstore_type& graph_store, Edges& edges, segment_manager_type *const segment_manager, const uint64_t chunk_size, const size_t edges_delete_ratio)
{
  int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt::havoqgt_env()->world_comm().size();

  havoqgt::havoqgt_env()->world_comm().barrier();
  if (mpi_rank == 0) std::cout << "-- Disp status of before generation --" << std::endl;
  if (mpi_rank == 0) print_usages(segment_manager);
  havoqgt::havoqgt_env()->world_comm().barrier();

  uint64_t loop_cnt = 0;
  auto edges_itr = edges.begin();
  auto edges_itr_end = edges.end();
  bool global_is_finished = false;
  request_vector_type update_request_vec = request_vector_type();

  size_t count_inserted = 0;
  size_t count_delete = 0;

  while (!global_is_finished) {
    if (mpi_rank == 0) std::cout << "\n[" << loop_cnt << "] : chunk_size =\t" << chunk_size << std::endl;

    update_request_vec.clear();
    generate_insertion_requests(edges_itr, edges_itr_end, chunk_size, update_request_vec, edges_delete_ratio);
    havoqgt::havoqgt_env()->world_comm().barrier();

    const double time_start = MPI_Wtime();
    unsigned char dummy = 0;
    for (auto request : update_request_vec) {
      auto edge = request.edge;
      if (request.is_delete) {
#if DEBUG_MODE
        ofs_edges << edge.first << " " << edge.second << " 1" << "\n";
#endif
        count_delete += graph_store.erase_edge(edge.first, edge.second);
      } else {
#if DEBUG_MODE
        ofs_edges << edge.first << " " << edge.second << " 0" << "\n";
#endif
        count_inserted += graph_store.insert_edge(edge.first, edge.second, dummy);
      }
    }
    havoqgt::havoqgt_env()->world_comm().barrier();
    const double time_end = MPI_Wtime();
    if (mpi_rank == 0) std::cout << "TIME: Execution time (sec.) =\t" << time_end - time_start << std::endl;
    if (mpi_rank == 0) print_usages(segment_manager);
    havoqgt::havoqgt_env()->world_comm().barrier();
#if VERBOSE
    for (int i = 0; i < mpi_size; ++i) {
      if (i == mpi_rank) {
        std::cout << "Proc no. " << mpi_rank << std::endl;
        graph_store.print_status();
      }
      havoqgt::havoqgt_env()->world_comm().barrier();
    }
#endif

    ++loop_cnt;

    bool local_is_finished = (edges_itr == edges_itr_end);
    MPI_Allreduce(&local_is_finished, &global_is_finished, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
  }
  if (mpi_rank == 0) print_usages(segment_manager);
  havoqgt::havoqgt_env()->world_comm().barrier();

  std::cout << "inserted edges : " << count_inserted << std::endl;
  std::cout << "deleted edges : " << count_delete << std::endl;
}

int main(int argc, char** argv) {

  havoqgt::havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt::havoqgt_env()->world_comm().size();
    havoqgt::get_environment();

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      std::cout << "CMD line:";
      for (int i=0; i<argc; ++i) {
        std::cout << " " << argv[i];
      }
      std::cout << std::endl;
      havoqgt::get_environment().print();
    }
    havoqgt::havoqgt_env()->world_comm().size();

    uint64_t num_vertices = 1;
    uint64_t vert_scale;
    uint64_t edge_factor;
    uint32_t delete_file;
    uint64_t chunk_size_exp;
    uint64_t segmentfile_init_size;
    uint64_t edges_delete_ratio = 0;
    std::string fname_output;
    std::vector<std::string> fname_edge_list;

    if (argc < 8) {
      std::cerr << "usage: <Scale> <Edge factor> <segmentfile name>"
      << " <segmentfile_init_size (exp)> <delete file on exit>"
      << " <chunk_size_exp> <edges_delete_ratio>"
      << " <edgelist file>"
      << " (argc:" << argc << " )." << std::endl;
      exit(-1);
    } else {
      int pos = 1;
      vert_scale            = boost::lexical_cast<uint64_t>(argv[pos++]);
      edge_factor           = boost::lexical_cast<uint64_t>(argv[pos++]);
      fname_output          = argv[pos++];
      segmentfile_init_size = boost::lexical_cast<uint64_t>(argv[pos++]);
      delete_file           = boost::lexical_cast<uint32_t>(argv[pos++]);
      chunk_size_exp        = boost::lexical_cast<uint64_t>(argv[pos++]);
      edges_delete_ratio    = boost::lexical_cast<uint64_t>(argv[pos++]);
      if (pos < argc) {
        std::string fname(argv[pos++]);
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
         if (num_local_files * mpi_rank <= count && count < num_local_files * (mpi_rank + 1)) {
            fname_edge_list.push_back(line);
          }
          ++count;
        }
      }
    }
    num_vertices <<= vert_scale;

    if (mpi_rank == 0) {
      std::cout << "Segment file name = " << fname_output << std::endl;
      std::cout << "Initialize segment filse size = " << segmentfile_init_size << std::endl;
      std::cout << "Delete on Exit = " << delete_file << std::endl;
      std::cout << "Chunk size exp = " << chunk_size_exp << std::endl;
      std::cout << "Edges Delete Ratio = " << edges_delete_ratio << std::endl;
      std::cout << "Midle-high degree threshold = " << midle_high_degree_threshold << std::endl;

      if (fname_edge_list.empty()) {
        std::cout << "Building RMAT graph Scale: " << vert_scale << std::endl;
        std::cout << "Building RMAT graph Edge factor: " << edge_factor << std::endl;
      } else {
        for (auto itr = fname_edge_list.begin(), itr_end = fname_edge_list.end(); itr != itr_end; ++itr)
          std::cout << mpi_rank << " : " << "Load edge list from " << *itr << std::endl;
      }

#if DEBUG_MODE
      {
      if (mpi_size > 1) {
          assert(false);
        }
        std::stringstream fname;
        fname << fname_output << ".debug_edges_raw";
        ofs_edges.open(fname.str());
      }
#endif
    }
    havoqgt::havoqgt_env()->world_comm().barrier();


    /// --- create a segument file --- ///
    std::stringstream fname;
    fname << fname_output << "_" << mpi_rank;
    boost::interprocess::file_mapping::remove(fname.str().c_str());
    if (mpi_rank == 0) {
      std::cout << "\n<<Construct segment>>" << std::endl;
      std::cout << "Create and map a segument file" << std::endl;
    }
    uint64_t graph_capacity = std::pow(2, segmentfile_init_size) / mpi_size;
    mapped_file_type asdf = mapped_file_type(boost::interprocess::create_only, fname.str().c_str(), graph_capacity);

#if 0
    boost::interprocess::mapped_region::advice_types advise = boost::interprocess::mapped_region::advice_types::advice_random;
    assert(asdf.advise(advise));
    std::cout << "Call adise_randam" << std::endl;
#endif
    if (mpi_rank == 0) {
      std::cout << "Call posix_fallocate\n";
    }
    fallocate(fname.str().c_str(), graph_capacity, asdf);
    havoqgt::havoqgt_env()->world_comm().barrier();


    /// --- create a segument --- ///
    segment_manager_type* segment_manager = asdf.get_segment_manager();
    if (mpi_rank == 0) print_usages(segment_manager);
    havoqgt::havoqgt_env()->world_comm().barrier();


    std::cout << "Allocate graphstore_rhh_matrix" << std::endl;
    graphstore_type graph_store(segment_manager);
    havoqgt::havoqgt_env()->world_comm().barrier();


    if (mpi_rank == 0) {
      std::cout << "\n<<Update edges>>" << std::endl;
    }
    if (fname_edge_list.empty()) {
      uint64_t num_edges_per_rank = num_vertices * edge_factor / mpi_size;
      havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
        vert_scale, num_edges_per_rank,
        0.57, 0.19, 0.19, 0.05, true, false);
      apply_edges_update_requests(graph_store,
                                  rmat,
                                  segment_manager,
                                  static_cast<uint64_t>(std::pow(2, chunk_size_exp)),
                                  edges_delete_ratio);
    } else {
      const double time_start = MPI_Wtime();
      havoqgt::parallel_edge_list_reader edgelist(fname_edge_list);
      if (mpi_rank == 0) {
        std::cout << "TIME: Initializing a edge list reader (sec.) =\t" << MPI_Wtime() - time_start << std::endl;
      }
      havoqgt::havoqgt_env()->world_comm().barrier();

      apply_edges_update_requests(graph_store,
                                  edgelist,
                                  segment_manager,
                                  static_cast<uint64_t>(std::pow(2, chunk_size_exp)),
                                  edges_delete_ratio);
    }


#if DEBUG_MODE
    {
      std::cout << "dumping all elements for debug" << std::endl;
      std::stringstream ofname;
      ofname << fname_output << ".debug_edges_graph";
      std::ofstream  of(ofname.str());
      graph_store.fprint_all_elements(of);
      std::cout << "done" << std::endl;
    }
#endif

    havoqgt::havoqgt_env()->world_comm().barrier();

  }
}
