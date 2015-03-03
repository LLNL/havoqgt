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

#include <havoqgt/construct_dynamicgraph.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/upper_triangle_edge_generator.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <iostream>
#include <assert.h>
#include <deque>
#include <utility>
#include <algorithm>
#include <functional>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#ifdef __linux__
#include <sys/sysinfo.h>
#endif

#define SORT_BY_CHUNK 0

using namespace havoqgt;
namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

typedef bip::managed_mapped_file mapped_t;
typedef mapped_t::segment_manager segment_manager_t;
typedef hmpi::construct_dynamicgraph<segment_manager_t> graph_type;
typedef boost::container::vector<EdgeUpdateRequest> request_vector_type;

bool get_system_memory_usages( size_t* const mem_unit, size_t* const totalram, size_t* const freeram, size_t* const usedram,
                               size_t* const bufferram, size_t* const totalswap, size_t* const freeswap )
{
  bool is_succeed = false;
#ifdef __linux__
  struct sysinfo info;
  if (sysinfo(&info) == 0) {
    *mem_unit = static_cast<size_t>(info.mem_unit);
    *totalram = static_cast<size_t>(info.totalram);
    *freeram = static_cast<size_t>(info.freeram);
    *usedram = static_cast<size_t>(info.totalram - info.freeram);
    *bufferram = static_cast<size_t>(info.bufferram);
    *totalswap = static_cast<size_t>(info.totalswap);
    *freeswap = static_cast<size_t>(info.freeswap);
    is_succeed = true;
  }
#endif
  return is_succeed;
}


bool get_my_memory_usages(size_t* const size, size_t* const resident, size_t* const share,
                          size_t* const text, size_t* const lib, size_t* const data, size_t* const dt)
{
  bool is_succeed = false;
#ifdef __linux__
  const char* statm_path = "/proc/self/statm";
  FILE *f = fopen(statm_path, "r");
  if (f) {
    if (7 == fscanf(f,"%ld %ld %ld %ld %ld %ld %ld", size, resident, share, text, lib, data, dt)) {
      is_succeed = true;
    }
  }
  fclose(f);
#endif
  return is_succeed;
}

void print_usages(segment_manager_t *const segment_manager)
{
  const size_t usages = segment_manager->get_size() - segment_manager->get_free_memory();
  std::cout << "Usage: segment size =\t"<< usages << "\n";
  std::cout << "----------------------------" << std::endl;
  size_t mem_unit, totalram, freeram, usedram, bufferram, totalswap, freeswap;
  if (get_system_memory_usages(&mem_unit, &totalram, &freeram, &usedram, &bufferram, &totalswap, &freeswap)) {
    std::cout << "Usage: mem_unit:\t" << mem_unit << "\n";
    std::cout << "Usage: totalram(GiB):\t" << static_cast<double>(totalram) / (1<<30ULL) << "\n";
    std::cout << "Usage: freeram(GiB):\t" << static_cast<double>(freeram) / (1<<30ULL) << "\n";
    std::cout << "Usage: usedram(GiB):\t" << static_cast<double>(usedram) / (1<<30ULL) << "\n";
    std::cout << "Usage: bufferram(GiB):\t" << static_cast<double>(bufferram) / (1<<30ULL) << "\n";
    std::cout << "Usage: totalswap(GiB):\t" << static_cast<double>(totalswap) / (1<<30ULL) << "\n";
    std::cout << "Usage: freeswap(GiB):\t" << static_cast<double>(freeswap) / (1<<30ULL) << "\n";
  }
  std::cout << "----------------------------" << std::endl;
  size_t size, resident, share, text, lib, data, dt;
  if (get_my_memory_usages(&size, &resident, &share, &text, &lib, &data, &dt)) {
    std::cout << "Usage: VmSize(GiB):\t" << static_cast<double>(size) / (1<<30ULL) << "\n";
    std::cout << "Usage: VmRSS(GiB):\t" << static_cast<double>(resident) / (1<<30ULL) << "\n";
    std::cout << "Usage: sharedpages(GiB):\t" << static_cast<double>(share) / (1<<30ULL) << "\n";
    std::cout << "Usage: text(GiB):\t" << static_cast<double>(text) / (1<<30ULL) << "\n";
    std::cout << "Usage: library(GiB):\t" << static_cast<double>(lib) / (1<<30ULL) << "\n";
    std::cout << "Usage: data+stack(GiB):\t" << static_cast<double>(data) / (1<<30ULL) << "\n";
    std::cout << "Usage: dirtypages(GiB):\t" << static_cast<double>(dt) / (1<<30ULL) << "\n";
  }
  std::cout << "----------------------------" << std::endl;
}

double sort_requests(request_vector_type& requests)
{
  const double time_start1 = MPI_Wtime();
  std::sort(requests.begin(), requests.end(), edgerequest_asc);
  return (MPI_Wtime() - time_start1);
}

template <typename Edges_itr>
void generate_insertion_requests(Edges_itr& edges_itr, const uint64_t chunk_size, request_vector_type& requests)
{
  assert(requests.size() == 0);
  requests.reserve(chunk_size);
  const double time_start = MPI_Wtime();
  for (uint64_t i = 0; i < chunk_size; ++i, ++edges_itr) {
    EdgeUpdateRequest request(*edges_itr, false);
    requests.push_back(request);
  }
  std::cout << "TIME: Generate edges into DRAM (sec.) =\t" << MPI_Wtime() - time_start << std::endl;
  std::cout << "Status: # insetion requests =\t"<< requests.size() << std::endl;
#if SORT_BY_CHUNK
  const double elapsed_time = sort_requests(requests);
  std::cout << "TIME: Sorting chunk (sec.) =\t" << elapsed_time << std::endl;
 #endif
}

template <typename Edges>
void apply_edges_update_requests(graph_type *const graph, Edges& edges, segment_manager_t *const segment_manager, const uint64_t chunk_size)
{
  std::cout << "-- Disp status of before generation --" << std::endl;
  print_usages(segment_manager);

  const uint64_t num_edges = edges.size();
  const uint64_t num_requests = num_edges;
  const uint64_t num_loops = num_requests / chunk_size;
  auto edges_itr = edges.begin();
  for (uint64_t i = 0; i < num_loops; ++i) {
    std::cout << "\n[" << i+1 << " / " << num_loops << "] : chunk_size =\t" << chunk_size << std::endl;
    request_vector_type update_request_vec = request_vector_type();
    generate_insertion_requests(edges_itr, chunk_size, update_request_vec);
    auto requests_itr = update_request_vec.begin();
    graph->add_edges_adjacency_matrix(requests_itr, chunk_size);
    print_usages(segment_manager);
  }
  const uint64_t remains = num_edges - num_loops*chunk_size;
  if (remains > 0) {
    std::cout << "\n[ * / " << num_loops << "] : remains =\t" << remains << std::endl;
    request_vector_type update_request_vec = request_vector_type();
    generate_insertion_requests(edges_itr, remains, update_request_vec);
    auto requests_itr = update_request_vec.begin();
    graph->add_edges_adjacency_matrix(requests_itr, remains);
  }
  print_usages(segment_manager);
  graph->print_profile();
}


int main(int argc, char** argv) {

  int mpi_rank(0), mpi_size(0);

  havoqgt_init(&argc, &argv);
  {
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();
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
    MPI_Barrier(MPI_COMM_WORLD);

    uint64_t num_vertices = 1;
    uint64_t vert_scale;
    uint64_t edge_factor;
    uint32_t delete_file;
    uint64_t chunk_size_exp;
    uint64_t low_degree_threshold;
    uint64_t edges_delete_ratio = 0;
    std::string fname_output;
    std::string data_structure_type;
    std::vector<std::string> fname_edge_list;

    if (argc < 9) {
      std::cerr << "usage: <Scale> <Edge factor> <file name>"
      << " <delete file on exit>"
      << " <chunk_size_exp> <VC_VC/MP_VC/RB_HS/RB_MP/RB_MX> <low_degree_threshold> <edges_delete_ratio>"
      << " <edgelist file>"
      << " (argc:" << argc << " )." << std::endl;
      exit(-1);
    } else {
      int pos = 1;
      vert_scale      = boost::lexical_cast<uint64_t>(argv[pos++]);
      edge_factor     = boost::lexical_cast<uint64_t>(argv[pos++]);
      fname_output    = argv[pos++];
      delete_file     = boost::lexical_cast<uint32_t>(argv[pos++]);
      chunk_size_exp  = boost::lexical_cast<uint64_t>(argv[pos++]);
      data_structure_type = argv[pos++];
      low_degree_threshold = boost::lexical_cast<uint64_t>(argv[pos++]);
      edges_delete_ratio = boost::lexical_cast<uint64_t>(argv[pos++]);
      if (pos < argc) {
        std::string fname(argv[pos]);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while(std::getline(fin, line)) {
          fname_edge_list.push_back(line);
        }
      }
    }
    num_vertices <<= vert_scale;

    if (mpi_rank == 0) {
      std::cout << "File name = " << fname_output << std::endl;
      std::cout << "Delete on Exit = " << delete_file << std::endl;
      std::cout << "Chunk size exp = " << chunk_size_exp << std::endl;
      std::cout << "Data structure type: " << data_structure_type << std::endl;
      std::cout << "Low Degree Threshold = " << low_degree_threshold << std::endl;
      std::cout << "Edges Delete Ratio = " << edges_delete_ratio << std::endl;
      if (fname_edge_list.empty()) {
        std::cout << "Building RMAT graph Scale: " << vert_scale << std::endl;
        std::cout << "Building RMAT graph Edge factor: " << edge_factor << std::endl;
      } else {
        for (auto itr = fname_edge_list.begin(), itr_end = fname_edge_list.end(); itr != itr_end; ++itr)
        std::cout << "Load edge list from " << *itr << std::endl;
      }
    }

    /// --- create a segument file --- ///
    std::stringstream fname;
    fname << fname_output << "_" << mpi_rank;
    if (std::ifstream(fname.str().c_str()).good())
    {
      HAVOQGT_ERROR_MSG("File already exists.");
    }


    std::cout << "\n<<Construct segment>>" << std::endl;
    uint64_t graph_capacity = std::pow(2, 39);
    mapped_t asdf = mapped_t(bip::create_only, fname.str().c_str(), graph_capacity);

#if 0
    boost::interprocess::mapped_region::advice_types advise = boost::interprocess::mapped_region::advice_types::advice_random;
    assert(asdf.advise(advise));
    std::cout << "Calling adise_randam" << std::endl;
#endif
    int fd  = open(fname.str().c_str(), O_RDWR);
    assert(fd != -1);
    int ret = posix_fallocate(fd,0,graph_capacity);
    assert(ret == 0);
    close(fd);
    asdf.flush();

    /// --- create a segument --- ///
    segment_manager_t* segment_manager = asdf.get_segment_manager();
    bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);
    print_usages(segment_manager);

    graph_type *graph;
    if (data_structure_type == "VC_VC") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (asdf, alloc_inst, graph_type::kUseVecVecMatrix, low_degree_threshold);
    } else if (data_structure_type == "MP_VC") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (asdf, alloc_inst, graph_type::kUseMapVecMatrix, low_degree_threshold);
    } else if (data_structure_type == "RH_AR") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (asdf, alloc_inst, graph_type::kUseRHHAsArray, low_degree_threshold);
    } if (data_structure_type == "RH_MX") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (asdf, alloc_inst, graph_type::kUseRHHAsMatrix, low_degree_threshold);
    } else if (data_structure_type == "HY_DA") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (asdf, alloc_inst, graph_type::kUseHybridDegreeAwareModel, low_degree_threshold);
    } else {
      std::cerr << "Unknown data structure type: " << data_structure_type << std::endl;
      exit(-1);
    }


    std::cout << "\n<<Update edges>>" << std::endl;
    if (fname_edge_list.empty()) {
      uint64_t num_edges_per_rank = num_vertices * edge_factor / mpi_size;
      havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
        vert_scale, num_edges_per_rank,
        0.57, 0.19, 0.19, 0.05, true, false);
      apply_edges_update_requests(graph,
                                  rmat,
                                  segment_manager,
                                  static_cast<uint64_t>(std::pow(2, chunk_size_exp)));
    } else {
      double time_start = MPI_Wtime();
      havoqgt::parallel_edge_list_reader edgelist(fname_edge_list);
      std::cout << "TIME: Initializing a edge list reader (sec.) =\t" << MPI_Wtime() - time_start << std::endl;
      apply_edges_update_requests(graph,
                                  edgelist,
                                  segment_manager,
                                  static_cast<uint64_t>(std::pow(2, chunk_size_exp)));
    }
    havoqgt_env()->world_comm().barrier();

    std::cout << "\n<<Profile>>" << std::endl;
    print_usages(segment_manager);

    havoqgt_env()->world_comm().barrier();

    std::cout << "\nDelete the segment" << std::endl;
    /// this function call graph_type's destructor
    segment_manager->destroy<graph_type>("graph_obj");
    if (delete_file) {
      std::cout << "Delete the mapped file" << std::endl;
      bip::file_mapping::remove(fname.str().c_str());
    }
   havoqgt_env()->world_comm().barrier();

  } //END Main MPI
  havoqgt_finalize();

  return 0;
}
