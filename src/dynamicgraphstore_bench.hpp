/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef DYNAMICGRAPHSTORE_BENCH_HPP
#define DYNAMICGRAPHSTORE_BENCH_HPP

#include <iostream>
#include <utility>
#include <algorithm>
#include <functional>
#include <fcntl.h>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/rhh/rhh_common.hpp>

#include <havoqgt/graphstore/graphstore_rhhda.hpp>


#define SORT_BY_CHUNK 0

using vertex_id_type        = uint64_t;
using vertex_meta_data_type = unsigned char;
using edge_weight_type      = unsigned char;

using mapped_file_type     = boost::interprocess::managed_mapped_file;
using segment_manager_type = graphstore::rhh::segment_manager_t;


struct EdgeUpdateRequest
{
  EdgeUpdateRequest(std::pair<vertex_id_type, vertex_id_type> _edge, bool _is_delete)
  {
    edge = _edge;
    is_delete = _is_delete;
  }

  std::pair<vertex_id_type, vertex_id_type> edge;
  bool is_delete;
};

bool edgerequest_asc( const EdgeUpdateRequest& left, const EdgeUpdateRequest& right ) {
  return left.edge.first < right.edge.first;
}
using request_vector_type = std::vector<EdgeUpdateRequest>;


void print_usages(segment_manager_type *const segment_manager)
{
  const size_t usages = segment_manager->get_size() - segment_manager->get_free_memory();
  std::cout << "Usage: segment size (GiB) =\t"<< static_cast<double>(usages) / (1ULL << 30) << "\n";
  std::cout << "----------------------------" << std::endl;
  size_t mem_unit, totalram, freeram, usedram, bufferram, totalswap, freeswap;
  if (graphstore::utility::get_system_memory_usages(&mem_unit, &totalram, &freeram, &usedram, &bufferram, &totalswap, &freeswap)) {
    std::cout << "Usage: mem_unit:\t" << mem_unit << "\n";
    std::cout << "Usage: totalram(GiB):\t" << static_cast<double>(totalram) / (1<<30ULL) << "\n";
    std::cout << "Usage: freeram(GiB):\t" << static_cast<double>(freeram) / (1<<30ULL) << "\n";
    std::cout << "Usage: usedram(GiB):\t" << static_cast<double>(usedram) / (1<<30ULL) << "\n";
    std::cout << "Usage: bufferram(GiB):\t" << static_cast<double>(bufferram) / (1<<30ULL) << "\n";
    std::cout << "Usage: totalswap(GiB):\t" << static_cast<double>(totalswap) / (1<<30ULL) << "\n";
    std::cout << "Usage: freeswap(GiB):\t" << static_cast<double>(freeswap) / (1<<30ULL) << "\n";
    std::cout << "----------------------------" << std::endl;
  }
  size_t size, resident, share, text, lib, data, dt;
  if (graphstore::utility::get_my_memory_usages(&size, &resident, &share, &text, &lib, &data, &dt)) {
    std::cout << "Usage: VmSize(GiB):\t" << static_cast<double>(size) / (1<<30ULL) << "\n";
    std::cout << "Usage: VmRSS(GiB):\t" << static_cast<double>(resident) / (1<<30ULL) << "\n";
    std::cout << "Usage: sharedpages(GiB):\t" << static_cast<double>(share) / (1<<30ULL) << "\n";
    std::cout << "Usage: text(GiB):\t" << static_cast<double>(text) / (1<<30ULL) << "\n";
    std::cout << "Usage: library(GiB):\t" << static_cast<double>(lib) / (1<<30ULL) << "\n";
    std::cout << "Usage: data+stack(GiB):\t" << static_cast<double>(data) / (1<<30ULL) << "\n";
    std::cout << "Usage: dirtypages(GiB):\t" << static_cast<double>(dt) / (1<<30ULL) << "\n";
    std::cout << "----------------------------" << std::endl;
  }
}


void fallocate(const char* const fname, size_t size, mapped_file_type& asdf)
{
#if _XOPEN_SOURCE >= 600 || _POSIX_C_SOURCE >= 200112L
    int fd  = open(fname, O_RDWR);
    assert(fd != -1);
    int ret = posix_fallocate(fd, 0, size);
    assert(ret == 0);
    close(fd);
    asdf.flush();
#endif
}


double sort_requests(request_vector_type& requests)
{
  const double time_start1 = MPI_Wtime();
  std::sort(requests.begin(), requests.end(), edgerequest_asc);
  return (MPI_Wtime() - time_start1);
}

template <typename Edges_itr>
void generate_insertion_requests(Edges_itr& edges_itr, Edges_itr& edges_itr_last, const size_t chunk_size, request_vector_type& requests, const size_t delete_ratio)
{
  const int mpi_rank = havoqgt::havoqgt_env()->world_comm().rank();

  assert(requests.size() == 0);
  requests.reserve(chunk_size);

  const double time_start = MPI_Wtime();
  for (size_t cnt = 0; edges_itr != edges_itr_last && cnt < chunk_size; ++edges_itr, ++cnt) {
    const bool is_delete = (rand() % 100 < delete_ratio);
    EdgeUpdateRequest request(*edges_itr, is_delete);
    requests.push_back(request);
    //    std::cerr << edges_itr->first << " " << edges_itr->second << "\n";
  }
  havoqgt::havoqgt_env()->world_comm().barrier();
  if (mpi_rank == 0) std::cout << "TIME: Generate edges into DRAM (sec.) =\t" << MPI_Wtime() - time_start << std::endl;
  std::cout << mpi_rank << ": Status: # generated insetion requests =\t"<< requests.size() << std::endl;
#if SORT_BY_CHUNK
  const double elapsed_time = sort_requests(requests);
  havoqgt::havoqgt_env()->world_comm().barrier();
  if (mpi_rank == 0) std::cout << "TIME: Sorting chunk (sec.) =\t" << elapsed_time << std::endl;
 #endif
}


#endif // DYNAMICGRAPHSTORE_BENCH_HPP

