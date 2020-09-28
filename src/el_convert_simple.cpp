// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <deque>
#include <experimental/filesystem>
#include <functional>
#include <havoqgt/detail/hash.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <string>
#include <ygm/mailbox_p2p_nrroute.hpp>
#include <ygm/mpi.hpp>

using namespace havoqgt;
using namespace std::experimental::filesystem;
using havoqgt::detail::hash32;

typedef double edge_data_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr << "Usage: -o <string> -d <int> [file ...]\n"
              << " -o <string>   - output graph base filename (required)\n"
              << "[file ...] - list of edge list files to ingest\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& output_filename,
                    std::vector<std::string>& input_filenames) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_output_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "o:h ")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 'o':
        found_output_filename = true;
        output_filename       = optarg;
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_output_filename) {
    usage();
    exit(-1);
  }

  for (int index = optind; index < argc; index++) {
    /// std::cout << "Input file = " << argv[index] << std::endl;
    input_filenames.push_back(argv[index]);
  }
}

int main(int argc, char** argv) {
  init(&argc, &argv);

  std::string              output_filename;
  std::vector<std::string> input_filenames;

  parse_cmd_line(argc, argv, output_filename, input_filenames);

  parallel_edge_list_reader<edge_data_type> pelr(input_filenames, false);
  bool has_edge_data = pelr.has_edge_data();

  if (comm_world().rank() == 0) {
    std::cout << "Ingesting graph from " << input_filenames.size()
              << " files.  has_edge_data = " << has_edge_data << std::endl;
  }

  using edge_type = std::pair<uint64_t, uint64_t>;
  std::deque<edge_type> partitioned_edges;
  uint64_t              local_parse_count(0);
  uint64_t              local_recv_count(0);
  uint64_t              local_set_count(0);
  auto recvr = [&partitioned_edges, &local_recv_count](auto* mail, bool bcast,
                                                       const edge_type& edge) {
    partitioned_edges.push_back(edge);
    ++local_recv_count;
  };

  {
    ygm::mailbox_p2p_nrroute<edge_type, decltype(recvr)> mailbox(recvr,
                                                                 1024 * 1024);

    for (auto itr = pelr.begin(); itr != pelr.end(); ++itr) {
      ++local_parse_count;
      auto edge_tuple = *itr;
      //
      // skip self loops
      if (std::get<0>(edge_tuple) == std::get<1>(edge_tuple)) continue;

      //
      // orient edge from low to high vertex id
      uint64_t v_min =
          std::min(std::get<0>(edge_tuple), std::get<1>(edge_tuple));
      uint64_t v_max =
          std::max(std::get<0>(edge_tuple), std::get<1>(edge_tuple));
      edge_type edge = edge_type(v_min, v_max);

      //
      // compute MPI owner by hash
      size_t owner =
          (hash32(edge.first) ^ hash32(edge.second)) % comm_world().size();

      //
      // send to owner
      mailbox.send(owner, edge);
    }
  }
  std::sort(partitioned_edges.begin(), partitioned_edges.end());
  partitioned_edges.erase(
      std::unique(partitioned_edges.begin(), partitioned_edges.end()),
      partitioned_edges.end());

  // std::cout << whoami()
  //           << " partitioned_edges.size() = " << partitioned_edges.size()
  //           << std::endl;
  // comm_world().barrier();
  // for (auto edge : partitioned_edges) {
  //   std::cout << whoami() << " (" << edge.first << "," << edge.second << ")"
  //             << std::endl;
  // }

  local_set_count           = partitioned_edges.size();
  uint64_t global_set_count = comm_world().all_reduce(local_set_count, MPI_SUM);
  uint64_t global_parse_count =
      comm_world().all_reduce(local_parse_count, MPI_SUM);
  uint64_t global_recv_count =
      comm_world().all_reduce(local_recv_count, MPI_SUM);

  if (comm_world().rank() == 0) {
    std::cout << "global_set_count = " << global_set_count << std::endl;
    std::cout << "global_parse_count = " << global_parse_count << std::endl;
    std::cout << "global_recv_count = " << global_recv_count << std::endl;
  }

  return 0;
}