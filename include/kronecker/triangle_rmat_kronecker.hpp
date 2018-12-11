#pragma once

#include <fstream>
#include <havoqgt/mpi.hpp>
#include <iostream>
#include <kronecker/triangle_kronecker_edge_generator.hpp>
#include <string>
#include <tuple>
#include <vector>

void read_graph_file(
    std::string                                            filename,
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& edge_list) {
  std::ifstream filestream(filename);

  if (filestream.is_open()) {
    std::string line;
    while (std::getline(filestream, line)) {
      std::istringstream iss(line);
      uint64_t           src, dest, tri;
      if (!(iss >> src >> dest >> tri)) {
        std::cerr << "Malformed line in input\n";
      } else {
        edge_list.push_back(std::make_tuple(src, dest, tri));
      }
    }
    filestream.close();
  } else {
    std::cerr << "Unable to open file " << filename << std::endl;
  }
}

triangle_kronecker_edge_generator<
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>, uint64_t>
rmat_kronecker(int scale) {
  std::string graph_prefix = "/p/lustre1/steil1/graphs/";

  int scale1, scale2;
  scale1 = scale / 2;
  scale2 = scale / 2 + scale % 2;

  std::string graph_filename1, graph_filename2;
  graph_filename1 = graph_prefix + "rmat_" + std::to_string(scale1) + "_a";
  graph_filename2 = graph_prefix + "rmat_" + std::to_string(scale2) + "_b";

  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> edge_list1, edge_list2;

  if (ygm::comm_world().rank() == 0) {
    std::cout << "Reading graphs" << std::endl;
    read_graph_file(graph_filename1, edge_list1);
    read_graph_file(graph_filename2, edge_list2);
  }

  size_t graph1_size, graph2_size;
  graph1_size = edge_list1.size();
  graph2_size = edge_list2.size();

  if (ygm::comm_world().rank() == 0) {
    std::cout << "Broadcasting edge lists" << std::endl;
  }
  havoqgt::mpi_bcast(graph1_size, 0, ygm::comm_world().mpi_comm());
  havoqgt::mpi_bcast(graph2_size, 0, ygm::comm_world().mpi_comm());
  edge_list1.resize(graph1_size);
  edge_list2.resize(graph2_size);
  MPI_Bcast(edge_list1.data(), sizeof(edge_list1[0]) * graph1_size, MPI_BYTE, 0,
            ygm::comm_world().mpi_comm());
  MPI_Bcast(edge_list2.data(), sizeof(edge_list2[0]) * graph2_size, MPI_BYTE, 0,
            ygm::comm_world().mpi_comm());

  if (ygm::comm_world().rank() == 0) {
    std::cout << "Constructing Kronecker product" << std::endl;
  }
  triangle_kronecker_edge_generator<
      std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
      std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>, uint64_t>
      kron(edge_list1, edge_list2, edge_list1.size(), edge_list2.size(), false);

  return kron;
}
