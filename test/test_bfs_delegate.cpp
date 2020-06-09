/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License. http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process disclosed, or represents that its use would not infringe
 * privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by the United States Government or Lawrence Livermore National Security, LLC.
 * The views and opinions of authors expressed herein do not necessarily state
 * or reflect those of the United States Government or Lawrence Livermore
 * National Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <bfs.hpp>
#include <havoqgt/mpi.hpp>
#include <util.hpp>

using namespace havoqgt::test;

/// \brief Validate BFS results on a star graph
void validate_bfs_result(const graph_type &          graph,
                         const std::size_t           max_vertex_id,
                         const bfs_level_data_type & level_data,
                         const bfs_parent_data_type &parent_data) {
  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  // Load BFS level table
  std::unordered_map<uint64_t, uint64_t> level_table;
  {
    std::ifstream ifs("./datasets/bfs_level_rmat_s17");
    if (!ifs) {
      std::cerr << "Cannot open" << std::endl;
      CHK_MPI(MPI_Abort(MPI_COMM_WORLD, -1));
    }

    uint64_t vid;
    uint64_t level;
    while (ifs >> vid >> level) {
      level_table[vid] = level;
    }
  }

  for (uint64_t i = 0; i <= max_vertex_id; ++i) {
    graph_type::vertex_locator vertex = graph.label_to_locator(i);
    if (vertex.owner() == mpi_rank) {
      size_t level = level_data[vertex];
      if (level_table[i] != level) {
        std::cout << "Unexpected BFS level" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    havoqgt::comm_world().barrier();
  }
}

int main(int argc, char **argv) {
  havoqgt::init(&argc, &argv);
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    havoqgt::cout_rank0() << "MPI initialized with " << mpi_size << " ranks."
                          << std::endl;

    havoqgt::cout_rank0() << "Data store path:\t"
                          << gen_test_dir_path(k_test_name) << std::endl;

    std::vector<std::string> f;
    f.emplace_back("./datasets/edge_list_rmat_s17-0");
    f.emplace_back("./datasets/edge_list_rmat_s17-1");
    f.emplace_back("./datasets/edge_list_rmat_s17-2");
    f.emplace_back("./datasets/edge_list_rmat_s17-3");
    f.emplace_back("./datasets/edge_list_rmat_s17-4");
    f.emplace_back("./datasets/edge_list_rmat_s17-5");
    f.emplace_back("./datasets/edge_list_rmat_s17-6");
    f.emplace_back("./datasets/edge_list_rmat_s17-7");
    ingest_edges(f);

    static constexpr uint64_t max_vertex_id = 131070;
    run_bfs(max_vertex_id, validate_bfs_result);

    havoqgt::comm_world().barrier();
    havoqgt::cout_rank0() << "Succeeded all BFS tests" << std::endl;
  }  // End of MPI

  return 0;
}