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

#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/new_triangle_count.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <havoqgt/distributed_db.hpp>
#include <assert.h>

#include <deque>
#include <string>
#include <utility>
#include <algorithm>
#include <functional>

#include <boost/interprocess/managed_heap_memory.hpp>

using namespace havoqgt;

int main(int argc, char** argv) {
  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
  typedef havoqgt::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  havoqgt::havoqgt_init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
    havoqgt::get_environment().print();
    //print_system_info(false);
  }
  MPI_Barrier(MPI_COMM_WORLD);


/*  std::string graph_input;
  std::string graph_backup

  if (argc < 2) {
    std::cerr << "usage: <graph input file name>"
      << " (argc:" << argc << " )." << std::endl;
    exit(-1);
  } else {
    int pos = 1;
    graph_input = argv[pos++];
  }
*/

  MPI_Barrier(MPI_COMM_WORLD);
  if(mpi_rank == 0) {
    std::cout << "Transfering " << argv[2] << " to " << argv[1] << std::endl;
    std::cout << "Outputing degree distributions to " << argv[3] << std::endl;
  }

  distributed_db::transfer(argv[2], argv[1]);
  havoqgt::distributed_db ddb(havoqgt::db_open(), argv[1]);

  graph_type *graph = ddb.get_segment_manager()->
    find<graph_type>("graph_obj").first;
  assert(graph != nullptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Loaded Ready." << std::endl;
  }
  //graph->print_graph_statistics();
  MPI_Barrier(MPI_COMM_WORLD);

//  for(int i=0; i<100; ++i) {
    uint64_t count = new_triangle_count(*graph, argv[3]);
    if(mpi_rank == 0) {
      std::cout << "Graph has " << count <<  " triangles." << std::endl;
    }
//  }

  }  // END Main MPI
  havoqgt::havoqgt_finalize();

  return 0;
}
