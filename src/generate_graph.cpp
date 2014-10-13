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

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
#include <havoqgt/upper_triangle_edge_generator.hpp>
#include <havoqgt/gen_preferential_attachment_edge_list.hpp>
#include <havoqgt/environment.hpp>
#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/distributed_db.hpp>
#include <iostream>
#include <assert.h>
#include <deque>
#include <utility>
#include <algorithm>
#include <functional>
#include <fstream>      // std::ifstream


// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

namespace hmpi = havoqgt::mpi;

using namespace havoqgt::mpi;

int main(int argc, char** argv) {

  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;

  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;

  int mpi_rank(0), mpi_size(0);

  CHK_MPI(MPI_Init(&argc, &argv));
  {

  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
  CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );
  havoqgt::get_environment();



  if (mpi_rank == 0) {

    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
    havoqgt::get_environment().print();

    //print_system_info(false);

  }
  MPI_Barrier(MPI_COMM_WORLD);


  uint64_t num_vertices = 1;
  uint64_t vert_scale;
  double   pa_beta;
  uint64_t hub_threshold;
  uint32_t load_from_disk;

  std::string type;
  std::string fname_output;
  if (argc < 7) {
    std::cerr << "usage: <RMAT/PA> <Scale> <PA-beta> <hub_threshold> <file name>"
              << " <load_from_disk> (argc:" << argc <<" )." << std::endl;
    exit(-1);
  } else {
    int pos = 1;
    type = argv[pos++];
    vert_scale    = boost::lexical_cast<uint64_t>(argv[pos++]);
    pa_beta       = boost::lexical_cast<double>(argv[pos++]);
    hub_threshold = boost::lexical_cast<uint64_t>(argv[pos++]);
    fname_output = argv[pos++];
    load_from_disk = boost::lexical_cast<uint32_t>(argv[pos++]);
  }
  num_vertices <<= vert_scale;
  if (mpi_rank == 0) {
    std::cout << "Building graph type: " << type << std::endl
      << "Building graph Scale: " << vert_scale << std::endl
      << "Hub threshold = " << hub_threshold << std::endl
      << "PA-beta = " << pa_beta << std::endl
      << "File name = " << fname_output << std::endl
      << "Load from disk = " << load_from_disk << std::endl
      << "Debuging = " << IS_DEBUGING << std::endl;
  }

/*
  std::stringstream fname;
  fname << fname_output << "_" << mpi_rank;

  bool file_exists;
  {
    std::ifstream fin(fname.str().c_str());
    file_exists = fin.good();
    fin.close();
  }

  uint64_t flash_capacity = std::pow(2,34) + std::pow(2,33) +  std::pow(2,32);
  assert (flash_capacity <= (751619276800.0/24.0));

  if (load_from_disk == 0) {
    if (file_exists) {
      std::cerr << "Graph file already exists." << std::endl;
      exit(-1);
    } else {
      //std::cout << "Creating new graph file:" << fname.str().c_str() << "." << std::endl;
    }
  }

  if (load_from_disk == 1) {
    if(file_exists) {
      //std::cout << "Loading from disk:" << fname.str().c_str() << "." << std::endl;
    } else {
      std::cerr << "Load From file specified, but file not found." << std::endl;
      exit(-1);
    }
  }
*/
  havoqgt::distributed_db ddb(havoqgt::db_create(), MPI_COMM_WORLD, fname_output.c_str());

  // boost::interprocess::mapped_region::advice_types advice;
  // advice = boost::interprocess::mapped_region::advice_types::advice_random;
  // bool assert_res = asdf.advise(advice);
  // assert(assert_res);



  segment_manager_t* segment_manager = ddb.get_segment_manager();
  bip::allocator<void, segment_manager_t> alloc_inst(segment_manager);

  graph_type *graph;

  if(type == "UPTRI") {
    uint64_t num_edges = num_vertices * 16;
    havoqgt::upper_triangle_edge_generator uptri(num_edges, mpi_rank, mpi_size,
         false);

    if (load_from_disk == 1) {
      if (mpi_rank == 0) {
        std::cout << "Loading graph from file." << std::endl;
      }

      graph = segment_manager->find<graph_type>("graph_obj").first;
      graph->complete_construction(alloc_inst, MPI_COMM_WORLD, uptri);
    } else {
      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }

      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (alloc_inst, MPI_COMM_WORLD, uptri, uptri.max_vertex_id(), hub_threshold);
    }


  } else if(type == "RMAT") {
    uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
    havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                      vert_scale, num_edges_per_rank,
                                      1.57, 0.19, 0.19, 0.05, true, true);


    if (load_from_disk == 1) {
      if (mpi_rank == 0) {
        std::cout << "Loading graph from file." << std::endl;
      }
      graph = segment_manager->find<graph_type>("graph_obj").first;


      // graph_type *graph_new;
      // graph_new = segment_manager->construct<graph_type>
      //   ("graph_obj_new")
      //   (alloc_inst, MPI_COMM_WORLD, rmat, rmat.max_vertex_id(), hub_threshold,
      //     graph_type::MetaDataGenerated);

      // assert(graph_new->compare(graph));

      graph->complete_construction(alloc_inst, MPI_COMM_WORLD, rmat);

    } else {
      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }
      graph = segment_manager->construct<graph_type>
        ("graph_obj")
        (alloc_inst, MPI_COMM_WORLD, rmat, rmat.max_vertex_id(), hub_threshold);
          // , graph_type::MetaDataGenerated);
        // graph->complete_construction(MPI_COMM_WORLD, rmat);
    }



  } else if(type == "PA") {
    std::vector< std::pair<uint64_t, uint64_t> > input_edges;

    gen_preferential_attachment_edge_list(input_edges, uint64_t(5489),
      vert_scale, vert_scale+4, pa_beta, 0.0, MPI_COMM_WORLD);

    if (load_from_disk == 1) {
      if (mpi_rank == 0) {
        std::cout << "Loading graph from file." << std::endl;
      }
      graph = segment_manager->find<graph_type>("graph_obj").first;
      graph->complete_construction(alloc_inst, MPI_COMM_WORLD, input_edges);
    } else {
      if (mpi_rank == 0) {
        std::cout << "Generating new graph." << std::endl;
      }
      graph = segment_manager->construct<graph_type>
        ("graph_obj")
        (alloc_inst, MPI_COMM_WORLD, input_edges, uint64_t(5489),
          hub_threshold);
    }

    {
      std::vector< std::pair<uint64_t, uint64_t> > empty(0);
      input_edges.swap(empty);
    }
  } else {
    std::cerr << "Unknown graph type: " << type << std::endl;  exit(-1);
  }



  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Graph Ready, Running Tests. (free/capacity) " << std::endl;
  }



  for (int i = 0; i < mpi_size; i++) {
    if (i == mpi_rank) {
      double percent = double(segment_manager->get_free_memory()) /
        double(segment_manager->get_size());
      std::cout << "[" << mpi_rank << "] " << segment_manager->get_free_memory()
      << "/" << segment_manager->get_size() << " = " << percent << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  graph->print_graph_statistics();


  //
  // Calculate max degree
  uint64_t max_degree(0);
  for (auto citr = graph->controller_begin(); citr != graph->controller_end(); ++citr) {
    max_degree = std::max(max_degree, graph->degree(*citr));
  }

  uint64_t global_max_degree = havoqgt::mpi::mpi_all_reduce(max_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);

  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));

  if (mpi_rank == 0) {
    std::cout << "Max Degree = " << global_max_degree << std::endl;
  }

//  asdf.flush();
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));

  if (mpi_rank == 0) {
    //print_system_info(false);
    //print_dmesg();
  }


  } //END Main MPI
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
  if (mpi_rank == 0)
    std::cout << "Pre Finalize." << std::endl;
  CHK_MPI(MPI_Finalize());
  if (mpi_rank == 0)
    std::cout << "After Bracket." << std::endl;

  std::cout << "FIN: " << mpi_rank << std::endl;
  return 0;
}
