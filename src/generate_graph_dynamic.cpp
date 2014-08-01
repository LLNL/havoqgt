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
#include <iostream>
#include <assert.h>
#include <deque>
#include <utility>
#include <algorithm>
#include <functional>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#ifdef PROFILE_DETAIL
 #warning PROFILE_DETAIL is enabled.
#endif


// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;
typedef bip::managed_mapped_file mapped_t;
typedef mapped_t::segment_manager segment_manager_t;
typedef hmpi::construct_dynamicgraph<segment_manager_t> graph_type;


template <typename Edges>
void add_edges_loop (graph_type *graph,
  mapped_t& asdf,
  bip::allocator<void, segment_manager_t>& alloc_inst,
  Edges& edges, uint64_t chunk_size) 
{

  const uint64_t num_edges = edges.size();
  chunk_size = std::min(chunk_size, num_edges);
  const uint64_t num_loop  = num_edges / chunk_size;

  auto edges_itr = edges.begin();

  for (uint64_t i = 0; i < num_loop; i++ ) {
    std::cout << "\n[" << i << " / " << num_loop << "]" << std::endl;
    boost::container::vector<std::pair<uint64_t, uint64_t>> onmemory_edges;
    const double time_start = MPI_Wtime();
    for (uint64_t j = 0; j < chunk_size && edges_itr != edges.end(); j++, edges_itr++) {
      onmemory_edges.push_back(*edges_itr);
      //std::cout << onmemory_edges[j].first << " " << onmemory_edges[j].second << std::endl;
    }
    const double time_end = MPI_Wtime();
    std::cout << "TIME: Generation edges into DRAM (sec.) =\t" << time_end - time_start << std::endl;
    graph->add_edges_adjacency_matrix(asdf, alloc_inst, onmemory_edges);

  }
  std::cout << "<< Results >>" << std::endl;
  graph->print_profile();

}



int main(int argc, char** argv) {


  CHK_MPI(MPI_Init(&argc, &argv));
  {
    int mpi_rank(0), mpi_size(0);
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
    }
    MPI_Barrier(MPI_COMM_WORLD);

  //typedef std::pair<uint64_t, uint64_t> edge_type;
  //std::deque<edge_type> tmp_edges;

  //havoqgt::mpi::edge_list<uint64_t> oned(MPI_COMM_WORLD);
  //havoqgt::mpi::edge_list<uint64_t> transpose_hubs(MPI_COMM_WORLD);
  //typedef extended_memory_arena arena_type;


    uint64_t num_vertices = 1;
    uint64_t vert_scale;
    double   pa_beta;
    uint64_t hub_threshold;
    uint32_t load_from_disk;
    uint32_t delete_file;
    uint64_t chunk_size_exp;
    std::string type;
    std::string fname_output;
    std::string fname_compare = "";
    std::string data_structure_type;

    if (argc < 10) {
      std::cerr << "usage: <RMAT/PA> <Scale> <PA-beta> <hub_threshold> <file name>"
      << " <load_from_disk> <delete file on exit>"
      << " <chunk_size_exp> <VC_VC/MP_VC>"
      << " <file to compare to>"
      << " (argc:" << argc << " )." << std::endl;
      exit(-1);
    } else {
      int pos = 1;
      type = argv[pos++];
      vert_scale      = boost::lexical_cast<uint64_t>(argv[pos++]);
      pa_beta         = boost::lexical_cast<double>(argv[pos++]);
      hub_threshold   = boost::lexical_cast<uint64_t>(argv[pos++]);
      fname_output    = argv[pos++];
      delete_file     = boost::lexical_cast<uint32_t>(argv[pos++]);
      load_from_disk  = boost::lexical_cast<uint32_t>(argv[pos++]);
      chunk_size_exp  = boost::lexical_cast<uint64_t>(argv[pos++]);
      data_structure_type = argv[pos++]; 
      if (pos < argc) {
        fname_compare = argv[pos++];
      }
    }
    num_vertices <<= vert_scale;

    if (mpi_rank == 0) {
      std::cout << "Building graph type: " << type << std::endl;
      std::cout << "Building graph Scale: " << vert_scale << std::endl;
      std::cout << "Hub threshold = " << hub_threshold << std::endl;
      std::cout << "PA-beta = " << pa_beta << std::endl;
      std::cout << "File name = " << fname_output << std::endl;
      std::cout << "Load from disk = " << load_from_disk << std::endl;
      std::cout << "Delete on Exit = " << delete_file << std::endl;
      std::cout << "Chunk size exp = " << chunk_size_exp << std::endl;
      std::cout << "Data structure type: " << data_structure_type << std::endl;
      if (fname_compare != "") {
        std::cout << "Comparing graph to " << fname_compare << std::endl;
      }
    }

    std::stringstream fname;
    //fname << "/l/ssd/"<< fname_output << "_" << mpi_rank;
    fname << fname_output << "_" << mpi_rank;

    if (load_from_disk  == 0) {
      remove(fname.str().c_str());
    }

  // TODO change to (2^34 + 2^33 + 2^32 )
    uint64_t graph_capacity = std::pow(2,34) + std::pow(2,33) +  std::pow(2,32);
    assert (graph_capacity <= (751619276800.0/24.0));
    mapped_t  asdf(bip::open_or_create, fname.str().c_str(),
      graph_capacity);
    segment_manager_t* segment_manager = asdf.get_segment_manager();
    bip::allocator<void,segment_manager_t> alloc_inst(segment_manager);

    graph_type *graph;

    if (data_structure_type == "VC_VC") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (alloc_inst, graph_type::kUseVecVecMatrix);
    } else if (data_structure_type == "MP_VC") {
      graph = segment_manager->construct<graph_type>
      ("graph_obj")
      (alloc_inst, graph_type::kUseMapVecMatrix);        
    } else {
      std::cerr << "Unknown data structure type: " << data_structure_type << std::endl;
      exit(-1);
    }

    if (load_from_disk  == 0) {

      if(type == "UPTRI") {
        uint64_t num_edges = num_vertices * 16;
        havoqgt::upper_triangle_edge_generator uptri(num_edges, mpi_rank, mpi_size,
         false);

        add_edges_loop(graph, asdf, alloc_inst, uptri, static_cast<uint64_t>(std::pow(2, chunk_size_exp)));

      } else if(type == "RMAT") {
        uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;

        havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
          vert_scale, num_edges_per_rank,
          0.57, 0.19, 0.19, 0.05, true, true);

        add_edges_loop(graph, asdf, alloc_inst, rmat, static_cast<uint64_t>(std::pow(2, chunk_size_exp)));

        // uint64_t num_edges = rmat.size();
        // uint64_t chunk_size = std::min((uint64_t)std::pow(2, chunk_size_exp), (uint64_t)num_edges);
        // uint64_t num_loop  = num_edges / chunk_size;

        // auto edges_itr = rmat.begin();

        // for (uint64_t i = 0; i < num_loop; i++ ) {
        //   boost::container::vector<std::pair<uint64_t, uint64_t>> onmemory_edges;
        //   double time_start = MPI_Wtime();
        //   for (uint64_t j = 0; j < chunk_size && edges_itr != rmat.end(); j++, edges_itr++) {
        //     onmemory_edges.push_back(*edges_itr);
        //     //std::cout << onmemory_edges[j].first << " " << onmemory_edges[j].second << std::endl;
        //   }
        //   double time_end = MPI_Wtime();
        //   std::cout << "TIME: Generation edges into DRAM (sec.) = " << time_end - time_start << std::endl;
        //   graph->add_edges_adjacency_matrix_map_vec_tor(asdf, alloc_inst, onmemory_edges);
        // }


      } else if(type == "PA") {
        std::vector< std::pair<uint64_t, uint64_t> > input_edges;
        gen_preferential_attachment_edge_list(input_edges, uint64_t(5489), vert_scale, vert_scale+4, pa_beta, 0.0, MPI_COMM_WORLD);

        add_edges_loop(graph, asdf, alloc_inst, input_edges, static_cast<uint64_t>(std::pow(2, chunk_size_exp)));
        {
          std::vector< std::pair<uint64_t, uint64_t> > empty(0);
          input_edges.swap(empty);
        }
      } else {
        std::cerr << "Unknown graph type: " << type << std::endl;  exit(-1);
      }
    } else {
      if (mpi_rank == 0) {
        std::cout << "Loading Graph from file." << std::endl;
      }
      graph = segment_manager->find<graph_type>("graph_obj").first;
    }

    MPI_Barrier(MPI_COMM_WORLD);
  // if (mpi_rank == 0) {
  //   std::cout << "Graph Ready, Running Tests. (free/capacity) " << std::endl;
  // }

    for (int i = 0; i < mpi_size; i++) {
      if (i == mpi_rank) {
        size_t usages = segment_manager->get_size() - segment_manager->get_free_memory();
        double percent = double(segment_manager->get_free_memory()) / double(segment_manager->get_size());
        std::cout << "[" << mpi_rank << "] " 
        << usages << " "
        << segment_manager->get_free_memory()
        << "/" << segment_manager->get_size() 
        << " = " << percent << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  //
  // Calculate max degree
  // uint64_t max_degree(0);
  // for (auto citr = graph->controller_begin(); citr != graph->controller_end(); ++citr) {
  //   max_degree = std::max(max_degree, graph->degree(*citr));
  // }

  // uint64_t global_max_degree = havoqgt::mpi::mpi_all_reduce(max_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
  // if (mpi_rank == 0) {
  //   std::cout << "Max Degree = " << global_max_degree << std::endl;
  // }

    // if (fname_compare != "") {
    //   std::stringstream fname2;
    //   fname2 << "/l/ssd/"<< fname_compare << "_" << mpi_rank;

    //   mapped_t  asdf2(bip::open_or_create, fname2.str().c_str(), 1024ULL*1024*1024*16);
    //   segment_manager_t* segment_manager2 = asdf2.get_segment_manager();

    //   graph_type *graph_other =
    //   segment_manager2->find<graph_type>("graph_obj").first;

    //   std::cout << "[" << mpi_rank << "]Comparing member variables of the two graphs";
    //   // if (graph != nullptr && graph_other != nullptr && *graph == *graph_other) {
    //   //   std::cout << "...they are equivelent" << std::endl;
    //   // } else {
    //   //   std::cout << "...they are different." << std::endl;
    //   // }
    // }

    if (delete_file) {
      bip::file_mapping::remove(fname.str().c_str());
    }

  } //END Main MPI
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
  CHK_MPI(MPI_Finalize());
  return 0;
}
