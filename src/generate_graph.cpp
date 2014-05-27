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

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/rmat_edge_generator.hpp>
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


// class heap_arena 
// {
// public:
//   heap_arena() { }
//   template <typename T>
//   struct allocator {
//     typedef typename std::allocator<T> type;
//   };

//   template<typename T>
//   std::allocator<T> make_allocator() {
//     return std::allocator<T>();
//   }
// private:
//   heap_arena(const heap_arena&);
// };

// class extended_memory_arena {
//   // typedef boost::interprocess::basic_managed_mapped_file 
//   //   <char
//   //   ,boost::interprocess::rbtree_best_fit<boost::interprocess::null_mutex_family, void*>
//   //   ,boost::interprocess::iset_index>  managed_mapped_file;
  
// public:
//   template <typename T>
//   struct allocator {
//     typedef typename boost::interprocess::managed_mapped_file::template allocator<T>::type type;
//   };

//   template <typename T>
//   typename allocator<T>::type make_allocator() {
//     return m_mapped_file->template get_allocator<T>();
//   }
//   extended_memory_arena(const char* fname) 
//     {//: m_mapped_file(boost::interprocess::create_only, fname, 1024*1024*128) {} 
//       m_mapped_file = new boost::interprocess::managed_mapped_file(boost::interprocess::create_only, fname, 1024*1024*128);
//   }
//   void print_info()
//   {
//     //std::cout << "free/size = " << m_mapped_file->get_free_memory() << " / " << m_mapped_file->get_size() << std::endl;
//     std::cout << "check_sanity() == " << m_mapped_file->check_sanity() << std::endl;
//   }
// private:
//   boost::interprocess::managed_mapped_file* m_mapped_file;
// };



// notes for how to setup a good test
// take rank * 100 and make edges between (all local)
// Make one vert per rank a hub.

namespace hmpi = havoqgt::mpi;
using namespace havoqgt::mpi;

int main(int argc, char** argv) {
  CHK_MPI(MPI_Init(&argc, &argv)); 
  {
  int mpi_rank(0), mpi_size(0);
  CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
  CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );
  havoqgt::get_environment();
  
  if(mpi_rank == 0){
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    std::cout << "CMD line:";
    for(int i=0; i<argc; ++i) {
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
  std::stringstream fname;
  fname << "/l/ssd/graph_test_" << mpi_rank;
  //arena_type arena(fname.str().c_str());

  typedef boost::interprocess::managed_mapped_file arena_type;
  remove(fname.str().c_str());
  boost::interprocess::managed_mapped_file  asdf(boost::interprocess::create_only, fname.str().c_str(), 1024ULL*1024*1024*16);


  uint64_t num_vertices = 1;
  uint64_t vert_scale;  
  double   pa_beta;
  uint64_t hub_threshold;
  std::string type;
  if(argc != 5) {
    std::cerr << "usage: <RMAT/PA> <Scale> <PA-beta> <hub_threshold>" << std::endl;
    exit(-1);
  } else {
    type = argv[1];
    vert_scale    = boost::lexical_cast<uint64_t>(argv[2]);
    pa_beta       = boost::lexical_cast<double>(argv[3]);
    hub_threshold = boost::lexical_cast<uint64_t>(argv[4]);
  }
  num_vertices <<= vert_scale;
  if(mpi_rank == 0) {
    std::cout << "Building graph type: " << type << std::endl;
    std::cout << "Building graph Scale: " << vert_scale << std::endl;
    std::cout << "Hub threshold = " << hub_threshold << std::endl;
    std::cout << "PA-beta = " << pa_beta << std::endl;
  }

  std::vector< std::pair<uint64_t, uint64_t> > input_edges;

  if(type == "RMAT") {
    uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
    havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                      vert_scale, num_edges_per_rank,
                                      0.57, 0.19, 0.19, 0.05, true, false);
    input_edges.resize(num_edges_per_rank); //times 2 because its undirected
    std::copy(rmat.begin(), rmat.end(), input_edges.begin());
  } else if(type == "PA") {
    gen_preferential_attachment_edge_list(input_edges, uint64_t(5489), vert_scale, vert_scale+4, pa_beta, 0.0, MPI_COMM_WORLD);
  } else {
    std::cerr << "Unknown graph type: " << type << std::endl;  exit(-1);
  }

  typedef hmpi::delegate_partitioned_graph<arena_type> graph_type;
  graph_type graph(asdf, MPI_COMM_WORLD, input_edges, hub_threshold);

  //arena.print_info();

  {
    std::vector< std::pair<uint64_t, uint64_t> > empty(0);
    input_edges.swap(empty);
  }



  //
  // Calculate max degree
  uint64_t max_degree(0);
  for(graph_type::vertex_iterator vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    max_degree = std::max(max_degree, graph.degree(*vitr));
  }
  for(graph_type::controller_iterator citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    max_degree = std::max(max_degree, graph.degree(*citr));
  }
  uint64_t global_max_degree = havoqgt::mpi::mpi_all_reduce(max_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
  if(mpi_rank == 0) {
    std::cout << "Max Degree = " << global_max_degree << std::endl;
  }

  

  } //END Main MPI
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));                                                                                      
  CHK_MPI(MPI_Finalize());
  return 0;
}
