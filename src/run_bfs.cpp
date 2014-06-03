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
#include <havoqgt/breadth_first_search.hpp>
#include <havoqgt/single_source_shortest_path.hpp>
#include <havoqgt/page_rank.hpp>
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

  typedef boost::interprocess::managed_mapped_file mapped_t;
  typedef mapped_t::segment_manager segment_manager_t;

  remove(fname.str().c_str());
  mapped_t  asdf(bip::create_only, fname.str().c_str(), 1024ULL*1024*1024*16);


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

  typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
  bip::allocator<void, segment_manager_t> alloc_inst (asdf.get_segment_manager());
  graph_type graph(alloc_inst, MPI_COMM_WORLD, input_edges, hub_threshold);

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

  //BFS Experiments
  {
    graph_type::vertex_data<uint8_t, segment_manager_t >* bfs_level_data
          = graph.create_vertex_data< uint8_t, segment_manager_t >(asdf.get_segment_manager());

    graph_type::vertex_data<uint64_t, segment_manager_t >* bfs_parent_data
          = graph.create_vertex_data< uint64_t, segment_manager_t >(asdf.get_segment_manager());

    //arena.print_info();
    //
    //  Run BFS experiments
    double time(0);
    int count(0);
    uint64_t isource=0;
    for(int nrbfs = 0; nrbfs < 16; ++nrbfs) {
      graph_type::vertex_locator source = graph.label_to_locator(isource);
      uint64_t global_degree(0);
      do {
        uint64_t local_degree = 0;
        source = graph.label_to_locator(isource++);
        if(source.is_delegate()) break;
        if(mpi_rank == source.owner()) {
          local_degree = graph.degree(source);
        }
        global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
      } while(global_degree == 0);
      if(mpi_rank == source.owner()) {
        std::cout << "Starting vertex = " << isource << std::endl;
        std::cout << "delegate? = " << source.is_delegate() << std::endl;
        std::cout << "local_id = " << source.local_id() << std::endl;
        std::cout << "degree = " << graph.degree(source) << std::endl;
      }

      bfs_level_data->reset(128);

      MPI_Barrier( MPI_COMM_WORLD );
      double time_start = MPI_Wtime();
      hmpi::breadth_first_search(graph, *bfs_level_data, *bfs_parent_data, source);
      MPI_Barrier( MPI_COMM_WORLD );
      double time_end = MPI_Wtime();

      uint64_t visited_total(0);
      for(uint64_t level = 0; level < 15; ++level) {
        uint64_t local_count(0);
        graph_type::vertex_iterator vitr;
        for(vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
          if((*bfs_level_data)[*vitr] == level) {
            ++local_count;
          }
        }

        //Count the controllers!
        graph_type::controller_iterator citr;
        for(citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
          if((*bfs_level_data)[*citr] == level) {
            ++local_count;
          }
        }
        uint64_t global_count = mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
        visited_total += global_count;
        if(mpi_rank == 0 && global_count > 0) {
          std::cout << "Level " << level << ": " << global_count << std::endl;
        }
        if(global_count == 0) break;
      }

      if(mpi_rank == 0) {
        if(visited_total > 1000) {
          std::cout << "Visited total = " << visited_total << ", percentage visited = " << double(visited_total) / double(num_vertices) * 100 << "%" << std::endl;
          std::cout << "BFS Time = " << time_end - time_start << std::endl;
          time += time_end - time_start;
          ++count;
        }
      }
    }
    if(mpi_rank == 0) {
      std::cout << "Count BFS = " << count << std::endl;
      std::cout << "AVERAGE BFS= " << time / double(count) << std::endl;
    }

    //arena.print_info();
  } //End BFS Test
  //arena.print_info();
  { //PageRank Test
    graph_type::vertex_data<double, segment_manager_t >* pr_data
          = graph.create_vertex_data< double, segment_manager_t >(asdf.get_segment_manager());

    MPI_Barrier( MPI_COMM_WORLD );
    double time_start = MPI_Wtime();
    for(int i=0; i<10; ++i) {
    pr_data->reset(0);

    hmpi::page_rank(graph, *pr_data);

    }
    MPI_Barrier( MPI_COMM_WORLD );
    double time_end = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "Average Page Rank Time = " << (time_end - time_start) / double(10) << std::endl;
    }

    // std::cout << "Here?" << std::endl;

    // graph_type::vertex_iterator vitr;
    // for(vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
    //   assert( (*pr_data)[*vitr] == graph.degree(*vitr) );
    //   if( (*pr_data)[*vitr] != graph.degree(*vitr) ) {
    //     std::cout << "ERROR Reg: " << (*pr_data)[*vitr] << " != " << graph.degree(*vitr) << std::endl;
    //     exit(-1);
    //   }
    // }
    // graph_type::controller_iterator citr;
    // for(citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
    //   assert( (*pr_data)[*citr] == graph.degree(*citr) );
    //   if( (*pr_data)[*citr] != graph.degree(*citr) ) {
    //     std::cout << "ERROR Controller" << std::endl;
    //     exit(-1);
    //   }
    // }
    std::cout << "Ending PageRank Test!!" << std::endl;
  } //End Pagerank Test

  std::cout << "Between Tests !!" << std::endl;
  { //SSSP Test
    graph_type::vertex_data<uint64_t, segment_manager_t >* sssp_data
          = graph.create_vertex_data< uint64_t, segment_manager_t >(asdf.get_segment_manager());
    graph_type::edge_data<uint32_t, segment_manager_t >* edge_data
          = graph.create_edge_data< uint32_t, segment_manager_t >(asdf.get_segment_manager());

    //Generate Random Edge Values
    boost::random::mt19937 rng(mpi_rank*13);
    boost::random::uniform_int_distribution<> rand_num(1,1024*1024*1024);
    typedef graph_type::edge_data<uint32_t, segment_manager_t >::iterator edge_iterator;
    for(edge_iterator itr = edge_data->owned_begin(); itr != edge_data->owned_end(); ++itr)
    {
      *itr = rand_num(rng);
    }
    for(edge_iterator itr = edge_data->delegate_begin(); itr != edge_data->delegate_end(); ++itr)
    {
      *itr = rand_num(rng);
    }


    MPI_Barrier( MPI_COMM_WORLD );
    double time_start = MPI_Wtime();
    uint64_t isource=0;
    for(int i=0; i<10; ++i) {
      graph_type::vertex_locator source = graph.label_to_locator(isource);
      uint64_t global_degree(0);
      do {
        uint64_t local_degree = 0;
        source = graph.label_to_locator(isource++);
        if(source.is_delegate()) break;
        if(mpi_rank == source.owner()) {
          local_degree = graph.degree(source);
        }
        global_degree = mpi_all_reduce(local_degree, std::greater<uint64_t>(), MPI_COMM_WORLD);
      } while(global_degree == 0);
      if(mpi_rank == source.owner()) {
        std::cout << "Starting vertex = " << isource << std::endl;
        std::cout << "delegate? = " << source.is_delegate() << std::endl;
        std::cout << "local_id = " << source.local_id() << std::endl;
        std::cout << "degree = " << graph.degree(source) << std::endl;
      }
      sssp_data->reset(std::numeric_limits<uint64_t>::max());
      MPI_Barrier( MPI_COMM_WORLD );
      double time_start = MPI_Wtime();
      hmpi::single_source_shortest_path(graph, *sssp_data, *edge_data, source);
      MPI_Barrier( MPI_COMM_WORLD );
      double time_end = MPI_Wtime();


      uint64_t local_count(0);
      graph_type::vertex_iterator vitr;
      for(vitr = graph.vertices_begin(); vitr != graph.vertices_end(); ++vitr) {
        if((*sssp_data)[*vitr] > 0) {
          ++local_count;
        }
      }

      //Count the controllers!
      graph_type::controller_iterator citr;
      for(citr = graph.controller_begin(); citr != graph.controller_end(); ++citr) {
        if((*sssp_data)[*citr] > 0) {
          ++local_count;
        }
      }
      uint64_t global_count = mpi_all_reduce(local_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      if(global_count == 0) break;

      if(mpi_rank == 0) {
        if(global_count > 100) {
          std::cout << "Visited total = " << global_count << ", percentage visited = " << double(global_count) / double(num_vertices) * 100 << "%" << std::endl;
          std::cout << "BFS Time = " << time_end - time_start << std::endl;
        }
      }
    }
  } //End SSSP Test

  } //END Main MPI
  CHK_MPI(MPI_Barrier(MPI_COMM_WORLD));
  CHK_MPI(MPI_Finalize());
  return 0;
}
