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
#include <havoqgt/rmat_edge_generator.hpp>
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

using namespace havoqgt;

void usage()  {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -s <int> -d <int> -o <string>\n"
         << " -s <int>    - RMAT graph Scale (default 17)\n"
         << " -o <string> - output graph base filename\n"
         << " -h          - print help and exit\n\n";

  }
}

void parse_cmd_line(int argc, char** argv, uint64_t& scale, std::string& output_filename) {
  if(havoqgt_env()->world_comm().rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_output_filename = false;
  scale = 17;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:o:h ")) != -1) {
     switch (c) {
       case 'h':
         prn_help = true;
         break;
      case 's':
         scale = atoll(optarg);
         break;
      case 'o':
         found_output_filename = true;
         output_filename = optarg;
         break;
      default:
         std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
         prn_help = true;
         break;
     }
   }
   if (prn_help || !found_output_filename) {
     usage();
     exit(-1);
   }
}

int main(int argc, char** argv) {

  int mpi_rank(0), mpi_size(0);

  havoqgt_init(&argc, &argv);
  {
    std::string                output_filename;
    int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_size = havoqgt_env()->world_comm().size();
    havoqgt::get_environment();

    if (mpi_rank == 0) {

      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
      havoqgt::get_environment().print();
    }
    havoqgt_env()->world_comm().barrier();

    uint64_t      num_vertices = 1;
    uint64_t      vert_scale;

    parse_cmd_line(argc, argv, vert_scale, output_filename);

    num_vertices <<= vert_scale;
    if (mpi_rank == 0) {
      std::cout << "Building Graph500"<< std::endl
        << "Building graph Scale: " << vert_scale << std::endl;
    }

    //Generate RMAT graph
    uint64_t num_edges_per_rank = num_vertices * 16 / mpi_size;
    havoqgt::rmat_edge_generator rmat(uint64_t(5489) + uint64_t(mpi_rank) * 3ULL,
                                      vert_scale, num_edges_per_rank,
                                      0.57, 0.19, 0.19, 0.05, true, false);


    //
    //
    uint64_t global_num_edges(0);
    uint64_t local_edges_generated(0);
    std::deque< std::pair< uint64_t, uint64_t >> local_edges;
    auto rmat_itr = rmat.begin();
    do {
      std::vector< std::vector< std::pair<uint64_t,uint64_t> > > ata_send(mpi_size);
      std::vector< std::vector< std::pair<uint64_t,uint64_t> > > ata_recv(mpi_size);
      uint64_t global_edges_to_gen = std::min(uint64_t(1024)*1024*mpi_size, num_vertices * 16 - global_num_edges);
      uint64_t local_edges_to_gen = global_edges_to_gen / mpi_size;
      if(global_edges_to_gen % mpi_size > mpi_rank) ++local_edges_to_gen;
      for(uint64_t i=0; i<local_edges_to_gen; ++i, ++rmat_itr) {
        std::pair<uint64_t, uint64_t> edge = *rmat_itr;
        local_edges_generated++;
        if(edge.first != edge.second) {
          if(edge.first > edge.second) std::swap(edge.first, edge.second);
          ata_send[(edge.first ^ edge.second) % mpi_size].push_back(edge);
        }
      }
      mpi_all_to_all(ata_send, ata_recv, MPI_COMM_WORLD);

      for(size_t r = 0; r < mpi_size; ++r) {
        for( auto edge : ata_recv[r]) {
          local_edges.push_back(edge);
        }
      }
      std::sort(local_edges.begin(), local_edges.end());
      local_edges.erase( std::unique(local_edges.begin(), local_edges.end() ), local_edges.end() );
      global_num_edges = mpi_all_reduce((uint64_t) local_edges.size(), std::plus<uint64_t>(), MPI_COMM_WORLD);
      if(mpi_rank == 0) {
        std::cout << "Total Edges gen = " << global_num_edges << std::endl;
      }
    } while(global_num_edges < num_vertices * 16);

    uint64_t global_edges_generated = mpi_all_reduce(local_edges_generated, std::plus<uint64_t>(), MPI_COMM_WORLD);
    if(mpi_rank == 0) {
      std::cout << "Total Edges gen = " << global_num_edges << std::endl;
      std::cout << "Edges should be = " << num_vertices * 16 << std::endl;
      std::cout << "global_edges_generated = " << global_edges_generated << std::endl;
      std::cout << "% duplicate = " << double(global_edges_generated) / double(global_num_edges) << std::endl;
      //std::cout << "havoqgt_env()->node_local_comm().size() = " << havoqgt_env()->node_local_comm().size() << std::endl;
    }

    //
    // Write out edges
    {
      std::vector<char> vec(16*1024*1024);
      std::stringstream fname;
      fname << output_filename << "_" << mpi_rank << "_of_" << mpi_size;
      std::ofstream ofs(fname.str().c_str());
      ofs.rdbuf()->pubsetbuf(&vec.front(), vec.size());
      //for(size_t i=0; i< havoqgt_env()->node_local_comm().size(); ++i) {
      //    if(i == havoqgt_env()->node_local_comm().rank()) {
      for(auto& edge : local_edges) {
        ofs << edge.first << " " << edge.second << "\n" << edge.second << " " << edge.first <<  "\n";
      }
      //    }

    //havoqgt_env()->world_comm().barrier();
    //  }
    }


    havoqgt_env()->world_comm().barrier();
  } //END Main MPI
  havoqgt_finalize();
  return 0;
}
