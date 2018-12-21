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

#include <iostream>
#include <fstream>
#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/mpi.hpp>
#include <havoqgt/parallel_edge_stream_reader.hpp>
#include <havoqgt/distributed_db.hpp>
#include <havoqgt/graph_colour_dynamic.hpp>
#include <dynamic_graph_store/dist_dynamic_graphstore.hpp>

using mapped_file_type     = boost::interprocess::managed_mapped_file;
using segment_manager_type = havoqgt::distributed_db::segment_manager_type;

using vertex_id_type        = uint64_t;
using edge_property_type    = int;
using vertex_property_type  = uint32_t;

/// --- Baseline model --- ///
#include <dynamic_graph_store/baseline/baseline.hpp>
using baseline_type         = graphstore::graphstore_baseline<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              segment_manager_type>;

/// --- DegAwareRHH model --- ///
#include <dynamic_graph_store/degawarerhh/degawarerhh.hpp>
enum : size_t {
  middle_high_degree_threshold = 2 // must be more or equal than 1
};
using degawarerhh_type      = graphstore::degawarerhh<vertex_id_type,
                                                      vertex_property_type,
                                                      edge_property_type,
                                                      segment_manager_type,
                                                      middle_high_degree_threshold>;

/// --- chooes graphstore type ---- ///
/// you don't need to change any code other than here to change the graphstore model, Baseline or DegAwareRHH.
using graphstore_type       = degawarerhh_type; // baseline_type or degawarerhh_type

/// --- wrapper class to easy fit into the visitor program --- ///
using dist_graphstore_type  = graphstore::dist_dynamic_graphstore<graphstore_type>;

void usage()  {
  if(havoqgt::mpi_comm_rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
         << " -c            - Prints out count of each colour.\n"
         << " -v            - Verifies correctness of result.\n"
         << " -g <string>   - output graph base filename (required)\n"
         << " -e <string>   - filename that has a list of edgelist files (required)\n"
         << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_file,
                    std::vector<std::string>& edgelist_files,
                    bool* verify, bool* col_count) {
  if(havoqgt::mpi_comm_rank() == 0) {
    std::cout << "CMD line:";
    for (int i=0; i<argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_graph_filename_ = false;
  bool found_edgelist_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "cvg:e:h")) != -1) {
     switch (c) {
       case 'c':
         *col_count = true;
         break;
       case 'v':
         *verify = true;
         break;
       case 'h':
         prn_help = true;
         break;
       case 'g':
         found_graph_filename_ = true;
         graph_file = optarg;
         break;
       case 'e':
       {
         found_edgelist_filename = true;
         std::string fname(optarg);
         std::ifstream fin(fname);
         std::string line;
         if (!fin.is_open()) {
           std::cerr << fname << std::endl;
           HAVOQGT_ERROR_MSG("Unable to open a file");
         }
         while (std::getline(fin, line)) {
           edgelist_files.push_back(line);
         }
         break;
       }
       default:
         std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
         prn_help = true;
         break;
     }
   }
   if (prn_help || !found_graph_filename_ || !found_edgelist_filename) {
     usage();
     exit(-1);
   }
}


int main(int argc, char** argv) {
  double graph_capacity_gb_per_rank = 4.0;
  int mpi_rank(0), mpi_size(0);
  bool col_count = false;
  bool verify = false;

  havoqgt::init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);


  /// --- parse argments ---- ///
  std::string graph_file;
  std::vector<std::string> edgelist_files;
  parse_cmd_line(argc, argv, graph_file, edgelist_files, &verify, &col_count);
  MPI_Barrier(MPI_COMM_WORLD);


  /// --- init segment file --- ///
  std::stringstream fname_local;
  fname_local << graph_file << "_" << mpi_rank;
  havoqgt::distributed_db ddb(havoqgt::db_create(), fname_local.str().c_str(), graph_capacity_gb_per_rank);

  /// --- get a segment_manager --- ///
  segment_manager_type* segment_manager = ddb.get_segment_manager();

  // Get edgelist.
  havoqgt::parallel_edge_stream_reader<> edgelist(edgelist_files, true);

  /// --- allocate a graphstore --- ///
  graphstore_type graphstore(segment_manager);
  dist_graphstore_type dist_graphstore(&graphstore);

  // Run algorithm.
  havoqgt::graph_colour_dynamic<dist_graphstore_type, havoqgt::parallel_edge_stream_reader<>, vertex_property_type>(&dist_graphstore, &edgelist, verify);

  if (mpi_rank == 0) {
    std::cout << "Collecting stats...\n";
  }

  // Collect stats.
  uint64_t visited_total = 0;
  bool finished = false;
  uint32_t num_colours = 0;
  const int AGGR_SIZE = 100;

  std::vector<uint32_t> colour_counts = *(new std::vector<uint32_t>());

  // Agregates 100 colours per all reduce.
  for (uint32_t colour = 0; finished != true; colour += AGGR_SIZE) {
    std::vector<uint64_t> local_count = *(new std::vector<uint64_t>());
    local_count.resize(AGGR_SIZE, 0);
    std::vector<uint64_t> global_count = *(new std::vector<uint64_t>());
    global_count.resize(AGGR_SIZE, 0);
    colour_counts.resize(colour_counts.size() + AGGR_SIZE, 0);

    for (auto v = graphstore.vertices_begin(); v != graphstore.vertices_end(); v++) {
      auto v_col = v.property_data();
      if (v_col >= colour && v_col < (colour + AGGR_SIZE)) {
        local_count[v_col % AGGR_SIZE]++;
      }
    }
    havoqgt::mpi_all_reduce(local_count, global_count, std::plus<uint64_t>(), MPI_COMM_WORLD);

    for (uint32_t i = 0; i < AGGR_SIZE; i++) {
      // Note: colour 0 (colour = 0, i = 0) will have 0 count.
      if (global_count[i] != 0 || (colour == 0 && i == 0)) {
        colour_counts[colour + i] = global_count[i];
        visited_total += global_count[i];

        // Print out per-colour count if desired.
        if (col_count && mpi_rank == 0) {
          std::cout << "Colour " << colour + i << ": " << global_count[i]
                    << std::endl;
        }
      } else {
        num_colours = colour + i - 1;  // One offset since 0 is not a colour.
        finished = true;
        break;
      }
    }
  }

  if (mpi_rank == 0 && visited_total > 1) {
    std::cout << "Number of Colours = " << num_colours <<  std::endl
              << "Visited total = " << visited_total << std::endl
              ;  // << "GC Time = " << time_end - time_start << std::endl;
  }

  // Print graph size
  if (mpi_rank == 0) std::cout << "Graph size (GB)" << std::endl;
  for (int i = 0; i < mpi_size; ++i) {
    if (mpi_rank == i) {
      const size_t usages = segment_manager->get_size() - segment_manager->get_free_memory();
      std::cout << " rank " << mpi_rank << " : " << static_cast<double>(usages) / (1ULL << 30) << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  }  // END Main MPI

  return 0;
}

