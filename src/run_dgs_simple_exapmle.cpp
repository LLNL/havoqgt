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

#include <havoqgt/parallel_edge_stream_reader.hpp>
#include <havoqgt/distributed_db.hpp>

using segment_manager_type = havoqgt::distributed_db::segment_manager_type;

using vertex_id_type        = uint64_t;
using edge_property_type    = int;
using vertex_property_type  = int;

#if 1
#include <dynamic_graph_store/baseline/baseline.hpp>
using graphstore_type       = graphstore::graphstore_baseline<vertex_id_type,
                                                              vertex_property_type,
                                                              edge_property_type,
                                                              segment_manager_type>;
#else
#include <dynamic_graph_store/degawarerhh/degawarerhh.hpp>
using graphstore_type       = graphstore::degawarerhh<vertex_id_type,
                                                      vertex_property_type,
                                                      edge_property_type,
                                                      segment_manager_type>;
#endif

using namespace havoqgt;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -g <string>   - output graph base filename (required)\n"
              << " -e <string>   - filename that has a list of edgelist files (required)\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char **argv, std::string &segmentfile_name, std::vector<std::string> &edgelist_files) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD line:";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool found_segmentfile_name_ = false;
  bool found_edgelist_filename = false;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:e:h")) != -1) {
    switch (c) {
      case 'h':prn_help = true;
        break;
      case 's':found_segmentfile_name_ = true;
        segmentfile_name = optarg;
        break;
      case 'e': {
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
      default:std::cerr << "Unrecognized option: " << c << ", ignore." << std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_segmentfile_name_ || !found_edgelist_filename) {
    usage();
    exit(-1);
  }
}

int main(int argc, char **argv) {
  int mpi_rank(0), mpi_size(0);

  havoqgt::init(&argc, &argv);
  {
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    if (mpi_rank == 0) {
      std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    /// --- parse argments ---- ///
    std::string segmentfile_name;
    std::vector<std::string> edgelist_files;
    parse_cmd_line(argc, argv, segmentfile_name, edgelist_files);
    MPI_Barrier(MPI_COMM_WORLD);


    /// --- create a segument file --- ///
    size_t graph_capacity_per_rank_gb = std::pow(2, 1);
    havoqgt::distributed_db ddb(havoqgt::db_create(), segmentfile_name.c_str(), graph_capacity_per_rank_gb);

    /// --- allocate a graphstore --- ///
    graphstore_type graphstore(ddb.get_segment_manager());


    /// ------- insert edges using parallel_edge_list_reader() ------- ///
    {
      /// --- setup a parallel edgelist reader --- ///
      havoqgt::parallel_edge_stream_reader<> edgelist(edgelist_files, true);

      for (const auto edge : edgelist) {
        vertex_id_type src = std::get<0>(edge);
        vertex_id_type dst = std::get<1>(edge);
        edge_property_type weight = 0;

        // uniquely insert a edge; return true if the edge is inserted (duplicated edge wasn't found)
        bool is_inserted = graphstore.insert_edge(src, dst, weight);

        graphstore.edge_property_data(src, dst) = weight + 1; // update edge's property
        // edge_property_data() return a reference

      }
    }



    /// ------- delete edges and update vertices' property data ------- ///
    {
      for (int i = 0; i < 10; ++i) {
        vertex_id_type src = i % 2;
        vertex_id_type dst = i;

        vertex_property_type v_prop = graphstore.vertex_property_data(src);  // get a vertex property data or
        graphstore.vertex_property_data(src) = v_prop;                       // update a vertex property data.
        // vertex_property_data() return a reference
        if (i % 2)
          bool is_erased_edge = graphstore.erase_edge(src, dst);             // return true if the edge is deleted
      }
    }


    /// ------- iterat over an adjacent-list ------- ///
    {
      vertex_id_type src_vrtx = 0;
      for (auto adj_edges = graphstore.adjacent_edge_begin(src_vrtx), end = graphstore.adjacent_edge_end(src_vrtx);
           adj_edges != end;
           ++adj_edges) {
        adj_edges.property_data() = 1; // update edge weight
        std::cout << "destination vertex: " << adj_edges.target_vertex() << ", weight: " << adj_edges.property_data()
                  << std::endl;
      }
    }

    std::cout << "END" << std::endl;
  }  // END Main MPI

  return 0;
}

