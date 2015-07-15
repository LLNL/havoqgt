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
#ifndef RHHDA_HPP
#define RHHDA_HPP

#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/rhhda/rhhda_defs.hpp>
#include <havoqgt/graphstore/rhhda/rhhda_common.hpp>
#include <havoqgt/graphstore/rhhda/rhh_container.hpp>
//#include <havoqgt/graphstore/rhhda/rhhda_child.h>
//#include <havoqgt/graphstore/rhhda/rhhda_parent.h>

namespace graphstore {
namespace rhhda {

/// Robin Hood Hashing Degree Aware
template <typename VertexIDType, typename vertex_attribute_type, typename edge_attribute_type>
class rhhda {

 public:
  using rhhda_child_controller_type = rhh_child_controller<VertexIDType, edge_attribute_type>;
  using rhhda_child_type = rhh_container<rhhda_child_controller_type>;
  using rhhda_parent_controller_type = rhh_parent_controller<VertexIDType, VertexIDType, vertex_attribute_type, edge_attribute_type, rhhda_child_type>;
  using rhhda_parent_type = rhh_container<rhhda_parent_controller_type>;

  explicit rhhda(segment_manager_t* segment_manager)
    : m_rhh(),
      m_rhh_controller(segment_manager),
      m_num_holding_edges(0)
  {
    disp_configuration();
    m_rhh->allocate_self_with_capacity(m_rhh_controller, kRHHDAMainInitialSize);
  }

  ~rhhda() {
    m_rhh.free(m_allocators);
  }

  /// --- Explicitly Deleted Copy and Move Functions -- ///
  rhhda(const rhhda&) =delete;
  rhhda& operator=(const rhhda&) =delete;
  rhhda(const rhhda&&) =delete;
  rhhda& operator=(rhhda&& old_obj) =delete;

  bool add_vertex(VertexIDType& vertex, vertex_attribute_type& vertex_attribute)
  {
  }

  bool add_edge(VertexIDType& src, VertexIDType& dst, edge_attribute_type& edge_attribute)
  {
    m_temp_value_block.value = dst;
    m_temp_value_block.attribute = edge_attribute;
    const bool is_succ =  m_rhh.insert_uniquely(m_allocators, src, m_temp_value_block, m_num_holding_edges);
    m_num_holding_edges += is_succ;
    return is_succ;
  }

//  bool erase_edge(VertexIDType& src, VertexIDType& dst)
//  {
//    const bool is_succ = m_rhh.delete_item(m_allocators, src, dst);
//    m_num_holding_edges -= is_succ;
//    return is_succ;
//  }

  size_t size()
  {
    return m_num_holding_edges;
  }

  void disp_status()
  {
//    m_rhh.disp_profileinfo();
  }

  void dump_detail_status(std::string prefix)
  {
//    {
//      std::stringstream fname;
//      fname << prefix << "rhhda_adjlist_length.log";
//      std::ofstream fout(fname);
//      m_rhh.fprint_value_lengths(fout);
//      fout.close();
//    }

//    {
//      std::stringstream fname;
//      fname << prefix << "rhhda_adjlists_prbdist.log";
//      std::ofstream fout(fname);
//      m_rhh.fprint_adjlists_prbdist(fout);
//      fout.close();
//    }

//    {
//      std::stringstream fname;
//      fname << prefix << "rhhda_adjlists_depth.log";
//      std::ofstream fout(fname);
//      m_rhh.fprint_adjlists_depth(fout);
//      fout.close();
//    }
  }

private:
  rhhda_parent_type m_rhh;
  rhhda_parent_controller_type m_rhh_controller;
  size_t m_num_holding_edges;
  m_rhh::rhh_controlloer_type::element::value_type m_temp_value_block;
};
}
}
#endif // RHHDA_HPP

