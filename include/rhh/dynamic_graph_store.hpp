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

#ifndef HAVOQGT_INCLUDE_RHH_DYNAMIC_GRAPH_STORE_HPP
#define HAVOQGT_INCLUDE_RHH_DYNAMIC_GRAPH_STORE_HPP

#include <memory>
#include <rhh/hash.hpp>
#include <rhh/rhh_map.hpp>

namespace rhh {

template <typename vertex_type, typename vertex_value_type,
          typename edge_value_type,
          typename allocator_type = std::allocator<std::byte>>
class dynamic_graph_store {
 private:
  using edge_table_t = rhh_map<vertex_type, edge_value_type, hash<vertex_type>,
                               std::equal_to<vertex_type>, allocator_type>;

  struct vertex_data_type {
    vertex_value_type value;
    edge_table_t     edges;
  };

  using vertex_table_t =
      rhh_map<vertex_type, vertex_data_type, hash<vertex_type>,
              std::equal_to<vertex_type>, allocator_type>;

 public:
  explicit dynamic_graph_store(const allocator_type& allocator = allocator_type())
      : m_vertex_table(allocator) {}

 private:
  vertex_table_t m_vertex_table;
};

}  // namespace rhh
#endif  // HAVOQGT_INCLUDE_RHH_DYNAMIC_GRAPH_STORE_HPP
