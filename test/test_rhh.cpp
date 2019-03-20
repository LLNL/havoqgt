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

#include <gtest/gtest.h>
#include <rhh/dynamic_robin_hood_hashing.hpp>

using map_type = rhh::detail::dynamic_robin_hood_hashing<uint32_t, uint32_t>;

namespace {
TEST(dynamic_robin_hood_hashing_test, init) {
  map_type map;

  ASSERT_EQ(map.capacity(), 1);
  ASSERT_EQ(map.size(), 0);
}


TEST(dynamic_robin_hood_hashing_test, insert) {
  map_type map;

  typename map_type::value_type val1{1, 1};
  ASSERT_EQ(map.at(map.insert(val1)), val1);
  ASSERT_EQ(map.size(), 1);

  typename map_type::value_type val2{1, 2};
  ASSERT_EQ(map.at(map.insert(val2)), val2);
  ASSERT_EQ(map.size(), 2);

  typename map_type::value_type val3{2, 1};
  ASSERT_EQ(map.at(map.insert(val3)), val3);
  ASSERT_EQ(map.size(), 3);
}

TEST(dynamic_robin_hood_hashing_test, erase) {
  map_type map;

  ASSERT_GE(map.insert(typename map_type::value_type{1, 1}), 0);
  ASSERT_GE(map.insert(typename map_type::value_type{1, 2}), 0);
  ASSERT_GE(map.insert(typename map_type::value_type{2, 3}), 0);

  ASSERT_EQ(map.erase(1), 2);
  ASSERT_EQ(map.size(), 1);

  ASSERT_EQ(map.erase(2), 1);
  ASSERT_EQ(map.size(), 0);
}

TEST(dynamic_robin_hood_hashing_test, find) {
  map_type map;

  typename map_type::value_type val1{1, 1};
  typename map_type::value_type val2{1, 2};
  typename map_type::value_type val3{2, 3};
  ASSERT_GE(map.insert(val1), 0);
  ASSERT_GE(map.insert(val2), 0);
  ASSERT_GE(map.insert(val3), 0);

  const auto pos1 = map.find(1);
  ASSERT_TRUE(map.at(pos1) == val1 || map.at(pos1) == val2);

  const auto pos2 = map.find(1, pos1 + 1);
  ASSERT_TRUE(map.at(pos2) == val1 || map.at(pos2) == val2);

  const auto pos3 = map.find(2);
  ASSERT_EQ(map.at(pos3), val3);


  if (map.at(pos1) == val1) {
    map.erase_at(pos1);
  } else {
    map.erase_at(pos2);
  }
  ASSERT_EQ(map.at(map.find(1)), val2);

  map.erase_at(pos3);
  ASSERT_EQ(map.find(2), map.capacity());
}

TEST(dynamic_robin_hood_hashing_test, clear) {
  map_type map;

  ASSERT_GE(map.insert(typename map_type::value_type{1, 1}), 0);
  ASSERT_GE(map.insert(typename map_type::value_type{1, 2}), 0);
  ASSERT_GE(map.insert(typename map_type::value_type{2, 3}), 0);

  map.clear();

  ASSERT_EQ(map.size(), 0);
  ASSERT_EQ(map.find(1), map.capacity());
  ASSERT_EQ(map.find(2), map.capacity());
}

TEST(dynamic_robin_hood_hashing_test, swap) {
  map_type map1;
  map_type map2;

  typename map_type::value_type val1{1, 1};
  typename map_type::value_type val2{1, 2};
  typename map_type::value_type val3{2, 3};

  ASSERT_GE(map1.insert(val1), 0);
  ASSERT_GE(map1.insert(val2), 0);
  ASSERT_GE(map1.insert(val3), 0);

  typename map_type::value_type val4{2, 4};
  ASSERT_GE(map2.insert(val4), 0);

  map1.swap(map2);

  ASSERT_TRUE(map2.at(map2.find(1)) == val1 || map2.at(map2.find(1)) == val2);
  ASSERT_TRUE(map2.at(map2.find(1, map2.find(1) + 1)) == val1 || map2.at(map2.find(1, map2.find(1) + 1)) == val2);
  ASSERT_EQ(map2.at(map2.find(2)), val3);

  ASSERT_EQ(map1.at(map1.find(2)), val4);
}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}