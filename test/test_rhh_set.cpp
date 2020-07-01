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
#include <rhh/rhh_set.hpp>

namespace {

TEST(rhh_set_test, init) {
  rhh::rhh_set<char> set;

  ASSERT_EQ(set.capacity(), 1);
  ASSERT_EQ(set.size(), 0);
}

TEST(rhh_set_test, insert) {
  rhh::rhh_set<int> set;

  const int val1(10);
  ASSERT_EQ(set.at(set.insert(val1)), val1);
  ASSERT_EQ(set.size(), 1);

  const int val2(20);
  ASSERT_EQ(set.at(set.insert(val2)), val2);
  ASSERT_EQ(set.size(), 2);

  const int val3(30);
  ASSERT_EQ(set.at(set.insert(val3)), val3);
  ASSERT_EQ(set.size(), 3);
}

TEST(rhh_set_test, insertString) {
  rhh::rhh_set<std::string> set;

  const std::string val1("10");
  ASSERT_EQ(set.at(set.insert(val1)), val1);
  ASSERT_EQ(set.size(), 1);

  const std::string val2("20");
  ASSERT_EQ(set.at(set.insert(val2)), val2);
  ASSERT_EQ(set.size(), 2);

  const std::string val3("30");
  ASSERT_EQ(set.at(set.insert(val3)), val3);
  ASSERT_EQ(set.size(), 3);
}

TEST(rhh_set_test, erase) {
  rhh::rhh_set<int> set;

  ASSERT_GE(set.insert(10), 0);
  ASSERT_GE(set.insert(20), 0);
  ASSERT_GE(set.insert(10), 0);

  ASSERT_EQ(set.erase(10), 2);
  ASSERT_EQ(set.size(), 1);

  ASSERT_EQ(set.erase(20), 1);
  ASSERT_EQ(set.size(), 0);
}

TEST(rhh_set_test, find) {
  rhh::rhh_set<int> set;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set.insert(val1), 0);
  ASSERT_GE(set.insert(val1), 0);
  ASSERT_GE(set.insert(val2), 0);

  const auto pos1 = set.find(val1);
  ASSERT_TRUE(set.at(pos1) == val1 || set.at(pos1) == val1);

  const auto pos1_2 = set.find(val1, pos1 + 1);
  ASSERT_TRUE(set.at(pos1_2) == val1 || set.at(pos1_2) == val1);

  ASSERT_EQ(set.find(val1, pos1 + 1), set.find_next(pos1 + 1));
  ASSERT_EQ(set.find(val2), set.find_next(set.find_next(pos1 + 1) + 1));

  const auto pos2 = set.find(val2);
  ASSERT_EQ(set.at(pos2), val2);


  if (set.at(pos1) == val1) {
    set.erase_at(pos1);
  } else {
    set.erase_at(pos1_2);
  }
  ASSERT_EQ(set.at(set.find(val1)), val1);

  set.erase_at(pos2);
  ASSERT_EQ(set.find(val2), set.capacity());
}

TEST(rhh_set_test, clear) {
  rhh::rhh_set<int> set;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set.insert(val1), 0);
  ASSERT_GE(set.insert(val1), 0);
  ASSERT_GE(set.insert(val2), 0);

  set.clear();

  ASSERT_EQ(set.size(), 0);
  ASSERT_EQ(set.find(val1), set.capacity());
  ASSERT_EQ(set.find(val2), set.capacity());
}

TEST(rhh_set_test, swap) {
  rhh::rhh_set<int> set1;
  rhh::rhh_set<int> set2;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set1.insert(val1), 0);
  ASSERT_GE(set1.insert(val1), 0);
  ASSERT_GE(set1.insert(val2), 0);

  const int val3(20);
  ASSERT_GE(set2.insert(val3), 0);

  set1.swap(set2);

  ASSERT_LT(set1.find(val3), set1.capacity());
  ASSERT_LT(set2.find(val1), set2.capacity());
  ASSERT_LT(set2.find(val2), set2.capacity());
}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}