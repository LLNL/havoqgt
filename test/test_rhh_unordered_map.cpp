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

#include <random>
#include <unordered_map>

#include <rhh/unordered_map.hpp>

TEST(rhh_unordered_map_test, insert) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  ASSERT_TRUE(map.insert(value_t{1, 1}).second);
  ASSERT_FALSE(map.insert(value_t{1, 1}).second);
  ASSERT_TRUE(map.insert(value_t{2, 2}).second);
  ASSERT_FALSE(map.insert(value_t{2, 2}).second);
}

TEST(rhh_unordered_map_test, erase_with_key) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  ASSERT_EQ(map.erase(1), 0);
  ASSERT_EQ(map.erase(2), 0);

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  ASSERT_EQ(map.erase(1), 1);
  ASSERT_EQ(map.erase(1), 0);
  ASSERT_EQ(map.size(), 1);

  ASSERT_EQ(map.erase(2), 1);
  ASSERT_EQ(map.erase(2), 0);
  ASSERT_EQ(map.size(), 0);
}

TEST(rhh_unordered_map_test, erase_with_iterator) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  auto itr1_next = map.find(1);
  ++itr1_next;
  ASSERT_EQ(map.erase(map.find(1)), itr1_next);
  ASSERT_EQ(map.size(), 1);

  auto itr2_next = map.find(2);
  ++itr2_next;
  ASSERT_EQ(map.erase(map.find(2)), itr2_next);
  ASSERT_EQ(map.size(), 0);
}

TEST(rhh_unordered_map_test, erase_with_range) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});
  map.insert(value_t{3, 3});
  map.insert(value_t{4, 4});

  auto itr1 = map.find(1);
  auto itr1_next = itr1;
  ++itr1_next;

  // Erase a value
  ASSERT_EQ(map.erase(itr1, itr1_next), itr1_next);
  ASSERT_EQ(map.size(), 3);

  // Erase all values
  ASSERT_EQ(map.erase(map.begin(), map.end()), map.end());
  ASSERT_EQ(map.size(), 0);
}

TEST(rhh_unordered_map_test, count) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  ASSERT_EQ(map.count(1), 0);
  map.insert(value_t{1, 1});
  map.insert(value_t{1, 1});
  ASSERT_EQ(map.count(1), 1);
  map.erase(1);
  ASSERT_EQ(map.count(1), 0);

  map.insert(value_t{1, 1});
  ASSERT_EQ(map.count(2), 0);
  map.insert(value_t{2, 2});
  ASSERT_EQ(map.count(2), 1);
  map.erase(2);
  ASSERT_EQ(map.count(2), 0);
}

TEST(rhh_unordered_map_test, find) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  const auto iterator1 = map.find(1);
  ASSERT_EQ(*iterator1, value_t(1, 1));

  const auto iterator2 = map.find(2);
  ASSERT_EQ(*iterator2, value_t(2, 2));
}

TEST(rhh_unordered_map_test, begin_end) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  for (auto &val : map) {
    ASSERT_TRUE(val == value_t(1, 1) || val == value_t(2, 2));
  }
}

TEST(rhh_unordered_map_test, random_insert_and_erase) {
  std::unordered_map<uint64_t, double> map_std;
  rhh::unordered_map<uint64_t, double> map_rhh;
  using value_t = std::pair<uint64_t, double>;

  unsigned seed = 1122;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<uint64_t> dis(0, 256);

  for (int i = 0; i < 4096; i++) {
    const uint64_t key = dis(gen);
    auto val = value_t{key, key};

    if (i % 2 == 0) {
      ASSERT_EQ(map_rhh.insert(val).second, map_std.insert(val).second) << "i, key = " << i << ", " << key;
    } else {
      ASSERT_EQ(map_rhh.erase(key), map_std.erase(key)) << "i, key = " << i << ", " << key;
    }

    ASSERT_EQ(map_rhh.count(key), map_std.count(key)) << "i, key = " << i << ", " << key;
    ASSERT_EQ(map_rhh.size(), map_std.size()) << "i, key = " << i << ", " << key;
  }

  for (auto val : map_std) {
    ASSERT_EQ(map_rhh.count(val.first), 1) << "key = " << val.first;
  }

  for (auto val : map_rhh) {
    ASSERT_EQ(map_std.count(val.first), 1) << "key = " << val.first;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}