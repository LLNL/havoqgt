// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include <rhh/rhh_map.hpp>

using map_type = rhh::rhh_map<uint32_t, uint32_t>;

namespace {

TEST(rhhMapTest, init) {
  map_type map;

  ASSERT_EQ(map.capacity(), 1ULL);
  ASSERT_EQ(map.size(), 0ULL);
}


TEST(rhhMapTest, insert) {
  map_type map;

  typename map_type::value_type val1{1, 1};
  ASSERT_EQ(map.at(map.insert(val1)), val1);
  ASSERT_EQ(map.size(), 1ULL);

  typename map_type::value_type val2{1, 2};
  ASSERT_EQ(map.at(map.insert(val2)), val2);
  ASSERT_EQ(map.size(), 2ULL);

  typename map_type::value_type val3{2, 1};
  ASSERT_EQ(map.at(map.insert(val3)), val3);
  ASSERT_EQ(map.size(), 3ULL);
}

TEST(rhhSetTest, insertString) {
  rhh::rhh_map<std::string, std::string> map;

  const std::tuple<std::string, std::string> val1("10", "100");
  ASSERT_EQ(map.at(map.insert(val1)), val1);
  ASSERT_EQ(map.size(), 1ULL);

  const std::tuple<std::string, std::string> val2("20", "200");
  ASSERT_EQ(map.at(map.insert(val2)), val2);
  ASSERT_EQ(map.size(), 2ULL);

  const std::tuple<std::string, std::string> val3("30", "300");
  ASSERT_EQ(map.at(map.insert(val3)), val3);
  ASSERT_EQ(map.size(), 3ULL);
}

TEST(rhhMapTest, erase) {
  map_type map;

  ASSERT_GE(map.insert(typename map_type::value_type{1, 1}), 0ULL);
  ASSERT_GE(map.insert(typename map_type::value_type{1, 2}), 0ULL);
  ASSERT_GE(map.insert(typename map_type::value_type{2, 3}), 0ULL);

  ASSERT_EQ(map.erase(1), 2ULL);
  ASSERT_EQ(map.size(), 1ULL);

  ASSERT_EQ(map.erase(2), 1ULL);
  ASSERT_EQ(map.size(), 0ULL);
}

TEST(rhhMapTest, find) {
  map_type map;

  typename map_type::value_type val1{1, 1};
  typename map_type::value_type val2{1, 2};
  typename map_type::value_type val3{2, 3};
  ASSERT_GE(map.insert(val1), 0ULL);
  ASSERT_GE(map.insert(val2), 0ULL);
  ASSERT_GE(map.insert(val3), 0ULL);

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
  ASSERT_EQ(map.find(2), map.npos);
}

TEST(rhhMapTest, clear) {
  map_type map;

  ASSERT_GE(map.insert(typename map_type::value_type{1, 1}), 0ULL);
  ASSERT_GE(map.insert(typename map_type::value_type{1, 2}), 0ULL);
  ASSERT_GE(map.insert(typename map_type::value_type{2, 3}), 0ULL);

  map.clear();

  ASSERT_EQ(map.size(), 0ULL);
  ASSERT_EQ(map.find(1), map.npos);
  ASSERT_EQ(map.find(2), map.npos);
}

TEST(rhhMapTest, swap) {
  map_type map1;
  map_type map2;

  typename map_type::value_type val1{1, 1};
  typename map_type::value_type val2{1, 2};
  typename map_type::value_type val3{2, 3};

  ASSERT_GE(map1.insert(val1), 0ULL);
  ASSERT_GE(map1.insert(val2), 0ULL);
  ASSERT_GE(map1.insert(val3), 0ULL);

  typename map_type::value_type val4{2, 4};
  ASSERT_GE(map2.insert(val4), 0ULL);

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