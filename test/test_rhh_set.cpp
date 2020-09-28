// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include <rhh/rhh_set.hpp>

namespace {

TEST(rhh_set_test, init) {
  rhh::rhh_set<char> set;

  ASSERT_EQ(set.capacity(), 1ULL);
  ASSERT_EQ(set.size(), 0ULL);
}

TEST(rhh_set_test, insert) {
  rhh::rhh_set<int> set;

  const int val1(10);
  ASSERT_EQ(set.at(set.insert(val1)), val1);
  ASSERT_EQ(set.size(), 1ULL);

  const int val2(20);
  ASSERT_EQ(set.at(set.insert(val2)), val2);
  ASSERT_EQ(set.size(), 2ULL);

  const int val3(30);
  ASSERT_EQ(set.at(set.insert(val3)), val3);
  ASSERT_EQ(set.size(), 3ULL);
}

TEST(rhh_set_test, insertString) {
  rhh::rhh_set<std::string> set;

  const std::string val1("10");
  ASSERT_EQ(set.at(set.insert(val1)), val1);
  ASSERT_EQ(set.size(), 1ULL);

  const std::string val2("20");
  ASSERT_EQ(set.at(set.insert(val2)), val2);
  ASSERT_EQ(set.size(), 2ULL);

  const std::string val3("30");
  ASSERT_EQ(set.at(set.insert(val3)), val3);
  ASSERT_EQ(set.size(), 3ULL);
}

TEST(rhh_set_test, erase) {
  rhh::rhh_set<int> set;

  ASSERT_GE(set.insert(10), 0ULL);
  ASSERT_GE(set.insert(20), 0ULL);
  ASSERT_GE(set.insert(10), 0ULL);

  ASSERT_EQ(set.erase(10), 2ULL);
  ASSERT_EQ(set.size(), 1ULL);

  ASSERT_EQ(set.erase(20), 1ULL);
  ASSERT_EQ(set.size(), 0ULL);
}

TEST(rhh_set_test, find) {
  rhh::rhh_set<int> set;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set.insert(val1), 0ULL);
  ASSERT_GE(set.insert(val1), 0ULL);
  ASSERT_GE(set.insert(val2), 0ULL);

  const auto pos1 = set.find(val1);
  ASSERT_EQ(set.at(pos1), val1);

  const auto pos1_2 = set.find(val1, pos1 + 1);
  ASSERT_EQ(set.at(pos1_2), val1);

  ASSERT_EQ(set.find(val1, pos1 + 1), set.find_any(pos1 + 1));
  ASSERT_EQ(set.find(val2), set.find_any(set.find_any(pos1 + 1) + 1));

  const auto pos2 = set.find(val2);
  ASSERT_EQ(set.at(pos2), val2);


  if (set.at(pos1) == val1) {
    set.erase_at(pos1);
  } else {
    set.erase_at(pos1_2);
  }
  ASSERT_EQ(set.at(set.find(val1)), val1);

  set.erase_at(pos2);
  ASSERT_EQ(set.find(val2), set.npos);
}

TEST(rhh_set_test, clear) {
  rhh::rhh_set<int> set;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set.insert(val1), 0ULL);
  ASSERT_GE(set.insert(val1), 0ULL);
  ASSERT_GE(set.insert(val2), 0ULL);

  set.clear();

  ASSERT_EQ(set.size(), 0ULL);
  ASSERT_EQ(set.find(val1), set.npos);
  ASSERT_EQ(set.find(val2), set.npos);
}

TEST(rhh_set_test, swap) {
  rhh::rhh_set<int> set1;
  rhh::rhh_set<int> set2;

  const int val1(10);
  const int val2(20);
  ASSERT_GE(set1.insert(val1), 0ULL);
  ASSERT_GE(set1.insert(val1), 0ULL);
  ASSERT_GE(set1.insert(val2), 0ULL);

  const int val3(20);
  ASSERT_GE(set2.insert(val3), 0ULL);

  set1.swap(set2);

  ASSERT_LT(set1.find(val3), set1.npos);
  ASSERT_LT(set2.find(val1), set2.npos);
  ASSERT_LT(set2.find(val2), set2.npos);
}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}