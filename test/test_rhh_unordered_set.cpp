// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>

#include <random>
#include <unordered_set>

#include <rhh/unordered_set.hpp>

TEST(rhh_unordered_set_test, insert) {
  rhh::unordered_set<uint64_t> set;

  ASSERT_TRUE(set.insert(1).second);
  ASSERT_FALSE(set.insert(1).second);
  ASSERT_TRUE(set.insert(2).second);
  ASSERT_FALSE(set.insert(2).second);
  ASSERT_TRUE(set.insert(3).second);
  ASSERT_FALSE(set.insert(3).second);
}

TEST(rhh_unordered_set_test, erase_with_key) {
  rhh::unordered_set<uint64_t> set;

  ASSERT_EQ(set.erase(1), 0ULL);
  ASSERT_EQ(set.erase(2), 0ULL);

  set.insert(1);
  set.insert(2);

  ASSERT_EQ(set.erase(1), 1ULL);
  ASSERT_EQ(set.erase(1), 0ULL);
  ASSERT_EQ(set.size(), 1ULL);

  ASSERT_EQ(set.erase(2), 1ULL);
  ASSERT_EQ(set.erase(2), 0ULL);
  ASSERT_EQ(set.size(), 0ULL);

}

TEST(rhh_unordered_set_test, erase_with_iterator) {
  rhh::unordered_set<uint64_t> set;

  set.insert(1);
  set.insert(2);

  auto itr1_next = set.find(1);
  ++itr1_next;
  ASSERT_EQ(set.erase(set.find(1)), itr1_next);
  ASSERT_EQ(set.size(), 1ULL);

  auto itr2_next = set.find(2);
  ++itr2_next;
  ASSERT_EQ(set.erase(set.find(2)), itr2_next);
  ASSERT_EQ(set.size(), 0ULL);

}

TEST(rhh_unordered_set_test, erase_with_range) {
  rhh::unordered_set<uint64_t> set;


  set.insert(1);
  set.insert(2);
  set.insert(3);
  set.insert(4);

  auto itr1 = set.find(1);
  auto itr1_next = itr1;
  ++itr1_next;

  // Erase a value
  ASSERT_EQ(set.erase(itr1, itr1_next), itr1_next);
  ASSERT_EQ(set.size(), 3ULL);

  // Erase all values
  ASSERT_EQ(set.erase(set.begin(), set.end()), set.end());
  ASSERT_EQ(set.size(), 0ULL);
}

TEST(rhh_unordered_set_test, count) {
  rhh::unordered_set<uint64_t> set;

  ASSERT_EQ(set.count(1), 0ULL);
  set.insert(1);
  set.insert(1);
  ASSERT_EQ(set.count(1), 1ULL);
  set.erase(1);
  ASSERT_EQ(set.count(1), 0ULL);

  set.insert(1);
  ASSERT_EQ(set.count(2), 0ULL);
  set.insert(2);
  ASSERT_EQ(set.count(2), 1ULL);
  set.erase(2);
  ASSERT_EQ(set.count(2), 0ULL);
}

TEST(rhh_unordered_set_test, find) {
  rhh::unordered_set<uint64_t> set;

  set.insert(1);
  set.insert(2);

  const auto iterator1 = set.find(1);
  ASSERT_EQ(*iterator1, 1ULL);

  const auto iterator2 = set.find(2);
  ASSERT_EQ(*iterator2, 2ULL);
}

TEST(rhh_unordered_set_test, begin_end) {
  rhh::unordered_set<uint64_t> set;

  set.insert(1);
  set.insert(2);

  for (auto &val : set) {
    ASSERT_TRUE(val == 1 || val == 2);
  }
}

TEST(rhh_unordered_set_test, random_insert_and_erase) {
  std::unordered_set<uint64_t> set_std;
  rhh::unordered_set<uint64_t> set_rhh;

  unsigned seed = 1122;
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<uint64_t> dis(0, 256);

  for (int i = 0; i < 4096; i++) {
    const uint64_t key = dis(gen);

    if (i % 2 == 0) {
      ASSERT_EQ(set_rhh.insert(key).second, set_std.insert(key).second) << "i, key = " << i << ", " << key;
    } else {
      ASSERT_EQ(set_rhh.erase(key), set_std.erase(key)) << "i, key = " << i << ", " << key;
    }

    ASSERT_EQ(set_rhh.count(key), set_std.count(key)) << "i, key = " << i << ", " << key;
    ASSERT_EQ(set_rhh.size(), set_std.size()) << "i, key = " << i << ", " << key;
  }

  for (auto key : set_std) {
    ASSERT_EQ(set_rhh.count(key), 1ULL) << "key = " << key;
  }

  for (auto key : set_rhh) {
    ASSERT_EQ(set_std.count(key), 1ULL) << "key = " << key;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}