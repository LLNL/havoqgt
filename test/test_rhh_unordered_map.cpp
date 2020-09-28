// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

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

  ASSERT_EQ(map.erase(1), 0ULL);
  ASSERT_EQ(map.erase(2), 0ULL);

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  ASSERT_EQ(map.erase(1), 1ULL);
  ASSERT_EQ(map.erase(1), 0ULL);
  ASSERT_EQ(map.size(), 1ULL);

  ASSERT_EQ(map.erase(2), 1ULL);
  ASSERT_EQ(map.erase(2), 0ULL);
  ASSERT_EQ(map.size(), 0ULL);
}

TEST(rhh_unordered_map_test, erase_with_iterator) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  map.insert(value_t{1, 1});
  map.insert(value_t{2, 2});

  auto itr1_next = map.find(1);
  ++itr1_next;
  ASSERT_EQ(map.erase(map.find(1)), itr1_next);
  ASSERT_EQ(map.size(), 1ULL);

  auto itr2_next = map.find(2);
  ++itr2_next;
  ASSERT_EQ(map.erase(map.find(2)), itr2_next);
  ASSERT_EQ(map.size(), 0ULL);
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
  ASSERT_EQ(map.size(), 3ULL);

  // Erase all values
  ASSERT_EQ(map.erase(map.begin(), map.end()), map.end());
  ASSERT_EQ(map.size(), 0ULL);
}

TEST(rhh_unordered_map_test, count) {
  rhh::unordered_map<uint64_t, double> map;
  using value_t = rhh::unordered_map<uint64_t, double>::value_type;

  ASSERT_EQ(map.count(1), 0ULL);
  map.insert(value_t{1, 1});
  map.insert(value_t{1, 1});
  ASSERT_EQ(map.count(1), 1ULL);
  map.erase(1);
  ASSERT_EQ(map.count(1), 0ULL);

  map.insert(value_t{1, 1});
  ASSERT_EQ(map.count(2), 0ULL);
  map.insert(value_t{2, 2});
  ASSERT_EQ(map.count(2), 1ULL);
  map.erase(2);
  ASSERT_EQ(map.count(2), 0ULL);
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
    ASSERT_EQ(map_rhh.count(val.first), 1ULL) << "key = " << val.first;
  }

  for (auto val : map_rhh) {
    ASSERT_EQ(map_std.count(val.first), 1ULL) << "key = " << val.first;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}