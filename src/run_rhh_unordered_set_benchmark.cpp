// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <iostream>
#include <assert.h>
#include <unordered_set>
#include <random>
#include <chrono>
#include <boost/functional/hash.hpp>
#include <havoqgt/detail/hash.hpp>
#include <rhh/unordered_set.hpp>

std::chrono::high_resolution_clock::time_point duration_time_misec() {
  return std::chrono::high_resolution_clock::now();
}

int64_t duration_time_misec(const std::chrono::high_resolution_clock::time_point &tic) {
  auto duration_time = std::chrono::high_resolution_clock::now() - tic;
  return std::chrono::duration_cast<std::chrono::microseconds>(duration_time).count();
}

template <typename value_t, typename set_type>
double insert_benchmark_kernel(std::vector<value_t> &dataset, set_type &set) {
  uint64_t cnt_inserted = 0;
  auto start = duration_time_misec();
  for (auto key : dataset) {
    const auto ret = set.insert(key);
    cnt_inserted += ret.second;
  }
  const double duration_time = duration_time_misec(start) / 1000000.0;
  // std::cout << "Partial time (sec) " << duration_time << std::endl;
  // std::cout << "Ratio of inserts : " << static_cast<double >(cnt_inserted) / dataset.size() << std::endl;

  return duration_time;
}

template <typename set_type>
void run_insert_benchmark(unsigned random_seed, const size_t length_dataset, const size_t length_chunk, set_type &set) {
  std::mt19937_64 gen(random_seed);
  std::uniform_int_distribution<uint64_t> dis(0, length_dataset - 1);
  std::vector<uint64_t> dataset(length_chunk);

  std::cout << "length_dataset : " << length_dataset << std::endl;
  std::cout << "length_chunk : " << length_chunk << std::endl;

  double total_time = 0;
  for (size_t i = 0; i < length_dataset; i += length_chunk) {
    // std::cout << "\nChunk# : " << i / length_chunk << std::endl;
    dataset.clear();
    for (size_t j = 0; j < length_chunk; ++j) {
      dataset.push_back(dis(gen));
    }
    total_time += insert_benchmark_kernel(dataset, set);
  }
  std::cout << "Total time (sec) " << total_time << std::endl;
}

template <typename set_type>
void run_insert_pair_benchmark(unsigned random_seed,
                               const size_t length_dataset,
                               const size_t length_chunk,
                               set_type &set) {
  std::mt19937_64 gen(random_seed);
  std::uniform_int_distribution<uint64_t> dis(0, length_dataset - 1);
  std::vector<std::pair<uint64_t, uint64_t>> dataset(length_chunk);

  std::cout << "length_dataset : " << length_dataset << std::endl;
  std::cout << "length_chunk : " << length_chunk << std::endl;

  double total_time = 0;
  for (size_t i = 0; i < length_dataset; i += length_chunk) {
    // std::cout << "\nChunk# : " << i / length_chunk << std::endl;
    dataset.clear();
    for (size_t j = 0; j < length_chunk; ++j) {
      dataset.emplace_back(std::make_pair(dis(gen), dis(gen)));
    }
    total_time += insert_benchmark_kernel(dataset, set);
  }
  std::cout << "Total time (sec) " << total_time << std::endl;
}

#include <havoqgt/rmat_edge_generator.hpp>

template <typename set_type>
void run_insert_rmat_edge_benchmark(unsigned random_seed, const size_t vertex_scale, const size_t num_edges,
                                    const size_t length_chunk, set_type &set) {
  havoqgt::rmat_edge_generator rmat(random_seed, vertex_scale,
                                    num_edges, 0.57, 0.19, 0.19, 0.05, true, true);

  std::vector<std::pair<uint64_t, uint64_t>> dataset(length_chunk);

  std::cout << "vertex_scale : " << vertex_scale << std::endl;
  std::cout << "num_edges : " << num_edges << std::endl;
  std::cout << "length_chunk : " << length_chunk << std::endl;

  double total_time = 0;
  std::size_t count_loop = 0;
  for (auto edge = rmat.begin(), end = rmat.end(); edge != end;) {
    // std::cout << "\nChunk# : " << count_loop / length_chunk << std::endl;
    dataset.clear();
    for (size_t j = 0; j < length_chunk; ++j) {
      assert(edge != rmat.end());
      dataset.emplace_back(*edge);
      ++edge;
    }
    total_time += insert_benchmark_kernel(dataset, set);
    ++count_loop;
  }
  std::cout << "Total time (sec) " << total_time << std::endl;
}

// User defined hash function
struct hash_pair {
  uint64_t operator()(const std::pair<uint64_t, uint64_t> &element) const {
    return boost::hash<std::pair<uint64_t, uint64_t>>()(std::make_pair<uint64_t, uint64_t>(havoqgt::detail::hash32(element.first),
                                                                                           havoqgt::detail::hash32(element.second)));
  }
};

int main() {

  const uint32_t random_seed = static_cast<uint32_t>(std::chrono::system_clock::now().time_since_epoch().count());
  std::cout << "random seed " << random_seed << std::endl;

  std::size_t num_elements_log2 = 25;

  {
    std::cout << "\n----- std::unordered<uint64_t> -----" << std::endl;
    std::unordered_set<uint64_t> set;
    run_insert_benchmark(random_seed, (1ULL << num_elements_log2), (1ULL << 20), set);
  }
  {
    std::cout << "\n----- rhh::unordered<uint64_t> -----" << std::endl;
    rhh::unordered_set<uint64_t> set;
    run_insert_benchmark(random_seed, (1ULL << num_elements_log2), (1ULL << 20), set);
    std::cout << "Average probe distance = " << set.load_factor() << std::endl;
  }

  {
    std::cout << "\n----- std::unordered<pair<uint64_t, uint64_t>> -----" << std::endl;
    std::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_pair_benchmark(random_seed, (1ULL << num_elements_log2), (1ULL << 20), set);
  }

  {
    std::cout << "\n----- rhh::unordered<pair<uint64_t, uint64_t>> -----" << std::endl;
    rhh::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_pair_benchmark(random_seed, (1ULL << num_elements_log2), (1ULL << 20), set);
    std::cout << "Average probe distance = " << set.load_factor() << std::endl;
  }

  {
    std::cout << "\n----- std::unordered<pair<uint64_t, uint64_t>> rmat edge -----" << std::endl;
    std::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_rmat_edge_benchmark(random_seed, num_elements_log2 - 4, (1ULL << num_elements_log2), (1ULL << 20), set);
  }

  {
    std::cout << "\n----- rhh::unordered<pair<uint64_t, uint64_t>> rmat edge -----" << std::endl;
    rhh::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_rmat_edge_benchmark(random_seed, num_elements_log2 - 4, (1ULL << num_elements_log2), (1ULL << 20), set);
    std::cout << "Average probe distance = " << set.load_factor() << std::endl;
  }

  return 0;
}