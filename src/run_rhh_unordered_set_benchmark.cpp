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

#include <iostream>
#include <assert.h>
#include <unordered_set>
#include <random>
#include <chrono>

#include <rhh/unordered_set.hpp>

std::chrono::high_resolution_clock::time_point duration_time_misec()
{
  return std::chrono::high_resolution_clock::now();
}

int64_t duration_time_misec(const std::chrono::high_resolution_clock::time_point& tic)
{
  auto duration_time = std::chrono::high_resolution_clock::now() - tic;
  return std::chrono::duration_cast<std::chrono::microseconds>(duration_time).count();
}


template <typename value_t, typename set_type>
double insert_benchmark_kernel(std::vector<value_t>& dataset, set_type& set)
{
  uint64_t cnt_inserted = 0;
  auto start = duration_time_misec();
  for (auto key : dataset) {
    const auto ret = set.insert(key);
    cnt_inserted += ret.second;
  }
  const double duration_time = duration_time_misec(start) / 1000000.0;
  std::cout << "Partial time (sec) " << duration_time << std::endl;
  std::cout << "Ratio of inserts : " << static_cast<double >(cnt_inserted) / dataset.size() << std::endl;

  return duration_time;
}

template <typename set_type>
void run_insert_benchmark(unsigned random_seed, const size_t length_dataset, const size_t length_chunk, set_type& set)
{
  std::mt19937_64 gen(random_seed);
  std::uniform_int_distribution<uint64_t> dis(0, length_dataset - 1);
  std::vector<uint64_t> dataset(length_chunk);

  std::cout << "length_dataset : " << length_dataset << std::endl;
  std::cout << "length_chunk : " << length_chunk << std::endl;

  double total_time = 0;
  for (size_t i = 0; i < length_dataset; i += length_chunk)
  {
    std::cout << "\nChunk# : " << i / length_chunk << std::endl;
    dataset.clear();
    for (size_t j = 0; j < length_chunk; ++j) {
      dataset.push_back(dis(gen));
    }
    total_time += insert_benchmark_kernel(dataset, set);
  }
  std::cout << "Total time (sec) " << total_time << std::endl;
}

template <typename set_type>
void run_insert_benchmark_pair(unsigned random_seed, const size_t length_dataset, const size_t length_chunk, set_type& set)
{
  std::mt19937_64 gen(random_seed);
  std::uniform_int_distribution<uint64_t> dis(0, length_dataset - 1);
  std::vector<std::pair<uint64_t, uint64_t>> dataset(length_chunk);

  std::cout << "length_dataset : " << length_dataset << std::endl;
  std::cout << "length_chunk : " << length_chunk << std::endl;

  double total_time = 0;
  for (size_t i = 0; i < length_dataset; i += length_chunk)
  {
    std::cout << "\nChunk# : " << i / length_chunk << std::endl;
    dataset.clear();
    for (size_t j = 0; j < length_chunk; ++j) {
      dataset.emplace_back(std::make_pair(dis(gen), dis(gen)));
    }
    total_time += insert_benchmark_kernel(dataset, set);
  }
  std::cout << "Total time (sec) " << total_time << std::endl;
}

// User defined hash function
struct hash_pair
{
  size_t operator()(const std::pair<uint64_t, uint64_t>& element) const
  {
    uint64_t seed = std::hash<uint64_t>()(element.first);
    seed ^=  std::hash<uint64_t>()(element.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;
  }
};


int main()
{

  const uint32_t random_seed = static_cast<uint32_t>(std::chrono::system_clock::now().time_since_epoch().count());
  std::cout << "random seed " << random_seed << std::endl;


  {
    std::cout << "\n----- std::unordered -----" << std::endl;
    std::unordered_set<uint64_t> set;
    run_insert_benchmark(random_seed, (1ULL << 25), (1ULL << 20), set);
  }
  {
    std::cout << "\n----- rhh::unordered -----" << std::endl;
    rhh::unordered_set<uint64_t> set;
    run_insert_benchmark(random_seed, (1ULL << 25), (1ULL << 20), set);
    std::cout << "Average probe distance = " << set.load_factor() << std::endl;
  }

  {
    std::cout << "\n----- std::unordered -----" << std::endl;
    std::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_benchmark_pair(random_seed, (1ULL << 25), (1ULL << 20), set);
  }

  {
    std::cout << "\n----- rhh::unordered -----" << std::endl;
    rhh::unordered_set<std::pair<uint64_t, uint64_t>, hash_pair> set;
    run_insert_benchmark_pair(random_seed, (1ULL << 25), (1ULL << 20), set);
    std::cout << "Average probe distance = " << set.load_factor() << std::endl;
  }

  return 0;
}