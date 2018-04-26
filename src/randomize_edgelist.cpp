//
// Created by Iwabuchi, Keita on 4/24/18.
//

//
// Created by Iwabuchi, Keita on 4/24/18.
//

// g++ -std=c++11 -O3 -fopenmp ../src/randomize_edgelist.cpp -o randomize_edgelist

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <strstream>
#include <random>
#include <sstream>
#include <iomanip>
#include <utility>
#include <string>
#include <chrono>

int main(int argc, char **argv) {

  std::string out_path_base(argv[1]);
  size_t num_out_files = std::stod(argv[2]);
  std::vector<std::string> edge_list_file;
  for (int i = 3; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  std::vector<std::ofstream*> out_file_list(num_out_files);
  std::vector<omp_lock_t> m_locks(num_out_files);
  for (size_t i = 0; i < num_out_files; ++i) {
    std::stringstream file_name;
    file_name << out_path_base << std::setfill('0') << std::setw(5) << i;
    out_file_list[i] = new std::ofstream(file_name.str());
    omp_init_lock(&m_locks[i]);
  }

  {
#pragma omp parallel for
    for (size_t i = 0; i < edge_list_file.size(); ++i) {
      const auto &f = edge_list_file[i];
      std::ifstream ifs(f);
      std::cout << "Open " << f << std::endl;
      if (!ifs.is_open()) {
        std::cerr << "Can not open: " << f << std::endl;
        std::abort();
      }

      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::mt19937_64 gen(seed);
      std::uniform_int_distribution<uint64_t> dis(0, num_out_files - 1);

      uint64_t src;
      uint64_t dst;
      while(ifs >> src >> dst) {
        const uint64_t target_file_no = dis(gen);
        while (!omp_test_lock(&m_locks[target_file_no]));
        *out_file_list[target_file_no] << src << " " << dst << "\n";
        omp_unset_lock(&m_locks[target_file_no]);
      }
    }
  }

  for (auto fp : out_file_list) {
    fp->close();
  }

  std::cout << "Done" << std::endl;

  return 0;
}