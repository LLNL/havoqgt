//
// Created by Iwabuchi, Keita on 4/24/18.
//

// g++ -std=c++11 -O3 -fopenmp ../src/degree_count.cpp -o degree_count

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <strstream>

constexpr size_t max_vid = 3563602788;

int main(int argc, char **argv) {

  std::string out_path(argv[1]);
  std::vector<std::string> edge_list_file;
  for (int i = 2; i < argc; ++i) {
    edge_list_file.emplace_back(argv[i]);
  }

  std::vector<uint64_t> degree_table(max_vid + 1, 0);
  {
#pragma omp parallel for
    for (size_t i = 0; i < edge_list_file.size(); ++i) {
      const auto &f = edge_list_file[i];
      std::ifstream ifs(f);
      std::cout << "Open " << f << std::endl;
      if (!ifs.is_open()) std::abort();

      for (std::string line; std::getline(ifs, line);) {
        std::istrstream is(line.c_str());

        uint64_t src;
        uint64_t dst;
        is >> src >> dst;

#pragma omp atomic
        ++degree_table[src];
      }
    }
  }

  std::ofstream ofs(out_path);
  for (uint64_t i = 0; i < degree_table.size(); ++i) {
    ofs << i << "\t" << degree_table[i] << "\n";
  }
  ofs.close();

  std::cout << "Done" << std::endl;

  return 0;
}