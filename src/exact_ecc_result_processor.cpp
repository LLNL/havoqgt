//
// Created by Iwabuchi, Keita on 2/22/18.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cassert>

using edge_list_t = std::vector<uint64_t>;
struct vertex_data_t
{
  uint16_t ecc;
  edge_list_t edge_list;
};
using graph_t = std::unordered_map<uint64_t, vertex_data_t>;

int main(void)
{
  std::string vid_file;
  std::string ecc_file;

  std::ifstream fin_vid(vid_file);
  std::ifstream fin_ecc(ecc_file);

  assert(fin_vid.is_open() && fin_ecc.is_open());

  graph_t graph;

  while (true) {
    uint64_t vid;
    uint16_t ecc;

    assert(fin_vid >> vid);
    assert(fin_ecc >> ecc);

    assert()
    if (fin_vid.eof() || fin_ecc.eof()) {
      assert(fin_vid.eof() && fin_ecc.eof());
      break;
    }
  }
}