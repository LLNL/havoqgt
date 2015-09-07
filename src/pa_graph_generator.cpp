#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <boost/lexical_cast.hpp>

using vertex_type = uint64_t;

int main(int argc, char** argv)
{
  size_t num_verticess;
  size_t num_min_degree;
  std::string edgefile_name;

  if (argc < 4) {
    std::cerr << "usage: <num_verticess> <num_min_degree> <edgefile_name>"
    << " (argc:" << argc << " )." << std::endl;
    exit(-1);
  } else {
    int pos = 1;
    num_verticess  = boost::lexical_cast<uint64_t>(argv[pos++]);
    num_min_degree = boost::lexical_cast<uint64_t>(argv[pos++]);
    edgefile_name  = argv[pos++];
  }

  std::vector<std::pair<vertex_type, vertex_type>> edge_vec;

  /// --- init random generator --- ///
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << " use seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);

  /// --- generate init graph --- ///
  edge_vec.push_back(std::make_pair(0, 1));
  edge_vec.push_back(std::make_pair(0, 2));
  edge_vec.push_back(std::make_pair(0, 3));

  edge_vec.push_back(std::make_pair(1, 0));
  edge_vec.push_back(std::make_pair(1, 2));
  edge_vec.push_back(std::make_pair(1, 3));

  edge_vec.push_back(std::make_pair(2, 0));
  edge_vec.push_back(std::make_pair(2, 1));
  edge_vec.push_back(std::make_pair(2, 3));

  edge_vec.push_back(std::make_pair(3, 0));
  edge_vec.push_back(std::make_pair(3, 1));
  edge_vec.push_back(std::make_pair(3, 2));


  /// --- generate edges --- ///
  for (uint64_t v = 4; v < num_verticess; ++v) {
    std::uniform_int_distribution<vertex_type> dis(0, (edge_vec.size() - 1));
    for (uint64_t i = 0; i < num_min_degree; ++i) {
      vertex_type edge_no = dis(gen);
      vertex_type src = edge_vec[edge_no].first;
      edge_vec.push_back(std::make_pair(v, src));
      edge_vec.push_back(std::make_pair(src, v));
    }
  }

  std::ofstream ofs(edgefile_name);
  for_each(edge_vec.begin(), edge_vec.end(),
                  [&ofs](std::pair<vertex_type, vertex_type> edge){
                  ofs << edge.first << " " << edge.second << std::endl;
                  }
          );
  std::cout << "generated graph: " << edge_vec.size() << std::endl;


  /// --- count degree ---- ///
#if 0
  std::vector<size_t> deg_vec(num_verticess);
  for_each(edge_vec.begin(), edge_vec.end(),
                  [&deg_vec](std::pair<vertex_type, vertex_type> edge){
                    ++deg_vec[edge.first];
                  }
          );

  for_each(deg_vec.begin(), deg_vec.end(),
                  [](uint64_t d){
                  std::cout << d << std::endl;
                  }
          );
#endif

  return 0;

}
