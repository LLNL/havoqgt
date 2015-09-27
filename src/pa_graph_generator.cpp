#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <boost/lexical_cast.hpp>

using vertex_type = uint64_t;

#define GENERATE_EDGE_BOTH_DIRECTION 1
#define COUNT_DEGREE_TABLE 0

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

  std::cout << "GENERATE_EDGE_BOTH_DIRECTION: " << GENERATE_EDGE_BOTH_DIRECTION << std::endl;

  std::vector<std::pair<vertex_type, vertex_type>> edge_vec;

  /// --- init random generator --- ///
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "use seed: " << seed << std::endl;
  std::mt19937_64 gen(seed);


  /// --- generate init graph (aomplete graph) --- ///
  /// k + 1 vertices are requred
  /// so that each vertex has at least k edges
  const size_t num_min_vertex = num_min_degree + 1;
  for (uint64_t v1 = 0; v1 < num_min_vertex; ++v1) {
    for (uint64_t v2 = 0; v2 < num_min_vertex; ++v2) {
      if (v1 == v2) continue;
      edge_vec.push_back(std::make_pair(v1, v2));
#if GENERATE_EDGE_BOTH_DIRECTION
      edge_vec.push_back(std::make_pair(v2, v1));
#endif
    }
  }


  /// --- generate edges --- ///
  /// \brief dis_2
#if !GENERATE_EDGE_BOTH_DIRECTION
  std::uniform_int_distribution<int> dis_2(0, 1);
#endif
  for (uint64_t v = num_min_vertex; v < num_verticess; ++v) {
    for (uint64_t i = 0; i < num_min_degree; ++i) {
      std::uniform_int_distribution<vertex_type> dis(0, (edge_vec.size() - 1));
      vertex_type edge_no = dis(gen);
#if GENERATE_EDGE_BOTH_DIRECTION
      const vertex_type src = edge_vec[edge_no].first;
      edge_vec.push_back(std::make_pair(src, v));
      edge_vec.push_back(std::make_pair(v, src));
#else
      const vertex_type src = (dis_2(gen)) ? edge_vec[edge_no].first : edge_vec[edge_no].second;
      edge_vec.push_back(std::make_pair(v, src));
#endif
    }
  }

  std::ofstream ofs(edgefile_name);
  for_each(edge_vec.begin(), edge_vec.end(),
                  [&ofs](std::pair<vertex_type, vertex_type> edge){
                  ofs << edge.first << " " << edge.second << std::endl;
                  }
          );
  ofs.close();
  std::cout << "GENERATE_EDGE_BOTH_DIRECTION: " << GENERATE_EDGE_BOTH_DIRECTION << std::endl;
  std::cout << "generated graph (undirected graph): " << edge_vec.size() << std::endl;


  /// --- count degree and make ---- ///
#if COUNT_DEGREE_TABLE
  std::vector<size_t> deg_vec(num_verticess, 0);
  for_each(edge_vec.begin(), edge_vec.end(),
                  [&deg_vec](std::pair<vertex_type, vertex_type> edge){
                    ++deg_vec[edge.first];
#if !GENERATE_EDGE_BOTH_DIRECTION
                    ++deg_vec[edge.second];
#endif
                  }
          );


  auto max_elem = std::max_element(deg_vec.begin(), deg_vec.end());
  std::cout << "max degree =\t" << *max_elem << std::endl;
  std::cout << "-- degree table --" << std::endl;
  std::vector<size_t> deg_tbl_vec(*max_elem+1, 0);
  for_each(deg_vec.begin(), deg_vec.end(),
                  [&deg_tbl_vec](uint64_t d){
                    ++deg_tbl_vec[d];
                  }
          );

 for (auto itr = deg_tbl_vec.begin(), end = deg_tbl_vec.end(); itr != end; ++itr) {
   if (*itr > 0)
    std::cout << itr - deg_tbl_vec.begin() << ":\t" << *itr << std::endl;
 }
#endif

  return 0;

}
