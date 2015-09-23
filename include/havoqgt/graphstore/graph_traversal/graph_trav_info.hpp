#ifndef GRAPH_TRAV_UTILITIES_HPP
#define GRAPH_TRAV_UTILITIES_HPP

namespace graph_trv_info {

template<typename vertex_type>
struct trv_inf {
  const size_t num_vertices;
  const size_t num_edges;
  bool* is_visited;
  vertex_type* tree;
  size_t count_visited_vertices;
  size_t count_visited_edges;

  trv_inf(size_t max_vertex_id, size_t _num_edges)
    : num_vertices(max_vertex_id+1),
      num_edges(_num_edges),
      is_visited(nullptr),
      tree(nullptr),
      count_visited_vertices(0),
      count_visited_edges(0)
  {
    init();
  }

  ~trv_inf()
  {
    delete[] is_visited;
    delete[] tree;
  }

  void init()
  {
    is_visited = new bool[num_vertices];
    std::cout << "Allocate is_visited:\t" <<  (double)num_vertices * sizeof(bool) / (1ULL<<30) << " GB" << std::endl;

    tree = new vertex_type[num_vertices];
    std::cout << "Allocate tree:\t" <<  (double)num_vertices * sizeof(vertex_type) / (1ULL<<30) << " GB" << std::endl;

    for (size_t i = 0 ; i < num_vertices; ++i) {
      is_visited[i] = false;
      tree[i] = std::numeric_limits<vertex_type>::max();
    }
    count_visited_vertices = 0;
    count_visited_edges = 0;
  }


};

}
#endif // GRAPH_TRAV_UTILITIES_HPP

