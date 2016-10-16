/**
 * A grig graph with three rows and five columns, 15 vertices and 
 * 44 edges in total.
 * If a delegate partitiond graph is created with the parameter -d 4, 
 * vertices 6, 7 and 8 are the delegates. The third value represents the 
 * weight associated with the corresponding edge.
 */
std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> 
  grid_graph_weighted_edges() { 
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> graph;
  graph.push_back(std::make_tuple(0, 1, 10));
  graph.push_back(std::make_tuple(0, 5, 5));  
  graph.push_back(std::make_tuple(1, 0, 21));
  graph.push_back(std::make_tuple(1, 2, 31));
  graph.push_back(std::make_tuple(1, 6, 51));
  graph.push_back(std::make_tuple(2, 1, 16));
  graph.push_back(std::make_tuple(2, 3, 19));
  graph.push_back(std::make_tuple(2, 7, 53));
  graph.push_back(std::make_tuple(3, 2, 12));
  graph.push_back(std::make_tuple(3, 4, 10));
  graph.push_back(std::make_tuple(3, 8, 57));
  graph.push_back(std::make_tuple(4, 3, 13));
  graph.push_back(std::make_tuple(4, 9, 55));
  graph.push_back(std::make_tuple(5, 0, 57));
  graph.push_back(std::make_tuple(5, 6, 18));
  graph.push_back(std::make_tuple(5, 10, 59));
  graph.push_back(std::make_tuple(6, 1, 51));
  graph.push_back(std::make_tuple(6, 5, 12));
  graph.push_back(std::make_tuple(6, 7, 14));
  graph.push_back(std::make_tuple(6, 11, 50));
  graph.push_back(std::make_tuple(7, 2, 55));
  graph.push_back(std::make_tuple(7, 6, 14));
  graph.push_back(std::make_tuple(7, 8, 13));
  graph.push_back(std::make_tuple(7, 12, 51));
  graph.push_back(std::make_tuple(8, 3, 52));
  graph.push_back(std::make_tuple(8, 7, 13));
  graph.push_back(std::make_tuple(8, 9, 14));
  graph.push_back(std::make_tuple(8, 13, 55));
  graph.push_back(std::make_tuple(9, 4, 55));
  graph.push_back(std::make_tuple(9, 8, 16));
  graph.push_back(std::make_tuple(9, 14, 57));
  graph.push_back(std::make_tuple(10, 5, 58));
  graph.push_back(std::make_tuple(10, 11, 19));
  graph.push_back(std::make_tuple(11, 6, 54));
  graph.push_back(std::make_tuple(11, 10, 13));
  graph.push_back(std::make_tuple(11, 12, 12));
  graph.push_back(std::make_tuple(12, 7, 51));
  graph.push_back(std::make_tuple(12, 11, 12));
  graph.push_back(std::make_tuple(12, 13, 13));
  graph.push_back(std::make_tuple(13, 8, 54));
  graph.push_back(std::make_tuple(13, 12, 15));
  graph.push_back(std::make_tuple(13, 14, 18));
  graph.push_back(std::make_tuple(14, 9, 53));
  graph.push_back(std::make_tuple(14, 13, 19));
  return graph;
}

std::vector<uint64_t> grid_graph_weighted_edges_degree = 
{2, 3, 3, 3, 2,
 3, 4, 4, 4, 3,
 2, 3, 3, 3, 2
};

std::vector<uint64_t> grid_graph_weighted_edges_offset =
{0, 2, 5, 8, 11,
 13, 16, 20, 24, 28,
 31, 33, 36, 39, 42, 44
};

void print_adjacency_list(
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>& graph, 
  std::vector<uint64_t>& offset) {

  for (size_t i = 0; i < offset.size() - 1; ++i) { 
    std::cout << i << ": ";
    for(size_t j = offset[i]; j < offset[i+1] ; ++j) {
      std::cout << std::get<1>(graph[j]) << ", ";   
    }  
    std::cout << std::endl;      
  }
}
