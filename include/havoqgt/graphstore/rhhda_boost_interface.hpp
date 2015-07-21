#ifndef RHHDA_BOOST_INTERFACE
#define RHHDA_BOOST_INTERFACE

/// http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/leda_conversion.html
///
///

template <typename GraphProperty, typename VertexProperty, typename EdgeProperty>
class rhhda_graphstore {
 public:
  typedef uint64_t vertex_descriptor;
  typedef std::pair <vertex_descriptor, vertex_descriptor> edge_descriptor;
  typedef boost::directed_tag directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  typedef adjacency_graph_tag traversal_category;

  typedef OutEdgeIteratorImpl out_edge_iterator;
  typedef out_edge_iterator adjacency_iterator;

  typedef size_t degree_size_type;
  typedef size_t vertices_size_type;
  typedef size_t edges_size_type;

  typedef void vertex_iterator;
  typedef void in_edge_iterator;
  typedef void edge_iterator;
};

#endif // RHHDA_BOOST_INTERFACE

