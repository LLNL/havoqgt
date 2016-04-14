#ifndef CSR_GRAPH_HPP
#define CSR_GRAPH_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <utility>
#include <vector>
#include <stack>
#include <queue>
#include <limits>
#include <ctime>
#include <cstdint>
#include <random>
#include <chrono>
#include <algorithm>
#include <utility>
#include <tuple>
#include <algorithm>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/container/vector.hpp>
#include <boost/container/set.hpp>
#include <boost/range/algorithm.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>

namespace csr_graph {

template <typename index_type, typename vertex_type>
class csr_graph_container {
 public:

  using index_array_element_type =  index_type;
  using index_array_type         = boost::container::vector<index_array_element_type>;

  using adjlist_array_element_type = vertex_type;
  using adjlist_array_type = boost::container::vector<adjlist_array_element_type>;

  csr_graph_container(size_t nv, size_t ne) :
    m_num_vertices(nv),
    m_num_edges(ne),
    m_index_array(),
    m_adjlist_array()
  {
    m_index_array.resize(m_num_vertices + 1, 0);
    m_adjlist_array.resize(m_num_edges, adjlist_array_element_type());
  }

  ~csr_graph_container()
  {
    m_num_vertices = 0;
    m_num_edges = 0;
  }

  template<typename edge_list_type>
  void construct(edge_list_type& edge_list)
  {
    std::cout << "Construct index-array" << std::endl;
    const auto si = graphstore::utility::duration_time();
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
         edge_itr != edge_itr_end;
         ++edge_itr) {
      const vertex_type src = edge_itr->first;
      ++m_index_array[src+1];
    }
    for (size_t i = 0; i < m_num_vertices; ++i) {
      m_index_array[i+1] += m_index_array[i];
    }
    assert(m_index_array[m_num_vertices] == m_num_edges);
    std::cout << "finished:\t" << graphstore::utility::duration_time_sec(si) << std::endl;

    std::cout << "Construct adjlist-array" << std::endl;
    const auto sa = graphstore::utility::duration_time();
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
         edge_itr != edge_itr_end;
         ++edge_itr) {
      const vertex_type src = edge_itr->first;
      const vertex_type dst = edge_itr->second;
      m_adjlist_array[m_index_array[src]++] = dst;
    }
    for (size_t i = m_num_vertices; i > 0; --i) {
      m_index_array[i] = m_index_array[i-1];
    }
    m_index_array[0] = 0;
    assert(m_index_array[m_num_vertices] == m_num_edges);
    std::cout << "finished:\t" << graphstore::utility::duration_time_sec(sa) << std::endl;
  }

  typename adjlist_array_type::iterator adjacencylist(const vertex_type& src)
  {
    auto itr = m_adjlist_array.begin();
    const index_type s = m_index_array[src];
    itr += s;
    return itr;
  }

  typename adjlist_array_type::iterator adjacencylist_end(const vertex_type& src)
  {
    auto itr = m_adjlist_array.begin();
    const index_type s = m_index_array[src+1];
    itr += s;
    return itr;
  }


private:
  uint64_t m_num_vertices;
  uint64_t m_num_edges;
  index_array_type  m_index_array;
  adjlist_array_type m_adjlist_array;
};


///
/// \brief The prop_csr_graph_container class
///   property csr graph data structure
///
template <typename index_type, typename vertex_property_type, typename vertex_type, typename edge_property_type>
class prop_csr_graph_container {
 public:

  using index_array_element_type   = std::pair<index_type, vertex_property_type>;
  using index_array_type           = boost::container::vector<index_array_element_type>;

  using adjlist_array_element_type = vertex_type;
  using adjlist_array_type         = boost::container::vector<adjlist_array_element_type>;

  using edge_property_element_type = edge_property_type;
  using edge_property_array_type   = boost::container::vector<edge_property_element_type>;


  enum index_array_element {
    index = 0,
    vrt_prop = 1
  };

  prop_csr_graph_container(size_t nv, size_t ne) :
    m_num_vertices(nv),
    m_num_edges(ne),
    m_index_array(),
    m_adjlist_array(),
    m_edge_property_array()
  {
    index_array_element_type zero_index(0, vertex_property_type());
    m_index_array.resize(m_num_vertices + 1, zero_index);
    m_adjlist_array.resize(m_num_edges, adjlist_array_element_type());
    m_edge_property_array.resize(m_num_edges, edge_property_element_type());

  }

  ~prop_csr_graph_container()
  {
    m_num_vertices = 0;
    m_num_edges = 0;
  }

  template<typename edge_list_type>
  void construct(edge_list_type& edge_list)
  {
    std::cout << "Construct index-array" << std::endl;
    const auto si = graphstore::utility::duration_time();
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
         edge_itr != edge_itr_end;
         ++edge_itr) {
      const vertex_type src = std::get<0>(*edge_itr);
      ++std::get<index_array_element::index>(m_index_array[src+1]);
    }
    for (size_t i = 0; i < m_num_vertices; ++i) {
      std::get<index_array_element::index>(m_index_array[i+1]) += std::get<index_array_element::index>(m_index_array[i]);
    }
    assert(std::get<index_array_element::index>(m_index_array[m_num_vertices]) == m_num_edges);
    std::cout << "finished:\t" << graphstore::utility::duration_time_sec(si) << std::endl;

    std::cout << "Construct adjlist-array" << std::endl;
    const auto sa = graphstore::utility::duration_time();
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
         edge_itr != edge_itr_end;
         ++edge_itr) {
      const vertex_type src = std::get<0>(*edge_itr);
      const vertex_type dst = std::get<1>(*edge_itr);
      m_adjlist_array[std::get<index_array_element::index>(m_index_array[src])++] = dst;
    }
    for (size_t i = m_num_vertices; i > 0; --i) {
      std::get<index_array_element::index>(m_index_array[i]) = std::get<index_array_element::index>(m_index_array[i-1]);
    }
    std::get<index_array_element::index>(m_index_array[0]) = 0;
    assert(std::get<index_array_element::index>(m_index_array[m_num_vertices]) == m_num_edges);
    std::cout << "finished:\t" << graphstore::utility::duration_time_sec(sa) << std::endl;
  }

  void init_vertex_property(vertex_property_type& init_val)
  {
    for (index_type i = 0; i < m_num_vertices; ++i) {
      std::get<index_array_element::vrt_prop>(m_index_array[i]) = init_val;
    }
  }

  vertex_property_type& vertex_property(const vertex_type& src)
  {
    return std::get<index_array_element::vrt_prop>(m_index_array[src]);
  }

  typename adjlist_array_type::iterator adjacencylist(const vertex_type& src)
  {
    auto itr = m_adjlist_array.begin();
    const index_type s = std::get<index_array_element::index>(m_index_array[src]);
    itr += s;
    return itr;
  }

  typename adjlist_array_type::iterator adjacencylist_end(const vertex_type& src)
  {
    auto itr = m_adjlist_array.begin();
    const index_type s = std::get<index_array_element::index>(m_index_array[src+1]);
    itr += s;
    return itr;
  }

 private:
  uint64_t m_num_vertices;
  uint64_t m_num_edges;
  index_array_type  m_index_array;
  adjlist_array_type m_adjlist_array;
  edge_property_array_type m_edge_property_array;
};


///
/// \brief The hyper_csr_graph_container class
///   hyper property csr graph data structure
///
template <typename index_type, typename vertex_property_type, typename vertex_type, typename edge_property_type>
class hyper_prop_csr_graph_container {
 public:

  enum index_array_element {
    index = 0,
    vrt_prop = 1
  };

  enum : index_type {
    kNotFound = std::numeric_limits<index_type>::max()
  };

  using vertex_array_element_type  = vertex_type;
  using vertex_array_type          = boost::container::vector<vertex_array_element_type>;

  using index_array_element_type   = std::pair<index_type, vertex_property_type>;
  using index_array_type           = boost::container::vector<index_array_element_type>;

  using adjlist_array_element_type = vertex_type;
  using adjlist_array_type         = boost::container::vector<adjlist_array_element_type>;

  using edge_property_element      = edge_property_type;
  using edge_property_array_type   = boost::container::vector<edge_property_element>;


  hyper_prop_csr_graph_container(size_t nv, size_t ne) :
    m_num_vertices(nv),
    m_num_edges(ne),
    m_vertex_array(),
    m_index_array(),
    m_adjlist_array(),
    m_edge_property_array()
  {
    m_vertex_array.reserve(m_num_vertices); /// this is only reserve
    index_array_element_type zero_index(0, vertex_property_type());
    m_index_array.resize(m_num_vertices + 1, zero_index);
    m_adjlist_array.resize(m_num_edges, adjlist_array_element_type());
    m_edge_property_array.resize(m_num_edges, edge_property_element());
  }

  ~hyper_prop_csr_graph_container()
  {
    m_num_vertices = 0;
    m_num_edges = 0;
  }

//  struct lessthan {
//    inline bool operator()(const index_array_element_type &x, const index_array_element_type &y) const {
//      return std::get<index_array_element::vertex>(x) < std::get<index_array_element::vertex>(y);
//    }
//  };

  template<typename edge_list_type>
  void construct(edge_list_type& edge_list)
  {
    std::cout << "Construct index-array" << std::endl;

    {
      std::cout << "Load vertices" << std::endl;
      const auto start = graphstore::utility::duration_time();
      boost::container::set<vertex_type> set;
      for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
           edge_itr != edge_itr_end;
           ++edge_itr) {
        const vertex_type src = edge_itr->first;
        if (set.insert(src).second) { /// new vertex
          m_vertex_array.push_back(src);
        }
      }
      std::cout << "finished:\t" << graphstore::utility::duration_time_sec(start) << std::endl;
    }

    {
      /// --- sort vertices to use binary search --- ///
      std::cout << "Sort" << std::endl;
      const auto start = graphstore::utility::duration_time();
      std::sort(m_vertex_array.begin(), m_vertex_array.end());
      std::cout << "finished:\t" << graphstore::utility::duration_time_sec(start) << std::endl;
    }

    {
      std::cout << "Count degree" << std::endl;
      for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
           edge_itr != edge_itr_end;
           ++edge_itr) {
        const vertex_type src = edge_itr->first;
        const index_type index = compute_index(src);
        ++(std::get<index_array_element::index>(m_index_array[index]));
      }
    }

    {
      std::cout << "Shift" << std::endl;
      const auto start = graphstore::utility::duration_time();
      for (size_t i = 0; i < m_num_vertices; ++i) {
        std::get<index_array_element::index>(m_index_array[i+1]) += std::get<index_array_element::index>(m_index_array[i]);
      }
      for (size_t i = m_num_vertices; i > 0; --i) {
        std::get<index_array_element::index>(m_index_array[i]) = std::get<index_array_element::index>(m_index_array[i-1]);
      }
      std::get<index_array_element::index>(m_index_array[0]) = 0;
      assert(std::get<index_array_element::index>(m_index_array[m_num_vertices]) == m_num_edges);
      std::cout << "finished:\t" << graphstore::utility::duration_time_sec(start) << std::endl;
    }

    {
      std::cout << "Construct adjlist-array" << std::endl;
      const auto start = graphstore::utility::duration_time();
      for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end();
           edge_itr != edge_itr_end;
           ++edge_itr) {
        const vertex_type src = edge_itr->first;
        const vertex_type dst = edge_itr->second;
        const index_type index = compute_index(src);
        m_adjlist_array[std::get<index_array_element::index>(m_index_array[index])++] = dst;
      }
      for (size_t i = m_num_vertices; i > 0; --i) {
        std::get<index_array_element::index>(m_index_array[i]) = std::get<index_array_element::index>(m_index_array[i-1]);
      }
      std::get<index_array_element::index>(m_index_array[0]) = 0;
      assert(std::get<index_array_element::index>(m_index_array[m_num_vertices]) == m_num_edges);
      std::cout << "finished:\t" << graphstore::utility::duration_time_sec(start) << std::endl;
    }

  }

  void init_vertex_property(vertex_property_type& init_val)
  {
    for (index_type i = 0; i < m_num_vertices; ++i) {
      std::get<index_array_element::vrt_prop>(m_index_array[i]) = init_val;
    }
  }

  vertex_property_type& vertex_property(const vertex_type& src)
  {
    const index_type index = compute_index(src);
    return std::get<index_array_element::vrt_prop>(m_index_array[index]);
  }


  typename adjlist_array_type::iterator adjacencylist(const vertex_type& src)
  {
    const index_type index = compute_index(src);
    const index_type s = std::get<index_array_element::index>(m_index_array[index]);
    auto itr = m_adjlist_array.begin();
    itr += s;
    return itr;
  }

  typename adjlist_array_type::iterator adjacencylist_end(const vertex_type& src)
  {
    const index_type index = compute_index(src);
    const index_type s = std::get<index_array_element::index>(m_index_array[index+1]);
    auto itr = m_adjlist_array.begin();
    itr += s;
    return itr;
  }

 private:

  index_type compute_index(const vertex_type& src) const
  {
    /// Since using existing find functions, e.g., std::find, allocate iterators and will cause overhead,
    /// we just use original binary search implementaion.
    return graphstore::utility::binary_search(m_vertex_array, m_vertex_array.size(), src);
  }

  uint64_t m_num_vertices;
  uint64_t m_num_edges;
  vertex_array_type m_vertex_array;
  index_array_type  m_index_array;
  adjlist_array_type m_adjlist_array;
  edge_property_array_type m_edge_property_array;
};


#if 0

template <typename index_t, typename vertex_t, typename segment_manager_t>
class csr_graph_container {
public:

  using graph_self_type = csr_graph_container<index_t, vertex_t, segment_manager_t>;
  using index_type = index_t;
  using vertex_type = vertex_t;
  using index_allocator = boost::interprocess::allocator<index_t, segment_manager_t>;
  using adjlist_allocator = boost::interprocess::allocator<vertex_t, segment_manager_t>;

  friend class VertexForwardIterator;
  template <typename Type>
  class VertexForwardIterator : public std::iterator<std::forward_iterator_tag, Type>
  {
    friend class csr_graph_container;
    using vertex_iterator_selftype = VertexForwardIterator<Type>;

  public:

    VertexForwardIterator(graph_self_type* graph) :
      m_graph(graph),
      m_crt_vrt(0)
    { }

    VertexForwardIterator(graph_self_type* graph, Type& vrt) :
      m_graph(graph),
      m_crt_vrt(vrt)
    { }

    VertexForwardIterator(const VertexForwardIterator& src)
    {
      m_graph   = src.m_graph;
      m_crt_vrt = src.m_crt_vrt;
    }

    void swap(VertexForwardIterator &other) noexcept
    {
      using std::swap;
      swap(m_graph, other.m_graph);
      swap(m_crt_vrt, other.m_crt_vrt);
    }

    VertexForwardIterator &operator++ () // Pre-increment
    {
      next_vertex();
      return *this;
    }

    VertexForwardIterator operator++ (int) // Post-increment
    {
      VertexForwardIterator tmp(*this);
      next_vertex();
      return tmp;
    }

    // two-way comparison: v.begin() == v.cbegin() and vice versa
    template<class OtherType>
    bool operator == (const VertexForwardIterator<OtherType> &graph) const
    {
      return is_equal(graph);
    }

    template<class OtherType>
    bool operator != (const VertexForwardIterator<OtherType> &graph) const
    {
      return !is_equal(graph);
    }

    Type& operator* () const
    {
      return m_graph->m_index_array[m_crt_vrt];
    }

    Type* operator-> () const
    {
      return &(m_graph->m_index_array[m_crt_vrt]);
    }

    // One way conversion: iterator -> const_iterator
    operator VertexForwardIterator<const Type>() const
    {
      return VertexForwardIterator<const Type>(m_graph, m_crt_vrt);
    }


    /// --- performance optimized methods --- ///
    inline bool is_end() const
    {
      return (m_crt_vrt == m_graph->num_vertices());
    }

  private:

    VertexForwardIterator()
    { }

    void next_vertex()
    {
      ++m_crt_vrt;
    }

    template<class OtherType>
    inline bool is_equal(const VertexForwardIterator<OtherType> &other) const
    {
      return (m_graph == other.m_graph) && (m_crt_vrt == other.m_crt_vrt);
    }

    graph_self_type* m_graph;
    Type m_crt_vrt;
  };
  // iteratorをtypedefする
  using vertex_iterator = VertexForwardIterator<vertex_type>;
  using const_vertex_iterator = VertexForwardIterator<const vertex_type>;


  friend class EdgeForwardIterator;
  template <typename Type>
  class EdgeForwardIterator : public std::iterator<std::forward_iterator_tag, Type>
  {
    friend class csr_graph_container;
    using edge_iterator_selftype = EdgeForwardIterator<Type>;

  public:

    EdgeForwardIterator(graph_self_type* graph, Type& src_vrt) :
      m_graph(graph),
      m_src_vrt(src_vrt),
      m_crt_off(0)
    { }

    EdgeForwardIterator(graph_self_type* graph, Type& vrt, size_t off) :
      m_graph(graph),
      m_src_vrt(vrt),
      m_crt_off(off)
    { }

    EdgeForwardIterator(const EdgeForwardIterator& src)
    {
      m_graph   = src.m_graph;
      m_src_vrt = src.m_src_vrt;
      m_crt_off = src.m_crt_off;
    }

    void swap(EdgeForwardIterator &other) noexcept
    {
      using std::swap;
      swap(m_graph, other.m_graph);
      swap(m_src_vrt, other.m_src_vrt);
      swap(m_crt_off, other.m_crt_off);
    }

    EdgeForwardIterator &operator++ () // Pre-increment
    {
      next_edge();
      return *this;
    }

    EdgeForwardIterator operator++ (int) // Post-increment
    {
      EdgeForwardIterator tmp(*this);
      next_edge();
      return tmp;
    }

    // two-way comparison: v.begin() == v.cbegin() and vice versa
    template<class OtherType>
    bool operator == (const EdgeForwardIterator<OtherType> &graph) const
    {
      return is_equal(graph);
    }

    template<class OtherType>
    bool operator != (const EdgeForwardIterator<OtherType> &graph) const
    {
      return !is_equal(graph);
    }

    Type& operator* () const
    {
      return m_graph->m_adj_list[m_graph->m_index_array[m_src_vrt] + m_crt_off];
    }

    Type* operator-> () const
    {
      return &(m_graph->m_adj_list[m_graph->m_index_array[m_src_vrt] + m_crt_off]);
    }

    // One way conversion: iterator -> const_iterator
    operator EdgeForwardIterator<const Type>() const
    {
      return EdgeForwardIterator<const Type>(m_graph, m_src_vrt, m_crt_off);
    }


    /// --- performance optimized methods --- ///
    inline bool is_end() const
    {
      return (m_crt_off == (m_graph->m_index_array[m_src_vrt + 1] - m_graph->m_index_array[m_src_vrt]));
    }

  private:

    EdgeForwardIterator()
    { }

    void next_edge()
    {
      ++m_crt_off;
    }

    template<class OtherType>
    inline bool is_equal(const EdgeForwardIterator<OtherType> &other) const
    {
      return (m_graph == other.m_graph) && (m_src_vrt == other.m_src_vrt) && (m_crt_off == other.m_crt_off);
    }

    graph_self_type* m_graph;
    Type m_src_vrt;
    size_t m_crt_off;
  };
  // iteratorをtypedefする
  using edge_iterator = EdgeForwardIterator<vertex_type>;
  using const_edge_iterator = EdgeForwardIterator<const vertex_type>;



  csr_graph_container(size_t max_vertex_id,
                      size_t num_edges,
                      segment_manager_t* segment_manager) :
    m_num_vertices(max_vertex_id + 1),
    m_num_edges(num_edges),
    m_index_array(nullptr),
    m_adj_list(nullptr),
    m_index_allocator(segment_manager),
    m_adjlist_allocator(segment_manager)
  {
    allocate();
  }

  csr_graph_container(std::string& prefix,
                      segment_manager_t* segment_manager) :
    m_num_vertices(0),
    m_num_edges(0),
    m_index_array(nullptr),
    m_adj_list(nullptr),
    m_index_allocator(segment_manager),
    m_adjlist_allocator(segment_manager)
  {
    /// Load a constructed graph from file
    load_graph_info_from_file(prefix);
    allocate();
    load_graph_from_file(prefix);
  }

  ~csr_graph_container()
  {
    deallocate();
  }

  inline size_t& num_vertices ()
  {
    return m_num_vertices;
  }

  inline size_t& num_edges ()
  {
    return m_num_edges;
  }

  inline index_type* index_array()
  {
    return m_index_array;
  }

  inline vertex_type*  adj_list()
  {
    return m_adj_list;
  }

  inline size_t degree(vertex_type src_vrt)
  {
    return m_index_array[src_vrt + 1] - m_index_array[src_vrt];
  }

  inline vertex_type target_vertex(vertex_type src_vrt, size_t offset)
  {
    return m_adj_list[m_index_array[src_vrt] + offset];
  }

  /// --- Lookup --- ///
  inline vertex_iterator begin()
  {
    return vertex_iterator(this);
  }

  inline vertex_iterator end()
  {
    return vertex_iterator(this, m_num_vertices);
  }

  inline edge_iterator begin_adjlist(vertex_type srt_vrt)
  {
    return edge_iterator(this, srt_vrt);
  }

  inline edge_iterator end_adjlist(vertex_type srt_vrt)
  {
    return edge_iterator(this, srt_vrt, degree(srt_vrt));
  }

private:

  void allocate()
  {
    std::cout << "Allocating index-array and adj-list" << std::endl;
    m_index_pointer = m_index_allocator.allocate(m_num_vertices + 1);
    m_adjlist_pointer = m_adjlist_allocator.allocate(m_num_edges);
    m_index_array = m_index_pointer.get();
    m_adj_list = m_adjlist_pointer.get();

    std::cout << "Allocate index_array:\t" <<  (double)(m_num_vertices + 1) * sizeof(index_type) / (1ULL<<30) << " GB" << std::endl;
    graphstore::utility::print_time();

    std::cout << "Zero clear index-array" << std::endl;
    std::memset(m_index_array, 0, (m_num_vertices + 1ULL) * sizeof(index_type));
    for (uint64_t i = 0; i < m_num_vertices; ++i) {
      m_index_array[i] = index_type();
    }
    graphstore::utility::print_time();

    std::cout << "Allocate adj_list:\t" <<  (double)(m_num_edges) * sizeof(vertex_type) / (1ULL<<30) << " GB" << std::endl;
    graphstore::utility::print_time();

    std::cout << "Touch adj_list" << std::endl;
    for (uint64_t p = 0; p < m_num_edges; p += 4096 / sizeof(vertex_type)) {
      m_adj_list[p] = vertex_type();
    }
    graphstore::utility::print_time();

  }

  void deallocate()
  {
    m_num_vertices = 0;
    m_num_edges = 0;
    m_index_allocator.deallocate(m_index_pointer, m_num_vertices + 1);
    m_adjlist_allocator.deallocate(m_adjlist_pointer, m_num_edges);
  }


  void load_graph_info_from_file(std::string& prefix)
  {
    std::cout << "Load info from a file" << std::endl;

    std::ifstream inf_file(prefix + "_inf");
    std::cout << prefix << "_inf" << std::endl;

    std::string line;
    assert(std::getline(inf_file, line));

    assert(std::getline(inf_file, line));
    std::stringstream ssline(line);
    ssline >> m_num_vertices >> m_num_edges;

    std::cout << "m_num_vertices:\t" << m_num_vertices << std::endl;
    std::cout << "#edges:\t" << m_num_edges << std::endl;

    inf_file.close();
  }

  void load_graph_from_file(std::string& prefix)
  {
    std::cout << "Load graph from files" << std::endl;

    std::ifstream idx_file(prefix + "_idx", std::ifstream::binary);
    std::ifstream adj_file(prefix + "_adj", std::ifstream::binary);

    assert(m_num_vertices > 0 && m_index_array);
    assert(m_num_edges    > 0 && m_adj_list);

    idx_file.read(reinterpret_cast<char*>(m_index_array),
                  sizeof(index_type) * (m_num_vertices+1));
    adj_file.read(reinterpret_cast<char*>(m_adj_list),
                  sizeof(vertex_type) * m_num_edges);

    std::cout << "Loading graph from files done" << std::endl;
    std::cout << m_index_array[m_num_vertices] << std::endl;

    adj_file.close();
    idx_file.close();
  }

  void dump_csr_graph(std::string& prefix)
  {
    std::ofstream idx_file(prefix + "_idx", std::ofstream::binary);
    std::ofstream adj_file(prefix + "_adj", std::ofstream::binary);
    std::ofstream inf_file(prefix + "_inf");

    inf_file << "#vertices #edges\n";
    inf_file << m_num_vertices << " " << m_num_edges << "\n";

    idx_file.write(reinterpret_cast<char*>(m_index_array),
                   sizeof(index_type) * (m_num_vertices+1));
    adj_file.write(reinterpret_cast<char*>(m_adj_list),
                   sizeof(vertex_type) * m_num_edges);

    inf_file.close();
    idx_file.close();
    adj_file.close();
  }


  size_t m_num_vertices;
  size_t m_num_edges;
  index_type*  m_index_array;
  vertex_type* m_adj_list;
  index_allocator m_index_allocator;
  adjlist_allocator m_adjlist_allocator;
  typename index_allocator::pointer m_index_pointer;
  typename adjlist_allocator::pointer m_adjlist_pointer;
};


template <typename edge_list_type, typename index_type, typename vertex_type, typename segment_manager_t>
void construct_csr(csr_graph_container<index_type, vertex_type, segment_manager_t>& csr,
                   edge_list_type& edge_list)
{
  std::cout << "Counting degree" << std::endl;
  for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
    const vertex_type src = edge_itr->first;
    ++csr.index_array()[src+1];
  }
  utilities::cumulate_array(csr.index_array(), csr.num_vertices());
  assert(csr.index_array()[csr.num_vertices()] == csr.num_edges());
  graphstore::utility::print_time();

  std::cout << "Constructing adj-list" << std::endl;
  for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
    const vertex_type src = edge_itr->first;
    const vertex_type dst = edge_itr->second;
    csr.adj_list()[csr.index_array()[src]++] = dst;
  }
  utilities::right_shift_aray(csr.index_array(), csr.num_vertices());
  csr.index_array()[0] = 0;
  assert(csr.index_array()[csr.num_vertices()] == csr.num_edges());
  graphstore::utility::print_time();

}

template <typename edge_list_type, typename index_type, typename index_prop_type, typename vertex_type, typename segment_manager_t>
void construct_csr_ip(csr_graph_container<std::pair<index_type, index_prop_type>, vertex_type, segment_manager_t>& csr,
                   edge_list_type& edge_list)
{
  std::cout << "Counting degree" << std::endl;
  for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
    const vertex_type src = edge_itr->first;
    ++csr.index_array()[src+1].first;
  }
  utilities::cumulate_array(csr.index_array(), csr.num_vertices());
  assert(csr.index_array()[csr.num_vertices()] == csr.num_edges());
  graphstore::utility::print_time();

  std::cout << "Constructing adj-list" << std::endl;
  for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
    const vertex_type src = edge_itr->first;
    const vertex_type dst = edge_itr->second;
    csr.adj_list()[csr.index_array()[src]++] = dst;
  }
  utilities::right_shift_aray(csr.index_array(), csr.num_vertices());
  csr.index_array()[0] = 0;
  assert(csr.index_array()[csr.num_vertices()] == csr.num_edges());
  graphstore::utility::print_time();

}

template <typename edge_list_type, typename index_type, typename vertex_type, typename segment_manager_t>
void construct_csr_from_sorted_edgelist(csr_graph_container<index_type, vertex_type, segment_manager_t>& csr,
                                        edge_list_type& edge_list)
{
  std::cout << "Constructing adj-list from sorted edge list" << std::endl;
  size_t cnt_edges = 0;
  vertex_type last_src_vertex = 0;
  for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
    const vertex_type src = edge_itr->first;
    const vertex_type dst = edge_itr->second;
    assert(last_src_vertex <= src); /// must be sorted
    assert(src < csr.num_vertices());   /// must be small than the max vertex ID

    ++csr.index_array()[src];
    csr.adj_list()[cnt_edges++] = dst;

    if (last_src_vertex < src)
      last_src_vertex = src;
  }
  std::cout << "last_src_vertex + 1 == csr.num_vertices() :" << last_src_vertex + 1 << " <= " << csr.num_vertices() << std::endl;
  assert(last_src_vertex + 1 <= csr.num_vertices());
  std::cout << "cnt_edges == csr.num_edges() :" << cnt_edges << " == " << csr.num_edges() << std::endl;
  assert(cnt_edges == csr.num_edges());
  graphstore::utility::print_time();

  std::cout << "Modificate index-array" << std::endl;
  utilities::cumulate_array(csr.index_array(), csr.num_vertices());
  utilities::right_shift_aray(csr.index_array(), csr.num_vertices());
  csr.index_array()[0] = 0;
  std::cout << "csr.index_array()[csr.num_vertices()] == csr.num_edges() :" << csr.index_array()[csr.num_vertices()] << " == " << csr.num_edges() << std::endl;
  // assert(csr.index_array()[csr.num_vertices()] == csr.num_edges());
  graphstore::utility::print_time();
}

#endif

}
#endif // CSR_GRAPH_HPP

