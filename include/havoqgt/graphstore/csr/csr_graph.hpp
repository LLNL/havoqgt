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

#include "graph_trav_info.hpp"
#include "utilities.hpp"

namespace csr_graph_struct {


template <typename edge_list_type, typename index_t, typename vertex_t>
class csr_graph {
 public:

  using graph_self_type = csr_graph<edge_list_type, index_t, vertex_t>;
  using index_type = vertex_t;
  using vertex_type = vertex_t;

  friend class VertexForwardIterator;
  template <typename Type>
  class VertexForwardIterator : public std::iterator<std::forward_iterator_tag, Type>
  {
    friend class csr_graph;
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
    friend class csr_graph;
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



  csr_graph(edge_list_type& edge_list) :
    m_num_vertices(edge_list.max_vertex_id() + 1),
    m_num_edges(edge_list.size()),
    m_index_array(nullptr),
    m_adj_list(nullptr)
  {
    allocate_graph();
#if 1
    construct_csr(edge_list);
#else
    construct_csr_from_sorted_edgelist(edge_list);
#endif
  }

  csr_graph(std::string& prefix) :
    m_num_vertices(0),
    m_num_edges(0),
    m_index_array(nullptr),
    m_adj_list(nullptr)
  {
    load_info(prefix);
    allocate_graph();
    load_graph(prefix);
  }

  ~csr_graph()
  {
    deallocate_graph();
  }

  inline size_t num_vertices ()
  {
    return m_num_vertices;
  }

  inline size_t num_edges ()
  {
    return m_num_edges;
  }

  inline index_type *index_array()
  {
    return m_index_array;
  }

  inline vertex_type *adj_list()
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

  void dump_csr_graph(std::string& prefix)
  {
    std::ofstream idx_file(prefix + "_idx", std::ofstream::binary);
    std::ofstream adj_file(prefix + "_adj", std::ofstream::binary);
    std::ofstream inf_file(prefix + "_inf");

    inf_file << "#vertices #edges\n";
    inf_file << m_num_vertices << " " << m_num_edges << "\n";

    idx_file.write(reinterpret_cast<char*>(m_index_array), sizeof(index_type)  * (m_num_vertices+1));
    adj_file.write(reinterpret_cast<char*>(m_adj_list),    sizeof(vertex_type) * m_num_edges);

    inf_file.close();
    idx_file.close();
    adj_file.close();
  }

 private:

  void allocate_graph()
  {
    std::cout << "Allocating index-array and adj-list" << std::endl;
    m_index_array = new index_type[m_num_vertices + 1];
    m_adj_list = new vertex_type[m_num_edges];

    std::cout << "Allocate index_array:\t" <<  (double)(m_num_vertices + 1) * sizeof(index_type) / (1ULL<<30) << " GB" << std::endl;
    graph_trv_utilities::print_time();

    std::cout << "Zero clear index-array" << std::endl;
    std::memset(m_index_array, 0, (m_num_vertices + 1ULL) * sizeof(index_type));
    graph_trv_utilities::print_time();

    std::cout << "Allocate adj_list:\t" <<  (double)(m_num_edges) * sizeof(vertex_type) / (1ULL<<30) << " GB" << std::endl;
    graph_trv_utilities::print_time();

    std::cout << "Touch adj_list" << std::endl;
    for (uint64_t p = 0; p < m_num_edges; p += 4096 / sizeof(vertex_type)) {
      m_adj_list[p] = 0;
    }
    graph_trv_utilities::print_time();

  }


  void deallocate_graph()
  {
    m_num_vertices = 0;
    m_num_edges = 0;
    delete[] m_index_array;
    delete[] m_adj_list;
  }

  void construct_csr(edge_list_type& edge_list)
  {
    std::cout << "Counting degree" << std::endl;
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
      const vertex_type src = edge_itr->first;
      ++m_index_array[src+1];
    }
    cumulate_array(m_index_array, m_num_vertices);
    assert(m_index_array[m_num_vertices] == m_num_edges);
    graph_trv_utilities::print_time();

    std::cout << "Constructing adj-list" << std::endl;
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
      const vertex_type src = edge_itr->first;
      const vertex_type dst = edge_itr->second;
      m_adj_list[m_index_array[src]++] = dst;
    }
    right_shift_aray(m_index_array, m_num_vertices);
    m_index_array[0] = 0;
    assert(m_index_array[m_num_vertices] == m_num_edges);
    graph_trv_utilities::print_time();

//    for (size_t i = 0; i < m_num_vertices+1; ++i) std::cout << m_index_array[i] << " ";
//    std::cout << std::endl;

  }

  void construct_csr_from_sorted_edgelist(edge_list_type& edge_list)
  {
    std::cout << "Constructing adj-list from sorted edge list" << std::endl;
    size_t cnt_edges = 0;
    vertex_type last_src_vertex = 0;
    for (auto edge_itr = edge_list.begin(), edge_itr_end = edge_list.end(); edge_itr != edge_itr_end; ++edge_itr) {
      const vertex_type src = edge_itr->first;
      const vertex_type dst = edge_itr->second;
      assert(last_src_vertex <= src); /// must be sorted
      assert(src < m_num_vertices);   /// must be small than the max vertex ID

      ++m_index_array[src];
      m_adj_list[cnt_edges++] = dst;

      if (last_src_vertex < src)
        last_src_vertex = src;
    }
    std::cout << "last_src_vertex + 1 == m_num_vertices :" << last_src_vertex + 1 << " <= " << m_num_vertices << std::endl;
    assert(last_src_vertex + 1 <= m_num_vertices);
    std::cout << "cnt_edges == m_num_edges :" << cnt_edges << " == " << m_num_edges << std::endl;
    assert(cnt_edges == m_num_edges);
    graph_trv_utilities::print_time();

    std::cout << "Modificate index-array" << std::endl;
    cumulate_array(m_index_array, m_num_vertices);
    right_shift_aray(m_index_array, m_num_vertices);
    m_index_array[0] = 0;
    std::cout << "m_index_array[m_num_vertices] == m_num_edges :" << m_index_array[m_num_vertices] << " == " << m_num_edges << std::endl;
    // assert(m_index_array[m_num_vertices] == m_num_edges);
    graph_trv_utilities::print_time();
  }

  void cumulate_array(index_type* array, size_t length)
  {
    for (size_t i = 0; i < length; ++i) {
      array[i+1] += array[i];
    }
  }

  void right_shift_aray(index_type* array, size_t length)
  {
    for (size_t i = length; i > 0; --i) {
      array[i] = array[i-1];
    }
  }

  void load_info(std::string& prefix)
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

  void load_graph(std::string& prefix)
  {
    std::cout << "Load graph from files" << std::endl;

    std::ifstream idx_file(prefix + "_idx", std::ifstream::binary);
    std::ifstream adj_file(prefix + "_adj", std::ifstream::binary);

    assert(m_num_vertices > 0 && m_index_array);
    assert(m_num_edges    > 0 && m_adj_list);

    idx_file.read(reinterpret_cast<char*>(m_index_array), sizeof(index_type)  * (m_num_vertices+1));
    adj_file.read(reinterpret_cast<char*>(m_adj_list),    sizeof(vertex_type) * m_num_edges);
    // std::cout << "error: only " << idx_file.gcount() << " could be read" << sizeof(index_type)  * (m_num_vertices+1);
    // assert(idx_file);
    // assert(adj_file);

    std::cout << "Loading graph from files done" << std::endl;
    std::cout << m_index_array[m_num_vertices] << std::endl;

    adj_file.close();
    idx_file.close();
  }

  size_t m_num_vertices;
  size_t m_num_edges;
  index_type*  m_index_array;
  vertex_type* m_adj_list;
};

template <typename edge_list_type, typename index_t, typename vertex_t>
void dump_unvisited_edges(csr_graph<edge_list_type, index_t, vertex_t>& graph,
                          graph_trv_info::trv_inf<vertex_t>& inf,
                          std::string fname,
                          std::string fname_wk)
{
  using csr_graph_type = csr_graph<edge_list_type, index_t, vertex_t>;
  std::vector<std::pair<vertex_t, vertex_t>> unvisited_edges;
  size_t visited_edges = inf.count_total_visited_edges();

  #if 0
  const index_type* index_array = graph.index_array();
  const vertex_type* adj_list = graph.adj_list();
  for (vertex_type src = 0; src < graph.num_vertices(); ++src) {
    for (index_type i = index_array[src]; i < index_array[src+1]; ++i) {
      if (!inf.global_edge_is_visited[i]) {
        vertex_type dst = adj_list[i];
        unvisited_edges.push_back(std::make_pair(src, dst));
      }
    }
  }
  std::cout << "#unvisited_edges:\t" << unvisited_edges.size() << std::endl;
  assert(graph.num_edges() == inf.count_total_visited_edges() + unvisited_edges.size());

  #else
  {
    std::cout << "Dump unvisited edges" << std::endl;
    // assert(graph.num_edges() == visited_edges + unvisited_edges.size());

    std::ofstream ofs_wk(fname_wk);
    size_t count = 0;
    size_t fname_count = 0;
    const typename csr_graph_type::index_type* index_array = graph.index_array();
    const typename csr_graph_type::vertex_type* adj_list = graph.adj_list();
    for (vertex_t src = 0; src < graph.num_vertices(); ++src) {
      for (index_t i = index_array[src]; i < index_array[src+1]; ++i) {
        if (!inf.global_edge_is_visited[i]) {
         // if (count % (1ULL << 28) == 0) {
         //    std::stringstream fname;
         //    fname << fname_prefix << std::setfill('0') << std::setw(5) << fname_count;
         //    std::cout << "open : " << fname.str().c_str() << std::endl;
         //    if (ofs_wk.is_open()) ofs_wk.close();
         //    ofs_wk.open(fname.str().c_str());
         //    ++fname_count;
         //  }
          vertex_t dst = adj_list[i];
          ofs_wk << src << " " << dst << "\n";
          ++count;
        }
      }
    }
    assert(graph.num_edges() == visited_edges + count);
    unvisited_edges.reserve(count);

    std::cout << "#unvisited_edges:\t" << unvisited_edges.size() << std::endl;
  }

  std::cout << "!!! Delete the graph !!!" << std::endl;
  graph.~csr_graph();
  std::cout << "!!! Delete the trv_inf !!!" << std::endl;
  inf.~trv_inf();

  {
    std::cout << "Reload unvisited edges" << std::endl;
    std::ifstream ifs_wk(fname_wk);
    std::string line;
    while (std::getline(ifs_wk, line)) {
      std::stringstream ssline(line);
      uint64_t src_vtx, dst_vtx;
      ssline >> src_vtx >> dst_vtx;
      unvisited_edges.push_back(std::make_pair(src_vtx, dst_vtx));
    }
  }

  {
    std::cout << "Truncate work file" << std::endl;
    std::ofstream ofs_wk(fname_wk, std::ofstream::trunc);
  }

  #endif

  unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "shuffle vector. seed: " << seed2 << std::endl;
  std::shuffle(unvisited_edges.begin(), unvisited_edges.end(), std::default_random_engine(seed2));

  std::cout << "Writing unvisited edges" << std::endl;
  std::ofstream ofs_crawling(fname);
  assert(ofs_crawling.good());
  for (auto itr = unvisited_edges.begin(), end = unvisited_edges.end();
      itr != end;
      ++itr)
  {
    ofs_crawling << itr->first << " " << itr->second << "\n";
  }

  ofs_crawling.close();
}

}
#endif // CSR_GRAPH_HPP

