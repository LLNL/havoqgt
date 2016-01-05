#ifndef dist_dynamic_graphstore_HPP
#define dist_dynamic_graphstore_HPP

#include <havoqgt/mpi.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

class no_delegate_vertex_locator;

namespace graphstore {

template <typename graphstore_type, typename vertex_locator_type = no_delegate_vertex_locator>
class dist_dynamic_graphstore
{
 public:
  using vertex_locator            = vertex_locator_type;

  using vertex_property_data_type = typename graphstore_type::vertex_property_data_type;
  using edge_property_data_type   = typename graphstore_type::edge_property_data_type;

  /// --- iterators --- ///
  class vertex_iterator;
  class adjacent_edge_iterator;

  dist_dynamic_graphstore(graphstore_type* graphstore) :
  m_graphstore(graphstore)
  {}

  /// Converts a vertex_locator to the vertex label
  uint64_t locator_to_label(vertex_locator locator) const
  {
    return locator.id();
  }

  /// Converts a vertex label to a vertex_locator
  vertex_locator label_to_locator(uint64_t label) const
  {
    return vertex_locator(label);
  }

  inline uint32_t master(const vertex_locator& locator) const
  {
    return locator.owner();
  }


  /// -------- Modifiers ------- ////
  /// insert a edge uniquely
  inline bool insert_edge(const vertex_locator& src, const vertex_locator& trg, const edge_property_data_type& weight)
  {
    return m_graphstore->insert_edge(src.id(), trg.id(), weight);
  }

  /// insert a edge
  inline bool insert_edge_dup(const vertex_locator& src, const vertex_locator& trg, const edge_property_data_type& weight)
  {
    return m_graphstore->insert_edge_dup(src.id(), trg.id(), weight);
  }

  /// erase a edge
  inline bool erase_edge(const vertex_locator& src, const vertex_locator& trg)
  {
    return m_graphstore->erase_edge(src.id(), trg.id());
  }

  /// erase edges
  inline bool erase_edge_dup(const vertex_locator& src, const vertex_locator& trg)
  {
    return m_graphstore->erase_edge_dup(src.id(), trg.id());
  }

  inline vertex_property_data_type& vertex_property_data(const vertex_locator& vertex)
  {
    return m_graphstore->vertex_property_data(vertex.id());
  }

  inline edge_property_data_type& edge_property_data(const vertex_locator& src, const vertex_locator& trg)
  {
    return m_graphstore->edge_property_data(src.id(), trg.id());
  }


  /// -------- Lookup -------- ///
  inline vertex_iterator vertices_begin()
  {
    return vertex_iterator(m_graphstore->vertices_begin());
  }

  inline vertex_iterator vertices_end()
  {
    return vertex_iterator(m_graphstore->vertices_end());
  }

  inline adjacent_edge_iterator adjacent_edge_begin(const vertex_locator& src)
  {
    return adjacent_edge_iterator(m_graphstore->adjacent_edge_begin(src.id()));
  }

  inline adjacent_edge_iterator adjacent_edge_end(const vertex_locator& src)
  {
    return adjacent_edge_iterator(m_graphstore->adjacent_edge_end(src.id()));
  }

  /// Returns the degree of a vertex
  inline uint64_t degree(vertex_locator locator) const
  {
    return m_graphstore->degree(locator.id());
  }


 private:
  graphstore_type* m_graphstore;
};

template <typename graphstore_type, typename vertex_locator_type>
class dist_dynamic_graphstore<graphstore_type, vertex_locator_type>::vertex_iterator
{
private:
 using self_type                 = dist_dynamic_graphstore<graphstore_type, vertex_locator_type>::vertex_iterator;
 using internal_vertex_iterator  = typename graphstore_type::vertex_iterator;

public:
 using vertex_locator            = vertex_locator_type;
 using vertex_property_data_type = typename graphstore_type::vertex_property_data_type;

 explicit vertex_iterator (internal_vertex_iterator&& iterator) :
   m_iterator(iterator)
 { }

 void swap(self_type &other) noexcept
 {
   using std::swap;
   swap(m_iterator, other.m_iterator);
 }

 self_type &operator++ () // Pre-increment
 {
   find_next_value();
   return *this;
 }

 self_type operator++ (int) // Post-increment
 {
   self_type tmp(*this);
   find_next_value();
   return tmp;
 }

 // two-way comparison: v.begin() == v.cbegin() and vice versa
 bool operator == (const self_type &rhs) const
 {
   return is_equal(rhs);
 }

 bool operator != (const self_type &rhs) const
 {
   return !is_equal(rhs);
 }

 vertex_locator_type source_vertex()
 {
   return vertex_locator_type(m_iterator.source_vertex());
 }

 vertex_property_data_type& property_data()
 {
   return m_iterator.property_data();
 }


private:

 inline bool is_equal(const self_type &rhs) const
 {
   return (m_iterator == rhs.m_iterator);
 }

 inline void find_next_value()
 {
   ++m_iterator;
 }

 internal_vertex_iterator m_iterator;
};

template <typename graphstore_type, typename vertex_locator_type>
class dist_dynamic_graphstore<graphstore_type, vertex_locator_type>::adjacent_edge_iterator
{
private:
 using self_type                        = dist_dynamic_graphstore<graphstore_type, vertex_locator_type>::adjacent_edge_iterator;
 using internal_adjacent_edge_iterator  = typename graphstore_type::adjacent_edge_iterator;

public:
 using vertex_locator                   = vertex_locator_type;
 using edge_property_data_type          = typename graphstore_type::edge_property_data_type;

 explicit adjacent_edge_iterator (internal_adjacent_edge_iterator&& iterator) :
   m_iterator(iterator)
 { }


 void swap(self_type &other) noexcept
 {
   using std::swap;
   swap(m_iterator, other.m_iterator);
 }

 self_type &operator++ () // Pre-increment
 {
   find_next_value();
   return *this;
 }

 self_type operator++ (int) // Post-increment
 {
   self_type tmp(*this);
   find_next_value();
   return tmp;
 }

 bool operator == (const self_type &rhs) const
 {
   return is_equal(rhs);
 }

 bool operator != (const self_type &rhs) const
 {
   return !is_equal(rhs);
 }

 vertex_locator target_vertex()
 {
   return vertex_locator_type(m_iterator.target_vertex());
 }

 const edge_property_data_type& property_data()
 {
   return m_iterator.property_data();
 }


private:

 inline bool is_equal(const self_type &rhs) const
 {
   return (m_iterator == rhs.m_iterator);
 }

 inline void find_next_value()
 {
   ++m_iterator;
 }

 internal_adjacent_edge_iterator m_iterator;
};

}


class no_delegate_vertex_locator {
  using vertex_type = uint64_t;
  using self_type   = no_delegate_vertex_locator;

 public:
  no_delegate_vertex_locator() {
    m_is_bcast     = 0;
    m_is_intercept = 0;
    m_id           = std::numeric_limits<vertex_type>::max();
  }

  explicit no_delegate_vertex_locator(vertex_type src) {
    m_is_bcast     = 0;
    m_is_intercept = 0;
    m_id           = src;
  }

  bool is_valid() const {
    uint64_t _m_id = std::numeric_limits<vertex_type>::max();
    return (m_id != _m_id);
  }

  // Computes the hash of a vertex based off of ID.
  inline uint64_t hash() const {
    std::hash<vertex_type> hasher;
    return hasher(id());
  }

  // With hash id, determine if vertex A has lesser priority over vertex B.
  static inline bool lesser_hash_priority(const self_type& a,
                                          const self_type& b) {
    size_t a_hash = a.hash();
    size_t b_hash = b.hash();
    if ((a_hash < b_hash) || (a_hash == b_hash && (a < b))) {
      return true;
    }
    return false;
  }

  bool is_delegate() const {  return false;  }

  uint64_t id() const {  return m_id;  }

  uint32_t owner() const {
    return hash() % havoqgt::havoqgt_env()->world_comm().size();
  }

  void set_dest(uint32_t dest) {
    std::cerr << "ERROR:  Tried to set dest in dynamic." << std::endl; exit(-1);
  }

  bool is_equal(const self_type x) const {
    return m_is_bcast     == x.m_is_bcast
        && m_is_intercept == x.m_is_intercept
        && m_id           == x.m_id;
  }

  uint32_t get_bcast() const {  return m_is_bcast;  }
  void set_bcast(uint32_t bcast) {  m_is_bcast = bcast;  }

  bool is_intercept() const {  return m_is_intercept == 1;  }
  void set_intercept(bool intercept) {  m_is_intercept = intercept;  }


  friend bool operator == (const self_type& x,
                           const self_type& y) { return x.is_equal(y); }
  friend bool operator < (const self_type& x,
                          const self_type& y) {
    return x.m_id < y.m_id;
  }
  friend bool operator >= (const self_type& x,
                           const self_type& y) { return !(x < y); }
  friend bool operator != (const self_type& x,
                           const self_type& y) { return !(x.is_equal(y)); }

 private:
  unsigned int m_is_bcast     : 1;
  unsigned int m_is_intercept : 1;
  uint64_t     m_id           : 61;
} __attribute__ ((packed)) ;


/// --- example of simple visitor class to construct a graph --- ///
template <typename graphstore_type>
class dg_visitor {
 public:
  enum visit_type { BAD, INI, ADD, CHK, DEL };
  using vertex_locator = typename graphstore_type::vertex_locator;

  // Default constructor.
  dg_visitor() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor.
  explicit dg_visitor(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), vis_type(INI) {  }

  // Who I am, who notified me, and what type of visit it is.
  dg_visitor(vertex_locator _vertex, vertex_locator _caller, visit_type _vistype) :
      vertex(_vertex), caller(_caller), vis_type(_vistype){  }


  bool pre_visit() const {
    // Perform an action dependant on the visit type.
    switch(vis_type) {
      case ADD:
        return true;
        break;
      default:
        std::cerr << "ERROR:  Bad visit type." << std::endl; exit(-1);
        break;
    }

    return false;
  }


  // TODO
  template<typename VisitorQueueHandle>
  bool visit(const graphstore_type& graph, VisitorQueueHandle vis_queue) const {
    switch (vis_type) {
      case ADD: {
        graph_ref()->insert_edge(vertex, caller, 0);
        return false;
      }
    }
  }

  // TODO
  friend inline bool operator > (const dg_visitor& v1,
                                 const dg_visitor& v2) {
    return false;
  }


  // Static graph reference.
  static void set_graph_ref(graphstore_type* _graph_ref) {
    graph_ref() = _graph_ref;
  }
  static graphstore_type*& graph_ref() {
    static graphstore_type* _graph_ref;
    return _graph_ref;
  }


  // Instance variables.
  vertex_locator vertex;
  vertex_locator caller;
  visit_type        vis_type;
} __attribute__((packed));

#endif // dist_dynamic_graphstore_HPP

