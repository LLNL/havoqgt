#ifndef dist_dynamic_graphstore_HPP
#define dist_dynamic_graphstore_HPP

#include <havoqgt/mpi.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>

template <typename graphstore_type>
class dist_dynamic_graphstore
{
 public:
  class vertex_locator;
  using vertex_iterator = typename graphstore_type::vertex_iterator;

  using vertex_property_data_type = unsigned char;
  using edge_property_data_type   = unsigned char;


  dist_dynamic_graphstore(graphstore_type* graphstore) :
  m_graphstore(graphstore)
  {}

  /// Converts a vertex_locator to the vertex label
  uint64_t locator_to_label(vertex_locator locator) const
  {
    return locator.m_id;
  }

  /// Converts a vertex label to a vertex_locator
  vertex_locator label_to_locator(uint64_t label) const
  {
    return vertex_locator(label);
  }

  inline bool insert_edge(const vertex_locator& src, const vertex_locator& trg, const edge_property_data_type& weight)
  {
    return m_graphstore->insert_edge(src.m_id, trg.m_id, weight);
  }

  inline vertex_property_data_type& vertex_property_data(const vertex_locator& vertex)
  {
    return m_graphstore->vertex_property_data(vertex.m_id);
  }

  /// Returns the degree of a vertex
  inline uint64_t degree(vertex_locator locator) const
  {
    m_graphstore->degree(locator.m_id);
  }

  inline edge_property_data_type& edge_property_data(const vertex_locator& src, const vertex_locator& trg)
  {
    return vertex_property_data.edge_property_data(src.m_id, trg.m_id);
  }

  uint32_t master(const vertex_locator& locator) const {
    assert(false);
    return 0;
  }

 private:
  graphstore_type* m_graphstore;
};



template <typename graphstore_type>
class dist_dynamic_graphstore<graphstore_type>::vertex_locator {
  using vertex_type = uint64_t;

  friend dist_dynamic_graphstore;

 public:
  vertex_locator() {
    m_is_bcast     = 0;
    m_is_intercept = 0;
    m_id           = std::numeric_limits<vertex_type>::max();
  }

  explicit vertex_locator(vertex_type src) {
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
  static inline bool lesser_hash_priority(const vertex_locator& a,
                                          const vertex_locator& b) {
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

  bool is_equal(const vertex_locator x) const {
    return m_is_bcast     == x.m_is_bcast
        && m_is_intercept == x.m_is_intercept
        && m_id           == x.m_id;
  }

  uint32_t get_bcast() const {  return m_is_bcast;  }
  void set_bcast(uint32_t bcast) {  m_is_bcast = bcast;  }

  bool is_intercept() const {  return m_is_intercept == 1;  }
  void set_intercept(bool intercept) {  m_is_intercept = intercept;  }


  friend bool operator == (const vertex_locator& x,
                           const vertex_locator& y) { return x.is_equal(y); }
  friend bool operator < (const vertex_locator& x,
                          const vertex_locator& y) {
    return x.m_id < y.m_id;
  }
  friend bool operator >= (const vertex_locator& x,
                           const vertex_locator& y) { return !(x < y); }
  friend bool operator != (const vertex_locator& x,
                           const vertex_locator& y) { return !(x.is_equal(y)); }

 private:
  unsigned int m_is_bcast     : 1;
  unsigned int m_is_intercept : 1;
  uint64_t     m_id           : 61;
} __attribute__ ((packed)) ;



template <typename graphstore_type>
class dg_visitor {
 public:
  enum visit_t { BAD, INI, ADD, CHK, DEL };
  using vertex_locator = typename graphstore_type::vertex_locator;

  // Default constructor.
  dg_visitor() :
      vertex(), caller(), vis_type(BAD) {  }

  // Baseline constructor.
  explicit dg_visitor(vertex_locator _vertex) :
      vertex(_vertex), caller(_vertex), vis_type(INI) {  }

  // Who I am, who notified me, and what type of visit it is.
  dg_visitor(vertex_locator _vertex, vertex_locator _caller, visit_t _vistype) :
      vertex(_vertex), caller(_caller), vis_type(_vistype){  }


  bool pre_visit() const {
    // Perform an action dependant on the visit type.
    switch(vis_type) {
      case ADD:
        graph_ref()->insert_edge(vertex, caller, 0);
        /*
        if (vertex.id() <= 50) {
          std::cout << havoqgt::havoqgt_env()->world_comm().rank() << ":" << vertex.id() << "," << caller.id() << " ";
        }
        */
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

    return false;
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
  visit_t        vis_type;
} __attribute__((packed));


#endif // dist_dynamic_graphstore_HPP

