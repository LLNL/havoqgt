/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef GRAPHSTORE_RHHDA_HPP
#define GRAPHSTORE_RHHDA_HPP

#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_container.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>


namespace graphstore {

template <typename vertex_id_type,
          typename vertex_property_type,
          typename edge_property_type,
          typename segment_manager_type,
          size_t middle_high_degree_threshold = 2>
class graphstore_rhhda
{
 private:
  using size_type                      = size_t;
  using low_degree_table_value_type    = utility::packed_tuple<vertex_property_type, vertex_id_type, edge_property_type>;
  using low_degree_table_type          = rhh_container<vertex_id_type, low_degree_table_value_type, size_type, segment_manager_type>;
  using mid_high_edge_chunk_type       = rhh_container<vertex_id_type, edge_property_type, size_type, segment_manager_type>;
  using mid_high_src_vertex_value_type = utility::packed_pair<vertex_property_type, mid_high_edge_chunk_type*>;
  using mid_high_degree_table_type     = rhh_container<vertex_id_type, mid_high_src_vertex_value_type, size_type, segment_manager_type>;
  using graphstore_rhhda_selftype      = graphstore_rhhda<vertex_id_type, vertex_property_type, edge_property_type,
                                                          segment_manager_type, middle_high_degree_threshold>;


 public:

  friend class adjlist_inputiterator;
  using adjlist_edge_type = std::pair<vertex_id_type, edge_property_type>;

  class adjlist_inputiterator : public std::iterator<std::input_iterator_tag, adjlist_edge_type>
  {
    friend class graphstore_rhhda;
    using edge_iterator_selftype           = adjlist_inputiterator;
    using graphstore_type                  = graphstore_rhhda_selftype;
    using low_deg_edge_iteratir_type       = typename low_degree_table_type::value_iterator;
    using mid_high_deg_edge_iteratir_type  = typename mid_high_edge_chunk_type::whole_iterator;


   public:

    adjlist_inputiterator() = delete;

    adjlist_inputiterator(graphstore_type* gstore, const vertex_id_type& src_vrt) :
      m_low_itr(gstore->find_low_edge(src_vrt)),
      m_mh_itr(gstore->find_mid_high_edge(src_vrt)),
      m_current_edge()
    {
      if (!m_low_itr.is_end()) {
        m_current_edge = adjlist_edge_type(m_low_itr->second, m_low_itr->third);
      } else if (!m_mh_itr.is_end()) {
        m_current_edge = adjlist_edge_type(m_mh_itr->key, m_mh_itr->value);
      }
    }

    adjlist_inputiterator(low_deg_edge_iteratir_type low_itr, mid_high_deg_edge_iteratir_type mh_itr) :
      m_low_itr(low_itr),
      m_mh_itr(mh_itr),
      m_current_edge()
    {
      if (!m_low_itr.is_end()) {
        m_current_edge = adjlist_edge_type(m_low_itr->second, m_low_itr->third);
      } else if (!m_mh_itr.is_end()) {
        m_current_edge = adjlist_edge_type(m_mh_itr->key, m_mh_itr->value);
      }
    }


    void swap(edge_iterator_selftype &other) noexcept
    {
      using std::swap;
      swap(m_low_itr, other.m_low_itr);
      swap(m_mh_itr, other.m_mh_itr);
      swap(m_current_edge, other.m_current_edge);
    }

    edge_iterator_selftype &operator++ () // Pre-increment
    {
      find_next_value();
      return *this;
    }

    edge_iterator_selftype operator++ (int) // Post-increment
    {
      edge_iterator_selftype tmp(*this);
      find_next_value();
      return tmp;
    }

    // two-way comparison: v.begin() == v.cbegin() and vice versa
    bool operator == (const edge_iterator_selftype &rhs) const
    {
      return is_equal(rhs);
    }

    bool operator != (const edge_iterator_selftype &rhs) const
    {
      return !is_equal(rhs);
    }

    const adjlist_edge_type& operator* () const
    {
      return m_current_edge;
    }

    const adjlist_edge_type* operator-> () const
    {
      return &(m_current_edge);
    }

   private:

    inline bool is_equal(const adjlist_inputiterator &rhs) const
    {
      return (m_low_itr == rhs.m_low_itr) && (m_mh_itr == rhs.m_mh_itr);
    }

    inline void find_next_value()
    {
      if (!m_low_itr.is_end()) {
        ++m_low_itr;
        if (!m_low_itr.is_end()) {
          m_current_edge = adjlist_edge_type(m_low_itr->second, m_low_itr->third);
        }
      } else if (!m_mh_itr.is_end()) {
        ++m_mh_itr;
        if (!m_mh_itr.is_end()) {
          m_current_edge = adjlist_edge_type(m_mh_itr->key, m_mh_itr->value);
        }
      }
    }

    low_deg_edge_iteratir_type m_low_itr;
    mid_high_deg_edge_iteratir_type m_mh_itr;
    adjlist_edge_type m_current_edge;
  };


  explicit graphstore_rhhda(segment_manager_type* segment_manager) {
    // -- init allocator -- //
    rhh::init_allocator<typename low_degree_table_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename mid_high_edge_chunk_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename mid_high_degree_table_type::allocator, segment_manager_type>(segment_manager);

    m_low_degree_table = low_degree_table_type::allocate(2);
    m_mid_high_degree_table = mid_high_degree_table_type::allocate(2);

    std::cout << "Element size: \n"
              << " low_degree_table = " << low_degree_table_type::kElementSize << "\n"
              << " mid_high_edge_chunk = " << mid_high_edge_chunk_type::kElementSize << "\n"
              << " mid_high_degree_table = " << mid_high_degree_table_type::kElementSize << std::endl;
    std::cout << "middle_high_degree_threshold = " << middle_high_degree_threshold << std::endl;

#ifdef RHH_DETAILED_ANALYSYS
    graphstore::utility::rhh_log_holder::instance().init();
#endif
  }

  ~graphstore_rhhda() {
    clear();
    low_degree_table_type::deallocate(m_low_degree_table);
    mid_high_degree_table_type::deallocate(m_mid_high_degree_table);
    rhh::destroy_allocator<typename low_degree_table_type::allocator>();
    rhh::destroy_allocator<typename mid_high_edge_chunk_type::allocator>();
    rhh::destroy_allocator<typename mid_high_degree_table_type::allocator>();
  }

  adjlist_inputiterator adjacencylist(const vertex_id_type& srt_vrtx)
  {
    return adjlist_inputiterator(this, srt_vrtx);
  }

  static adjlist_inputiterator adjacencylist_end(const vertex_id_type& srt_vrtx)
  {
    return adjlist_inputiterator(low_degree_table_type::find_end(), mid_high_edge_chunk_type::end());
  }


  void opt()
  {
    rehash_low_table();
  }

  void shrink_to_fit_low_table()
  {
    rhh::shrink_to_fit(&m_low_degree_table);
  }

  void shrink_to_fit_mid_high_table()
  {
    rhh::shrink_to_fit(&m_mid_high_degree_table);
  }

  void rehash_low_table()
  {
    std::cout << "rehash_low_table()" << std::endl;
    m_low_degree_table->rehash();
  }

  void rehash_mid_high_table()
  {
    std::cout << "rehash_mid_high_table()" << std::endl;
    m_mid_high_degree_table->rehash();
  }

  ///
  /// \brief insert_edge
  ///   inert a edge uniquely
  /// \param src
  /// \param trg
  /// \param weight
  /// \return
  ///   true: if inserted
  ///   false: if a duplicated edge is found
  ///
  bool insert_edge(const vertex_id_type& src, const vertex_id_type& trg, const edge_property_type& weight)
  {

    /// -- count the degree of the source vertex in the low degree table -- ///
    size_type count_in_single = 0;
    for (auto itr_single = m_low_degree_table->find(src); !itr_single.is_end(); ++itr_single) {
      if ((*itr_single).second == trg) {
        return false;
      }
      ++count_in_single;
    }

    if (count_in_single > 0) { /// -- the low table has the source vertex -- ///
      if (count_in_single + 1 < middle_high_degree_threshold) {
        /// --- insert into the low table --- ///
        low_degree_table_value_type value(vertex_property_type(), trg, weight);
        rhh::insert(&m_low_degree_table, src, value);
      } else {

        /// --- move the elements from low table to high-mid table --- ///
        mid_high_edge_chunk_type* adj_list = mid_high_edge_chunk_type::allocate(middle_high_degree_threshold * 2);
        auto itr_single = m_low_degree_table->find(src);
        mid_high_src_vertex_value_type value((*itr_single).first, nullptr);
        for (; !itr_single.is_end(); ++itr_single) {
          rhh::insert(&adj_list, (*itr_single).second, (*itr_single).third);
          m_low_degree_table->erase(itr_single);
        }
        rhh::insert(&adj_list, trg, weight);
        value.second = adj_list;
        rhh::insert(&m_mid_high_degree_table, src, value);
        /// rhh::shrink_to_fit(&m_low_degree_table);
      }

    } else {
      auto itr_src = m_mid_high_degree_table->find(src);
      if (itr_src.is_end()) {
        /// --- since the high-mid table dosen't have the vertex, insert into the low table (new vertex) --- ///
        low_degree_table_value_type value(vertex_property_type(), trg, weight);
        rhh::insert(&m_low_degree_table, src, value);
      } else {
        /// --- the high-mid table has source vertex --- ///
        mid_high_edge_chunk_type* adj_list = itr_src->second;
        auto itr_trg = adj_list->find(trg);

        if (itr_trg.is_end()) {
          /// --- insert the edge --- ///
          rhh::insert(&adj_list, trg, weight);
          itr_src->second = adj_list;
        } else {
          /// --- if the same edge is found, do nothing --- ///
          return false;
        }
      }

    }

    /// TODO: insert the target vertex into the source vertex list
    return true;
  }

  inline bool insert_vertex(const vertex_id_type& vertex, const vertex_property_type& property_data)
  {
    auto itr_single = m_low_degree_table->find(vertex);
    if (itr_single.is_end()) {
      auto itr = m_mid_high_degree_table->find(vertex);
      if (itr.is_end()) {
        rhh::insert(m_low_degree_table,
                    vertex,
                    low_degree_table_value_type(property_data, vertex_id_type(), edge_property_type()));
        return true;
      }
    }
    return false;
  }


  ///
  /// \brief erase_edge
  ///         erase edges. this function can delete duplicated edges.
  /// \param src
  /// \param trg
  /// \return
  ///         the number of edges erased
  size_type erase_edge(const vertex_id_type& src, const vertex_id_type& trg)
  {
    size_type count = 0;
    for (auto itr = m_low_degree_table->find(src); !itr.is_end(); ++itr) {
      if ((*itr).second == trg) {
        m_low_degree_table->erase(itr);
        ++count;
      }
    }
    if (count > 0) {
      rhh::shrink_to_fit(&m_low_degree_table, 2.0);
      return count;
    }

    auto itr_matrix = m_mid_high_degree_table->find(src);
    /// has source vertex ?
    if (itr_matrix.is_end()) return 0;
    mid_high_edge_chunk_type* adj_list = itr_matrix->second;
    for (auto itr = adj_list->find(trg); !itr.is_end(); ++itr) {
      adj_list->erase(itr);
      ++count;
    }

    if (count > 0) {
      if (adj_list->size() < middle_high_degree_threshold) {
        const vertex_property_type& property_data = itr_matrix->first;
        for (auto itr = adj_list->begin(); !itr.is_end(); ++itr) {
          low_degree_table_value_type value(property_data, itr->key, itr->value);
          rhh::insert(&m_low_degree_table, src, value);
          adj_list->erase(itr);
        }
        mid_high_edge_chunk_type::deallocate(adj_list);
        m_mid_high_degree_table->erase(itr_matrix);
        rhh::shrink_to_fit(&m_mid_high_degree_table);
      } else {
        rhh::shrink_to_fit(&adj_list);
        itr_matrix->second = adj_list;
      }
    }

    return count;
  }

//  inline size_type erase_vertex(const vertex_id_type& vertex)
//  {
//    return m_mid_high_degree_table->erase(vertex);
//  }

  inline vertex_property_type& vertex_property_data(const vertex_id_type& vertex)
  {
    auto itr = m_low_degree_table->find(vertex);
    if (!itr.is_end()) {
      return itr->first;
    }
    auto itr_matrix = m_mid_high_degree_table->find(vertex);
    return itr_matrix->first;
  }

  void clear()
  {
    // for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr)
    for (const auto& itr : *m_mid_high_degree_table) {
      mid_high_edge_chunk_type* const adj_list = itr.value.second;
      adj_list->clear();
    }
    m_low_degree_table->clear();
  }

  typename low_degree_table_type::whole_iterator begin_low_edges()
  {
    return m_low_degree_table->begin();
  }

  typename mid_high_degree_table_type::whole_iterator begin_mid_high_edges()
  {
    return m_mid_high_degree_table->begin();
  }

  typename low_degree_table_type::value_iterator find_low_edge (const vertex_id_type& src_vrt)
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename mid_high_edge_chunk_type::whole_iterator find_mid_high_edge (const vertex_id_type& src_vrt)
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    if (!itr_matrix.is_end()) {
      mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
      return adj_list->begin();
    } else {
      return mid_high_edge_chunk_type::end();
    }
  }


  ///
  /// \brief print_status
  ///   Note: this function accesses entier data of the rhhda containers to compute statuses
  ///         thus, this function would affect pagecache and cause I/Os
  void print_status(const int level) const
  {
    std::cout << "<low degree table>:"
              << "\n size, capacity, rate: " << m_low_degree_table->size()
              << ", " << m_low_degree_table->capacity() * m_low_degree_table->depth()
              << ", " << (double)(m_low_degree_table->size()) / (m_low_degree_table->capacity() * m_low_degree_table->depth())
              << "\n chaine depth: " << m_low_degree_table->depth()
              << "\n capacity*element_size(GB): "
              << (double)m_low_degree_table->capacity() * m_low_degree_table->depth() * low_degree_table_type::kElementSize / (1ULL<<30) << std::endl;

    std::cout << "<high-middle degree table>: "
              << "\n size, capacity, rate: " << m_mid_high_degree_table->size()
              << ", " << m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth()
              << ", " << (double)(m_mid_high_degree_table->size()) / (m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth())
              << "\n chaine depth : " << m_mid_high_degree_table->depth()
              << "\n capacity*element_size(GB): "
              << (double)m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth() * mid_high_degree_table_type::kElementSize  / (1ULL<<30) << std::endl;
    if (level == 0) return;

    std::cout << "<low degree table>:"
              << "\n average probedistance: " << m_low_degree_table->load_factor()
              << std::endl;
    {
      size_type histgram_prbdist[low_degree_table_type::property_program::kLongProbedistanceThreshold] = {0};
      std::cout << "probedistance: ";
      m_low_degree_table->histgram_load_factor(histgram_prbdist);
      for (int i = 0; i < utility::array_length(histgram_prbdist); ++i) {
        std::cout << histgram_prbdist[i] << " ";
      }
      std::cout << std::endl;
    }

    std::cout << "<high-middle degree table>:"
              << "\n average probedistance: " << m_mid_high_degree_table->load_factor()
              << std::endl;
    {
      size_type histgram_ave_prbdist[mid_high_degree_table_type::property_program::kLongProbedistanceThreshold] = {0};
      size_type histgram_cap[50] = {0};
      size_type histgram_size[50] = {0};
      size_type histgram_dept[50] = {0};
      size_type size_sum = 0;
      size_type capacity_sum = 0;
      for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
        auto adj_list = itr->value.second;

        /// --- size ---- ///
        size_sum += adj_list->size();
        const size_type sz_log2 = std::min(std::log2(adj_list->size()),
                                        static_cast<double>(utility::array_length(histgram_size) - 1));
        ++histgram_size[sz_log2];

        /// --- average probe distance (laod factor) ---- ///
        assert(adj_list->load_factor() < utility::array_length(histgram_ave_prbdist));
        ++histgram_ave_prbdist[static_cast<size_t>(adj_list->load_factor())];

        /// --- capacity --- ///
        capacity_sum += adj_list->capacity() * adj_list->depth();
        const size_type cap_log2 = std::min(std::log2(adj_list->capacity()),
                                        static_cast<double>(utility::array_length(histgram_cap) - 1));
        ++histgram_cap[cap_log2];

        /// --- depth --- ///
        const size_type depth = std::min(adj_list->depth(),
                                      utility::array_length(histgram_dept) - 1);
        ++histgram_dept[depth];
      }

      std::cout << "<high-middle edge chunks>: "
                << "\n size: " << size_sum
                << "\n capacity: " << capacity_sum
                << "\n rate: " << (double)(size_sum) / capacity_sum
                << "\n capacity*element_size(GB): " << (double)capacity_sum * mid_high_edge_chunk_type::kElementSize  / (1ULL<<30) << std::endl;

      std::cout << "average probedistance: ";
      for (int i = 0; i < utility::array_length(histgram_ave_prbdist); ++i) {
        std::cout << histgram_ave_prbdist[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "capacity (log2): ";
      for (int i = 0; i < utility::array_length(histgram_cap); ++i) {
        std::cout << histgram_cap[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "size (log2): ";
      for (int i = 0; i < utility::array_length(histgram_size); ++i) {
        std::cout << histgram_size[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "depth: ";
      for (int i = 1; i < utility::array_length(histgram_dept); ++i) {
        std::cout << histgram_dept[i] << " ";
      }
      std::cout << std::endl;

    }
  }

  void fprint_all_elements(std::ofstream& of)
  {
    for (auto itr = m_low_degree_table->begin(); !itr.is_end(); ++itr) {
      of << itr->key << " " << (itr->value).second << "\n";
    }

    for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        of << itr->key << " " << itr2->key << "\n";
      }
    }
  }

  void print_all_elements_low()
  {
    std::cout << "------------------" << std::endl;
    for (auto itr = m_low_degree_table->begin(); !itr.is_end(); ++itr) {
      std::cout << itr->key << " " << (itr->value).second << "\n";
    }
  }

  void print_all_elements_mh()
  {
    std::cout << "------------------" << std::endl;
    for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        std::cout << itr->key << " " << itr2->key << "\n";
      }
    }
  }


 private:
  low_degree_table_type* m_low_degree_table;
  mid_high_degree_table_type* m_mid_high_degree_table;

};


///
/// \brief The graphstore_rhhda<vertex_id_type, vertex_property_type, edge_property_type, 1> class
///   partial speciallization class when middle_high_degree_threshold is 1
template <typename vertex_id_type,
          typename vertex_property_type,
          typename edge_property_type,
          typename segment_manager_type>
class graphstore_rhhda <vertex_id_type, vertex_property_type, edge_property_type, segment_manager_type, 1>
{
 private:
  using size_type = size_t;
  using low_degree_table_value_type    = utility::packed_tuple<vertex_property_type, vertex_id_type, edge_property_type>;
  using low_degree_table_type          = rhh_container<vertex_id_type, low_degree_table_value_type, size_type, segment_manager_type>;
  using mid_high_edge_chunk_type       = rhh_container<vertex_id_type, edge_property_type, size_type, segment_manager_type>;
  using mid_high_src_vertex_value_type = utility::packed_pair<vertex_property_type, mid_high_edge_chunk_type*>;
  using mid_high_degree_table_type     = rhh_container<vertex_id_type, mid_high_src_vertex_value_type, size_type, segment_manager_type>;


 public:

  explicit graphstore_rhhda(segment_manager_type* segment_manager) {
    // -- init allocator -- //
    rhh::init_allocator<typename low_degree_table_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename mid_high_edge_chunk_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename mid_high_degree_table_type::allocator, segment_manager_type>(segment_manager);

    m_low_degree_table = low_degree_table_type::allocate(2);
    m_mid_high_degree_table = mid_high_degree_table_type::allocate(2);

    std::cout << "Element size: \n"
              << " low_degree_table = " << low_degree_table_type::kElementSize << "\n"
              << " mid_high_edge_chunk = " << mid_high_edge_chunk_type::kElementSize << "\n"
              << " mid_high_degree_table = " << mid_high_degree_table_type::kElementSize << std::endl;
    std::cout << "middle_high_degree_threshold (using only m_h table) = " << 0 << std::endl;
  }

  ~graphstore_rhhda() {
    clear();
    low_degree_table_type::deallocate(m_low_degree_table);
    mid_high_degree_table_type::deallocate(m_mid_high_degree_table);
    rhh::destroy_allocator<typename mid_high_edge_chunk_type::allocator>();
    rhh::destroy_allocator<typename mid_high_degree_table_type::allocator>();
  }

  void opt()
  {
  }

  void shrink_to_fit_low_table()
  {
    rhh::shrink_to_fit(&m_low_degree_table);
  }

  void shrink_to_fit_mid_high_table()
  {
    rhh::shrink_to_fit(&m_mid_high_degree_table);
  }

  void rehash_low_table()
  {
    std::cout << "rehash_low_table()" << std::endl;
    m_low_degree_table->rehash();
  }

  void rehash_mid_high_table()
  {
    std::cout << "rehash_mid_high_table()" << std::endl;
    m_mid_high_degree_table->rehash();
  }

  ///
  /// \brief insert_edge
  ///   inert a edge uniquely
  /// \param src
  /// \param trg
  /// \param weight
  /// \return
  ///   true: if inserted
  ///   false: if a duplicated edge is found
  ///
  bool insert_edge(const vertex_id_type& src, const vertex_id_type& trg, const edge_property_type& weight)
  {


    auto itr_src = m_mid_high_degree_table->find(src);
    if (itr_src.is_end()) {
      /// --- new vertex --- ///
      mid_high_edge_chunk_type* adj_list = mid_high_edge_chunk_type::allocate(2);
      rhh::insert(&adj_list, trg, weight);
      mid_high_src_vertex_value_type value(vertex_property_type(), adj_list);
      rhh::insert(&m_mid_high_degree_table, src, value);
    } else {
      /// --- high-mid table has source vertex --- ///
      mid_high_edge_chunk_type* adj_list = itr_src->second;
      auto itr_trg = adj_list->find(trg);

      if (itr_trg.is_end()) {
        /// --- insert the edge --- ///
        rhh::insert(&adj_list, trg, weight);
        itr_src->second = adj_list;
      } else {
        /// --- if the same edge is found, do nothing --- ///
        return false;
      }

    }

    return true;
  }


  inline bool insert_vertex(const vertex_id_type& vertex, const vertex_property_type& property_data)
  {
    auto itr = m_mid_high_degree_table->find(vertex);
    if (itr.is_end()) {
      mid_high_src_vertex_value_type value(property_data, nullptr);
      rhh::insert(&m_mid_high_degree_table, vertex, value);
      return true;
    }
    return false;
  }


  ///
  /// \brief erase_edge
  ///         erase edges. this function can delete duplicated edges.
  /// \param src
  /// \param trg
  /// \return
  ///         the number of edges erased
  size_type erase_edge(const vertex_id_type& src, const vertex_id_type& trg)
  {
    size_type count = 0;

    auto itr_matrix = m_mid_high_degree_table->find(src);
    /// has source vertex ?
    if (itr_matrix.is_end()) return false;
    mid_high_edge_chunk_type* adj_list = itr_matrix->second;

    for (auto itr = adj_list->find(trg); !itr.is_end(); ++itr) {
      adj_list->erase(itr);
      ++count;
    }

    if (count > 0) {
      if (adj_list->size() == 0) {
        mid_high_edge_chunk_type::deallocate(adj_list);
        m_mid_high_degree_table->erase(itr_matrix);
        /// rhh::shrink_to_fit(&m_mid_high_degree_table);
      }
    }

    return count;
  }

//  size_type erase_vertex(const vertex_id_type& vertex)
//  {
//    return m_mid_high_degree_table->erase(vertex);
//  }

  vertex_property_type& vertex_property_data(const vertex_id_type& vertex)
  {
    auto itr_matrix = m_mid_high_degree_table->find(vertex);
    return itr_matrix->first;
  }

  void clear()
  {
    // for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr)
    for (const auto& itr : *m_mid_high_degree_table) {
      mid_high_edge_chunk_type* const adj_list = itr.value.second;
      adj_list->clear();
    }
  }

  typename low_degree_table_type::whole_iterator begin_low_edges()
  {
    return low_degree_table_type::end();
  }

  typename mid_high_degree_table_type::whole_iterator begin_mid_high_edges()
  {
    return m_mid_high_degree_table->begin();
  }

  typename low_degree_table_type::value_iterator find_low_edge (const vertex_id_type& src_vrt)
  {
    return low_degree_table_type::find_end();
  }

  typename low_degree_table_type::const_value_iterator find_low_edge (const vertex_id_type& src_vrt) const
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename mid_high_edge_chunk_type::whole_iterator find_mid_high_edge (const vertex_id_type& src_vrt)
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->begin();
  }

  typename mid_high_edge_chunk_type::const_whole_iterator find_mid_high_edge (const vertex_id_type& src_vrt) const
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    const mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->cbegin();
  }


  ///
  /// \brief print_status
  ///   Note: this function accesses entier data of the rhhda containers to compute statuses
  ///         thus, this function would affect pagecache and cause I/Os
  void print_status(const int level) const
  {

    std::cout << "<high-middle degree table>: "
              << "\n size, capacity, rate: " << m_mid_high_degree_table->size()
              << ", " << m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth()
              << ", " << (double)(m_mid_high_degree_table->size()) / (m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth())
              << "\n chaine depth : " << m_mid_high_degree_table->depth()
              << "\n capacity*element_size(GB): "
              << (double)m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth() * mid_high_degree_table_type::kElementSize  / (1ULL<<30) << std::endl;
    if (level == 0) return;

    std::cout << "<high-middle degree table>:"
              << "\n average probedistance: " << m_mid_high_degree_table->load_factor()
              << std::endl;
    {
      size_type histgram_ave_prbdist[mid_high_degree_table_type::property_program::kLongProbedistanceThreshold] = {0};
      size_type histgram_cap[50] = {0};
      size_type histgram_size[50] = {0};
      size_type histgram_dept[50] = {0};
      size_type size_sum = 0;
      size_type capacity_sum = 0;
      for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
        auto adj_list = itr->value.second;

        /// --- size ---- ///
        size_sum += adj_list->size();
        const size_type sz_log2 = std::min(std::log2(adj_list->size()),
                                        static_cast<double>(utility::array_length(histgram_size) - 1));
        ++histgram_size[sz_log2];

        /// --- average probe distance (laod factor) ---- ///
        assert(adj_list->load_factor() < utility::array_length(histgram_ave_prbdist));
        ++histgram_ave_prbdist[static_cast<size_t>(adj_list->load_factor())];

        /// --- capacity --- ///
        capacity_sum += adj_list->capacity() * adj_list->depth();
        const size_type cap_log2 = std::min(std::log2(adj_list->capacity()),
                                        static_cast<double>(utility::array_length(histgram_cap) - 1));
        ++histgram_cap[cap_log2];

        /// --- depth --- ///
        const size_type depth = std::min(adj_list->depth(),
                                      utility::array_length(histgram_dept) - 1);
        ++histgram_dept[depth];
      }

      std::cout << "<high-middle edge chunks>: "
                << "\n size: " << size_sum
                << "\n capacity: " << capacity_sum
                << "\n rate: " << (double)(size_sum) / capacity_sum
                << "\n capacity*element_size(GB): " << (double)capacity_sum * mid_high_edge_chunk_type::kElementSize  / (1ULL<<30) << std::endl;

      std::cout << "average probedistance: ";
      for (int i = 0; i < utility::array_length(histgram_ave_prbdist); ++i) {
        std::cout << histgram_ave_prbdist[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "capacity (log2): ";
      for (int i = 0; i < utility::array_length(histgram_cap); ++i) {
        std::cout << histgram_cap[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "size (log2): ";
      for (int i = 0; i < utility::array_length(histgram_size); ++i) {
        std::cout << histgram_size[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "depth: ";
      for (int i = 1; i < utility::array_length(histgram_dept); ++i) {
        std::cout << histgram_dept[i] << " ";
      }
      std::cout << std::endl;

    }
  }

  void fprint_all_elements(std::ofstream& of)
  {
    for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        of << itr->key << " " << itr2->key << "\n";
      }
    }
  }

  void print_all_elements_low()
  {
    std::cout << "------------------" << std::endl;
  }

  void print_all_elements_mh()
  {
    std::cout << "------------------" << std::endl;
    for (auto itr = m_mid_high_degree_table->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        std::cout << itr->key << " " << itr2->key << "\n";
      }
    }
  }


 private:
  low_degree_table_type* m_low_degree_table;
  mid_high_degree_table_type* m_mid_high_degree_table;
};

}
#endif // GRAPHSTORE_RHHDA_HPP

