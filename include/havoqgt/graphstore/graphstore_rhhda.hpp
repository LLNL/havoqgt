/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef GRAPHSTORE_RHHDA_HPP
#define GRAPHSTORE_RHHDA_HPP

#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>

#include <havoqgt/graphstore/rhh/rhh_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_container.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>


namespace graphstore {

template <typename vertex_id_type,
          typename vertex_meta_data_type,
          typename edge_weight_type,
          typename segment_manager_type,
          size_t middle_high_degree_threshold>
class graphstore_rhhda
{
 private:
  using size_type                      = size_t;
  using low_degree_table_value_type    = utility::packed_tuple<vertex_meta_data_type, vertex_id_type, edge_weight_type>;
  using low_degree_table_type          = rhh_container<vertex_id_type, low_degree_table_value_type, size_type, segment_manager_type>;
  using mid_high_edge_chunk_type       = rhh_container<vertex_id_type, edge_weight_type, size_type, segment_manager_type>;
  using mid_high_src_vertex_value_type = utility::packed_pair<vertex_meta_data_type, mid_high_edge_chunk_type*>;
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
  bool insert_edge(const vertex_id_type& src, const vertex_id_type& trg, const edge_weight_type& weight)
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
        low_degree_table_value_type value(vertex_meta_data_type(), trg, weight);
        rhh::insert(&m_low_degree_table, src, value);
      } else {

        /// --- move the elements from low table to high-mid table --- ///
        mid_high_edge_chunk_type* adj_list = mid_high_edge_chunk_type::allocate(middle_high_degree_threshold);
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
        low_degree_table_value_type value(vertex_meta_data_type(), trg, weight);
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

  inline bool insert_vertex(const vertex_id_type& vertex, const vertex_meta_data_type& meta_data)
  {
    auto itr_single = m_low_degree_table->find(vertex);
    if (itr_single.is_end()) {
      auto itr = m_mid_high_degree_table->find(vertex);
      if (itr.is_end()) {
        rhh::insert(m_low_degree_table,
                    vertex,
                    low_degree_table_value_type(meta_data, vertex_id_type(), edge_weight_type()));
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
  size_type erase_edge(vertex_id_type& src, vertex_id_type& trg)
  {
    size_type count = 0;
    for (auto itr = m_low_degree_table->find(src); !itr.is_end(); ++itr) {
      if ((*itr).second == trg) {
        m_low_degree_table->erase(itr);
        ++count;
      }
    }

    if (count > 0) {
      /// rhh::shrink_to_fit(&m_low_degree_table);
      return count;
    }

    auto itr_matrix = m_mid_high_degree_table->find(src);
    /// has source vertex ?
    if (itr_matrix.is_end()) return false;
    mid_high_edge_chunk_type* adj_list = itr_matrix->second;

    for (auto itr = adj_list->find(trg); !itr.is_end(); ++itr) {
      adj_list->erase(itr);
      ++count;
    }

    if (count > 0) {
      if (adj_list->size() < middle_high_degree_threshold) {
        const vertex_meta_data_type& meta_data = itr_matrix->first;
        for (auto itr = adj_list->begin(); !itr.is_end(); ++itr) {
          low_degree_table_value_type value(meta_data, itr->key, itr->value);
          rhh::insert(&m_low_degree_table, src, value);
          adj_list->erase(itr);
        }
        mid_high_edge_chunk_type::deallocate(adj_list);
        m_mid_high_degree_table->erase(itr_matrix);
        /// rhh::shrink_to_fit(&m_mid_high_degree_table);
      } else {
        /// rhh::shrink_to_fit(&adj_list);
        itr_matrix->second = adj_list;
      }
    }

    return count;
  }

//  inline size_type erase_vertex(vertex_id_type& vertex)
//  {
//    return m_mid_high_degree_table->erase(vertex);
//  }

  inline vertex_meta_data_type& vertex_meta_data(const vertex_id_type& vertex)
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

  typename low_degree_table_type::value_iterator find_low_edge (vertex_id_type& src_vrt)
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename low_degree_table_type::const_value_iterator find_low_edge (vertex_id_type& src_vrt) const
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename mid_high_edge_chunk_type::whole_iterator find_mid_high_edge (vertex_id_type& src_vrt)
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->begin();
  }

  typename mid_high_edge_chunk_type::const_whole_iterator find_mid_high_edge (vertex_id_type& src_vrt) const
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
/// \brief The graphstore_rhhda<vertex_id_type, vertex_meta_data_type, edge_weight_type, 0> class
///   partial speciallization class when middle_high_degree_threshold is 0
template <typename vertex_id_type,
          typename vertex_meta_data_type,
          typename edge_weight_type,
          typename segment_manager_type>
class graphstore_rhhda <vertex_id_type, vertex_meta_data_type, edge_weight_type, segment_manager_type, 0>
{
 private:
  using size_type = size_t;
  using low_degree_table_value_type    = utility::packed_tuple<vertex_meta_data_type, vertex_id_type, edge_weight_type>;
  using low_degree_table_type          = rhh_container<vertex_id_type, low_degree_table_value_type, size_type, segment_manager_type>;
  using mid_high_edge_chunk_type       = rhh_container<vertex_id_type, edge_weight_type, size_type, segment_manager_type>;
  using mid_high_src_vertex_value_type = utility::packed_pair<vertex_meta_data_type, mid_high_edge_chunk_type*>;
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

  void shrink_to_fit_low_table()
  {
    rhh::shrink_to_fit(&m_low_degree_table);
  }

  void shrink_to_fit_mid_high_table()
  {
    rhh::shrink_to_fit(&m_mid_high_degree_table);
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
  bool insert_edge(vertex_id_type& src, vertex_id_type& trg, edge_weight_type& weight)
  {


    auto itr_src = m_mid_high_degree_table->find(src);
    if (itr_src.is_end()) {
      /// --- new vertex --- ///
      mid_high_edge_chunk_type* adj_list = mid_high_edge_chunk_type::allocate(2);
      rhh::insert(&adj_list, trg, weight);
      mid_high_src_vertex_value_type value(vertex_meta_data_type(), adj_list);
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


  inline bool insert_vertex(vertex_id_type& vertex, vertex_meta_data_type& meta_data)
  {
    auto itr = m_mid_high_degree_table->find(vertex);
    if (itr.is_end()) {
      mid_high_src_vertex_value_type value(meta_data, nullptr);
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
  size_type erase_edge(vertex_id_type& src, vertex_id_type& trg)
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

  size_type erase_vertex(vertex_id_type& vertex)
  {
    return m_mid_high_degree_table->erase(vertex);
  }

  vertex_meta_data_type& vertex_meta_data(const vertex_id_type& vertex)
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
    return m_low_degree_table->end();
  }

  typename mid_high_degree_table_type::whole_iterator begin_mid_high_edges()
  {
    return m_mid_high_degree_table->begin();
  }

  typename low_degree_table_type::value_iterator find_low_edge (vertex_id_type& src_vrt)
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename low_degree_table_type::const_value_iterator find_low_edge (vertex_id_type& src_vrt) const
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename mid_high_edge_chunk_type::whole_iterator find_mid_high_edge (vertex_id_type& src_vrt)
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->begin();
  }

  typename mid_high_edge_chunk_type::const_whole_iterator find_mid_high_edge (vertex_id_type& src_vrt) const
  {
    const auto itr_matrix = m_mid_high_degree_table->find(src_vrt);
    const mid_high_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->cbegin();
  }


  ///
  /// \brief print_status
  ///   Note: this function accesses entier data of the rhhda containers to compute statuses
  ///         thus, this function would affect pagecache and cause I/Os
  void print_status()
  {

    std::cout << "<high-middle degree table>: "
              << "\n size, capacity, rate: " << m_mid_high_degree_table->size() << ", " << m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth()
              << ", " << (double)(m_mid_high_degree_table->size()) / (m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth())
              << "\n chaine depth : " << m_mid_high_degree_table->depth()
              << "\n average probedistance: " << m_mid_high_degree_table->load_factor()
              << "\n capacity*element_size(GB): "
                << (double)m_mid_high_degree_table->capacity() * m_mid_high_degree_table->depth() * mid_high_degree_table_type::kElementSize  / (1ULL<<30) << std::endl;
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
        ++histgram_ave_prbdist[adj_list->load_factor()];

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

