/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef GRAPHSTORE_RHHDA_HPP
#define GRAPHSTORE_RHHDA_HPP

#include <tuple>

#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.h>
#include <havoqgt/graphstore/rhh/rhh_container.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

namespace graphstore {

template <typename vertex_id_type, typename vertex_meta_data_type, typename edge_weight_type, size_t middle_high_degree_threshold>
class graphstore_rhhda
{
private:
  using size_type = size_t;
  using low_degree_table_value_type    = utility::packed_tuple<vertex_meta_data_type, vertex_id_type, edge_weight_type>;
  using low_degree_table_type          = rhh_container_base<vertex_id_type, low_degree_table_value_type, size_type>;
  using high_mid_edge_chunk_type       = rhh_container_base<vertex_id_type, edge_weight_type, size_type>;
  using high_mid_src_vertex_value_type = utility::packed_pair<vertex_meta_data_type, high_mid_edge_chunk_type*>;
  using high_mid_degree_table_type     = rhh_container_base<vertex_id_type, high_mid_src_vertex_value_type, size_type>;
  using segment_manager_type           = rhh::segment_manager_t;


public:

  explicit graphstore_rhhda(segment_manager_type* segment_manager) {
    // -- init allocator -- //
    rhh::init_allocator<typename low_degree_table_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename high_mid_edge_chunk_type::allocator, segment_manager_type>(segment_manager);
    rhh::init_allocator<typename high_mid_degree_table_type::allocator, segment_manager_type>(segment_manager);

    m_low_degree_table = low_degree_table_type::allocate(2);
    m_high_mid_degree_table = high_mid_degree_table_type::allocate(2);

    std::cout << "Element size: \n"
              << " low_degree_table " << low_degree_table_type::kElementSize << "\n"
              << " high_mid_edge_chunk " << high_mid_edge_chunk_type::kElementSize << "\n"
              << " high_mid_degree_table " << high_mid_degree_table_type::kElementSize << std::endl;
  }

  ~graphstore_rhhda() {
    clear();
    low_degree_table_type::deallocate(m_low_degree_table);
    high_mid_degree_table_type::deallocate(m_high_mid_degree_table);
    rhh::destroy_allocator<typename low_degree_table_type::allocator>();
    rhh::destroy_allocator<typename high_mid_edge_chunk_type::allocator>();
    rhh::destroy_allocator<typename high_mid_degree_table_type::allocator>();
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

    // count degree of the source vertex in low degree table
    size_t count_in_single = 0;
    for (auto itr_single = m_low_degree_table->find(src); !itr_single.is_end(); ++itr_single) {
      if ((*itr_single).second == trg) {
        return false;
      }
      ++count_in_single;
    }

    if (count_in_single > 0) {
      // insert into a low table or move to middle-high one
      if (count_in_single < middle_high_degree_threshold - 1) {
        low_degree_table_value_type value(vertex_meta_data_type(), trg, weight);
        rhh_container_utility::insert(&m_low_degree_table, src, value);
      } else {
        high_mid_edge_chunk_type* adj_list = high_mid_edge_chunk_type::allocate(middle_high_degree_threshold);
        auto itr_single2 = m_low_degree_table->find(src);
        high_mid_src_vertex_value_type value((*itr_single2).first, nullptr);
        for (; !itr_single2.is_end(); ++itr_single2) {
          rhh_container_utility::insert(&adj_list, (*itr_single2).second, (*itr_single2).third);
          m_low_degree_table->erase(itr_single2);
        }
        rhh_container_utility::insert(&adj_list, trg, weight);
        value.second = adj_list;
        rhh_container_utility::insert(&m_high_mid_degree_table, src, value);
        rhh_container_utility::shrink_to_fit(&m_low_degree_table);
      }

    } else {
      auto itr_src = m_high_mid_degree_table->find(src);
      if (itr_src.is_end()) {
        low_degree_table_value_type value(vertex_meta_data_type(), trg, weight);
        rhh_container_utility::insert(&m_low_degree_table, src, value);
      } else {
        /// has source vertex
        high_mid_edge_chunk_type* adj_list = itr_src->second;
        auto itr_trg = adj_list->find(trg);

        /// if same edge is found do nothing
        if (itr_trg.is_end()) {
          /// actually insert the edge (insert the target vertex into the adjacency list)
          rhh_container_utility::insert(&adj_list, trg, weight);
          itr_src->second = adj_list;
        } else {
          return false;
        }
      }

    }

    /// TODO: insert the target vertex into the source vertex list
EDGE_INSERTED:
    return true;
  }

  inline bool insert_vertex(vertex_id_type& vertex, vertex_meta_data_type& meta_data)
  {
    auto itr_single = m_low_degree_table->find(vertex);
    if (itr_single.is_end()) {
      auto itr = m_high_mid_degree_table->find(vertex);
      if (itr.is_end()) {
        rhh_container_utility::insert(m_low_degree_table, vertex, low_degree_table_value_type(meta_data, vertex_id_type(), edge_weight_type()));
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
  size_t erase_edge(vertex_id_type& src, vertex_id_type& trg)
  {
    size_t count = 0;
    for (auto itr = m_low_degree_table->find(src); !itr.is_end(); ++itr) {
      if ((*itr).second == trg) {
        m_low_degree_table->erase(itr);
        ++count;
      }
    }

    if (count > 0) {
      rhh_container_utility::shrink_to_fit(&m_low_degree_table);
      return count;
    }

    auto itr_matrix = m_high_mid_degree_table->find(src);
    /// has source vertex ?
    if (itr_matrix.is_end()) return false;
    high_mid_edge_chunk_type* adj_list = itr_matrix->second;

    for (auto itr = adj_list->find(trg); !itr.is_end(); ++itr) {
      adj_list->erase(itr);
      ++count;
    }

    if (count > 0) {
      if (adj_list->size() < middle_high_degree_threshold) {
        const vertex_meta_data_type& meta_data = itr_matrix->first;
        for (auto itr = adj_list->begin(); !itr.is_end(); ++itr) {
          low_degree_table_value_type value(meta_data, itr->key, itr->value);
          rhh_container_utility::insert(&m_low_degree_table, src, value);
          adj_list->erase(itr);
        }
        high_mid_edge_chunk_type::deallocate(adj_list);
        m_high_mid_degree_table->erase(itr_matrix);
        rhh_container_utility::shrink_to_fit(&m_high_mid_degree_table);
      } else {
        rhh_container_utility::shrink_to_fit(&adj_list);
        itr_matrix->second = adj_list;
      }
    }

    return count;
  }

  size_t erase_vertex(vertex_id_type& vertex)
  {
    return m_high_mid_degree_table->erase(vertex);
  }

  void clear()
  {
    for (auto itr = m_high_mid_degree_table->begin(); !itr.is_end(); ++itr) {
      high_mid_edge_chunk_type* adj_list = itr->value.second;
      adj_list->clear();
    }
    m_low_degree_table->clear();
  }


  typename low_degree_table_type::value_iterator find_low_edge (vertex_id_type& src_vrt)
  {
    return m_low_degree_table->find(src_vrt);
  }

  typename high_mid_edge_chunk_type::whole_iterator find_mid_high_edge (vertex_id_type& src_vrt)
  {
    const auto itr_matrix = m_high_mid_degree_table->find(src_vrt);
    high_mid_edge_chunk_type* const adj_list = itr_matrix->second;
    return adj_list->begin();
  }

  ///
  /// \brief print_status
  ///   Note: this function accesses entier data of the rhhda containers to compute statuses
  ///         thus, this function would affect pagecache and cause I/Os
  void print_status()
  {
    std::cout << "<low degree table>:"
              << "\n size, capacity, rate: " << m_low_degree_table->size() << ", " << m_low_degree_table->capacity() * m_low_degree_table->depth()
              << ", " << (double)(m_low_degree_table->size()) / (m_low_degree_table->capacity() * m_low_degree_table->depth())
              << "\n chaine depth: " << m_low_degree_table->depth()
              << "\n average probedistance: " << m_low_degree_table->load_factor()
              << "\n capacity*element_size(GB): " << (double)m_low_degree_table->capacity() * m_low_degree_table->depth() * low_degree_table_type::kElementSize / (1ULL<<30) << std::endl;

    std::cout << "<high-middle degree table>: "
              << "\n size, capacity, rate: " << m_high_mid_degree_table->size() << ", " << m_high_mid_degree_table->capacity() * m_high_mid_degree_table->depth()
              << ", " << (double)(m_high_mid_degree_table->size()) / (m_high_mid_degree_table->capacity() * m_high_mid_degree_table->depth())
              << "\n chaine depth : " << m_high_mid_degree_table->depth()
              << "\n average probedistance: " << m_high_mid_degree_table->load_factor()
              << "\n capacity*element_size(GB): " << (double)m_high_mid_degree_table->capacity() * m_high_mid_degree_table->depth() * high_mid_degree_table_type::kElementSize  / (1ULL<<30) << std::endl;
    {
      size_t histgram_load_factor[high_mid_degree_table_type::property_program::kLongProbedistanceThreshold] = {0};
      size_t histgram_cap_log2[50] = {0};
      size_t histgram_dept[30] = {0};
      size_t size_sum = 0;
      size_t capacity_sum = 0;
      for (auto itr = m_high_mid_degree_table->begin(); !itr.is_end(); ++itr) {
        auto adj_list = itr->value.second;

        size_sum += adj_list->size();

        assert(adj_list->load_factor() < utility::array_length(histgram_load_factor));
        ++histgram_load_factor[adj_list->load_factor()];

        capacity_sum += adj_list->capacity() * adj_list->depth();
        size_t cap_log2 = std::log2l(adj_list->capacity() * adj_list->depth());
        if (cap_log2 >= utility::array_length(histgram_cap_log2))
          cap_log2 = utility::array_length(histgram_cap_log2) - 1;
        ++histgram_cap_log2[cap_log2];

        size_t depth = adj_list->depth();
        if (depth >= utility::array_length(histgram_dept))
          depth = utility::array_length(histgram_dept);
        ++histgram_dept[depth];
      }

      std::cout << "<high-middle edge chunks>: "
                << "\n size: " << size_sum
                << "\n capacity: " << capacity_sum
                << "\n rate: " << (double)(size_sum) / capacity_sum
                << "\n capacity*element_size(GB): " << (double)capacity_sum * high_mid_edge_chunk_type::kElementSize  / (1ULL<<30) << std::endl;

      std::cout << "average probedistance: ";
      for (int i = 0; i < utility::array_length(histgram_load_factor); ++i) {
        std::cout << histgram_load_factor[i] << " ";
      }
      std::cout << std::endl;

      std::cout << "capacity (log2): ";
      for (int i = 0; i < utility::array_length(histgram_cap_log2); ++i) {
        std::cout << histgram_cap_log2[i] << " ";
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

    for (auto itr = m_high_mid_degree_table->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        of << itr->key << " " << itr2->key << "\n";
      }
    }
  }

 private:
  low_degree_table_type* m_low_degree_table;
  high_mid_degree_table_type* m_high_mid_degree_table;

};

}
#endif // GRAPHSTORE_RHHDA_HPP

