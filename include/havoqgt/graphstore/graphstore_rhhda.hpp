/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef GRAPHSTORE_RHHDA_HPP
#define GRAPHSTORE_RHHDA_HPP

#include <tuple>

#include <havoqgt/graphstore/rhhda/rhhda_defs.hpp>
#include <havoqgt/graphstore/graphstore_common.hpp>
#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/rhhda/rhh_utilities.h>
#include <havoqgt/graphstore/rhhda/rhh_container.hpp>
#include <havoqgt/graphstore/rhhda/rhhda_allocator_holder.hpp>

namespace graphstore {

template <typename vertex_id_type, typename vertex_meta_data_type, typename edge_weight_type, size_t midle_high_degree_threshold>
class graphstore_rhhda
{
private:
  using size_type = size_t;
  using ld_singlelist_value_type   = std::tuple<vertex_meta_data_type, vertex_id_type, edge_weight_type>;
  using ld_singlelist_type         = rhh_container_base<vertex_id_type, ld_singlelist_value_type, size_type>;
  using hd_trg_vertex_adjlist_type = rhh_container_base<vertex_id_type, edge_weight_type, size_type>;
  using hd_src_vertex_value_type   = std::pair<vertex_meta_data_type, hd_trg_vertex_adjlist_type*>;
  using hd_adj_matrix_type         = rhh_container_base<vertex_id_type, hd_src_vertex_value_type, size_type>;
  using segment_manager_type       = rhhda::segment_manager_t;


public:

  explicit graphstore_rhhda(segment_manager_type* segment_manager) {
    // -- init allocator -- //
    rhhda::init_allocator<typename ld_singlelist_type::allocator, segment_manager_type>(segment_manager);
    rhhda::init_allocator<typename hd_trg_vertex_adjlist_type::allocator, segment_manager_type>(segment_manager);
    rhhda::init_allocator<typename hd_adj_matrix_type::allocator, segment_manager_type>(segment_manager);

    m_ld_singlelist = ld_singlelist_type::allocate(2);
    m_hd_adj_matrix = hd_adj_matrix_type::allocate(2);
  }

  ~graphstore_rhhda() {
    clear();
    ld_singlelist_type::deallocate(m_ld_singlelist);
    hd_adj_matrix_type::deallocate(m_hd_adj_matrix);
    rhhda::destroy_allocator<typename ld_singlelist_type::allocator>();
    rhhda::destroy_allocator<typename hd_trg_vertex_adjlist_type::allocator>();
    rhhda::destroy_allocator<typename hd_adj_matrix_type::allocator>();
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

    size_t count_in_single = 0;
    for (auto itr_single = m_ld_singlelist->find(src); !itr_single.is_end(); ++itr_single) {
      if (std::get<1>(*itr_single) == trg) {
        goto SKIP_EDGE_INSERTION;
      }
      ++count_in_single;
    }

    if (count_in_single > 0) {
      // insert into a low container or move to midle-high one
      if (count_in_single < midle_high_degree_threshold - 1) {
        ld_singlelist_value_type value(vertex_meta_data_type(), trg, weight);
        /// TODO: uniquly insertion
        rhh_container_utility::insert(&m_ld_singlelist, src, value);
      } else {
        hd_trg_vertex_adjlist_type* adj_list = hd_trg_vertex_adjlist_type::allocate(midle_high_degree_threshold);
        auto itr_single2 = m_ld_singlelist->find(src);
        hd_src_vertex_value_type value(std::get<0>(*itr_single2), nullptr);
        for (; !itr_single2.is_end(); ++itr_single2) {
          rhh_container_utility::insert(&adj_list, std::get<1>(*itr_single2), std::get<2>(*itr_single2));
          m_ld_singlelist->erase(itr_single2);
        }
        rhh_container_utility::insert(&adj_list, trg, weight);
        value.second = adj_list;
        rhh_container_utility::insert(&m_hd_adj_matrix, src, value);
      }

    } else {
      auto itr_src = m_hd_adj_matrix->find(src);
      if (itr_src.is_end()) {
        ld_singlelist_value_type value(vertex_meta_data_type(), trg, weight);
        rhh_container_utility::insert(&m_ld_singlelist, src, value);
      } else {
        /// has source vertex
        hd_trg_vertex_adjlist_type* adj_list = itr_src->second;
        auto itr_trg = adj_list->find(trg);

        /// if same edge is found, overwrite the edge weight
        if (!itr_trg.is_end()) {
          /// (*itr_trg) = weight;
          goto SKIP_EDGE_INSERTION;
        } else {
          /// actually insert the edge (insert the target vertex into the adjacency list)
          rhh_container_utility::insert(&adj_list, trg, weight);
          itr_src->second = adj_list;
        }
      }

    }

SKIP_EDGE_INSERTION:
    /// insert the target vertex into the source vertex list
    /// insert_vertex(trg, init_vertex_meta_data);

    return true;
  }

  inline bool insert_vertex(vertex_id_type& vertex, vertex_meta_data_type& meta_data)
  {
    auto itr_single = m_ld_singlelist->find(vertex);
    if (itr_single.is_end()) {
      auto itr = m_hd_adj_matrix->find(vertex);
      if (itr.is_end()) {
        rhh_container_utility::insert(m_ld_singlelist, vertex, ld_singlelist_value_type(meta_data, vertex_id_type(), edge_weight_type()));
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
    for (auto itr = m_ld_singlelist->find(src); !itr.is_end(); ++itr) {
      if (std::get<1>(*itr) == trg) {
        m_ld_singlelist->erase(itr);
        ++count;
      }
    }

    if (count > 0) {
      return count;
    }

    auto itr_matrix = m_hd_adj_matrix->find(src);
    /// has source vertex ?
    if (itr_matrix.is_end()) return false;
    hd_trg_vertex_adjlist_type* adj_list = itr_matrix->second;

    for (auto itr = adj_list->find(trg); !itr.is_end(); ++itr) {
      adj_list->erase(itr);
      ++count;
    }

    /// if the adj_list has no edges, deallocate the adj_list
    if (count > 0 && adj_list->size() < midle_high_degree_threshold) {
      const vertex_meta_data_type& meta_data = itr_matrix->first;
      for (auto itr = adj_list->begin(); !itr.is_end(); ++itr) {
        ld_singlelist_value_type value(meta_data, itr->key, itr->value);
        rhh_container_utility::insert(&m_ld_singlelist, src, value);
        adj_list->erase(itr);
      }
      hd_trg_vertex_adjlist_type::deallocate(adj_list);
      m_hd_adj_matrix->erase(itr_matrix);
    }

    return count;
  }

  size_t erase_vertex(vertex_id_type& vertex)
  {
    return m_hd_adj_matrix->erase(vertex);
  }

  void clear()
  {
    for (auto itr = m_hd_adj_matrix->begin(); !itr.is_end(); ++itr) {
      hd_trg_vertex_adjlist_type* adj_list = itr->value.second;
      adj_list->clear();
    }
    m_ld_singlelist->clear();
  }


  void print_status()
  {
    std::cout << "low degree table: "
              << m_ld_singlelist->size() << " / " << m_ld_singlelist->capacity() * m_ld_singlelist->depth()
              << " : " << m_ld_singlelist->load_factor() << std::endl;
    std::cout << "high-midle degree table: "
              << m_hd_adj_matrix->size() << " / " << m_hd_adj_matrix->capacity() * m_hd_adj_matrix->depth()
              << " : " << m_hd_adj_matrix->load_factor() << std::endl;
  }

  void fprint_all_elements(std::ofstream& of)
  {
    for (auto itr = m_ld_singlelist->begin(); !itr.is_end(); ++itr) {
      of << itr->key << " " << std::get<1>(itr->value) << std::endl;
    }

    for (auto itr = m_hd_adj_matrix->begin(); !itr.is_end(); ++itr) {
      auto adj_list = itr->value.second;
      for (auto itr2 = adj_list->begin(); !itr2.is_end(); ++itr2) {
        of << itr->key << "\t" << itr2->key << "\n";
      }
    }

  }

 private:
  ld_singlelist_type* m_ld_singlelist;
  hd_adj_matrix_type* m_hd_adj_matrix;

};

}
#endif // GRAPHSTORE_RHHDA_HPP

