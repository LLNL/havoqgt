/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHH_CONTAINER_HPP_INCLUDED
#define RHH_CONTAINER_HPP_INCLUDED

#include <limits>
#include <cassert>      // assert
#include <iterator>     // iterator
#include <utility>      // swap

#include <havoqgt/graphstore/graphstore_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_element_base.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

namespace graphstore {

template<typename _key_type,
         typename _value_type,
         typename _size_type,
         typename _segment_manager_type,
         typename _key_hash_func = rhh::key_hash_func_64bit_to_64bit<_key_type, _size_type>,
         typename _property_program = rhh::rhh_property_program_base<unsigned char>>
class rhh_container {
 public:
  using property_program     = _property_program;
  using key_type             = _key_type;
  using key_hash_func        = _key_hash_func;
  using value_type           = _value_type;
  using size_type            = _size_type;
  using segment_manager_type = _segment_manager_type;
  using property_type        = typename property_program::property_type;
  using probedistance_type   = typename property_program::probedistance_type;
  using rhh_contatiner_selftype = rhh_container<key_type, value_type, size_type,
                                                segment_manager_type,
                                                key_hash_func, property_program>;

  /// TODO: specializetion for no value case
  #pragma pack(1)
  typedef struct element
  {
    property_type property;
    key_type      key;
    value_type    value;
  } packed_element;
  using element_type = packed_element;

  enum : size_type {
    kKeyNotFound   = std::numeric_limits<size_type>::max(),
    kElementSize   = sizeof(element_type),
#if RHH_DETAILED_ANALYSYS
    kTableBaseSize = sizeof(size_type) + sizeof(size_type) + sizeof(void*) + sizeof(size_t) * 3
#else
    kTableBaseSize = sizeof(size_type) + sizeof(size_type) + sizeof(void*)
#endif
  };

  using allocator = graphstore::rhh::allocator_holder_sglt<segment_manager_type,
                                                           kElementSize,
                                                           kTableBaseSize>;

  class whole_iterator;
  class value_iterator;


#if RHH_DETAILED_ANALYSYS
  size_t total_probed_distance;
  size_t max_probed_distance;
  size_t total_num_insert_into_body_called;
#endif

  /// explicitly delete to prevent unexpected behaivers
  rhh_container()  = delete;
  ~rhh_container() = delete;
  rhh_container(const rhh_contatiner_selftype& src) = delete;
  rhh_container& operator=(const rhh_container&) = delete;
  rhh_container(const rhh_container&&) = delete;
  rhh_container& operator=(rhh_container&& old_obj) = delete;


  /// --- Capacity ---- ///
  inline size_t size() const
  {
    return m_num_elems;
  }

  inline size_t capacity() const
  {
    return m_capacity;
  }

  inline size_t table_mem_size() const
  {
    return kElementSize * capacity() + kTableBaseSize;
  }

  /// --- Lookup --- ///
  inline value_iterator find(const key_type& key)
  {
    return value_iterator(this, key);
  }

  static inline value_iterator find_end()
  {
    return value_iterator::end();
  }

  inline whole_iterator begin()
  {
    return whole_iterator(this);
  }

  static inline whole_iterator end()
  {
    return whole_iterator::end();
  }

  /// ---- Modifiers ---- ///
  inline bool insert(key_type&& key, value_type&& value, key_type& key_long_prbdst, value_type& val_long_prbdst)
  {
    const bool is_success = insert_into_body(std::move(key), std::move(value), key_long_prbdst, val_long_prbdst);
    m_num_elems += is_success;
    return is_success;
  }

  void clear()
  {
    const size_type cap = m_capacity;
    for (size_type i = 0; i < cap; ++i) {
      erase_element_at(i);
    }
    m_num_elems = 0;

    if (m_next) {
      m_next->clear();
    }
  }

  /// TODO: move to non-member function
  inline void erase(const whole_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_pos);
    --m_num_elems;
  }

  /// TODO: move to non-member function ?
  inline void erase(const value_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_pos);
    --m_num_elems;
  }

  inline void rehash()
  {
    rehash_elements();
  }

  /// --- Hash policy --- ///
  /// \brief load_factor
  /// \return
  ///     average probe distance
  double load_factor() const
  {
    const size_type total = sum_probedistane();

    if (m_num_elems > 0)
      return static_cast<double>(total) / static_cast<double>(m_num_elems);
    else
      return 0.0;
  }

  inline size_type depth() const
  {
    if (!m_next) return 1;
    return m_next->depth() + 1;
  }


  /// --- allocator & deallocator --- ///
  static rhh_contatiner_selftype* allocate(const size_t capacity)
  {
    rhh_contatiner_selftype* rhh = reinterpret_cast<rhh_contatiner_selftype*>(allocator::instance().allocate(capacity));
    rhh->m_num_elems = 0;
    rhh->m_capacity = capacity; /// XXX: if capacity is not power of 2, this is not actual capacity
    rhh->m_next = nullptr;

    for (size_type i = 0; i < capacity; ++i) {
      property_program::empty(rhh->m_body[i].property);
    }
    return rhh;
  }

  static void deallocate(rhh_contatiner_selftype* rhh)
  {
    while (rhh != nullptr) {
      // std::cout << rhh->m_num_elems << " " << rhh->m_capacity << " " << rhh->m_next << std::endl;
      rhh_contatiner_selftype* next_rhh = rhh->m_next;
      allocator::instance().deallocate(reinterpret_cast<void*>(rhh), rhh->m_capacity);
      rhh = next_rhh;
    }
  }

  inline rhh_contatiner_selftype* chained_rhh() const
  {
    return m_next;
  }

  inline void assign_to_chained_rhh(rhh_contatiner_selftype* next)
  {
    m_num_elems = next->m_num_elems;
    m_capacity = next->m_capacity;
    m_next = next;
  }

  ///
  /// \brief make_with_source_rhh
  ///         allocate new rhh with new_capacity and move all elements from
  ///         a source_rhh with handling long probe distance case
  /// \param source_rhh
  /// \param new_capacity
  /// \return
  ///         pointer to the new rhh
  static rhh_contatiner_selftype* make_with_source_rhh(rhh_contatiner_selftype* source_rhh, const size_type new_capacity)
  {
    assert(source_rhh->size() <= new_capacity);

    const size_type src_capacity = source_rhh->m_capacity;
    key_type wk_key;
    value_type wk_val;

    rhh_contatiner_selftype* new_rhh = allocate(new_capacity);

    while (source_rhh != nullptr) {
      for (size_type i = 0; i < src_capacity; ++i) {
        const property_type property = source_rhh->m_body[i].property;
        if (!property_program::is_empty(property) && !property_program::is_tombstone(property)) {
          bool is_successed = new_rhh->insert_into_body(std::move(source_rhh->m_body[i].key), std::move(source_rhh->m_body[i].value), wk_key, wk_val);
          if (!is_successed) {
#if RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE
            rhh::resize(&new_rhh, new_capacity * rhh::kCapacityGrowingFactor);
#else
            rhh::chain(&new_rhh);
#endif
            is_successed = new_rhh->insert_into_body(std::move(wk_key), std::move(wk_val), wk_key, wk_val);
            if (!is_successed) {
              rhh::chain(&new_rhh);
              assert(new_rhh->insert_into_body(std::move(wk_key), std::move(wk_val), wk_key, wk_val));
            }
          }
          ++new_rhh->m_num_elems;
        }
      }
      source_rhh = source_rhh->m_next;
    }

    return new_rhh;
  }


  /// --- debug --- ///
  void print_all_element()
  {
    for (size_type i = 0; i < m_capacity; ++i) {
      std::cout << "property " << m_body[i].property << ", key " << m_body[i].key << std::endl;
    }
    if (m_next != nullptr) {
      m_next->print_all_element();
    }
  }

  ///
  /// \brief histgram_load_factor
  /// \param histgram
  ///   An array hold load factors
  ///
  template <size_t size>
  void histgram_load_factor(size_type (&histgram)[size])
  {
    for (size_type i = 0; i < m_capacity; ++i) {
      if (!property_program::is_empty(m_body[i].property) &&
          !property_program::is_tombstone(m_body[i].property)) {
        const size_type d = property_program::extract_probedistance(m_body[i].property);
        if (d >= size) exit(1);
        ++histgram[d];
      }
    }
    if (m_next != nullptr) {
      m_next->histgram_load_factor(histgram);
    }
  }

#if RHH_DETAILED_ANALYSYS
  void init_detailed_analysis()
  {
    total_probed_distance = 0;
    max_probed_distance   = 0;
    total_probed_distance = 0;
  }

  void print_detailed_analysis()
  {
    std::cout << "total_probed_distance: " << total_probed_distance << std::endl;
    std::cout << "max_probed_distance: " << max_probed_distance << std::endl;
    std::cout << "total_num_insert_into_body_called: " << total_num_insert_into_body_called << std::endl;
    std::cout << "average probed distance: " << (double)total_probed_distance / total_num_insert_into_body_called << std::endl;
  }
#endif

 /// ---------------------------------------------------------- ///
 ///                         private
 /// ---------------------------------------------------------- ///
 private:

  inline void erase_element_at(size_type pos)
  {
    property_program::tombstone(m_body[pos].property);
    m_body[pos].key.~key_type();
    m_body[pos].value.~value_type();
  }

  inline size_type cal_desired_position(size_type hash) const
  {
    return hash & (m_capacity - 1);
  }

  static inline void internal_locate(const key_type& key,
                                     const rhh_contatiner_selftype** body_ptr,
                                     size_type& pos,
                                     probedistance_type& cur_prb_dist)
  {
    const size_type hash = key_hash_func::hash(key);
    pos = (*body_ptr)->cal_desired_position(hash);
    cur_prb_dist = 0;
    internal_locate_with_hint(key, body_ptr, pos, cur_prb_dist);
  }

  static void internal_locate_with_hint(const key_type& key,
                                        const rhh_contatiner_selftype** body_ptr,
                                        size_type& pos,
                                        probedistance_type& prb_dist)
  {
    const size_type mask = ((*body_ptr)->m_capacity - 1);

    while (true) {
      property_type exist_property = (*body_ptr)->m_body[pos].property;
      if (property_program::is_empty(exist_property)) {
        break;
      } else if (prb_dist > property_program::extract_probedistance(exist_property)) {
        break;
      } else if (!property_program::is_tombstone(exist_property) && key == (*body_ptr)->m_body[pos].key) {
        return;
      }
      pos = (pos+1) & mask;
      ++prb_dist;
    }

    *body_ptr = (*body_ptr)->m_next;
    if ((*body_ptr) != nullptr) {
      internal_locate(key, body_ptr, pos, prb_dist);
      return;
    }

    pos = kKeyNotFound;
  }

  inline void construct(const size_type pos, const probedistance_type dist,
                        key_type&& key, value_type&& value)
  {
    property_program::init_property(m_body[pos].property, dist);
    m_body[pos].key   = std::move(key);
    m_body[pos].value = std::move(value);
    if (dist >= m_capacity) {
      rehash_elements();
    }
  }


  ///
  /// \brief insert_into_body
  /// \param key
  /// \param value
  /// \param key_long_prbdst: be set a value if probe ditance exceed a threshhold
  /// \param val_long_prbdst: be set a value if probe ditance exceed a threshhold
  /// \return
  ///        true: successed insert
  ///        false: probe distance of a element exceeded long probe distance thleshold
  bool insert_into_body(key_type&& key, value_type&& value,
                        key_type& key_long_prbdst, value_type& val_long_prbdst)
  {
#if RHH_DETAILED_ANALYSYS
      ++total_num_insert_into_body_called;
    size_t num_probe = 0;
#endif

    probedistance_type prb_dist = 0;
    const size_type hash = key_hash_func::hash(key);
    size_type pos = cal_desired_position(hash);
    const size_type mask = (m_capacity - 1);


    while(true) {
      const property_type exist_property = m_body[pos].property;

      if(property_program::is_empty(exist_property))
      {
#if RHH_DETAILED_ANALYSYS
          max_probed_distance = std::max(max_probed_distance, num_probe);
#endif
        if (property_program::is_long_probedistance(prb_dist)) {
          key_long_prbdst = std::move(key);
          val_long_prbdst = std::move(value);
          return false;
        }
        construct(pos, prb_dist, std::move(key), std::move(value));
        return true;
      }

      /// If the existing elem has probed equal or less than new, then swap places with existing
      /// elem, and keep going to find another slot for that elem.
      if (property_program::extract_probedistance(exist_property) <= prb_dist)
      {
        if (property_program::is_long_probedistance(prb_dist)) {
          key_long_prbdst = std::move(key);
          val_long_prbdst = std::move(value);
#if RHH_DETAILED_ANALYSYS
          max_probed_distance = std::max(max_probed_distance, num_probe);
#endif
          return false;
        }
        if(property_program::is_tombstone(exist_property))
        {
          construct(pos, prb_dist, std::move(key), std::move(value));
#if RHH_DETAILED_ANALYSYS
          max_probed_distance = std::max(max_probed_distance, num_probe);
#endif
          return true;
        }
        m_body[pos].property = prb_dist;
        using std::swap;
        swap(key, m_body[pos].key);
        swap(value, m_body[pos].value);
        prb_dist = property_program::extract_probedistance(exist_property);
      }
      pos = (pos+1) & mask;
      ++prb_dist;

#if RHH_DETAILED_ANALYSYS
      ++total_probed_distance;
      ++num_probe;
#endif
    }


    return false;
  }


  /// --- Hash policy --- ///
  size_type sum_probedistane() const
  {
    size_type total = 0;
    for (size_type pos = 0; pos < m_capacity; ++pos) {
      if (!property_program::is_empty(m_body[pos].property) &&
          !property_program::is_tombstone(m_body[pos].property))
        total += property_program::extract_probedistance(m_body[pos].property);
    }

    if (m_next != nullptr) {
      total += m_next->sum_probedistane();
    }

    return total;
  }


  /// --- optimization --- ///
  void rehash_elements()
  {
    key_type wk_key;
    value_type wk_val;

    rhh_contatiner_selftype* const tmp_rhh = allocate(m_capacity);

    for (size_type i = 0; i < m_capacity; ++i) {
      const property_type property = m_body[i].property;
      if (!property_program::is_empty(property) && !property_program::is_tombstone(property)) {
        /// should check the return value ?
        tmp_rhh->insert_into_body(std::move(m_body[i].key), std::move(m_body[i].value), wk_key, wk_val);
        ++tmp_rhh->m_num_elems;
      }
    }
    std::memcpy(&(m_body[0]), &(tmp_rhh->m_body[0]), m_capacity * kElementSize);

    deallocate(tmp_rhh);
  }


  /// --- private valiable --- ///
  size_type m_num_elems; /// including chaining rhh tables
  size_type m_capacity;  /// this table's size
  rhh_contatiner_selftype* m_next;
  packed_element m_body[0];

};

namespace rhh {
template<typename key_type,
         typename value_type,
         typename size_type,
         typename segment_manager_type,
         typename key_hash_func,
         typename property_program>
bool shrink_to_fit(rhh_container<key_type, value_type, size_type,
                                 segment_manager_type,
                                 key_hash_func,
                                 property_program>** rhh,
                   const double lazy_factor = 1.0)
{
  using rhh_type = rhh_container<key_type, value_type, size_type,
                                 segment_manager_type,
                                 key_hash_func,
                                 property_program>;

  const typename rhh_type::size_type cur_size = (*rhh)->size();
  typename rhh_type::size_type cur_capacity = (*rhh)->capacity() /* * (*rhh)->depth()*/;

  if (cur_size > static_cast<double>(cur_capacity / kCapacityGrowingFactor / lazy_factor) * kFullCapacitFactor) {
    /// --- current capacity is fit to current size, do nothing --- ///
    return false;
  }

  typename rhh_type::size_type new_capacity = 1;
  while ( cur_size > static_cast<double>(new_capacity) * kFullCapacitFactor ) {
    new_capacity *= kCapacityGrowingFactor;
  }

  resize(rhh, new_capacity);

  return true;
}
} /// namespace rhh

} /// namespace graphstore

#include <havoqgt/graphstore/rhh/rhh_whole_iterator.hpp>
#include <havoqgt/graphstore/rhh/rhh_value_iterator.hpp>

#endif /// RHH_HPP
