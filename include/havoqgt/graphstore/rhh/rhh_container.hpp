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

/// #define RHH_DETAILED_ANALYSYS

#define GROW_WHEN_LONG_PROBE_DISTANCE 0

namespace graphstore {
namespace rhh {

/// ---- Utility functions ---- ///
template <typename rhh_type>
inline bool has_key(rhh_type* const rhh, typename rhh_type::key_type& key)
{
  return (rhh->find(key) != rhh_type::kKeyNotFound);
}


template <typename rhh_type>
inline static void resize(rhh_type** rhh, const typename rhh_type::size_type new_capacity)
{
  rhh_type* new_rhh = rhh_type::make_with_source_rhh(*rhh, new_capacity);
  rhh_type::deallocate(*rhh);
  *rhh = new_rhh;
}

template <typename rhh_type>
inline static void grow(rhh_type** rhh)
{
#if 0
  if ((*rhh)->table_mem_size() > 1024) {
    chain(rhh);
  } else {
    resize(rhh, (*rhh)->capacity() * kCapacityGrowingFactor);
  }
#else
  resize(rhh, (*rhh)->capacity() * kCapacityGrowingFactor);
#endif
}

template <typename rhh_type>
inline static void chain(rhh_type** rhh)
{
  rhh_type* const new_rhh = rhh_type::allocate((*rhh)->capacity());
  new_rhh->assign_to_chained_rhh(*rhh);
  *rhh = new_rhh;
}


///
/// \brief insert
///   insert a element with checking the capacity and a probde distance.
///   if the probe distance exceed a threshold, allocate a chainged table
///   Note: this function causes copy of a key and a value !!
/// \param rhh
/// \param key
/// \param value
template <typename rhh_type, typename key_type, typename value_type>
void insert(rhh_type** rhh, key_type key, value_type value)
{
  /// --- check capacity --- ///
  if ((*rhh)->size() + 1 >= static_cast<size_t>(static_cast<double>((*rhh)->capacity()) * kFullCapacitFactor)) {
    grow(rhh);
  }

  /// --- Consider long probe distance --- ///
  while (!(*rhh)->insert(key, value, key, value)) {
#if GROW_WHEN_LONG_PROBE_DISTANCE
    grow(rhh);
#else
    chain(rhh);
#endif
  }
}

///
/// \brief insert_sizeup
///   insert a element with checking the capacity and a probde distance.
///   if the probe distance exceed a threshold, grow the table (not alocating chained table)
///   Note: this function causes copy of a key and a value !!
/// \param rhh
/// \param key
/// \param value
template <typename rhh_type, typename key_type, typename value_type>
void insert_sizeup(rhh_type** rhh, key_type key, value_type value)
{
  /// --- check capacity --- ///
  if ((*rhh)->size() + 1 >= static_cast<double>((*rhh)->capacity()) * kFullCapacitFactor) {
    grow(rhh);
  }

  /// --- Consider long probe distance --- ///
  while (!(*rhh)->insert(key, value, key, value)) {
    grow(rhh);
  }
}


template <typename rhh_type>
bool shrink_to_fit(rhh_type** rhh, const double lazy_factor = 1.0)
{
  const typename rhh_type::size_type cur_size = (*rhh)->size();
  typename rhh_type::size_type cur_capacity = (*rhh)->capacity() * (*rhh)->depth();

//  std::cout << cur_size << " " << cur_capacity << std::endl;
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
    kTableBaseSize = sizeof(size_type) + sizeof(size_type) + sizeof(void*)
  };

  using allocator = graphstore::rhh::allocator_holder_sglt<segment_manager_type,
                                                           kElementSize,
                                                           kTableBaseSize>;

  class whole_iterator;
  class value_iterator;


  /// explicitly delete due to prevent unexpected behaivers
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
    return value_iterator();
  }

  inline whole_iterator begin()
  {
    return whole_iterator(this);
  }

  static inline whole_iterator end()
  {
    return whole_iterator(nullptr, kKeyNotFound);
  }


  /// ---- Modifiers ---- ///
  inline bool insert(key_type& key, value_type& value, key_type& key_long_prbdst, value_type& val_long_prbdst)
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
    size_type total = 0;
    for (size_type pos = 0; pos < m_capacity; ++pos) {
      if (!property_program::is_empty(m_body[pos].property) &&
          !property_program::is_scratched(m_body[pos].property))
        total += property_program::extract_probedistance(m_body[pos].property);
    }

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
        property_type property = source_rhh->m_body[i].property;
        if (!property_program::is_empty(property) && !property_program::is_scratched(property)) {
          bool is_successed = new_rhh->insert_into_body(std::move(source_rhh->m_body[i].key), std::move(source_rhh->m_body[i].value), wk_key, wk_val);
          new_rhh->m_num_elems += is_successed;
          while (!is_successed) {
#if GROW_WHEN_LONG_PROBE_DISTANCE
            rhh_contatiner_selftype* old_rhh = new_rhh;
            new_rhh = make_with_source_rhh(old_rhh, new_capacity * rhh::kCapacityGrowingFactor);
            deallocate(old_rhh);
#else
            rhh_contatiner_selftype* chained_rhh = new_rhh;
            new_rhh = allocate(chained_rhh->m_capacity);
            new_rhh->assign_to_chained_rhh(chained_rhh);
#endif
            is_successed = new_rhh->insert_into_body(std::move(wk_key), std::move(wk_val), wk_key, wk_val);
            new_rhh->m_num_elems += is_successed;
          }
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
      const size_type d = property_program::extract_probedistance(m_body[i].property);
      if (d >= size && !property_program::is_empty(m_body[i].property)) {
        exit(1);
      }
      ++histgram[d];
    }
    if (m_next != nullptr) {
      m_next->histgram_load_factor(histgram);
    }
  }


 /// ---------------------------------------------------------- ///
 ///                         private
 /// ---------------------------------------------------------- ///
 private:

  inline void erase_element_at(size_type pos)
  {
    property_program::scratch(m_body[pos].property);
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
      } else if (!property_program::is_scratched(exist_property) && key == (*body_ptr)->m_body[pos].key) {
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
    probedistance_type prb_dist = 0;
    const size_type hash = key_hash_func::hash(key);
    size_type pos = cal_desired_position(hash);
    const size_type mask = (m_capacity - 1);


    while(true) {
      property_type exist_property = m_body[pos].property;

      if(property_program::is_empty(exist_property))
      {
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
          return false;
        }
        if(property_program::is_scratched(exist_property))
        {
          construct(pos, prb_dist, std::move(key), std::move(value));
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

#ifdef RHH_DETAILED_ANALYSYS
      ++(rhh::rhh_log_holder::instance().max_distance_to_moved);
#endif
    }


    return false;
  }

  void rehash_elements()
  {
    key_type wk_key;
    value_type wk_val;

    rhh_contatiner_selftype* const tmp_rhh = allocate(m_capacity);

    for (size_type i = 0; i < m_capacity; ++i) {
      const property_type property = m_body[i].property;
      if (!property_program::is_empty(property) && !property_program::is_scratched(property)) {
        /// should check the return value ?
        tmp_rhh->insert_into_body(std::move(m_body[i].key), std::move(m_body[i].value), wk_key, wk_val);
        ++tmp_rhh->m_num_elems;
      }
    }
    std::memcpy(&(m_body[0]), &(tmp_rhh->m_body[0]), m_capacity * kElementSize);

    deallocate(tmp_rhh);
  }


  /// --- private valiable --- ///
  size_type m_num_elems; /// including chained rhhs
  size_type m_capacity;
  rhh_contatiner_selftype* m_next;
  packed_element m_body[0];

};

} /// namespace graphstore

#include <havoqgt/graphstore/rhh/rhh_whole_iterator.hpp>
#include <havoqgt/graphstore/rhh/rhh_value_iterator.hpp>

#endif /// RHH_HPP
