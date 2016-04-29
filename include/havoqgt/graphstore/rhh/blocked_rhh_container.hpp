/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef BLOCKED_RHH_CONTAINER_HPP_INCLUDED
#define BLOCKED_RHH_CONTAINER_HPP_INCLUDED

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
class blocked_rhh_container {
 public:
  using property_program     = _property_program;
  using key_type             = _key_type;
  using key_hash_func        = _key_hash_func;
  using value_type           = _value_type;
  using size_type            = _size_type;
  using block_size_type      = uint8_t;
  using segment_manager_type = _segment_manager_type;
  using property_type        = typename property_program::property_type;
  using probedistance_type   = typename property_program::probedistance_type;
  using rhh_contatiner_selftype = blocked_rhh_container<key_type, value_type, size_type,
                                                        segment_manager_type,
                                                        key_hash_func, property_program>;

  enum : block_size_type {
    kBlockCapacity = 8,
    kKeyNotFoundInBlock = std::numeric_limits<block_size_type>::max()
  };

  /// TODO: specializetion for no value case
  #pragma pack(1)
  struct packed_element
  {
    property_type property;
    key_type      key;
    value_type    value;
  };
  using element_type = packed_element;
  using block_type   = std::array<element_type, kBlockCapacity>;

  enum : size_type {
    kKeyNotFound   = std::numeric_limits<size_type>::max(),
    kElementSize   = sizeof(element_type),
    kBlockSize     = sizeof(block_type),
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
  blocked_rhh_container()  = delete;
  ~blocked_rhh_container() = delete;
  blocked_rhh_container(const rhh_contatiner_selftype& src) = delete;
  blocked_rhh_container& operator=(const blocked_rhh_container&) = delete;
  blocked_rhh_container(const blocked_rhh_container&&) = delete;
  blocked_rhh_container& operator=(blocked_rhh_container&& old_obj) = delete;


  /// --- Capacity ---- ///
  inline size_t size() const
  {
    return m_num_elems;
  }

  inline size_t capacity() const
  {
    return m_num_block * kBlockCapacity;
  }

  inline size_t table_mem_size() const
  {
    return kBlockSize * capacity() + kTableBaseSize;
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
    for (size_type i = 0; i < m_num_block; ++i) {
      for (size_type j = 0; j < kBlockCapacity; ++j) {
        erase_element_at(i, j);
      }
    }
    m_num_elems = 0;

    if (m_next) {
      m_next->clear();
    }
  }

  /// TODO: move to non-member function
  inline void erase(const whole_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_block_pos, itr.m_elem_pos);
    --m_num_elems;
  }

  /// TODO: move to non-member function ?
  inline void erase(const value_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_block_pos, itr.m_elem_pos);
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
//    size_type total = 0;
//    for (size_type pos = 0; pos < m_num_block; ++pos) {
//      if (!property_program::is_empty(m_body[pos].property) &&
//          !property_program::is_tombstone(m_body[pos].property))
//        total += property_program::extract_probedistance(m_body[pos].property);
//    }

//    if (m_num_elems > 0)
//      return static_cast<double>(total) / static_cast<double>(m_num_elems);
//    else
//      return 0.0;

  }

  inline size_type depth() const
  {
    if (!m_next) return 1;
    return m_next->depth() + 1;
  }


  /// --- allocator & deallocator --- ///
  static rhh_contatiner_selftype* allocate(const size_t capacity)
  {
    const size_type num_block = (capacity + kBlockCapacity - 1) / kBlockCapacity;
    rhh_contatiner_selftype* rhh = reinterpret_cast<rhh_contatiner_selftype*>(allocator::instance().allocate(num_block * kBlockCapacity));
    rhh->m_num_elems = 0;
    rhh->m_num_block = num_block;
    rhh->m_next = nullptr;

    for (size_type i = 0; i < num_block; ++i) {
      block_type& block = rhh->m_body[i];
      for (size_type j = 0; j < kBlockCapacity; ++j) {
        property_program::empty(block[j].property);
      }
    }
    return rhh;
  }

  static void deallocate(rhh_contatiner_selftype* rhh)
  {
    while (rhh != nullptr) {
      // std::cout << rhh->m_num_elems << " " << rhh->m_num_block << " " << rhh->m_next << std::endl;
      rhh_contatiner_selftype* next_rhh = rhh->m_next;
      allocator::instance().deallocate(reinterpret_cast<void*>(rhh), rhh->m_num_block);
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
    m_num_block = next->m_num_block;
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

    const size_type src_num_block = source_rhh->m_num_block;
    key_type wk_key;
    value_type wk_val;

    rhh_contatiner_selftype* new_rhh = allocate(new_capacity);

    while (source_rhh != nullptr) {
      for (size_type i = 0; i < src_num_block; ++i) {
        block_type& block = source_rhh->m_body[i];
        for (size_type j = 0; j < kBlockCapacity; ++j) {
          const property_type property = block[j].property;
          if (property_program::is_empty(property)) {
            break;
          } else if (!property_program::is_tombstone(property)) {
            bool is_successed = new_rhh->insert_into_body(std::move(block[j].key), std::move(block[j].value), wk_key, wk_val);
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
      }
      source_rhh = source_rhh->m_next;
    }

    return new_rhh;
  }


  /// --- debug --- ///
  void print_all_element()
  {
    for (size_type i = 0; i < m_num_block; ++i) {
      for (block_size_type j = 0; j < kBlockCapacity; ++j) {
        std::cout << "property " << (uint64_t)m_body[i][j].property << ", key " << m_body[i][j].key << std::endl;
      }
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
//    for (size_type i = 0; i < m_num_block; ++i) {
//      const size_type d = property_program::extract_probedistance(m_body[i].property);
//      if (d >= size && !property_program::is_empty(m_body[i].property)) {
//        exit(1);
//      }
//      ++histgram[d];
//    }
//    if (m_next != nullptr) {
//      m_next->histgram_load_factor(histgram);
//    }
  }


 /// ---------------------------------------------------------- ///
 ///                         private
 /// ---------------------------------------------------------- ///
 private:

  /// --- utilities --- ///
  inline size_type cal_desired_block(size_type hash) const
  {
    return hash & (m_num_block - 1);
  }


  /// --- find ---- ///
  inline rhh_contatiner_selftype* internal_locate(const key_type& key,
                                                  size_type& block_pos,
                                                  block_size_type& elem_pos,
                                                  probedistance_type& cur_prbdist,
                                                  bool& has_reached_last_block)
  {
    const size_type hash = key_hash_func::hash(key);
    block_pos = cal_desired_block(hash);
    elem_pos = 0;
    cur_prbdist = 0;
    has_reached_last_block = false;
    return internal_locate_with_hint(key, block_pos, elem_pos, cur_prbdist, has_reached_last_block);
  }

  rhh_contatiner_selftype* internal_locate_with_hint(const key_type& key,
                                                     size_type& block_pos,
                                                     block_size_type& elem_pos,
                                                     probedistance_type& prbdist,
                                                     bool& has_reached_last_block)
  {
    const size_type mask = (m_num_block - 1);

    while (true) {
      auto ret = internal_locate_in_block_with_hint(m_body[block_pos], key, elem_pos, prbdist);
      elem_pos = ret.first;
      has_reached_last_block = ret.second;
      if (elem_pos == kKeyNotFoundInBlock) {
        if (has_reached_last_block) {
          break;
        }
      } else {
        return this;
      }

      /// Check next block
      block_pos = (block_pos+1) & mask;
      elem_pos = 0;
      ++prbdist;
    }

    if (m_next != nullptr) { /// Check chained table
      return m_next->internal_locate(key, block_pos, elem_pos, prbdist, has_reached_last_block);
    }

    block_pos = kKeyNotFound;
    elem_pos = kKeyNotFoundInBlock;
    return nullptr;
  }

  std::pair<block_size_type, bool>
  internal_locate_in_block_with_hint(const block_type& block,
                                     const key_type& key,
                                     const block_size_type init_elem_pos,
                                     const probedistance_type prbdist) const
  {
    bool has_reached_last_block = false;

    for (block_size_type elem_pos = init_elem_pos; elem_pos < kBlockCapacity; ++elem_pos) {
      const property_type exist_property = block[elem_pos].property;
      if (property_program::is_empty(exist_property)) { /// Must check if empty, first
        return std::make_pair(kKeyNotFoundInBlock, true);
      } else if (!property_program::is_tombstone(exist_property) && key == block[elem_pos].key) {
        /// Found
        return std::make_pair(elem_pos, has_reached_last_block);
      }
      if (prbdist > property_program::extract_probedistance(exist_property)) {
        has_reached_last_block = true;
      }
    }

    return std::make_pair(kKeyNotFoundInBlock, has_reached_last_block);
  }


  /// --- Insertion --- ///
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
    probedistance_type prbdist = 0;
    const size_type hash = key_hash_func::hash(key);
    size_type block_pos = cal_desired_block(hash);
    const size_type mask = (m_num_block - 1);

    while (!property_program::is_long_probedistance(prbdist)) {

      block_type& block = m_body[block_pos];
      const auto ret = find_avairable_space_in_block(block, prbdist);
      element_type& target_elem = block[ret.second];
      if (ret.first) { /// An empty space or tombstone found
        construct(target_elem, prbdist, std::move(key), std::move(value));
        return true;
      } else {
        const probedistance_type exist_prbdist = property_program::extract_probedistance(target_elem.property);
        if (prbdist > exist_prbdist) {
          property_program::init_property(target_elem.property, prbdist);
          using std::swap;
          swap(key, target_elem.key);
          swap(value, target_elem.value);
          prbdist = exist_prbdist;
        }
      }

      block_pos = (block_pos+1) & mask;
      ++prbdist;

    }

    key_long_prbdst = std::move(key);
    val_long_prbdst = std::move(value);
    return false;
  }

  ///
  /// \brief find_space_in_block
  /// \param block
  /// \return
  ///         first element (bool): whether an avairable space found
  ///         second element (block_size_type): the position of a space avairable or
  ///                                  minimum probe distance
  std::pair<bool, block_size_type> find_avairable_space_in_block(block_type& block, const probedistance_type prbdist) const
  {

    probedistance_type min_prbdist = std::numeric_limits<probedistance_type>::max();
    block_size_type min_prbdist_pos = kBlockCapacity;

    for (block_size_type elem_pos = 0; elem_pos < kBlockCapacity; ++elem_pos) {

      const property_type exist_property = block[elem_pos].property;

      /// Empty space found
      if (property_program::is_empty(exist_property)) {
        return std::make_pair(true, elem_pos);
      }

      /// Tombstone found
      const probedistance_type exist_prbdist = property_program::extract_probedistance(exist_property);
      if (property_program::is_tombstone(exist_property) && exist_prbdist <= prbdist) {
          return std::make_pair(true, elem_pos);
      }

      /// Update a minimum probe distance element, a candidate of swaped out.
      if (min_prbdist > exist_prbdist) {
        min_prbdist = exist_prbdist;
        min_prbdist_pos = elem_pos;
      }
    }

    return std::make_pair(false, min_prbdist_pos);;
  }

  inline void construct(element_type& space, const probedistance_type dist,
                        key_type&& key, value_type&& value)
  {
    property_program::init_property(space.property, dist);
    space.key   = std::move(key);
    space.value = std::move(value);
    if (dist >= m_num_block) {
      rehash_elements();
    }
  }


  /// --- deletion --- ///
  inline void erase_element_at(const size_type block_pos, const size_type elem_pos)
  {
    property_program::tombstone(m_body[block_pos][elem_pos].property);
    m_body[block_pos][elem_pos].key.~key_type();
    m_body[block_pos][elem_pos].value.~value_type();
  }


  /// --- optimization --- ///
  void rehash_elements()
  {
//    key_type wk_key;
//    value_type wk_val;

//    rhh_contatiner_selftype* const tmp_rhh = allocate(m_num_block);

//    for (size_type i = 0; i < m_num_block; ++i) {
//      const property_type property = m_body[i].property;
//      if (!property_program::is_empty(property) && !property_program::is_tombstone(property)) {
//        /// should check the return value ?
//        tmp_rhh->insert_into_body(std::move(m_body[i].key), std::move(m_body[i].value), wk_key, wk_val);
//        ++tmp_rhh->m_num_elems;
//      }
//    }
//    std::memcpy(&(m_body[0]), &(tmp_rhh->m_body[0]), m_num_block * kElementSize);

//    deallocate(tmp_rhh);
  }


  /// --- private valiable --- ///
  size_type m_num_elems; /// including chaining rhh tables
  size_type m_num_block;  /// this table's size
  rhh_contatiner_selftype* m_next;
  block_type m_body[0];

};

} /// namespace graphstore

#include <havoqgt/graphstore/rhh/blocked_rhh_whole_iterator.hpp>
#include <havoqgt/graphstore/rhh/blocked_rhh_value_iterator.hpp>

#endif /// RHH_HPP
