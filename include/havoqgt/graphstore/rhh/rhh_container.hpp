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
#include <havoqgt/graphstore/rhh/rhh_common.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.h>
#include <havoqgt/graphstore/rhh/rhh_element_base.h>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

namespace graphstore {

namespace rhh_container_utility {

/// ---- Utility functions ---- ///
template <typename rhh_type>
inline bool has_key(rhh_type* const rhh, typename rhh_type::key_type& key)
{
  return (rhh->find(key) != rhh_type::kKeyNotFound);
}

/// Note: this function causes copy of a key and a value !!
template <typename rhh_type, typename key_type, typename value_type>
inline void insert(rhh_type** rhh, key_type key, value_type value)
{
  if ((*rhh)->size() + 1 >= static_cast<size_t>(static_cast<double>((*rhh)->capacity()) * graphstore::rhh::kFullCapacitFactor)) {
    (*rhh) = rhh_type::resize((*rhh), (*rhh)->capacity() * graphstore::rhh::kCapacityGrowingFactor);
  }
  // --- Consider long probe distance --- //
  while (!(*rhh)->insert(key, value, key, value)) {
    rhh_type* new_rhh = rhh_type::allocate((*rhh)->capacity());
    new_rhh->assign_to_chained_rhh((*rhh));
    (*rhh) = new_rhh;
  }
}

template <typename rhh_type>
inline void shrink_to_fit(rhh_type** rhh)
{
  const typename rhh_type::size_type cur_size = (*rhh)->size();
  typename rhh_type::size_type new_capacity = (*rhh)->capacity();
  while ( cur_size <
            static_cast<double>(new_capacity / graphstore::rhh::kCapacityGrowingFactor) * graphstore::rhh::kFullCapacitFactor ) {
    new_capacity /= graphstore::rhh::kCapacityGrowingFactor;
  }

  if ((*rhh)->capacity() > new_capacity) {
    (*rhh) = rhh_type::resize((*rhh), new_capacity);
  }

}

} /// namespace rhh_container_utility


template<typename _key_type,
         typename _value_type,
         typename _size_type,
         typename _key_hash_func = rhh::key_hash_func_64bit_to_64bit<_key_type, _size_type>,
         typename _property_program = rhh_container_utility::rhh_property_program_base<unsigned char>>
class rhh_container_base {

public:
  using property_program   = _property_program;
  using key_type           = _key_type;
  using key_hash_func      = _key_hash_func;
  using value_type         = _value_type;
  using size_type          = _size_type;
  using property_type      = typename property_program::property_type;
  using probedistance_type = typename property_program::probedistance_type;

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
    kKeyNotFound = std::numeric_limits<size_type>::max(),
    kElementSize = sizeof(element_type)
  };

  using rhh_contatiner_selftype = rhh_container_base<key_type, value_type, size_type, key_hash_func, property_program>;
  using allocator               = graphstore::rhh::allocator_holder_sglt<kElementSize,
                                                                           sizeof(size_type) + sizeof(size_type) + sizeof(void*)>;

  enum : size_t{
    kCapacityGrowingFactor = 2ULL
  };


  friend class WholeForwardIterator;
  template <typename Type>
  class WholeForwardIterator : public std::iterator<std::forward_iterator_tag, Type>
  {
    friend class rhh_container_base;
    using whole_iterator_selftype = WholeForwardIterator<Type>;
    using rhh_type = rhh_contatiner_selftype;

   public:

    WholeForwardIterator(rhh_type* rhh) :
      m_rhh_ptr(rhh),
      m_pos(-1) /// Note: next_valid_element increment m_pos in the first line
    {
      next_valid_element();
    }

    WholeForwardIterator(rhh_type* rhh, rhh_type::size_type pos) :
      m_rhh_ptr(rhh),
      m_pos(pos)
    { }

    WholeForwardIterator(const WholeForwardIterator& src)
    {
      m_rhh_ptr = src.m_rhh_ptr;
      m_pos      = src.m_pos;
    }

    void swap(WholeForwardIterator &other) noexcept
    {
      using std::swap;
      swap(m_rhh_ptr, other.m_rhh_ptr);
      swap(m_pos, other.m_pos);
    }

    WholeForwardIterator &operator++ () // Pre-increment
    {
      next_valid_element();
      return *this;
    }

    WholeForwardIterator operator++ (int) // Post-increment
    {
      WholeForwardIterator tmp(*this);
      next_valid_element();
      return tmp;
    }

    // two-way comparison: v.begin() == v.cbegin() and vice versa
    template<class OtherType>
    bool operator == (const WholeForwardIterator<OtherType> &rhs) const
    {
      return is_equal(rhs);
    }

    template<class OtherType>
    bool operator != (const WholeForwardIterator<OtherType> &rhs) const
    {
      return !is_equal(rhs);
    }

    Type& operator* () const
    {
      return m_rhh_ptr->m_body[m_pos];
    }

    Type* operator-> () const
    {
      return &(m_rhh_ptr->m_body[m_pos]);
    }

    // One way conversion: iterator -> const_iterator
    operator WholeForwardIterator<const Type>() const
    {
      return WholeForwardIterator<const Type>(m_rhh_ptr, m_pos);
    }


    /// --- performance optimized methods --- ///
    inline bool is_end() const
    {
      return (m_pos == rhh_type::kKeyNotFound);
    }

   private:

    WholeForwardIterator()
    { }

    void next_valid_element()
    {
      ++m_pos;
      while(m_rhh_ptr) {
        for (; m_pos < m_rhh_ptr->capacity(); ++m_pos) {
          if (!property_program::is_empty(m_rhh_ptr->m_body[m_pos].property) && !property_program::is_scratched(m_rhh_ptr->m_body[m_pos].property))
            return;
        }
        m_pos = 0;
        m_rhh_ptr = m_rhh_ptr->chained_rhh();
      }
      m_pos = rhh_type::kKeyNotFound;
    }

    template<class OtherType>
    inline bool is_equal(const WholeForwardIterator<OtherType> &rhs) const
    {
      return (m_rhh_ptr == rhs.m_rhh_ptr) && (m_pos == rhs.m_pos);
    }

    rhh_type* m_rhh_ptr;
    rhh_type::size_type m_pos;
  };
  // iteratorをtypedefする
  using whole_iterator = WholeForwardIterator<element_type>;
  using const_whole_iterator = WholeForwardIterator<const element_type>;



  friend class ValueForwardIterator;
  template <typename Type>
  class ValueForwardIterator : public std::iterator<std::forward_iterator_tag, Type>
  {
    friend class rhh_container_base;
    using value_iterator_selftype = ValueForwardIterator<Type>;
    using rhh_type                = rhh_contatiner_selftype;

   public:

    ValueForwardIterator(rhh_type* rhh, const key_type& key) :
      m_rhh_ptr(rhh),
      m_key(key),
      m_pos(0),
      m_prb_dist(0)
    {
      m_rhh_ptr->internal_locate(m_key, &m_rhh_ptr, m_pos, m_prb_dist);
    }

    ValueForwardIterator(rhh_type* rhh, const key_type& key, rhh_type::size_type pos,  rhh_type::probedistance_type prb_dist) :
      m_rhh_ptr(rhh),
      m_key(key),
      m_pos(pos),
      m_prb_dist(prb_dist)
    { }

    void swap(value_iterator_selftype &other) noexcept
    {
      using std::swap;
      swap(m_rhh_ptr, other.m_rhh_ptr);
      swap(m_pos, other.m_pos);
      swap(m_prb_dist, other.m_prb_dist);
    }

    value_iterator_selftype &operator++ () // Pre-increment
    {
      find_next_value();
      return *this;
    }

    value_iterator_selftype operator++ (int) // Post-increment
    {
      value_iterator_selftype tmp(*this);
      find_next_value();
      return tmp;
    }

    // two-way comparison: v.begin() == v.cbegin() and vice versa
    template<class OtherType>
    bool operator == (const ValueForwardIterator<OtherType> &rhs) const
    {
      return is_equal(rhs);
    }

    template<class OtherType>
    bool operator != (const ValueForwardIterator<OtherType> &rhs) const
    {
      return !is_equal(rhs);
    }

    Type& operator* () const
    {
      return m_rhh_ptr->m_body[m_pos].value;
    }

    Type* operator-> () const
    {
      return &(m_rhh_ptr->m_body[m_pos].value);
    }

    // One way conversion: iterator -> const_iterator
    operator ValueForwardIterator<const Type>() const
    {
      return ValueForwardIterator<const Type>(m_rhh_ptr, m_key, m_pos);
    }


    /// --- performance optimized methods --- ///
    inline bool is_end() const
    {
      return (m_pos == rhh_type::kKeyNotFound);
    }


   private:

    ValueForwardIterator()
    { }

    template<class OtherType>
    inline bool is_equal(const ValueForwardIterator<OtherType> &rhs) const
    {
      return (m_rhh_ptr->m_body == rhs.m_rhh_ptr->m_body) && (m_pos == rhs.m_pos);
    }

    inline void find_next_value()
    {
      const rhh_type::size_type start_pos = (m_pos + 1) & (m_rhh_ptr->capacity() - 1);
      m_rhh_ptr->internal_locate_with_hint(m_key, start_pos, &m_rhh_ptr, m_pos, m_prb_dist);
    }

    rhh_type* m_rhh_ptr;
    const key_type m_key;
    rhh_type::size_type m_pos;
    rhh_type::probedistance_type m_prb_dist;
  };
  // iteratorをtypedefする
  using value_iterator       = ValueForwardIterator<value_type>;
  using const_value_iterator = ValueForwardIterator<const value_type>;


  /// --- Explicitly Delete -- ///
  rhh_container_base()  =delete;
  ~rhh_container_base() =delete;
  rhh_container_base(const rhh_contatiner_selftype& src) =delete;
  rhh_container_base& operator=(const rhh_container_base&) =delete;
  rhh_container_base(const rhh_container_base&&) =delete;
  rhh_container_base& operator=(rhh_container_base&& old_obj) =delete;

  /// --- operator ---- ///
  //  inline element_type& operator[](size_type pos) {
  //    return m_body[pos];
  //  }

  /// --- Capacity ---- ///
  inline size_t size() const
  {
    return m_num_elems;
  }

  inline size_t capacity() const
  {
    return m_capacity;
  }


  /// --- Lookup --- ///
  inline value_iterator find(const key_type& key)
  {
    return value_iterator(this, key);
  }

  inline std::pair<value_iterator, value_iterator> equal_range(const key_type& key)
  {
    return std::make_pair(fin(key), value_iterator(nullptr, kKeyNotFound));
  }

  whole_iterator begin()
  {
    return whole_iterator(this);
  }

  whole_iterator end()
  {
    return whole_iterator(this, kKeyNotFound);
  }


  /// ---- Modifiers ---- ///
  inline bool insert(key_type& key, value_type& value, key_type& key_long_prbdst, value_type& val_long_prbdst)
  {
    const bool is_success = insert_into_body(std::move(key), std::move(value), key_long_prbdst, val_long_prbdst);
    m_num_elems += is_success;
    return is_success;
  }

  inline void clear()
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

  inline void erase(whole_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_pos);
    --m_num_elems;
  }

  inline void erase(value_iterator& itr)
  {
    itr.m_rhh_ptr->erase_element_at(itr.m_pos);
    --m_num_elems;
  }

  /// --- Hash policy --- ///
  /// \brief load_factor
  /// \return
  ///     average probe distance
  inline size_type load_factor() const
  {
    size_type total = 0;
    for (size_type pos = 0; pos < m_capacity; ++pos) {
      if (!property_program::is_empty(m_body[pos].property) &&
          !property_program::is_scratched(m_body[pos].property))
        total += property_program::extract_probedistance(m_body[pos].property);
    }

    if (m_num_elems > 0)
      return total / m_num_elems;
    else
      return 0;

  }

  inline size_type depth() const
  {
    if (!m_next) return 1;
    return m_next->depth() + 1;
  }


  /// --- allocator & deallocator --- ///
  inline static rhh_contatiner_selftype* allocate(size_t capacity)
  {
    rhh_contatiner_selftype* rhh = reinterpret_cast<rhh_contatiner_selftype*>(allocator::instance().allocate(capacity));
    rhh->m_num_elems = 0;
    rhh->m_capacity = capacity; /// XXX: if capacity is not power of 2, this is not actual capacity
    rhh->m_next = nullptr;

    const size_type cap = rhh->m_capacity;
    for (size_type i = 0; i < cap; ++i) {
      property_program::empty(rhh->m_body[i].property);
    }
    return rhh;
  }

  inline static void deallocate(rhh_contatiner_selftype* rhh)
  {
    while (rhh != nullptr) {
      // std::cout << rhh->m_num_elems << " " << rhh->m_capacity << " " << rhh->m_next << std::endl;
      rhh_contatiner_selftype* next_rhh = rhh->m_next;
      allocator::instance().deallocate(reinterpret_cast<void*>(rhh), rhh->m_capacity);
      rhh = next_rhh;
    }
  }

  inline static rhh_contatiner_selftype* resize(rhh_contatiner_selftype* rhh, size_type new_capacity)
  {
    rhh_contatiner_selftype* old_rhh = rhh;
    rhh_contatiner_selftype* new_rhh = make_with_source_rhh(old_rhh, new_capacity);
    deallocate(old_rhh);
    return new_rhh;
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

  inline void internal_locate(const key_type& key,
                              rhh_contatiner_selftype** body_ptr,
                              size_type& pos,
                              probedistance_type& cur_prb_dist)
  {
    const size_type hash = key_hash_func::hash(key);
    size_type start_pos = cal_desired_position(hash);
    cur_prb_dist = 0;
    internal_locate_with_hint(key, start_pos, body_ptr, pos, cur_prb_dist);
  }

  void internal_locate_with_hint(const key_type& key,
                            const size_type start_pos,
                            rhh_contatiner_selftype** body_ptr,
                            size_type& found_pos,
                            probedistance_type& prb_dist)
  {
    size_type pos = start_pos;
//    probedistance_type prb_dist = cur_prb_dist;
    const size_type mask = (m_capacity - 1);

    while (true) {
      property_type exist_property = m_body[pos].property;
      if (property_program::is_empty(exist_property)) {
        break;
      } else if (prb_dist > property_program::extract_probedistance(exist_property)) {
        break;
      } else if (!property_program::is_scratched(exist_property) && key == m_body[pos].key) {
        *body_ptr = this;
        found_pos = pos;
        return;
      }
      pos = (pos+1) & mask;
      ++prb_dist;
    }

    if (m_next != nullptr) {
      m_next->internal_locate(key, body_ptr, found_pos, prb_dist);
    } else {
      *body_ptr = this;
      found_pos = kKeyNotFound;
    }
  }

  inline void construct(const size_type pos, const probedistance_type dist,
                        key_type&& key, value_type&& value)
  {
    property_program::init_property(m_body[pos].property, dist);
    m_body[pos].key   = std::move(key);
    m_body[pos].value = std::move(value);
  }


  ///
  /// \brief insert_into_body
  /// \param key
  /// \param value
  /// \param key_long_prbdst: be set a value if probe ditance exceed a threshhold
  /// \param val_long_prbdst: be set a value if probe ditance exceed a threshhold
  /// \return
  ///        true: successed
  ///        false: long probe distance
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
        break;
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
          break;
        }
        m_body[pos].property = prb_dist;
        using std::swap;
        swap(key, m_body[pos].key);
        swap(value, m_body[pos].value);
        prb_dist = property_program::extract_probedistance(exist_property);
      }
      pos = (pos+1) & mask;
      ++prb_dist;
    }
    return true;
  }

  static rhh_contatiner_selftype* make_with_source_rhh(rhh_contatiner_selftype* source_rhh, size_type new_capacity)
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
            rhh_contatiner_selftype* chained_rhh = new_rhh;
            new_rhh = allocate(chained_rhh->m_capacity);
            new_rhh->assign_to_chained_rhh(chained_rhh);
            is_successed = new_rhh->insert_into_body(std::move(wk_key), std::move(wk_val), wk_key, wk_val);
            new_rhh->m_num_elems += is_successed;
          }
        }
      }
      source_rhh = source_rhh->m_next;
    }

    return new_rhh;
  }


  /// --- private valiable --- ///
  size_type m_num_elems; /// including chained rhhs
  size_type m_capacity;
  rhh_contatiner_selftype* m_next;
  packed_element m_body[0];

};

} /// namespace graphstore
#endif /// RHH_HPP
