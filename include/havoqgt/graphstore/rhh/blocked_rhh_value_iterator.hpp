#ifndef BLOCKED_RHH_VALUE_ITERATOR_HPP_INCLUDED
#define BLOCKED_RHH_VALUE_ITERATOR_HPP_INCLUDED

#include <havoqgt/graphstore/rhh/blocked_rhh_container.hpp>

namespace graphstore {

template<typename _key_type,
         typename _value_type,
         typename _size_type,
         typename _segment_manager_type,
         typename _key_hash_func,
         typename _property_program>
class blocked_rhh_container<_key_type, _value_type, _size_type,
                            _segment_manager_type,
                            _key_hash_func,
                            _property_program>::value_iterator
{

 private:
  using rhh_type = blocked_rhh_container<_key_type, _value_type, _size_type,
                                         _segment_manager_type,
                                         _key_hash_func,
                                         _property_program>;
  friend blocked_rhh_container;

 public:

  /// ---- Constructors ----
  /// initialize to an 'end iterator'
  value_iterator() :
    m_rhh_ptr(nullptr),
    m_key(),
    m_block_pos(kKeyNotFound),
    m_elem_pos(kKeyNotFoundInBlock),
    m_prb_dist()
  { }

  value_iterator(rhh_type* const rhh, const key_type& key) :
    m_rhh_ptr(rhh),
    m_key(key),
    m_block_pos(0),
    m_elem_pos(0),
    m_prb_dist(0),
    m_has_reached_last_block(false)
  {
    m_rhh_ptr = m_rhh_ptr->internal_locate(m_key, m_block_pos, m_elem_pos, m_prb_dist, m_has_reached_last_block);
  }

  /// Copy constructor
  value_iterator(const value_iterator& other) :
    m_rhh_ptr(other.m_rhh_ptr),
    m_key(other.m_key),
    m_block_pos(other.m_block_pos),
    m_elem_pos(other.m_elem_pos),
    m_prb_dist(other.m_prb_dist),
    m_has_reached_last_block(other.m_has_reached_last_block)
  { }

  /// Move constructor
  value_iterator (value_iterator&& other) :
    m_rhh_ptr(std::move(other.m_rhh_ptr)),
    m_key(std::move(other.m_key)),
    m_block_pos(std::move(other.m_block_pos)),
    m_elem_pos(std::move(other.m_elem_pos)),
    m_prb_dist(std::move(other.m_prb_dist)),
    m_has_reached_last_block(std::move(other.m_has_reached_last_block))
  { }

  /// Assignment operators
  value_iterator& operator=(const value_iterator& other)
  {
    m_rhh_ptr   = other.m_rhh_ptr;
    m_key       = other.m_key;
    m_block_pos = other.m_block_pos;
    m_elem_pos  = other.m_elem_pos;
    m_prb_dist  = other.m_prb_dist;
    m_has_reached_last_block  = other.m_has_reached_last_block;
    return *this;
  }

  /// Move operators
  value_iterator& operator=(value_iterator&& other)
  {
    m_rhh_ptr   = std::move(other.m_rhh_ptr);
    m_key       = std::move(other.m_key);
    m_block_pos = std::move(other.m_block_pos);
    m_elem_pos  = std::move(other.m_elem_pos);
    m_prb_dist  = std::move(other.m_prb_dist);
    m_has_reached_last_block  = std::move(other.m_has_reached_last_block);
    return *this;
  }


  void swap(value_iterator &other) noexcept
  {
    using std::swap;
    swap(m_rhh_ptr, other.m_rhh_ptr);
    swap(m_key, other.m_key);
    swap(m_block_pos, other.m_block_pos);
    swap(m_elem_pos, other.m_elem_pos);
    swap(m_prb_dist, other.m_prb_dist);
    swap(m_has_reached_last_block, other.m_has_reached_last_block);
  }

  // Pre-increment
  value_iterator &operator++ ()
  {
    find_next_value();
    return *this;
  }

  // Post-increment
  value_iterator operator++ (int)
  {
    value_iterator tmp(*this);
    find_next_value();
    return tmp;
  }

  bool operator == (const value_iterator &rhs) const
  {
    return is_equal(rhs);
  }

  bool operator != (const value_iterator &rhs) const
  {
    return !is_equal(rhs);
  }

  typename rhh_type::value_type& operator* () const
  {
    return m_rhh_ptr->m_body[m_block_pos][m_elem_pos].value;
  }

  typename rhh_type::value_type* operator-> () const
  {
    return &(m_rhh_ptr->m_body[m_block_pos][m_elem_pos].value);
  }

  static value_iterator end()
  {
    return value_iterator();
  }

  /// --- performance optimized methods --- ///
  inline bool is_end() const
  {
    return ((m_rhh_ptr == nullptr) && (m_block_pos == rhh_type::kKeyNotFound));
  }


 private:

  inline bool is_equal(const value_iterator &rhs) const
  {
    return (m_rhh_ptr == rhs.m_rhh_ptr) && (m_block_pos == rhs.m_block_pos) && (m_elem_pos == rhs.m_elem_pos);
  }

  inline void find_next_value()
  {
    ++m_elem_pos;
    m_rhh_ptr = m_rhh_ptr->internal_locate_with_hint(m_key, m_block_pos, m_elem_pos, m_prb_dist, m_has_reached_last_block);
  }

  rhh_type* m_rhh_ptr;
  key_type m_key;
  size_type m_block_pos;
  block_size_type m_elem_pos;
  probedistance_type m_prb_dist;
  bool m_has_reached_last_block;
};

}

#endif // RHH_VALUE_ITERATOR_HPP

