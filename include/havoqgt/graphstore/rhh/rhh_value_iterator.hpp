#ifndef RHH_VALUE_ITERATOR_HPP
#define RHH_VALUE_ITERATOR_HPP

#include <havoqgt/graphstore/rhh/rhh_container.hpp>

namespace graphstore {

template<typename _key_type,
         typename _value_type,
         typename _size_type,
         typename _segment_manager_type,
         typename _key_hash_func = rhh::key_hash_func_64bit_to_64bit<_key_type, _size_type>,
         typename _property_program = rhh::rhh_property_program_base<unsigned char>>
class rhh_container<_key_type, _value_type, _size_type,
                    _segment_manager_type,
                    _key_hash_func,
                    _property_program>::value_iterator
{

 private:
  using rhh_type = rhh_container<_key_type, _value_type, _size_type,
                                  _segment_manager_type,
                                  _key_hash_func,
                                  _property_program>;
  friend rhh_container;

 public:

  value_iterator() :
    m_rhh_ptr(nullptr),
    m_key(),
    m_pos(rhh_type::kKeyNotFound),
    m_prb_dist()
  { }

  value_iterator(rhh_type* const rhh, const key_type& key) :
    m_rhh_ptr(rhh),
    m_key(key),
    m_pos(0),
    m_prb_dist(0)
  {
    internal_locate(m_key, const_cast<const rhh_type**>(&m_rhh_ptr), m_pos, m_prb_dist);
  }

  void swap(value_iterator &other) noexcept
  {
    using std::swap;
    swap(m_rhh_ptr, other.m_rhh_ptr);
    swap(m_pos, other.m_pos);
    swap(m_prb_dist, other.m_prb_dist);
  }

  value_iterator &operator++ () // Pre-increment
  {
    find_next_value();
    return *this;
  }

  value_iterator operator++ (int) // Post-increment
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
    return m_rhh_ptr->m_body[m_pos].value;
  }

  typename rhh_type::value_type* operator-> () const
  {
    return &(m_rhh_ptr->m_body[m_pos].value);
  }

  /// --- performance optimized methods --- ///
  inline bool is_end() const
  {
    return ((m_rhh_ptr == nullptr) && (m_pos == rhh_type::kKeyNotFound));
  }


 private:

  inline bool is_equal(const value_iterator &rhs) const
  {
    return (m_rhh_ptr == rhs.m_rhh_ptr) && (m_pos == rhs.m_pos);
  }

  inline void find_next_value()
  {
    m_pos = (m_pos + 1) & (m_rhh_ptr->capacity() - 1);
    ++m_prb_dist;
    internal_locate_with_hint(m_key, const_cast<const rhh_type**>(&m_rhh_ptr), m_pos, m_prb_dist);
  }

  rhh_type* m_rhh_ptr;
  const typename rhh_type::key_type m_key;
  typename rhh_type::size_type m_pos;
  typename rhh_type::probedistance_type m_prb_dist;
};

}

#endif // RHH_VALUE_ITERATOR_HPP

