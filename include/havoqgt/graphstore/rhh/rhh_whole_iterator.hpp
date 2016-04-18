#ifndef RHH_WHOLE_ITERATOR_HPP
#define RHH_WHOLE_ITERATOR_HPP

#include <havoqgt/graphstore/rhh/rhh_container.hpp>

namespace graphstore {

template<typename _key_type,
         typename _value_type,
         typename _size_type,
         typename _segment_manager_type,
         typename _key_hash_func,
         typename _property_program>
class rhh_container<_key_type, _value_type, _size_type,
                    _segment_manager_type,
                    _key_hash_func,
                    _property_program>::whole_iterator
{
private:
  using rhh_type = rhh_container<_key_type, _value_type, _size_type,
                                  _segment_manager_type,
                                  _key_hash_func,
                                  _property_program>;
  friend rhh_container;

public:

  /// initialize to an 'end iterator'
  whole_iterator() :
    m_rhh_ptr(nullptr),
    m_pos(kKeyNotFound)
  { }

  explicit whole_iterator(rhh_type* const rhh) :
    m_rhh_ptr(rhh),
    m_pos(-1) /// Note: next_valid_element increment m_pos in the first line
  {
    next_valid_element();
  }

  /// Copy constructor
  whole_iterator(const whole_iterator& other) :
    m_rhh_ptr(other.m_rhh_ptr),
    m_pos(other.m_pos)
  { }

  /// Move constructor
  whole_iterator (whole_iterator&& other) :
    m_rhh_ptr(std::move(other.m_rhh_ptr)),
    m_pos(std::move(other.m_pos))
  { }

  /// Copy assignment operators
  whole_iterator& operator=(const whole_iterator& other)
  {
    m_rhh_ptr = other.m_rhh_ptr;
    m_pos     = other.m_pos;
    return *this;
  }

  /// Move assignment operators
  whole_iterator& operator=(whole_iterator&& other)
  {
    m_rhh_ptr = std::move(other.m_rhh_ptr);
    m_pos     = std::move(other.m_pos);
    return *this;
  }

  void swap(whole_iterator &other) noexcept
  {
    using std::swap;
    swap(m_rhh_ptr, other.m_rhh_ptr);
    swap(m_pos, other.m_pos);
  }

  whole_iterator &operator++ () // Pre-increment
  {
    next_valid_element();
    return *this;
  }

  whole_iterator operator++ (int) // Post-increment
  {
    whole_iterator tmp(*this);
    next_valid_element();
    return tmp;
  }

  bool operator == (const whole_iterator &rhs) const
  {
    return is_equal(rhs);
  }

  bool operator != (const whole_iterator &rhs) const
  {
    return !is_equal(rhs);
  }

  typename rhh_type::element_type& operator* () const
  {
    return m_rhh_ptr->m_body[m_pos];
  }

  typename rhh_type::element_type* operator-> () const
  {
    return &(m_rhh_ptr->m_body[m_pos]);
  }

  static whole_iterator end()
  {
    return whole_iterator();
  }

  /// --- performance optimized methods --- ///
  inline bool is_end() const
  {
    return (m_pos == rhh_type::kKeyNotFound);
  }


 private:

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

  inline bool is_equal(const whole_iterator &rhs) const
  {
    return (m_rhh_ptr == rhs.m_rhh_ptr) && (m_pos == rhs.m_pos);
  }

  rhh_type* m_rhh_ptr;
  typename rhh_type::size_type m_pos;
};

}
#endif // RHH_WHOLE_ITERATOR_HPP

