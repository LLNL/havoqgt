/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHH_UTILITIES_HPP_INCLUDED
#define RHH_UTILITIES_HPP_INCLUDED

#include <havoqgt/detail/hash.hpp>

namespace graphstore {
namespace rhh {

/// ---- Utility functions ---- ///
template <typename rhh_type>
inline static void resize(rhh_type** rhh, const typename rhh_type::size_type new_capacity)
{
  rhh_type* new_rhh = rhh_type::make_with_source_rhh(*rhh, new_capacity);
  rhh_type::deallocate(*rhh);
  *rhh = new_rhh;
}

template <typename rhh_type>
inline static void chain(rhh_type** rhh)
{
  rhh_type* const new_rhh = rhh_type::allocate((*rhh)->capacity());
  new_rhh->assign_to_chained_rhh(*rhh);
  *rhh = new_rhh;
}

template <typename rhh_type>
inline static void grow(rhh_type** rhh)
{
#if RHH_CHAIN_AT_LARGE_TABLE_SIZE
  if ((*rhh)->table_mem_size() > 4096) {
    chain(rhh);
  } else {
    resize(rhh, (*rhh)->capacity() * kCapacityGrowingFactor);
  }
#else
  resize(rhh, (*rhh)->capacity() * kCapacityGrowingFactor);
#endif
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

  /// --- dealing with a long probe distance problem --- ///
  if (!(*rhh)->insert(std::move(key), std::move(value), key, value)) {
#if RHH_ATTEMPT_GROWING_TO_SOLVE_LONG_PROBE_DISTANCE
    grow(rhh);
#else
    chain(rhh);
#endif
    if (!(*rhh)->insert(std::move(key), std::move(value), key, value)) {
      /// This program enter this point when a long probe distance problem is occurred
      /// and growing capacity strategy can't solve the problem.
      /// When there are a lot of element whoes key is exactoly same,
      /// chaining operation is the only strategy that can slove the problem.
      chain(rhh);
      assert((*rhh)->insert(std::move(key), std::move(value), key, value));
    }
  }
}


/// --- A template helper caluculate the size of a element --- ///
///     if type is empy_type return 0;
///     the other cases, simply return sizeof(value)
struct empy_type {};

template<typename type>
struct sizeof_{
  enum {
    size = sizeof(type)
  };
};
template<>
struct sizeof_<empy_type> {
  enum {
    size = 0
  };
};

// A template helper used to return the element size.
template <typename property, typename key, typename value>
struct element_size
{
  enum : size_t {
    size = sizeof_<property>::size + sizeof_<key>::size + sizeof_<value>::size
  };
};


// ---- A template helper used to select A or B based on a condition --- ///
template<bool cond, typename A, typename B>
struct if_{
  using type = A;
};
template<typename A, typename B>
struct if_<false, A, B> {
  using type = B;
};


template<typename key_type, typename hash_type>
class key_hash_func_64bit_to_64bit {
 public:
  inline static hash_type hash(const key_type& key)
  {
#if 1
    return static_cast<hash_type>(static_cast<uint64_t>(havoqgt::detail::hash32(key>>32ULL)) << 32ULL
                                                              | havoqgt::detail::hash32(key));
#else
    return static_cast<hash_type>(key);
#endif
  }
};

}
}
#endif // RHH_UTILITY
