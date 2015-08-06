/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHH_UTILITIES_HPP_INCLUDED
#define RHH_UTILITIES_HPP_INCLUDED

namespace graphstore {
namespace rhh {


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
  typedef A type;
};
template<typename A, typename B>
struct if_<false, A, B> {
  typedef B type;
};



template<typename key_type, typename hash_type>
class key_hash_func_64bit_to_64bit {
 public:
  inline static hash_type hash(const key_type& key)
  {
    #include <havoqgt/detail/hash.hpp>
    return static_cast<hash_type>(static_cast<uint64_t>(havoqgt::detail::hash32(key>>32ULL)) << 32ULL
                                                              | havoqgt::detail::hash32(key));
    //    return static_cast<hash_type>(key);
  }
};

}
}
#endif // RHH_UTILITY
