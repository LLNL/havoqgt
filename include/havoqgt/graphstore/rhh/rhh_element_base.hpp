/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef rhh_property_program_base_HPP_INCLUDED
#define rhh_property_program_base_HPP_INCLUDED

#include <havoqgt/detail/hash.hpp>
#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.hpp>

namespace graphstore {
namespace rhh {

template<typename _property_type>
class rhh_property_program_base {
 public:

  using property_type      = _property_type;
  using probedistance_type = _property_type;

  enum : probedistance_type {
    kLongProbedistanceThreshold = std::numeric_limits<probedistance_type>::max() / static_cast<property_type>(4) - 1
  };

  /// --- Explicitly Deleted Functions -- ///
  rhh_property_program_base() =delete;
  ~rhh_property_program_base() =delete;
  rhh_property_program_base(const rhh_property_program_base&) =delete;
  rhh_property_program_base& operator=(const rhh_property_program_base&) =delete;
  rhh_property_program_base(const rhh_property_program_base&&) =delete;
  rhh_property_program_base& operator=(rhh_property_program_base&& old_obj) =delete;

  inline static void empty(property_type& prop)
  {
    prop = kEmptyProperyValue;
  }

  inline static void scratch(property_type& prop)
  {
    prop |= kTomstone;
  }

  inline static bool is_empty(const property_type prop)
  {
    return (prop == kEmptyProperyValue);
  }

  inline static bool is_scratched(const property_type prop)
  {
    return ((prop & kTomstoneExtractMask) == kTomstone);
  }

  inline static probedistance_type extract_probedistance(const property_type prop)
  {
    return (prop & kClearTombstoneMask);
  }

  /// if dist is equal or more than kTomstone, this function won't work collectoly
  inline static void init_property(property_type& prop, const probedistance_type dist)
  {
    prop = dist;
  }

  inline static bool is_long_probedistance(const probedistance_type dist)
  {
    return (dist >= kLongProbedistanceThreshold);
  }

 private:
  enum : property_type {
    kTomstone = static_cast<property_type>(1) << (sizeof(property_type) * static_cast<property_type>(8) - 1),
    kEmptyProperyValue = kTomstone - 1,
    kTomstoneExtractMask = kTomstone,
    kClearTombstoneMask  = kTomstone - 1,
  };
};


}}
#endif // rhh_property_program_base

