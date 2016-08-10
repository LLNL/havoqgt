/*
 * Written by Keita Iwabuchi.
 * LLNL / TokyoTech
 */

#ifndef RHH_ALLOCATOR_HOLDER_QUARTZ_HPP_INCLUDED
#define RHH_ALLOCATOR_HOLDER_QUARTZ_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/allocators/node_allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

#include "pmalloc.h"

namespace graphstore {
namespace rhh {

struct ALLOC_IN_QUARTZ {};
///
/// \brief allocator_holder_sglt class
///  singleton pattern
///  specialized class for in core allocation
template<size_t element_size, size_t extra_size>
class allocator_holder_sglt <ALLOC_IN_QUARTZ, element_size, extra_size>
{

 public:

  using segment_manager_type = ALLOC_IN_QUARTZ;
  using self_type = allocator_holder_sglt<ALLOC_IN_QUARTZ, element_size, extra_size>;

  static self_type& instance() {
    static self_type _instance;
    return _instance;
  }

  template <typename U1, size_t U2, size_t U3>
  struct rebind
  {
      using other = allocator_holder_sglt<U1, U2, U3>;
  };

  /// --- allocate & deallocate functions --- ///
  void* allocate(const size_t capacity)
  {
    return ::pmalloc(capacity * element_size + extra_size);
  }

  void deallocate(void *ptr, const size_t capacity)
  {
    ::pfree(ptr, capacity * element_size + extra_size);
  }


 private:
  allocator_holder_sglt() {}
  allocator_holder_sglt(const self_type &)  = delete;
  allocator_holder_sglt(const self_type &&) = delete;
  self_type &operator=(const self_type &)   = delete;
  self_type &operator=(const self_type &&)  = delete;

};

}}

#endif
