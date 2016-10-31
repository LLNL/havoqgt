#ifndef NUMA_ALLOCATOR_HPP
#define NUMA_ALLOCATOR_HPP

#include <memory>
#include <cstdlib>
#include <numa.h>

template<class T>
class numa_allocator : public std::allocator<T> {
 public:
  // typedefs
  using pointer = typename std::allocator<T>::pointer;
  using size_type = typename std::allocator<T>::size_type;
  using const_pointer = typename std::allocator<T>::const_pointer;

  template<class U>
  struct rebind
  {
    using other = numa_allocator<U>;
  };

  numa_allocator() = default;
  numa_allocator(const numa_allocator&) = default;

  template<class U>
  numa_allocator(const numa_allocator<U>& x) {}

  pointer allocate(size_type n, const_pointer hint = 0) {
      return reinterpret_cast<pointer>(numa_alloc_local(n * sizeof(T)));
  }

  void deallocate(pointer ptr, size_type size) {
      numa_free(ptr, size);
  }
};
#endif // NUMA_ALLOCATOR_HPP
