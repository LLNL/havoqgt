#ifndef GRAPHSTORE_UTILITIES_HPP
#define GRAPHSTORE_UTILITIES_HPP

#include <memory>
#include <stdlib.h>
namespace graphstore {
  void* aligned_alloc(size_t length, size_t align_size)
  {
    const size_t headspace_size = sizeof(void *);
    const size_t buffer_size = headspace_size + align_size + length;

#if 0
    /// std::align is not implemented in gcc ?
    void* actual_buffer = ::operator new(buffer_size);
    void* aligned_ptr = std::align(align_size, length, actual_buffer+headspace_size, buffer_size);
    ((void **)aligned_ptr)[-1] = actual_buffer;
#else
    void* actual_buffer;
    int result = posix_memalign(&actual_buffer, align_size, buffer_size);
    void* aligned_ptr = actual_buffer;
#endif
    return aligned_ptr;
  }

  void aligned_free(void* ptr)
  {
    if (!ptr) return;
#if 0
     void* actual_buffer = ((void **)ptr)[-1];
     ::operator delete(actual_buffer);
#else
    free(ptr);
#endif
  }

}
#endif // GRAPHSTORE_UTILITIES_HPP

