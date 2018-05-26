#pragma once

// using namespace std;
// Based upon
// https://stackoverflow.com/questions/21028299/is-this-behavior-of-
// vectorresizesize-type-n-under-c11-and-boost-container/21028912#21028912

template <class T>
class default_init_alloc : public std::allocator<T> {
 public:
  using std::allocator<T>::allocator;

  template <class U, class... Args>
  void construct(U*, Args&&...) {}
};

// template <typename T, typename A=std::allocator<T>>
//   class default_init_alloc : public A {

//   typedef std::allocator_traits<A> a_t;
// public:
//   template <typename U> struct rebind {
//     using other =
//       default_init_alloc<U, typename a_t::template rebind_alloc<U>>;
//   };

//   using A::A;

//   template <typename U>
//   void construct(U* ptr)
//     noexcept(std::is_nothrow_default_constructible<U>::value) {
//     ::new(static_cast<void*>(ptr)) U;
//   }
//   template <typename U, typename...Args>
//   void construct(U* ptr, Args&&... args) {
//     a_t::construct(static_cast<A&>(*this), ptr, std::forward<Args>(args)...);
//   }
// };

// Usage:
// std::vector<int, no_init_alloc<int>> std_vec;