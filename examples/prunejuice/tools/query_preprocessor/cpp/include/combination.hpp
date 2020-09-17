#pragma once

namespace prunejuice { namespace pattern {

//template <typename T>
class combination {
private:  
  /**
   * Generate combinations
   */
  //template <typename N, typename K, typename NCK>
  template <typename NCollection, typename ItemType, typename O,
    typename CombinationBuffer, typename S>
  static void generate_combinations_recursive(size_t s, size_t n, size_t k, size_t r,
    CombinationBuffer& buffer, O& combinations, NCollection& n_items,
    S& combinations_items) {
    if (r == k) {

      //for (auto& t : buffer) {
      //  std::cout << t;
      //}
      //std::cout << std::endl;

      combinations.push_back(buffer);

      std::unordered_set<ItemType> buffer_set(0);
      for (auto i : buffer) {
        buffer_set.insert(n_items[i]);
      }
      assert(buffer.size() == buffer_set.size());
      combinations_items.push_back(buffer_set);

      return;
    }

    for (size_t i = s; i < n ; i++) {
      buffer[r] = i;
      generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
        (i + 1, n, k, r + 1, buffer, combinations, n_items, combinations_items);
    } // for
  }

public:
  //template <typename N, typename K, typename NCK>
  template <typename NCollection, typename ItemType, typename O, typename S>
  static void generate_combinations_recursive(size_t n, size_t k,
    NCollection& n_items, O& combinations, S& combinations_items) {

    typedef std::vector<size_t> CombinationBuffer;
    CombinationBuffer buffer(k); // size_t - array index type

    generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
      (0, n, k, 0, buffer, combinations, n_items, combinations_items);
      // size_t s, size_t n, size_t k, size_t r, ...
  }

  combination() {}
  ~combination() {}  
 
};

}} // end namespace prunejuice::pattern
