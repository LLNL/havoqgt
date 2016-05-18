//
// Created by Poudel, Suraj on 4/12/16.
//

#ifndef INTERVAL_TREE_HPP
#define INTERVAL_TREE_HPP

#include <iostream>
#include <vector>
#include <unordered_set>

//TODO : shrink to condensed tree pending will make it faster by reducing
//the height of the tree
template<typename key, typename index_type, typename segment_data>
class interval_tree {
public:
  class Node {
  public:
    key k;
    std::vector<index_type> left;
    std::vector<index_type> right;

    Node() : k(0) { }

    explicit Node(key _k) : k(_k) { }

    Node &operator=(const Node &node) {
      k = node.k;
      left.clear();
      right.clear();
      left.resize(node.left.size());
      right.resize(node.right.size());
      std::copy(node.left.begin(), node.left.end(), left.begin());
      std::copy(node.right.begin(), node.right.end(), right.begin());
      return *this;
    }

  };

  segment_data *seg_data;
    
  interval_tree() { }
  
  explicit interval_tree(segment_data *_seg_data) : seg_data(_seg_data) { }

  void data_ref(segment_data *_seg_data) {
    seg_data = _seg_data;
  }

  void add_keys(key key1) {
    unique_keys.insert(key1);
  }

  template<typename index_iterator>
  void create_tree(index_iterator idx_start, index_iterator idx_end) {

    keys.resize(unique_keys.size());
    std::copy(unique_keys.begin(), unique_keys.end(), keys.begin());
    std::sort(keys.begin(), keys.end());
    tree.resize(keys.size() + 1);

    create_tree(1, keys.begin(), keys.end());

    for (auto itr = idx_start; itr != idx_end; ++itr) {
      register_edge(itr);
    }

    for (auto &node: tree) {
      if (node.left.size() == 0) continue;
      std::sort(node.left.begin(), node.left.end(), [this](const index_type &a, const index_type &b) {
	  return (*(this->seg_data))[a].start_time() < (*(this->seg_data))[b].start_time();
	});

      std::sort(node.right.begin(), node.right.end(), [this](const index_type &a, const index_type &b) {
	  return (*(this->seg_data))[a].end_time() < (*(this->seg_data))[b].end_time();
	});
    }
  }

  void print_tree() {
    for (auto &node: tree) {
      std::cout << node.k << " : left -> ( ";
      for (auto &idx : node.left)
	std::cout << idx << " ";
      std::cout << ") , right -> ( ";
      for (auto &idx : node.right) {
	std::cout << idx << " ";
      }
      std::cout << ")\n";
    }
  }


  /**
   * typename Operator:
   *  --> Function Object : (interval pointer iterator begin, interval pointer iterator end) -> bool {}
   *  --> Object with bool operator()(interval pointer iterator begin, interval pointer iterator end) overloading
   *
   *  interval pointer iterator accesses all intervals within [k_start, k_end] at this tree node
   */
  template<typename Operator>
  void query(std::size_t i, key k_start, key k_end, Operator &op) {
    if (i >= tree.size()) return;
    if (tree[i].k >= k_start && tree[i].k <= k_end) {
      bool next = op(tree[i].left.begin(), tree[i].left.end());
      if( next ) {
	query(i << 1, k_start, k_end, op);
	query((i << 1) + 1, k_start, k_end, op);
      }
    } else if (tree[i].k < k_start) {
      bool next = true;
      if (tree[i].right.size() > 0) {
	auto itr = std::lower_bound(tree[i].right.begin(), tree[i].right.end(), k_start,
				    [this](const index_type &a, const key &b) -> bool {
				      return (*(this->seg_data))[a].end_time() < b;
				    });
	next = op(itr, tree[i].right.end());
      }
      if(next) query((i << 1) + 1, k_start, k_end, op);
    } else if (tree[i].k > k_end) {
      bool next = true;
      if (tree[i].left.size() > 0) {
	auto itr = std::upper_bound(tree[i].left.begin(), tree[i].left.end(), k_end,
				    [this](const key &a, const index_type &b) -> bool {
				      return a < (*(this->seg_data))[b].start_time();
				    });
	next = op(tree[i].left.begin(), itr);
      }
      if(next) query((i << 1), k_start, k_end, op);
    }
  }

  void query_ith( key k_start, key k_end, std::size_t index_i, index_type &index) {
    auto op=[&index_i, &index](typename std::vector<index_type>::iterator begin, typename std::vector<index_type>::iterator end) -> bool {
      if( (end - begin) > index_i ) {
	index = *(begin + index_i);
	return false;
      }else {
	index_i = index_i - (end - begin);
	return true;
      }
    };
    query(1, k_start, k_end, op);
  }

  void query_size( key k_start, key k_end, std::size_t& size) {
    auto op = [&size](typename std::vector<index_type>::iterator begin, typename std::vector<index_type>::iterator end) -> bool{
      size += (end - begin);
      return true;
    };
    query(1, k_start, k_end, op);
  }

  static bool find_intersection(std::pair<key, key> interval1, std::pair<key, key> interval2,
				std::pair<key, key> &intersect_interval) {
    intersect_interval = std::make_pair(std::max(interval1.first, interval2.first),
					std::min(interval1.second, interval2.second));
    return intersect_interval.first <= intersect_interval.second;
  }
private:
  void register_edge(index_type index) {
    std::size_t i = 1;
    while (i < tree.size()) {
      if (tree[i].k >= (*seg_data)[index].start_time() && tree[i].k <= (*seg_data)[index].end_time()) {
	tree[i].left.push_back(index);
	tree[i].right.push_back(index);
	return;
      } else if (tree[i].k < (*seg_data)[index].start_time()) { //go right
	i = (i << 1) + 1;
      } else if (tree[i].k > (*seg_data)[index].end_time()) { // go left
	i = (i << 1);
      } else {
	std::cerr << "Logic error!!!\n";
      }
    }
  }

  template<typename key_iterator>
  void create_tree(std::size_t i, key_iterator key_begin, key_iterator key_end) {
    if (i >= tree.size()) return;
    auto key_mid = key_begin + ((key_end - key_begin) >> 1);
    tree[i] = Node(*key_mid);
    create_tree(i << 1, key_begin, key_mid);
    create_tree((i << 1) + 1, key_mid + 1, key_end);
  }

  std::unordered_set<key> unique_keys;
  std::vector<key> keys;
  std::vector<Node> tree;
};

#endif //INTERVAL_TREE_HPP
