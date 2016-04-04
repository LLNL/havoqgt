#ifndef _SEGMENT_TREE_HPP
#define _SEGMENT_TREE_HPP

#include <fstream>
#include <vector>
#include <limits>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <havoqgt/distributed_db.hpp>

using clock_type = std::chrono::high_resolution_clock;
using time_point_t = std::chrono::time_point<clock_type>;

template<typename T> T NegInf() { return std::numeric_limits<T>::min(); }
template<typename T> T PosInf() { return std::numeric_limits<T>::max(); }

template<typename T>
using interval_type = std::pair< T, T>; 

template<typename key, typename index_type>
  struct Node {
  typedef interval_type<key> interval;

  std::size_t size;
  interval range;
  std::vector<index_type> indexes;
  
  Node() : size(0), range( std::make_pair(PosInf<key>(), PosInf<key>())) { }
  
  Node(key start, key end, std::size_t _size) : size(_size), range( std::make_pair(start, end)) { }

  Node& operator=(const Node& node) {
    size = node.size;
    range = std::make_pair( node.range.first, node.range.second );
    indexes.clear();
    for( auto index: node.indexes)
      indexes.push_back(index);
    return *this;
  }

  };

template <typename key, typename index_type>
class segment_tree {
public:
  typedef interval_type<key> interval; 
  
  segment_tree() : size(0) {}

  segment_tree(std::size_t _size) : size(_size) { }
 
  template<typename sort_predicate>
  void balanced_tree(sort_predicate sort, std::vector<key>& temp) {
    //O(n logn )
    //std::cout << "Inside the balanced tree" << std::endl;
    std::sort( temp.begin(), temp.end(), sort );
    std::vector<key> uniq_temp;
    std::size_t curr = 0;
    uniq_temp.push_back(temp[curr]);
    for(std::size_t c = 1; c < temp.size(); c++) {
      if( uniq_temp[curr] != temp[c] ) {
	uniq_temp.push_back( temp[c] );
	curr++;
      }
    }

#if 0
    std::copy( temp.begin(), temp.end(), std::ostream_iterator<key>(std::cout, " "));
    std::cout << std::endl;
#endif

#if 0
    std::copy( uniq_temp.begin(), uniq_temp.end(), std::ostream_iterator<key>(std::cout, " "));
    std::cout << std::endl;
#endif

    std::size_t n = uniq_temp.size();
    uint32_t count = 0;
    while( n != 0 ) {
      count++;
      n = n >> 1;
    }
    n = ((std::size_t)1) << count;
    size = 4 * n;
    arr.resize( 4 * n );
    
    //std::cout << "Prepopulating elementary intervals" << std::endl;
    //O(n)
    std::size_t i, k;
    key prev = NegInf<key>();
    for(i = size / 2, k = 0; k < uniq_temp.size(); k++,i+=2 ) {
      arr[i]     = Node<key, index_type>( prev        , uniq_temp[k], 0);
      arr[i + 1] = Node<key, index_type>( uniq_temp[k], uniq_temp[k], 0);
      prev = uniq_temp[k];
    }

    //std::cout << "Filling up the tree" << std::endl;
    //O(n)
    std::size_t arr_size = 2 * uniq_temp.size();
    for( i = size/4 ; i >= 1; i /= 2 ) { // i >= 1
      for(k = 0; k < arr_size; k+=2 ) {
	arr[ i + k/2 ] = Node<key, index_type>( arr[2 * i + k].range.first, arr[ 2 * i + k + 1].range.second, 0);
      }
      arr_size = ( arr_size + 1 ) >> 1;
    }
  }

  //O(log n)
  void add( std::size_t i, interval range, index_type index) {
    auto &a = arr[i];
    if( range.first <= (a.range.first + 1) && range.second >= a.range.second ) {
      a.indexes.push_back(index);
      a.size++;
      return;
    } else {
      std::size_t d_i = 2 * i ;
      if( d_i >= arr.size() ) return;
      key c = arr[d_i].range.second;
      if( c >= range.second ) {
	//go left
	add( 2 * i , range, index);
      } else if( range.first > c ) {
	// go right
	add(2 * i + 1, range, index);
      } else {
	//split
	//add( 2 * i, std::make_pair(range.first, c ), index);
	//add( 2 * i + 1, std::make_pair( arr[d_i + 1].range.first, range.second ), index ); 
	add( d_i, range, index); add( d_i + 1, range, index);
      }
      
    }
  }

  void add( interval range, index_type index) {
    add(1, range, index);
  }

  //O( logn + sizeof(indexes) ) // Output Sensitive
  void query( key k, std::size_t i, std::vector<index_type>& indexes) {
    if( i >= arr.size() ) return;
    if( k <= arr[i].range.first || k > arr[i].range.second) return;
    else {
      for( auto index : arr[i].indexes ) indexes.push_back(index);
      query(k, 2 * i, indexes );
      query(k, 2 * i + 1, indexes);
    }
  }

  //O( log n )
  void query( key k, std::size_t i, uint64_t& size) {
    //    std::cout << arr.size() << " vs " << i << std::endl;
    if( i >= arr.size() ) return;
    if( k <= arr[i].range.first || k > arr[i].range.second) return;
    else {
      size += arr[i].size;
      query(k, 2*i, size);
      query(k, 2*i + 1, size);
    }
  }

  //O (log n)
  bool query_ith(key k, std::size_t itr, std::size_t i_th, index_type& index) {
    if( itr >= arr.size() ) return false;
    if( k <= arr[itr].range.first || k> arr[itr].range.second) return false;
    else {
      if( i_th < arr[itr].size ) {
	index = arr[itr].indexes[i_th];
	return true;
      } else {
	 bool found = false;
	 found =  query_ith( k, 2 * itr, i_th - arr[itr].size, index );
	 if( !found )
	  found = query_ith( k, 2 * itr + 1, i_th - arr[itr].size, index );
	 return found;
      }
    }
  }

  void query( key k, uint64_t& size) {
    query(k , 1, size);
  }

  void query( key k, std::vector<index_type>& indexes) {
    query(k, 1, indexes );
  }


  bool query_ith( key k, uint64_t i_th, index_type& index) {
    return query_ith(k, 1, i_th, index );
  }

  void print_tree(std::ofstream &ofs) {
    for( std::size_t i =  1; i < arr.size(); i *= 2 ) {
      for( std::size_t k = i; k < 2 * i; k++) {
	ofs << " ( "<< arr[k].range.first << "," << arr[k].range.second << ") [ ";
		std::copy( arr[k].indexes.begin(), arr[k].indexes.end(), std::ostream_iterator<index_type>(ofs, "|" ) );
		ofs << "]" << std::endl;
      }
    }
  }
  
  segment_tree& operator=(const segment_tree& tree) {
    size = tree.size;
    arr.clear();
    for(auto node : tree.arr) {
      arr.push_back(node);
    }
    return *this;
  }
private:
  std::size_t size;
  std::vector<Node<key, index_type>> arr;
};

#endif
