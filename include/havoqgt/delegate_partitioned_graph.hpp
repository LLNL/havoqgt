/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED
#define HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED

#include <havoqgt/mpi.hpp>
//#include <havoqgt/distributed_edge_list.hpp>
#include <havoqgt/detail/iterator.hpp>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/map.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <stdint.h>
#include <utility>
#include <limits>

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/vector.hpp>


namespace havoqgt {
namespace mpi {

namespace bip = boost::interprocess;

/**
 * Delegate partitioned graph using MPI for communication.
 *
 * @todo Test using simple deterministic patterns.
 * @todo Add edge_data
 * @todo Make vertex_iterator a random access iterator
 * @todo Add invalid bit or state to vertex_locator
 * @todo Verify low-degree CSR creation:  ipp line 167
 * @todo Boostify controller locator
 */
template <typename SegementManager>
class delegate_partitioned_graph {
 public:

 	template<typename T>
 	using SegmentAllocator = bip::allocator<T, SegementManager>;

  /// Object that uniquely locates a vertex, both MPI rank and local offset.
  class vertex_locator {
<<<<<<< HEAD
  public:
=======
   public:
>>>>>>> develop
    vertex_locator() {
      m_is_delegate  = 0;
      m_is_bcast     = 0;
      m_is_intercept = 0;
<<<<<<< HEAD
      m_owner_dest   = std::numeric_limits<uint64_t>::max();
      m_local_id     = std::numeric_limits<uint64_t>::max();
=======
      m_owner_dest   = std::numeric_limits<uint32_t>::max();
      m_local_id     = std::numeric_limits<uint64_t>::max();
    }

    bool is_valid() {
    	struct Temp {
    		uint64_t local_id : 39;
    		uint32_t owner_dest : 20;
    	};

    	Temp conv;
    	conv.local_id = std::numeric_limits<uint64_t>::max();
    	conv.owner_dest = std::numeric_limits<uint64_t>::max();


    	return (m_local_id != conv.local_id || m_owner_dest != conv.owner_dest);
>>>>>>> develop
    }

    bool is_delegate() const { return m_is_delegate == 1;}
    uint32_t owner() const { return m_owner_dest; }
    void set_dest(uint32_t dest) { m_owner_dest = dest; assert(m_owner_dest == dest);}
    uint64_t local_id() const { return m_local_id;}
    bool is_equal(const vertex_locator x) const;
    uint32_t get_bcast() const { return m_is_bcast ; }
    void set_bcast(uint32_t bcast) { m_is_bcast = bcast; }
    bool is_intercept() const { return m_is_intercept == 1;}
    void set_intercept(bool intercept) { m_is_intercept = intercept; }

    friend bool operator==(const vertex_locator& x,
                           const vertex_locator& y) {return x.is_equal(y); }
    friend bool operator<(const vertex_locator& x,
                           const vertex_locator& y) {
    	if (x.m_is_delegate == y.m_is_delegate) {
    		if (x.m_owner_dest == y.m_owner_dest) {
	    		return x.m_local_id < y.m_local_id;
	    	}
	    	else {
	    		return x.m_owner_dest < y.m_owner_dest;
	    	}

    	} else {
    		return x.m_is_delegate < y.m_is_delegate;
    	}

<<<<<<< HEAD
=======
    }

>>>>>>> develop
    friend bool operator!=(const vertex_locator& x,
                           const vertex_locator& y) {return !(x.is_equal(y)); }
   private:
    friend class delegate_partitioned_graph;
    unsigned int m_is_delegate  : 1;
<<<<<<< HEAD
    unsigned int m_is_bcast     : 3;
    unsigned int m_is_intercept : 1;
    unsigned int m_owner_dest   : 20;
    uint64_t     m_local_id     : 39;
=======
		unsigned int m_is_bcast     : 3;
		unsigned int m_is_intercept : 1;
	  uint32_t m_owner_dest   : 20;
    uint64_t m_local_id     : 39;

>>>>>>> develop
    vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest);
  };

  /// Edge Iterator class for delegate partitioned graph
  class edge_iterator {
   public:
    edge_iterator()
    	: m_ptr_graph(NULL) {};

    edge_iterator& operator++();
    edge_iterator operator++(int);

    bool is_equal(const edge_iterator& x) const;

    friend bool operator==(const edge_iterator& x,
                           const edge_iterator& y) {return x.is_equal(y); }

    friend bool operator!=(const edge_iterator& x,
                           const edge_iterator& y) {return !(x.is_equal(y)); }

    vertex_locator source() const { return m_source; }
    vertex_locator target() const;
<<<<<<< HEAD
  protected:
=======

   protected:
>>>>>>> develop
    friend class delegate_partitioned_graph;
    template <typename T1, typename T2> friend class edge_data;
    edge_iterator(vertex_locator source, uint64_t edge_offset,
                  const delegate_partitioned_graph* const pgraph);

    vertex_locator                          m_source;
    uint64_t                                m_edge_offset;
    const delegate_partitioned_graph* const m_ptr_graph;
  };

  /// Vertex Iterator class for delegate partitioned graph
  class vertex_iterator
    : public std::iterator<std::input_iterator_tag, vertex_locator, ptrdiff_t,
                           const vertex_locator* const, const vertex_locator&> {
   public:
    vertex_iterator()
    	: m_ptr_graph(NULL) {};
    vertex_iterator& operator++();
    vertex_iterator operator++(int);

    bool is_equal(const vertex_iterator& x) const;

    friend bool operator==(const vertex_iterator& x,
                           const vertex_iterator& y) { return x.is_equal(y); }

    friend bool operator!=(const vertex_iterator& x,
                           const vertex_iterator& y) { return !(x.is_equal(y)); }


    const vertex_locator& operator*()        const { return m_locator; }
    const vertex_locator* const operator->() const { return &m_locator; }

   private:
    friend class delegate_partitioned_graph;

    vertex_iterator(uint64_t index, const delegate_partitioned_graph*  pgraph);
    void update_locator();

    const delegate_partitioned_graph*  m_ptr_graph;
    uint64_t                                m_owned_vert_index;
    vertex_locator                          m_locator;
  };

  /// Vertex Data storage
<<<<<<< HEAD
  template <typename T, typename SegmentManager>
=======
  template <typename T, typename SegManagerOther>
>>>>>>> develop
  class vertex_data {
   public:
   	typedef T value_type;
    vertex_data() {}

    T&       operator[] (const vertex_locator& locator);
    const T& operator[] (const vertex_locator& locator) const;

    void reset(const T& r) {
      for(size_t i=0; i<m_owned_vert_data.size(); ++i) {
        m_owned_vert_data[i] = r;
      }
      for(size_t i=0; i<m_delegate_data.size(); ++i) {
        m_delegate_data[i] = r;
      }
    }

    void all_reduce() {
      std::vector<T> tmp_in(m_delegate_data.begin(), m_delegate_data.end());
      std::vector<T> tmp_out(tmp_in.size(), 0);
      mpi_all_reduce(tmp_in, tmp_out, std::plus<T>(), MPI_COMM_WORLD);
      std::copy(tmp_out.begin(), tmp_out.end(), m_delegate_data.begin());
    }

  //private:
  //  friend class delegate_partitioned_graph;
    vertex_data(uint64_t owned_data_size, uint64_t delegate_size,
    	SegManagerOther* segment_manager);
    vertex_data(uint64_t owned_data_size, uint64_t delegate_size,
    	const T& init, SegManagerOther* segment_manager);

   private:
    bip::vector<T, bip::allocator<T, SegManagerOther>  > m_owned_vert_data;
    bip::vector<T, bip::allocator<T, SegManagerOther>  > m_delegate_data;
  };

  /// Edge Data storage
<<<<<<< HEAD
  template <typename T, typename SegmentManager>
=======
  template <typename T, typename SegManagerOther>
>>>>>>> develop
  class edge_data {
   public:
   	typedef typename bip::vector< T, bip::allocator<T, SegManagerOther> >
   			::iterator iterator;
    typedef T value_type;

    edge_data() {}

    T&       operator[] (const edge_iterator& itr);
    const T& operator[] (const edge_iterator& itr) const;

    void reset(const T& r) {
      for(size_t i=0; i<m_owned_edge_data.size(); ++i) {
        m_owned_edge_data[i] = r;
      }
      for(size_t i=0; i<m_delegate_edge_data.size(); ++i) {
        m_delegate_edge_data[i] = r;
      }
    }

    iterator delegate_begin() { return m_delegate_edge_data.begin(); }
    iterator delegate_end()   { return m_delegate_edge_data.end(); }
    iterator owned_begin()    { return m_owned_edge_data.begin(); }
    iterator owned_end()      { return m_owned_edge_data.end(); }

  //private:
  //  friend class delegate_partitioned_graph;
    edge_data(uint64_t owned_size, uint64_t delegate_size,
    		SegManagerOther* sm);
    edge_data(uint64_t owned_size, uint64_t delegate_size, const T& init,
    		SegManagerOther* sm);

   private:
    bip::vector< T, bip::allocator<T, SegManagerOther> >
    		m_owned_edge_data;
    bip::vector< T, bip::allocator<T, SegManagerOther> >
    		m_delegate_edge_data;
  };

<<<<<<< HEAD
  typedef typename std::vector<vertex_locator>::const_iterator controller_iterator;

  /// Constructor that initializes given and unsorted sequence of edges
  template <typename Container>
  delegate_partitioned_graph(Arena& arena,
=======
  typedef typename bip::vector<vertex_locator, SegmentAllocator<vertex_locator> >
  		::const_iterator controller_iterator;

  /// Constructor that initializes given and unsorted sequence of edges
  template <typename Container>
  delegate_partitioned_graph(const SegmentAllocator<void>& seg_allocator,
>>>>>>> develop
                             MPI_Comm mpi_comm,
                             Container& edges,
                             uint64_t delegate_degree_threshold);






  /// Converts a vertex_locator to the vertex label
  uint64_t locator_to_label(vertex_locator locator) const;

  /// Converts a vertex label to a vertex_locator
  vertex_locator label_to_locator(uint64_t label) const;

  /// Returns a begin iterator for edges of a vertex
  edge_iterator edges_begin(vertex_locator locator) const;

  /// Returns an end iterator for edges of a vertex
  edge_iterator edges_end(vertex_locator locator) const;

  /// Returns the degree of a vertex
  uint64_t degree(vertex_locator locator) const;

  /// Returns the local degree of a vertex
  uint64_t local_degree(vertex_locator locator) const;

  /// Returns a begin iterator for all local vertices
  vertex_iterator vertices_begin() const;

  /// Returns an end iterator for all local vertices
  vertex_iterator vertices_end() const;

  /// Tests if vertex label is a delegate
  bool is_label_delegate(uint64_t label) const;

  /// Creates vertex_data of type T
  template <typename T, typename SegManagerOther>
  vertex_data<T, SegManagerOther>* create_vertex_data(
  		SegManagerOther*,
  		const char *obj_name = nullptr) const;

  /// Creates vertex_data of type T, with initial value
  template <typename T, typename SegManagerOther>
  vertex_data<T, SegManagerOther>* create_vertex_data(
  		const T& init, SegManagerOther*,
  		const char *obj_name = nullptr) const;

  /// Creates edge_data of type T
  template <typename T, typename SegManagerOther>
  edge_data<T, SegManagerOther>* create_edge_data(
  		SegManagerOther*,
  		const char *obj_name = nullptr) const;

  /// Creates edge_data of type T, with initial value
  template <typename T, typename SegManagerOther>
  edge_data<T, SegManagerOther>* create_edge_data(
  		const T& init, SegManagerOther*,
  		const char *obj_name = nullptr) const;

  size_t num_local_vertices() const {
  	return m_owned_info.size();
  }

  size_t num_delegates() const { return m_delegate_degree.size(); }

  uint32_t master(const vertex_locator& locator) const {
  	return locator.m_local_id % m_mpi_size;
  }

  controller_iterator controller_begin() const {
  	return m_controller_locators.begin();
  }

  controller_iterator controller_end()   const {
  	return m_controller_locators.end();
  }

//private:
  /// Synchronizes hub set amongst all processes.
  void sync_global_hub_set(const boost::unordered_set<uint64_t>& local_hubs,
                           boost::unordered_set<uint64_t>& global_hubs,
                           bool local_change, MPI_Comm mpi_comm);


  template <typename InputIterator>
	void count_high_degree_transpose(MPI_Comm mpi_comm,
               InputIterator unsorted_itr,
               InputIterator unsorted_itr_end,
               boost::unordered_set<uint64_t>& global_hub_set,
               std::vector<uint64_t>& high_count_per_rank);


	template <typename InputIterator>
	void count_low_degree(MPI_Comm mpi_comm,
               InputIterator unsorted_itr,
               InputIterator unsorted_itr_end,
               boost::unordered_set<uint64_t>& global_hub_set,
               std::vector<uint64_t>& vertex_low_degree_count,
               uint64_t delegate_degree_threshold);

	template <typename InputIterator>
	void partition_low_degree(MPI_Comm mpi_comm,
               InputIterator unsorted_itr,
               InputIterator unsorted_itr_end,
               boost::unordered_set<uint64_t>& global_hub_set,
               std::deque<std::pair<uint64_t, uint64_t> >& edges_low);

	template <typename InputIterator>
	void partition_high_degree(MPI_Comm mpi_comm,
               InputIterator unsorted_itr,
               InputIterator unsorted_itr_end,
               boost::unordered_set<uint64_t>& global_hub_set,
               std::deque<std::pair<uint64_t, uint64_t> >& edges_high,
               std::deque<std::pair<uint64_t, uint64_t> >& edges_high_overflow,
               std::map<int, uint64_t>& overflow_schedule);


  /// Stores information about owned vertices
  class vert_info {
  public:
    vert_info(bool in_is_delegate, uint64_t in_delegate_id,
              uint64_t in_low_csr_idx);

    uint32_t is_delegate :  1;
    uint32_t delegate_id : 24;
    uint64_t low_csr_idx : 39;

    inline bool operator==(const vert_info& other){
    	return (is_delegate == other.is_delegate && delegate_id == other.delegate_id
    		&& low_csr_idx == other.low_csr_idx);
    }

    inline bool operator!=(const vert_info& other){
    	return !(*this == other);
    }
  };

  int m_mpi_size;
  int m_mpi_rank;

  uint64_t m_max_vertex;

	bip::vector<vert_info, SegmentAllocator<vert_info>> m_owned_info;
  bip::vector<vertex_locator, SegmentAllocator<vertex_locator>> m_owned_targets;

	// Delegate Storage
  uint64_t m_delegate_degree_threshold;
<<<<<<< HEAD
  bip::vector<uint64_t, typename Arena::template allocator<uint64_t>::type > m_delegate_info;
  bip::vector<uint64_t, typename Arena::template allocator<uint64_t>::type > m_delegate_degree;
  bip::vector<uint64_t, typename Arena::template allocator<uint64_t>::type > m_delegate_label;
  boost::unordered_map<uint64_t, vertex_locator, boost::hash<uint64_t>,
          std::equal_to<uint64_t>, typename Arena::template allocator<std::pair<uint64_t,vertex_locator> >::type > m_map_delegate_locator;
  bip::vector<vertex_locator, typename Arena::template allocator<vertex_locator>::type > m_delegate_targets;
  std::vector<vertex_locator> m_controller_locators;
=======
  bip::vector< uint64_t, SegmentAllocator<uint64_t> > m_delegate_info;
  bip::vector< uint64_t, SegmentAllocator<uint64_t> > m_delegate_degree;
  bip::vector< uint64_t, SegmentAllocator<uint64_t> > m_delegate_label;
  bip::vector< vertex_locator, SegmentAllocator<vertex_locator> >
  		m_delegate_targets;

  //Note: BIP only contains a map, not an unordered_map object.
  boost::unordered_map<
  		uint64_t, vertex_locator, boost::hash<uint64_t>, std::equal_to<uint64_t>,
  		SegmentAllocator< std::pair<uint64_t,vertex_locator> >
     > m_map_delegate_locator;

  bip::vector<vertex_locator, SegmentAllocator<vertex_locator> >
  	m_controller_locators;




  inline bool operator==(const delegate_partitioned_graph<SegementManager>&
  			other)
  {
 	  if (m_mpi_size != other.m_mpi_size) {
  		return false;
  	} else if (m_mpi_rank != other.m_mpi_rank) {
  		return false;
  	} else if (m_delegate_degree_threshold != other.m_delegate_degree_threshold) {
  		return false;
  	}

  	{
  		size_t size = m_owned_info.size();
  		if (size != other.m_owned_info.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_owned_info[i] != other.m_owned_info[i]) {
  					return false;
  				}
  			}
  		}
		}

		{
  		size_t size = m_owned_targets.size();
  		if (size != other.m_owned_targets.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_owned_targets[i] != other.m_owned_targets[i]) {
  					return false;
  				}
  			}
  		}
		}


		{
  		size_t size = m_delegate_info.size();
  		if (size != other.m_delegate_info.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_delegate_info[i] != other.m_delegate_info[i]) {
  					return false;
  				}
  			}
  		}
		}



		{
  		size_t size = m_delegate_degree.size();
  		if (size != other.m_delegate_degree.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_delegate_degree[i] != other.m_delegate_degree[i]) {
  					return false;
  				}
  			}
  		}
		}


		{
  		size_t size = m_delegate_label.size();
  		if (size != other.m_delegate_label.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_delegate_label[i] != other.m_delegate_label[i]) {
  					return false;
  				}
  			}
  		}
		}

		{
  		size_t size = m_delegate_targets.size();
  		if (size != other.m_delegate_targets.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_delegate_targets[i] != other.m_delegate_targets[i]) {
  					return false;
  				}
  			}
  		}
		}


		{
  		size_t size = m_controller_locators.size();
  		if (size != other.m_controller_locators.size()) {
  			return false;
  		} else {
  			for (size_t i = 0; i < size; i++) {
  				if (m_controller_locators[i] != other.m_controller_locators[i]) {
  					return false;
  				}
  			}
  		}
		}

		{
			if (m_map_delegate_locator != other.m_map_delegate_locator) {
				return false;
			}
		}

		return true;
	}

	inline bool operator!=(const delegate_partitioned_graph<SegementManager>&
				other){
		return !(*this == other);
	}
>>>>>>> develop
};

/// Frees the container of edges
template <typename Container>
void free_edge_container(Container &edges) {};
<<<<<<< HEAD
=======

template<>
void free_edge_container<std::vector<std::pair<uint64_t, uint64_t> > >(std::vector<std::pair<uint64_t, uint64_t> > &edges){
	std::vector< std::pair<uint64_t, uint64_t> >empty(0);
  edges.swap(empty);
};
>>>>>>> develop

template<>
void free_edge_container<std::vector<std::pair<uint64_t, uint64_t> > >(std::vector<std::pair<uint64_t, uint64_t> > &edges){
	std::vector< std::pair<uint64_t, uint64_t> >empty(0);
  edges.swap(empty);
};

} // namespace mpi
} // namespace havoqgt

#include <havoqgt/impl/delegate_partitioned_graph.ipp>

#endif //HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED
