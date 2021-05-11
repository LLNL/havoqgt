// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
#define HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED

#include <algorithm>
#include <boost/container/deque.hpp>
#include <fstream>
#include <havoqgt/visitor_queue.hpp>
#include <utility>
//#include <boost/container/flat_map.hpp>
#include <map>
#include <unordered_map>

namespace havoqgt {

struct dogr_edge {
  uint32_t target_degree       = 0;
  uint32_t edge_triangle_count = 0;
  bool     is_j                = false;
};

template <typename Visitor>
class lifo_queue {
 protected:
  std::vector<Visitor> m_data;

 public:
  lifo_queue() {}

  bool push(Visitor const& task) {
    m_data.push_back(task);
    return true;
  }

  void pop() { m_data.pop_back(); }

  Visitor const& top()  // const
  {
    return m_data.back();
  }

  size_t size() const {
    return m_data.size();
    ;
  }

  bool empty() const { return m_data.empty(); }

  void clear() { m_data.clear(); }
};

template <typename Graph>
class core2_visitor {
 public:
  typedef core2_visitor<Graph>           my_type;
  typedef typename Graph::vertex_locator vertex_locator;

  core2_visitor() : vertex(), init(true) {}

  core2_visitor(vertex_locator v) : vertex(v), init(true) {}

  core2_visitor(vertex_locator v, bool _init) : vertex(v), init(_init) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if (vertex.is_delegate()) {
      if (!vertex.is_delegate_master()) {
        return true;
      }
    }
    if (std::get<1>(alg_data)[vertex]) {
      --(std::get<0>(alg_data)[vertex]);

      if (std::get<0>(alg_data)[vertex] < 2) {
        // remove from 2 core
        std::get<1>(alg_data)[vertex] = false;
        std::get<0>(alg_data)[vertex] = 0;
        return true;
      }
    }

    // Retrun true if alive
    // return std::get<1>(alg_data)[vertex];
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (std::get<1>(alg_data)[vertex]) {
      //   --(std::get<0>(alg_data)[vertex]);
      if (std::get<0>(alg_data)[vertex] < 2) {
        // remove from 2 core
        std::get<1>(alg_data)[vertex] = false;
        std::get<0>(alg_data)[vertex] = 0;
        for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
             ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type        new_visitor(neighbor, false);
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return true;
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (init) {
      if (std::get<0>(alg_data)[vertex] < 2) {
        // remove from 2 core
        std::get<1>(alg_data)[vertex] = false;
        std::get<0>(alg_data)[vertex] = 0;
        for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
             ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type        new_visitor(neighbor, false);
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return true;
    }

    if (std::get<1>(alg_data)[vertex]) {
      std::cerr << "LOGIC ERROR" << std::endl;
      exit(-1);
    }
    // if(std::get<1>(alg_data)[vertex]) {
    //  --(std::get<0>(alg_data)[vertex]);
    //  if(std::get<0>(alg_data)[vertex] < 2) {
    //    //remove from 2 core
    //    std::get<1>(alg_data)[vertex] = false;
    //    std::get<0>(alg_data)[vertex] = 0;
    for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
         ++eitr) {
      vertex_locator neighbor = eitr.target();
      my_type        new_visitor(neighbor, false);
      vis_queue->queue_visitor(new_visitor);
    }
    //  }
    //  return true;
    //}
    // return false;
    return true;
  }

  friend inline bool operator>(const core2_visitor& v1,
                               const core2_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const core2_visitor& v1,
                               const core2_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  bool           init;
};

template <typename Graph>
class directed_core2 {
 public:
  typedef directed_core2<Graph>          my_type;
  typedef typename Graph::vertex_locator vertex_locator;

  directed_core2() : vertex(), init(true) {}

  directed_core2(vertex_locator v) : vertex(v), init(true) {}

  directed_core2(vertex_locator v, vertex_locator _from, uint32_t _from_degree)
      : vertex(v), from_label(_from), from_degree(_from_degree), init(false) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    // if(std::get<0>(alg_data)[vertex] >= 2) {
    if (from_degree >= std::get<2>(alg_data).degree(
                           vertex) /*std::get<0>(alg_data)[vertex]*/) {
      // previously returned true, but changing here --- return true;
      if (vertex.is_delegate()) {
        if (!vertex.is_delegate_master()) {
          return true;
        }
      }
      if (/*std::get<0>(alg_data)[vertex]*/ std::get<2>(alg_data).degree(
              vertex) < 2)
        return false;
      // only here should be low-degree & masters
      if (from_degree > std::get<2>(alg_data).degree(
                            vertex) /*std::get<0>(alg_data)[vertex]*/
          || (from_degree ==
                  std::get<2>(alg_data).degree(
                      vertex) /*std::get<0>(alg_data)[vertex]*/
              && vertex < from_label)) {
        auto vv = vertex;
        vv.set_bcast(0);
        vv.set_intercept(0);
        auto fl = from_label;
        fl.set_bcast(0);
        fl.set_intercept(0);

        std::get<1>(alg_data)[vv][fl].target_degree = from_degree;
        // std::get<1>(alg_data)[vv].add(fl, dogr_edge{from_degree, 0, 0});
        // attempting to replicate dogr graph
        if (vertex.is_delegate()) {
          if (vertex.is_delegate_master()) {
            return true;
          }
        }
      }
    }
    //}
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (/*std::get<0>(alg_data)[vertex]*/ g.degree(vertex) >= 2) {
      // if in 2core, send degree to neighbors
      uint32_t my_degree = /*std::get<0>(alg_data)[vertex];*/ g.degree(vertex);
      for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
           ++eitr) {
        vertex_locator neighbor = eitr.target();
        if (neighbor == vertex) continue;
        my_type new_visitor(neighbor, vertex, my_degree);
        vis_queue->queue_visitor(new_visitor);
      }
      return true;
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (init) {
      if (/*std::get<0>(alg_data)[vertex]*/ g.degree(vertex) >= 2) {
        // if in 2core, send degree to neighbors
        uint32_t my_degree =
            /*std::get<0>(alg_data)[vertex];*/ g.degree(vertex);
        for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
             ++eitr) {
          vertex_locator neighbor = eitr.target();
          if (neighbor == vertex) continue;
          my_type new_visitor(neighbor, vertex, my_degree);
          vis_queue->queue_visitor(new_visitor);
        }
        return true;
      } else {
        return false;
      }
    } else {
      auto vv = vertex;
      vv.set_bcast(0);
      vv.set_intercept(0);
      auto fl = from_label;
      fl.set_bcast(0);
      fl.set_intercept(0);

      std::get<1>(alg_data)[vv][fl].target_degree = from_degree;
      // std::get<1>(alg_data)[vv].add(fl, dogr_edge{from_degree, 0, 0});
      return true;
    }

    //    if(from_degree > g.degree(vertex) /*std::get<0>(alg_data)[vertex]*/ ||
    //      (from_degree == g.degree(vertex)/*std::get<0>(alg_data)[vertex]*/ &&
    //      vertex < from_label )) {
    //
    //      auto vv = vertex;
    //      vv.set_bcast(0);
    //      vv.set_intercept(0);
    //      auto fl = from_label;
    //      fl.set_bcast(0);
    //      fl.set_intercept(0);
    //
    //      std::get<1>(alg_data)[vv][fl] = from_degree;
    //    }
    // return false;
  }

  friend inline bool operator>(const directed_core2& v1,
                               const directed_core2& v2) {
    return false;
  }

  friend inline bool operator<(const directed_core2& v1,
                               const directed_core2& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator from_label;
  uint32_t       from_degree;
  bool           init;
};

template <typename Graph>
class core2_wedges {
 public:
  typedef core2_wedges<Graph>            my_type;
  typedef typename Graph::vertex_locator vertex_locator;

  core2_wedges() : vertex() {}

  core2_wedges(vertex_locator v) : vertex(v) {}

  core2_wedges(vertex_locator v, vertex_locator _check_close,
               vertex_locator _from_vertex, bool _do_check_close)
      : vertex(v),
        check_close(_check_close),
        from_vertex(_from_vertex),
        do_check_close(_do_check_close) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if (vertex.is_delegate()) {
      if (vertex.is_delegate_master()) {
        if (do_check_close) {
          ++std::get<2>(alg_data);
          ++std::get<3>(alg_data);
          // std::get<3>(alg_data) <<
          // std::get<4>(alg_data).locator_to_label(vertex) << " " <<
          // check_close
          // << std::endl;
          // std::get<4>(alg_data)[vertex]++;
          if (std::get<0>(alg_data)[vertex].count(check_close) > 0) {
            ++std::get<1>(alg_data);
            // std::get<5>(alg_data)[vertex]++;
            // std::cout << "FOUND" << std::endl;
            std::get<0>(alg_data)[vertex][check_close].edge_triangle_count++;
            std::get<0>(alg_data)[vertex][check_close].is_j = true;
            return true;  /// now returning true to count
          }
        } else {  // already found, edge counting
          if (std::get<0>(alg_data)[vertex].count(check_close) == 0 ||
              std::get<0>(alg_data)[vertex].count(from_vertex) == 0) {
            std::cerr << "Error in edge counting" << std::endl;
          }
        }
        return false;
      } else {
        if (do_check_close) {
          if (std::get<0>(alg_data)[vertex].count(check_close) == 0) {
            ++std::get<4>(alg_data);
            // std::cout << "Skipped via replication" << std::endl;
            // exit(-1);
            return false;
          }
        }
        return true;
      }
    }
    if (do_check_close) {
      ++std::get<2>(alg_data);
      // std::get<4>(alg_data)[vertex]++;
      // std::get<3>(alg_data) << std::get<4>(alg_data).locator_to_label(vertex)
      // << " " << check_close << std::endl;
      if (std::get<0>(alg_data)[vertex].count(check_close) > 0) {
        ++std::get<1>(alg_data);
        // std::get<5>(alg_data)[vertex]++;
        // std::cout << "FOUND" << std::endl;
        std::get<0>(alg_data)[vertex][check_close].edge_triangle_count++;
        std::get<0>(alg_data)[vertex][check_close].is_j = true;
        return true;
      }
    } else {  // already found, edge counting
      if (std::get<0>(alg_data)[vertex].count(check_close) == 0 ||
          std::get<0>(alg_data)[vertex].count(from_vertex) == 0) {
        std::cerr << "Error in edge counting" << std::endl;
      } else {
        std::get<0>(alg_data)[vertex][check_close].edge_triangle_count++;
        std::get<0>(alg_data)[vertex][from_vertex].edge_triangle_count++;
      }
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (std::get<0>(alg_data)[vertex].size() > 1) {
      for (const auto& pair_a : std::get<0>(alg_data)[vertex]) {
        for (const auto& pair_b : std::get<0>(alg_data)[vertex]) {
          if (pair_a.second.target_degree < pair_b.second.target_degree ||
              (pair_a.second.target_degree == pair_b.second.target_degree &&
               pair_a.first < pair_b.first)) {
            my_type new_visitor(pair_a.first, pair_b.first, vertex, true);
            vis_queue->queue_visitor(new_visitor);
            // std::get<3>(alg_data)[vertex]++;
            // fake ++std::get<1>(alg_data);  //fake counting here
          }
        }
      }
      /*for(auto itr_i = std::get<0>(alg_data)[vertex].begin(); itr_i !=
      std::get<0>(alg_data)[vertex].end(); ++itr_i) {
        for(auto itr_j = itr_i+1; itr_j != std::get<0>(alg_data)[vertex].end();
      ++itr_j) {
          auto pair_a = *itr_i;
          auto pair_b = *itr_j;
          if(pair_a.second < pair_b.second ||
            (pair_a.second == pair_b.second && pair_a.first < pair_b.first)) {
              my_type new_visitor(g.label_to_locator(pair_a.first),
      pair_b.first);
              //vis_queue->queue_visitor(new_visitor);
              ++std::get<1>(alg_data);  //fake counting here
          } else {
            my_type new_visitor(g.label_to_locator(pair_b.first), pair_a.first);
            //vis_queue->queue_visitor(new_visitor);
             ++std::get<1>(alg_data);  //fake counting here
          }
        }
      }*/

      //        ++std::get<1>(alg_data);  //fake counting here
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (do_check_close) {
      my_type new_visitor(from_vertex, check_close, vertex, false);
      vis_queue->queue_visitor(new_visitor);
      return false;
    }
    std::cout << "Shoudn't be here" << std::endl;
    exit(-1);
    /*    ++std::get<2>(alg_data);
        std::get<3>(alg_data) << g.locator_to_label(vertex) << " " <<
       check_close << std::endl;
        if(std::get<0>(alg_data)[vertex].count(check_close) > 0) {
          ++std::get<1>(alg_data);
          //std::cout << "FOUND" << std::endl;
        }
        return false;*/
  }

  friend inline bool operator>(const core2_wedges& v1, const core2_wedges& v2) {
    return false;
  }

  friend inline bool operator<(const core2_wedges& v1, const core2_wedges& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator check_close;
  vertex_locator from_vertex;
  bool           do_check_close;
};

void output_degree_distribution(
    const std::map<uint64_t, uint64_t>& local_degree_count, const char* fname) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> send_p_vec(mpi_size);
  std::vector<std::vector<std::pair<uint64_t, uint64_t>>> recv_p_vec(mpi_size);
  for (auto deg_cnt : local_degree_count) {
    send_p_vec[deg_cnt.first % uint64_t(mpi_size)].push_back(deg_cnt);
  }

  mpi_all_to_all(send_p_vec, recv_p_vec, MPI_COMM_WORLD);

  std::map<uint64_t, uint64_t> partitioned_deg_count;
  for (auto& vec_deg_cnt : recv_p_vec) {
    for (auto deg_cnt : vec_deg_cnt) {
      partitioned_deg_count[deg_cnt.first] += deg_cnt.second;
    }
  }

  //
  // Send all to Rank 0 -- this is not efficient way
  send_p_vec.clear();
  send_p_vec.resize(mpi_size);
  recv_p_vec.clear();
  recv_p_vec.resize(mpi_size);
  for (auto vec_deg_cnt : partitioned_deg_count) {
    send_p_vec[0].push_back(vec_deg_cnt);
  }
  mpi_all_to_all(send_p_vec, recv_p_vec, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::vector<std::pair<uint64_t, uint64_t>> all_sorted;
    for (auto& vec_deg_cnt : recv_p_vec) {
      for (auto deg_cnt : vec_deg_cnt) {
        all_sorted.push_back(deg_cnt);
      }
    }
    std::sort(all_sorted.begin(), all_sorted.end());
    std::ofstream outfile(fname);
    for (auto deg_cnt : all_sorted) {
      outfile << deg_cnt.first << "\t" << deg_cnt.second << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
struct vertex_locator_hash {
  template <typename T>
  std::size_t operator()(const T& k) const {
    return k.hash();
  }
};

template <typename K, typename V>
class lazy_flat_map {
  using data = std::pair<K, V>;
  // struct data {
  //   K key;
  //   V value;
  // };  // __attribute__((packed));

 public:
  lazy_flat_map() { m_sorted = true; }

  typename std::vector<data>::iterator begin() {
    do_sort();
    return m_data.begin();
  }
  typename std::vector<data>::iterator end() { return m_data.end(); }

  size_t size() const { return m_data.size(); }

  void add(const K& key, const V& value) {
    m_data.push_back(data{key, value});
    m_sorted = false;
  }

  V& operator[](const K& inkey) {
    do_sort();
    data tmp_key{inkey, V{}};
    auto lb = std::lower_bound(m_data.begin(), m_data.end(), tmp_key,
                               [](const data& left, const data& right) {
                                 return left.first < right.first;
                               });
    if ((*lb).first == inkey) {
      return (*lb).second;
    } else {
      std::cerr << "NOT FOUND" << std::endl;
      exit(-1);
    }
  }

  size_t count(const K& inkey) {
    if (m_data.empty()) return 0;
    do_sort();
    data tmp_key{inkey, V{}};
    auto lb = std::lower_bound(m_data.begin(), m_data.end(), tmp_key,
                               [](const data& left, const data& right) {
                                 return left.first < right.first;
                               });
    if ((*lb).first == inkey) {
      return 1;
    } else {
      return 0;
    }
  }

  void do_sort() {
    if (!m_sorted) {
      m_sorted = true;
      std::sort(m_data.begin(), m_data.end(),
                [](const data& left, const data& right) {
                  return left.first < right.first;
                });
      auto last = std::unique(m_data.begin(), m_data.end(),
                              [](const data& left, const data& right) {
                                return left.first == right.first;
                              });
      m_data.erase(last, m_data.end());
    }
  }

 private:
  bool              m_sorted;
  std::vector<data> m_data;
};

template <typename TGraph>
uint64_t new_triangle_count(TGraph& g, const char* deg_output_fname) {
  typedef TGraph                          graph_type;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  // std::stringstream fname_sv_partition, fname_wc_partition;
  // fname_sv_partition << "/p/lscratchf/havoqgtu/wedge_check_experiment/DOD_"
  // << mpi_rank << "_of_" << mpi_size;
  // fname_wc_partition <<
  // "/p/lscratchf/havoqgtu/wedge_check_experiment/wedges_" << mpi_rank <<
  // "_of_" << mpi_size;
  //
  // std::ofstream ofs_sv_partition(fname_sv_partition.str().c_str());
  // std::ofstream ofs_wc_partition(fname_wc_partition.str().c_str());

  //
  // 1)  Calculate 2core degree.   0 degree means that vertex is not in 2core
  typename graph_type::template vertex_data<
      /*lazy_flat_map*/ std::map<vertex_locator, dogr_edge>,
      std::allocator</*lazy_flat_map*/ std::map<vertex_locator, dogr_edge>>>
      core2_directed(g);
  {
    typename graph_type::template vertex_data<uint32_t,
                                              std::allocator<uint32_t>>
        core2_degree(g);
    {
      typename graph_type::template vertex_data<bool, std::allocator<bool>>
          core2_alive(g);
      core2_alive.reset(true);
      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        core2_degree[*vitr] = g.degree(*vitr);
      }
      for (auto ditr = g.delegate_vertices_begin();
           ditr != g.delegate_vertices_end(); ++ditr) {
        core2_degree[*ditr] = g.degree(*ditr);
      }

      // // /// This computes the 2core
      // double start_time = MPI_Wtime();
      // {
      //   auto alg_data = std::forward_as_tuple(core2_degree, core2_alive);
      //   auto vq = create_visitor_queue<core2_visitor<graph_type>,
      //   lifo_queue>(
      //       &g, alg_data);
      //   vq.init_visitor_traversal();
      // }
      // double end_time = MPI_Wtime();
      // if (mpi_rank == 0) {
      //   std::cout << "2Core time = " << end_time - start_time << std::endl;
      // }

      // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr)
      // {
      //   std::cout << "Vertex " << g.locator_to_label(*vitr) << " has
      //   core2_degree = " << core2_degree[*vitr] << std::endl;
      // }
      // for (auto ditr = g.delegate_vertices_begin(); ditr !=
      // g.delegate_vertices_end(); ++ditr) {
      //   std::cout << "Delegate " << g.locator_to_label(*ditr) << " has
      //   core2_degree = " <<  core2_degree[*ditr] << std::endl;
      // }
    }

    //
    // 2)  Calculate directed 2core edges
    // typename graph_type::template
    // vertex_data<boost::container::flat_map<vertex_locator,uint32_t>,
    // std::allocator<boost::container::flat_map<vertex_locator,uint32_t>> >
    // core2_directed(g);
    // typename graph_type::template
    // vertex_data<std::unordered_map<vertex_locator,uint32_t,vertex_locator_hash>,
    // std::allocator<std::unordered_map<vertex_locator,uint32_t,vertex_locator_hash>
    // > >   core2_directed(g);
    double start_time = MPI_Wtime();
    {
      auto alg_data = std::forward_as_tuple(core2_degree, core2_directed, g);
      auto vq = create_visitor_queue<directed_core2<graph_type>, lifo_queue>(
          &g, alg_data);
      vq.init_visitor_traversal();
    }
    double end_time = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Directed 2Core time = " << end_time - start_time
                << std::endl;
    }

    uint64_t local_core2_directed_edge_count(0);
    for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
      local_core2_directed_edge_count += core2_directed[*vitr].size();
    }

    uint64_t global_core2_directed_edge_count =
        comm_world().all_reduce(local_core2_directed_edge_count, MPI_SUM);
    if (comm_world().rank() == 0) {
      std::cout << "global_core2_directed_edge_count = "
                << global_core2_directed_edge_count << std::endl;
    }

    // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    //    ofs_sv_partition << g.locator_to_label(*vitr);
    //    for(auto pair : core2_directed[*vitr]) {
    //      ofs_sv_partition << " " << pair.first;
    //    }
    //    ofs_sv_partition << std::endl;
    //}

    // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    //   std::cout << "Vertex " << g.locator_to_label(*vitr) << " :: ";
    //     for(auto p : core2_directed[*vitr]) {
    //       std::cout << p.first << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Sort Directed 2 core -- lazy
    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t local_max_dod(0), local_max_deg(0);
    start_time = MPI_Wtime();
    {
      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        // core2_directed[*vitr].do_sort();
        local_max_dod =
            std::max(local_max_dod, uint64_t(core2_directed[*vitr].size()));
        local_max_deg = std::max(local_max_deg, uint64_t(g.degree(*vitr)));
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
        // core2_directed[*citr].do_sort();
        local_max_dod =
            std::max(local_max_dod, uint64_t(core2_directed[*citr].size()));
        local_max_deg = std::max(local_max_deg, uint64_t(g.degree(*citr)));
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    uint64_t global_max_dod =
        mpi_all_reduce(local_max_dod, std::greater<uint64_t>(), MPI_COMM_WORLD);
    uint64_t global_max_deg =
        mpi_all_reduce(local_max_deg, std::greater<uint64_t>(), MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Largest DOD out degree = " << global_max_dod << std::endl;
      std::cout << "Largest orig degree = " << global_max_deg << std::endl;
    }
    {
      //
      // 4)  Compute distributions
      std::map<uint64_t, uint64_t> local_orig_degree;
      std::map<uint64_t, uint64_t> local_2core_degree;
      std::map<uint64_t, uint64_t> local_2core_out_degree;

      uint64_t local_edge_count(0), local_dod_edge_count(0),
          local_in_zero_count(0), local_in_zero_edges_count(0);

      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        local_orig_degree[g.degree(*vitr)]++;
        local_2core_degree[core2_degree[*vitr]]++;
        local_2core_out_degree[core2_directed[*vitr].size()]++;
        local_edge_count += g.degree(*vitr);
        local_dod_edge_count += core2_directed[*vitr].size();
        if (core2_degree[*vitr] == core2_directed[*vitr].size()) {
          ++local_in_zero_count;
          local_in_zero_edges_count += core2_directed[*vitr].size();
        }
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
        local_orig_degree[g.degree(*citr)]++;
        local_2core_degree[core2_degree[*citr]]++;
        local_2core_out_degree[core2_directed[*citr].size()]++;
        local_edge_count += g.degree(*citr);
        local_dod_edge_count += core2_directed[*citr].size();
        if (core2_degree[*citr] == core2_directed[*citr].size()) {
          ++local_in_zero_count;
          local_in_zero_edges_count += core2_directed[*citr].size();
        }
      }

      uint64_t global_edge_count = mpi_all_reduce(
          local_edge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t global_dod_edge_count = mpi_all_reduce(
          local_dod_edge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t global_in_zero_count = mpi_all_reduce(
          local_in_zero_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t global_in_zero_edge_count = mpi_all_reduce(
          local_in_zero_edges_count, std::plus<uint64_t>(), MPI_COMM_WORLD);

      if (mpi_rank == 0) {
        std::cout << "global_edge_count = " << global_edge_count << std::endl;
        std::cout << "global_dod_edge_count = " << global_dod_edge_count
                  << std::endl;
        std::cout << "global_in_zero_count = " << global_in_zero_count
                  << std::endl;
        std::cout << "global_in_zero_edge_count = " << global_in_zero_edge_count
                  << std::endl;
      }

      // std::stringstream orig_degree, tcore_degree, tcore_out_degree;
      // orig_degree << deg_output_fname << "_orig_degree.txt";
      // tcore_degree << deg_output_fname << "_2core_degree.txt";
      // tcore_out_degree << deg_output_fname << "_2core_out_degree.txt";

      // output_degree_distribution(local_orig_degree,
      // orig_degree.str().c_str());
      // output_degree_distribution(local_2core_degree,
      // tcore_degree.str().c_str());
      // output_degree_distribution(local_2core_out_degree,
      // tcore_out_degree.str().c_str());
    }
  }

  //
  // 3)  Build wedges & count
  uint64_t local_triangle_count(0), local_wedge_count(0),
      local_master_wedge_count(0), local_delegate_skip(0);
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();
  {
    auto alg_data = std::forward_as_tuple(
        core2_directed, local_triangle_count, local_wedge_count,
        local_master_wedge_count,
        local_delegate_skip /*, count_wedges_created, count_wedges_checked, count_triangles_found*/);
    auto vq = create_visitor_queue<core2_wedges<graph_type>, lifo_queue>(
        &g, alg_data);
    vq.init_visitor_traversal();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double   end_time = MPI_Wtime();
  uint64_t global_wedge_count =
      mpi_all_reduce(local_wedge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
  uint64_t global_master_wedge_count =
      comm_world().all_reduce(local_master_wedge_count, MPI_SUM);
  uint64_t global_delegate_skip =
      comm_world().all_reduce(local_delegate_skip, MPI_SUM);
  if (mpi_rank == 0) {
    std::cout << "TC on directed 2core time = " << end_time - start_time
              << std::endl;
    std::cout << "Total wedges checked = " << global_wedge_count << std::endl;
    std::cout << "global_master_wedge_count = " << global_master_wedge_count
              << std::endl;
    std::cout << "global_delegate_skip = " << global_delegate_skip << std::endl;
  }

  uint64_t local_edge_triangle_count(0);
  uint64_t local_j_count(0);
  uint64_t local_edge_count(0);
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    for (const auto& edge_info : core2_directed[*vitr]) {
      local_edge_triangle_count += edge_info.second.edge_triangle_count;
      ++local_edge_count;
      if (edge_info.second.is_j) {
        ++local_j_count;
      }
    }
  }
  uint64_t local_controller_dogr_degree(0);
  for (auto citr = g.controller_begin(); citr != g.controller_end(); ++citr) {
    // std::cout << "Controller  Info: " << whoami()
    //           << " orig deg = " << g.degree(*citr)
    //           << " dogr degree = " << core2_directed[*citr].size() <<
    //           std::endl;
    local_controller_dogr_degree += core2_directed[*citr].size();
    for (const auto& edge_info : core2_directed[*citr]) {
      local_edge_triangle_count += edge_info.second.edge_triangle_count;
      ++local_edge_count;
      if (edge_info.second.is_j) {
        ++local_j_count;
      }
    }
  }
  uint64_t global_edge_triangle_count =
      comm_world().all_reduce(local_edge_triangle_count, MPI_SUM);

  uint64_t global_j_count = comm_world().all_reduce(local_j_count, MPI_SUM);

  uint64_t global_edge_count =
      comm_world().all_reduce(local_edge_count, MPI_SUM);
  uint64_t global_controller_dogr_degree =
      comm_world().all_reduce(local_controller_dogr_degree, MPI_SUM);

  if (comm_world().rank() == 0) {
    std::cout << "global_controller_dogr_degree = "
              << global_controller_dogr_degree << std::endl;
    std::cout << "global_edge_triangle_count / 3 = "
              << global_edge_triangle_count / 3 << std::endl;
    std::cout << "global_j_count = " << global_j_count << " ratio = "
              << double(global_j_count) / double(global_edge_count)
              << std::endl;
  }

  uint64_t global_triangle_count = mpi_all_reduce(
      local_triangle_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
  return global_triangle_count;
}

}  // end namespace havoqgt

#endif  // HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
