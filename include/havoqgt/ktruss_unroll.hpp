// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma include once

#include <boost/container/deque.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <fstream>
#include <havoqgt/visitor_queue.hpp>
#include <map>
#include <set>

namespace havoqgt {

struct dod_graph_edge {
  uint32_t target_degree = 0;
};

struct dod_graph_truss_edge {
  dod_graph_truss_edge()
      : target_degree(0),
        edge_triangle_count(0),
        mark_for_deletion(0),
        fully_deleted(0) {}
  dod_graph_truss_edge(dod_graph_edge dge)
      : target_degree(dge.target_degree),
        edge_triangle_count(0),
        mark_for_deletion(0),
        fully_deleted(0) {}
  void decrement_edge() {
    if (edge_triangle_count > 0) --edge_triangle_count;
  }
  void increment_edge() {
    if (edge_triangle_count ==
        std::numeric_limits<decltype(edge_triangle_count)>::max()) {
      std::cout << "WARNING: edge_triangle_count overflow!" << std::endl;
    }
    ++edge_triangle_count;
  }

  uint32_t target_degree : 30;
  uint32_t mark_for_deletion : 1;
  uint32_t fully_deleted : 1;
  uint32_t edge_triangle_count;
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

template <typename VLOC>
inline bool edge_order_gt(uint32_t deg_a, uint32_t deg_b, VLOC v_a, VLOC v_b) {
  // auto max_deg = std::max(deg_a, deg_b);
  // auto min_deg = std::min(deg_a, deg_b);
  // if (max_deg / (1.5f) > min_deg) {
  //   return v_a.hash() > v_b.hash();
  // }
  if (deg_a > deg_b) {
    return true;
  }
  if (deg_a < deg_b) {
    return false;
  } else {
    return v_b < v_a;
  }
  // if (from_degree > /*std::get<2>(alg_data).degree(
  //                       vertex)*/ std::get<0>(alg_data)[vertex] ||
  //     (from_degree ==
  //          /*std::get<2>(alg_data).degree(
  //                 vertex)*/ std::get<0>(alg_data)[vertex] &&
  //      vertex.hash() < from_vertex.hash())) {
}

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
      : vertex(v), from_vertex(_from), from_degree(_from_degree), init(false) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    // if(std::get<0>(alg_data)[vertex] >= 2) {
    if (from_degree < std::get<2>(alg_data).degree(vertex)) return false;
    // previously returned true, but changing here --- return true;
    if (vertex.is_delegate()) {
      if (!vertex.is_delegate_master()) {
        return true;
      }
    }
    if (std::get<0>(alg_data)[vertex] < 2) return false;
    // only here should be low-degree & masters
    if (edge_order_gt(from_degree, std::get<2>(alg_data).degree(vertex),
                      from_vertex, vertex)) {
      // if (from_degree > /*std::get<2>(alg_data).degree(
      //                       vertex)*/ std::get<0>(alg_data)[vertex] ||
      //     (from_degree ==
      //          /*std::get<2>(alg_data).degree(
      //                 vertex)*/ std::get<0>(alg_data)[vertex] &&
      //      vertex.hash() < from_vertex.hash())) {
      auto vv = vertex;
      vv.set_bcast(0);
      vv.set_intercept(0);
      auto fl = from_vertex;
      fl.set_bcast(0);
      fl.set_intercept(0);

      std::get<1>(alg_data)[vv][fl].target_degree = from_degree;
      // std::get<1>(alg_data)[vv].add(fl, from_degree);
    }
    //}
    //}
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (std::get<0>(alg_data)[vertex] >= 2) {
      // if in 2core, send degree to neighbors
      uint32_t my_degree = g.degree(vertex);
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
      if (std::get<0>(alg_data)[vertex] >= 2) {
        // if in 2core, send degree to neighbors
        uint32_t my_degree = g.degree(vertex);
        for (auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex);
             ++eitr) {
          vertex_locator neighbor = eitr.target();
          if (neighbor == vertex) continue;
          my_type new_visitor(neighbor, vertex, my_degree);
          vis_queue->queue_visitor(new_visitor);
        }
        return true;
      }
    }
    return false;
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
  vertex_locator from_vertex;
  uint32_t       from_degree;
  bool           init;
};

template <typename Graph, bool Decompose>
class core2_wedges {
 public:
  typedef core2_wedges<Graph, Decompose> my_type;
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
      if (!vertex.is_delegate_master()) {
        return true;
      }
    }
    if (do_check_close) {       // Checking for closing edge (j,k)
      ++std::get<2>(alg_data);  // counts wedge checks
      if (std::get<0>(alg_data)[vertex].count(check_close) > 0) {
        if (std::get<0>(alg_data)[vertex][check_close].fully_deleted) {
          return false;
        }
        ++std::get<1>(alg_data);  // counts total triangles
        // counts edge partipation in triangles
        if (!Decompose) {
          std::get<0>(alg_data)[vertex][check_close].increment_edge();
          // std::get<5>(alg_data)[vertex].insert(from_vertex);  // i
          std::get<5>(alg_data)[vertex][from_vertex]++;
        } else {
          std::get<0>(alg_data)[vertex][check_close].decrement_edge();
          if (std::get<5>(alg_data)[vertex].count(from_vertex) == 0) {
            std::cout << "NOT EXPECTED" << std::endl;
          }
          if (--(std::get<5>(alg_data)[vertex][from_vertex]) == 0) {
            std::get<5>(alg_data)[vertex].erase(from_vertex);
          }
          // if iedges was map, could decrement here
        }
        return true;
      }
    } else {  // already found, not counting
      if (std::get<0>(alg_data)[vertex].count(check_close) == 0 ||
          std::get<0>(alg_data)[vertex].count(from_vertex) == 0) {
        std::cerr << "Error in edge counting" << std::endl;
      } else {
        if (!Decompose) {
          std::get<0>(alg_data)[vertex][check_close].increment_edge();
          std::get<0>(alg_data)[vertex][from_vertex].increment_edge();
        } else {
          std::get<0>(alg_data)[vertex][check_close].decrement_edge();
          std::get<0>(alg_data)[vertex][from_vertex].decrement_edge();
        }
      }
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    std::set<vertex_locator> to_delete;
    if (Decompose) {
      int k = std::get<3>(alg_data);
      for (auto& pair_a : std::get<0>(alg_data)[vertex]) {
        if (pair_a.second.fully_deleted) continue;
        if (!pair_a.second.mark_for_deletion &&
            pair_a.second.edge_triangle_count < k - 2) {
          pair_a.second.mark_for_deletion = true;
          to_delete.insert(pair_a.first);
          std::get<4>(alg_data)++;
        }
      }
      if (to_delete.size() == 0) return false;  // nothing to decompose
    }
    if (std::get<0>(alg_data)[vertex].size() > 1) {
      for (const auto& pair_a : std::get<0>(alg_data)[vertex]) {
        if (pair_a.second.fully_deleted) continue;
        for (const auto& pair_b : std::get<0>(alg_data)[vertex]) {
          if (pair_a.first == pair_b.first) continue;
          if (pair_b.second.fully_deleted) continue;
          if (edge_order_gt(pair_b.second.target_degree,
                            pair_a.second.target_degree, pair_b.first,
                            pair_a.first)) {
            if (!Decompose) {
              my_type new_visitor(pair_a.first, pair_b.first, vertex, true);
              vis_queue->queue_visitor(new_visitor);
            } else {
              if (to_delete.count(pair_a.first) > 0 ||
                  to_delete.count(pair_b.first) > 0) {
                // if eather edge is ready to delete
                // send delete request!
                my_type new_visitor(pair_a.first, pair_b.first, vertex, true);
                vis_queue->queue_visitor(new_visitor);
              }
            }
          }
        }
      }
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

template <typename Graph>
class unroll_iedges_visitor {
 public:
  typedef unroll_iedges_visitor<Graph>   my_type;
  typedef typename Graph::vertex_locator vertex_locator;

  unroll_iedges_visitor() : vertex() {}

  unroll_iedges_visitor(vertex_locator v) : vertex(v) {}

  unroll_iedges_visitor(vertex_locator v, vertex_locator j, vertex_locator k)
      : vertex(v), jvertex(j), kvertex(k) {}

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if (vertex.is_delegate()) {
      if (!vertex.is_delegate_master()) {
        return true;
      }
    }
    // if (std::get<0>(alg_data)[vertex].count(kvertex) > 0) {
    //   if(!std::get<0>(alg_data)[vertex][jvertex].fully_deleted || )
    // if (!std::get<0>(alg_data)[vertex][jvertex].fully_deleted &&
    //     !std::get<0>(alg_data)[vertex][kvertex].fully_deleted) {
    // kvertex exists, decrement j & k edges.
    // if (!std::get<0>(alg_data)[vertex][jvertex].mark_for_deletion) {
    if (std::get<0>(alg_data)[vertex].count(jvertex) > 0 &&
        std::get<0>(alg_data)[vertex].count(kvertex) > 0) {
      if (std::get<0>(alg_data)[vertex][jvertex].mark_for_deletion ||
          std::get<0>(alg_data)[vertex][kvertex].mark_for_deletion ||
          std::get<0>(alg_data)[vertex][jvertex].fully_deleted ||
          std::get<0>(alg_data)[vertex][kvertex].fully_deleted) {
        return false;
      }
      std::get<0>(alg_data)[vertex][jvertex].decrement_edge();
      std::get<0>(alg_data)[vertex][kvertex].decrement_edge();
    }
    //}
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    for (auto& adje : std::get<0>(alg_data)[vertex]) {
      if (adje.second.fully_deleted) continue;
      if (adje.second.mark_for_deletion &&
          !std::get<1>(alg_data)[vertex].empty()) {
        // send notice to all i's
        vertex_locator j = vertex;
        vertex_locator k = adje.first;
        std::get<2>(alg_data)++;
        for (auto& ic : std::get<1>(alg_data)[vertex]) {
          my_type new_visitor(ic.first, j, k);
          vis_queue->queue_visitor(new_visitor);
        }
        // std::get<1>(alg_data)[vertex].clear();  // only do this once
      }
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    std::cout << "Shoudn't be here" << std::endl;
    exit(-1);
  }

  friend inline bool operator>(const unroll_iedges_visitor& v1,
                               const unroll_iedges_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const unroll_iedges_visitor& v1,
                               const unroll_iedges_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator jvertex;
  vertex_locator kvertex;
};

template <typename TGraph, typename DOGR, typename IEDGES>
void count_all_triangles_from_scratch(TGraph& g, DOGR& dogr, IEDGES& iedges) {
  // typedef TGraph                          graph_type;
  // typedef typename TGraph::vertex_locator vertex_locator;
  // typename graph_type::template vertex_data<
  //     std::set<vertex_locator>, std::allocator<std::set<vertex_locator>>>
  //     iedges_tmp(g);
  //
  // reset all edges
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    for (auto& kvp : dogr[*vitr]) {
      kvp.second.edge_triangle_count = 0;
      kvp.second.mark_for_deletion   = 0;
    }
  }
  for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
    for (auto& kvp : dogr[*vitr]) {
      kvp.second.edge_triangle_count = 0;
      kvp.second.mark_for_deletion   = 0;
    }
  }

  uint64_t local_triangle_count(0), local_wedge_count(0);
  // MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();
  {
    int  k         = -1;  // dummy val
    int  cut_count = 0;   // dummy
    auto alg_data  = std::forward_as_tuple(
        dogr, local_triangle_count, local_wedge_count, k, cut_count, iedges);
    auto vq = create_visitor_queue<core2_wedges<TGraph, false>, lifo_queue>(
        &g, alg_data);
    vq.init_visitor_traversal();
  }

  // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
  //   iedges[*vitr].insert(iedges[*vitr].end(), iedges_tmp[*vitr].begin(),
  //                        iedges_tmp[*vitr].end());
  // }
  // for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr)
  // {
  //   iedges[*vitr].insert(iedges[*vitr].end(), iedges_tmp[*vitr].begin(),
  //                        iedges_tmp[*vitr].end());
  // }
  // MPI_Barrier(MPI_COMM_WORLD);
  // double   end_time           = MPI_Wtime();
  // uint64_t global_wedge_count = mpi_all_reduce(
  //     local_wedge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
  // if (mpi_rank == 0) {
  //   std::cout << "TC on directed 2core time = " << end_time -
  //   start_time
  //             << std::endl;
  //   std::cout << "Total wedges checked = " << global_wedge_count
  //             << std::endl;
  // }
}

template <typename TGraph, typename DOGR, typename IEDGES>
uint64_t unroll_jk(TGraph& g, DOGR& dogr, IEDGES& iedges) {
  uint64_t local_jk_cut_count(0);
  {
    auto alg_data = std::forward_as_tuple(dogr, iedges, local_jk_cut_count);
    auto vq = create_visitor_queue<unroll_iedges_visitor<TGraph>, lifo_queue>(
        &g, alg_data);
    vq.init_visitor_traversal();
  }
  return comm_world().all_reduce(local_jk_cut_count, MPI_SUM);
}

template <typename TGraph, typename DOGR, typename IEDGES>
void decompose_truss(TGraph& g, DOGR& dogr, IEDGES& iedges, int k,
                     uint64_t& global_edges_remaining,
                     uint64_t& global_vertices_remaining,
                     uint64_t& global_triangles_remaining) {
  uint64_t global_cut_count(0);
  uint64_t local_total_removed_count(0), local_edges_remaining(0),
      local_vertices_remaining(0), local_triangles_remaining(0);
  do {
    local_total_removed_count = 0;
    local_edges_remaining     = 0;
    local_vertices_remaining  = 0;
    local_triangles_remaining = 0;
    uint64_t local_cut_count  = 0;
    uint64_t local_triangle_count(0), local_wedge_count(0);

    double start_time = MPI_Wtime();
    {
      auto alg_data =
          std::forward_as_tuple(dogr, local_triangle_count, local_wedge_count,
                                k, local_cut_count, iedges);
      auto vq = create_visitor_queue<core2_wedges<TGraph, true>, lifo_queue>(
          &g, alg_data);
      vq.init_visitor_traversal();
    }
    global_cut_count = comm_world().all_reduce(local_cut_count, MPI_SUM);
    // if (global_cut_count > 0) {
    global_cut_count += unroll_jk(g, dogr, iedges);

    //
    // Remove marked edges.
    for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
      size_t vert_alive_degree(0);
      for (auto itr = dogr[*vitr].begin(); itr != dogr[*vitr].end(); ++itr) {
        if (itr->second.fully_deleted) {
          continue;
        }
        if (itr->second.mark_for_deletion) {
          ++local_total_removed_count;
          itr->second.fully_deleted = true;
        } else {
          ++local_edges_remaining;
          ++vert_alive_degree;
          local_triangles_remaining += itr->second.edge_triangle_count;
        }
      }
      if (vert_alive_degree == 0) {
        dogr[*vitr].clear();
      } else {
        ++local_vertices_remaining;
      }
    }
    for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
      size_t vert_alive_degree(0);
      for (auto itr = dogr[*vitr].begin(); itr != dogr[*vitr].end(); ++itr) {
        if (itr->second.fully_deleted) {
          continue;
        }
        if (itr->second.mark_for_deletion) {
          ++local_total_removed_count;
          itr->second.fully_deleted = true;
        } else {
          ++local_edges_remaining;
          ++vert_alive_degree;
          local_triangles_remaining += itr->second.edge_triangle_count;
        }
      }
      if (vert_alive_degree == 0) {
        dogr[*vitr].clear();
      } else {
        ++local_vertices_remaining;
      }
    }
  } while (comm_world().all_reduce(local_total_removed_count, MPI_SUM) > 0);

  global_edges_remaining =
      comm_world().all_reduce(local_edges_remaining, MPI_SUM);
  global_vertices_remaining =
      comm_world().all_reduce(local_vertices_remaining, MPI_SUM);
  global_triangles_remaining =
      comm_world().all_reduce(local_triangles_remaining, MPI_SUM) / 3;
}

template <typename TGraph, typename DODgraph>
void construct_dod_graph(TGraph& g, DODgraph& dod_graph_truss) {
  typedef TGraph                          graph_type;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  typename graph_type::template vertex_data<
      std::map<vertex_locator, dod_graph_edge>,
      std::allocator<std::map<vertex_locator, dod_graph_edge>>>
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
    }

    //
    // 2)  Calculate directed 2core edges
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

    // (was) Sort Directed 2 core
    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t local_max_dod(0), local_max_deg(0), local_max_dod_orig_deg(0);
    start_time = MPI_Wtime();
    {
      std::stringstream ss;
      // ss << "degrees_" << comm_world().rank();
      // std::ofstream ofs_degrees(ss.str().c_str());
      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        // ofs_degrees << core2_directed[*vitr].size() << " " << g.degree(*vitr)
        //            << std::endl;
        if (core2_directed[*vitr].size() > local_max_dod) {
          local_max_dod_orig_deg = uint64_t(g.degree(*vitr));
        }
        local_max_dod =
            std::max(local_max_dod, uint64_t(core2_directed[*vitr].size()));
        local_max_deg = std::max(local_max_deg, uint64_t(g.degree(*vitr)));
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
        // ofs_degrees << core2_directed[*citr].size() << " " << g.degree(*citr)
        //            << std::endl;
        if (core2_directed[*citr].size() > local_max_dod) {
          local_max_dod_orig_deg = uint64_t(g.degree(*citr));
        }
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
    // std::cout << whoami() << " max_local_dod_deg = " << local_max_dod
    //           << ", orig deg was " << local_max_dod_orig_deg << std::endl;
    {  // 4)  Compute distributions
      //

      uint64_t local_edge_count(0), local_dod_edge_count(0),
          local_in_zero_count(0), local_in_zero_edges_count(0);

      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        local_edge_count += g.degree(*vitr);
        local_dod_edge_count += core2_directed[*vitr].size();
        if (core2_degree[*vitr] == core2_directed[*vitr].size()) {
          ++local_in_zero_count;
          local_in_zero_edges_count += core2_directed[*vitr].size();
        }
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
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
    }
  }

  // Copy edges into dod_graph_truss
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    dod_graph_truss[*vitr].insert(core2_directed[*vitr].begin(),
                                  core2_directed[*vitr].end());
  }
  for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
  }
}

template <typename TGraph>
uint64_t ktruss_unroll(TGraph& g) {
  typedef TGraph                          graph_type;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  typename graph_type::template vertex_data<
      boost::container::flat_map<vertex_locator, dod_graph_truss_edge>,
      std::allocator<
          boost::container::flat_map<vertex_locator, dod_graph_truss_edge>>>
      dod_graph_truss(g);

  construct_dod_graph(g, dod_graph_truss);

  MPI_Barrier(MPI_COMM_WORLD);
  double watchj_start_time = MPI_Wtime();
  typename graph_type::template vertex_data<
      std::map<vertex_locator, uint32_t>,
      std::allocator<std::map<vertex_locator, uint32_t>>>
           iedges(g);
  int      k = 3;
  uint64_t global_edges_remaining(0);
  do {
    uint64_t global_vertices_remaining(0), global_triangles_remaining(0);
    double   this_truss_start = MPI_Wtime();
    if (k == 3) {
      count_all_triangles_from_scratch(g, dod_graph_truss, iedges);
    }
    decompose_truss(g, dod_graph_truss, iedges, k, global_edges_remaining,
                    global_vertices_remaining, global_triangles_remaining);
    double this_truss_end = MPI_Wtime();
    if (comm_world().rank() == 0) {
      std::cout << "K= " << k << " edges_remaining= " << global_edges_remaining
                << " vertices_remaining= " << global_vertices_remaining
                << " triangles_remaining= " << global_triangles_remaining
                << " time= " << this_truss_end - this_truss_start << std::endl;
    }
    ++k;
  } while (global_edges_remaining > 0);
  MPI_Barrier(MPI_COMM_WORLD);
  double watchj_end_time = MPI_Wtime();
  if (comm_world().rank() == 0) {
    std::cout << "TOTAL KTRUSS TIME = " << watchj_end_time - watchj_start_time
              << std::endl;
  }

  return 0;
}

}  // end namespace havoqgt
