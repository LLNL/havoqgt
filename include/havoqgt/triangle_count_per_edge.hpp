// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma include once

#include <boost/container/deque.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <deque>
#include <fstream>
#include <havoqgt/visitor_queue.hpp>
#include <map>

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
  std::deque<Visitor> m_data;

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
  } else if (deg_a < deg_b) {
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
    if (from_degree < std::get<1>(alg_data).degree(vertex)) return false;

    if (vertex.is_delegate()) {
      if (!vertex.is_delegate_master()) {
        return true;
      }
    }

    // only here should be low-degree & masters
    if (edge_order_gt(from_degree, std::get<1>(alg_data).degree(vertex),
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

      std::get<0>(alg_data)[vv][fl].target_degree = from_degree;
    }
    //}
    //}
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (g.degree(vertex) >= 2) {
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
      if (g.degree(vertex) >= 2) {
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
      } else {
        std::cout << "My degreess < 2??" << std::endl;
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
        std::get<0>(alg_data)[vertex][check_close].increment_edge();
        // std::get<5>(alg_data)[vertex].insert(from_vertex);  // i
        return true;
      }
    } else {  // already found, not counting
      if (std::get<0>(alg_data)[vertex].count(check_close) == 0 ||
          std::get<0>(alg_data)[vertex].count(from_vertex) == 0) {
        std::cerr << "Error in edge counting" << std::endl;
      } else {
        std::get<0>(alg_data)[vertex][check_close].increment_edge();
        std::get<0>(alg_data)[vertex][from_vertex].increment_edge();
      }
    }
    return false;
  }

  template <typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
                  AlgData& alg_data) const {
    if (std::get<0>(alg_data)[vertex].size() > 1) {
      for (const auto& pair_a : std::get<0>(alg_data)[vertex]) {
        if (pair_a.second.fully_deleted) continue;
        for (const auto& pair_b : std::get<0>(alg_data)[vertex]) {
          if (pair_a.first == pair_b.first) continue;
          if (pair_b.second.fully_deleted) continue;
          if (edge_order_gt(pair_b.second.target_degree,
                            pair_a.second.target_degree, pair_b.first,
                            pair_a.first)) {
            my_type new_visitor(pair_a.first, pair_b.first, vertex, true);
            vis_queue->queue_visitor(new_visitor);
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

template <typename TGraph, typename DOGR>
void count_all_triangles_from_scratch(TGraph& g, DOGR& dogr) {
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
    auto alg_data  = std::forward_as_tuple(dogr, local_triangle_count,
                                          local_wedge_count, k, cut_count);

    auto vq =
        create_visitor_queue<core2_wedges<TGraph>, lifo_queue>(&g, alg_data);
    vq.init_visitor_traversal();
  }
  uint64_t global_triangle_count =
      comm_world().all_reduce(local_triangle_count, MPI_SUM);
  uint64_t global_wedge_count =
      comm_world().all_reduce(local_wedge_count, MPI_SUM);
  if (comm_world().rank() == 0) {
    std::cout << "Triangle Count = " << global_triangle_count << std::endl
              << "Number of Wedge Checks = " << global_wedge_count << std::endl;
  }
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
    //
    // 2)  Calculate directed 2core edges
    double start_time = MPI_Wtime();
    {
      auto alg_data = std::forward_as_tuple(core2_directed, g);
      auto vq = create_visitor_queue<directed_core2<graph_type>, lifo_queue>(
          &g, alg_data);
      vq.init_visitor_traversal();
    }
    double end_time = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Time to construct DODGraph (seconds) = "
                << end_time - start_time << std::endl;
    }

    uint64_t local_core2_directed_edge_count(0);
    for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
      local_core2_directed_edge_count += core2_directed[*vitr].size();
    }

    uint64_t global_core2_directed_edge_count =
        comm_world().all_reduce(local_core2_directed_edge_count, MPI_SUM);

    // (was) Sort Directed 2 core
    MPI_Barrier(MPI_COMM_WORLD);
    uint64_t local_max_dod(0), local_max_deg(0), local_max_dod_orig_deg(0);
    start_time = MPI_Wtime();
    {
      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        if (core2_directed[*vitr].size() > local_max_dod) {
          local_max_dod_orig_deg = uint64_t(g.degree(*vitr));
        }
        local_max_dod =
            std::max(local_max_dod, uint64_t(core2_directed[*vitr].size()));
        local_max_deg = std::max(local_max_deg, uint64_t(g.degree(*vitr)));
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
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
      std::cout << "Largest original degree = " << global_max_deg << std::endl;
      std::cout << "Largest DODGraph out degree = " << global_max_dod
                << std::endl;
    }

    {  // 4)  Compute distributions
      //
      uint64_t local_edge_count(0), local_dod_edge_count(0);

      for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
        local_edge_count += g.degree(*vitr);
        local_dod_edge_count += core2_directed[*vitr].size();
      }
      for (auto citr = g.controller_begin(); citr != g.controller_end();
           ++citr) {
        local_edge_count += g.degree(*citr);
        local_dod_edge_count += core2_directed[*citr].size();
      }

      uint64_t global_edge_count = mpi_all_reduce(
          local_edge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);
      uint64_t global_dod_edge_count = mpi_all_reduce(
          local_dod_edge_count, std::plus<uint64_t>(), MPI_COMM_WORLD);

      if (mpi_rank == 0) {
        std::cout << "Count of nonzeros in original graph = "
                  << global_edge_count << std::endl;
        std::cout << "Count of directed edges in DODGraph = "
                  << global_dod_edge_count << std::endl;
      }
    }
  }

  // Copy edges into dod_graph_truss
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    dod_graph_truss[*vitr].insert(core2_directed[*vitr].begin(),
                                  core2_directed[*vitr].end());
    core2_directed[*vitr].clear();
  }
  for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
    dod_graph_truss[*vitr].insert(core2_directed[*vitr].begin(),
                                  core2_directed[*vitr].end());
    core2_directed[*vitr].clear();
  }
}

template <typename TGraph>
uint64_t triangle_count_per_edge(TGraph& g, const std::string& output_tc) {
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

  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();
  construct_dod_graph(g, dod_graph_truss);

  count_all_triangles_from_scratch(g, dod_graph_truss);
  double end_time = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "Total Triangle Count Time (seconds) = "
              << end_time - start_time << std::endl;
  }

  //
  // Print output if file path is provided
  if (output_tc.length() > 0) {
    std::stringstream ssoutfname;
    ssoutfname << output_tc << "_" << mpi_rank << "_of_" << mpi_size;
    std::ofstream ofs(ssoutfname.str().c_str());
    for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
      for (auto& kvp : dod_graph_truss[*vitr]) {
        ofs << g.locator_to_label(*vitr) << " " << g.locator_to_label(kvp.first)
            << " " << kvp.second.edge_triangle_count << std::endl;
      }
    }
    for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
      for (auto& kvp : dod_graph_truss[*vitr]) {
        ofs << g.locator_to_label(*vitr) << " " << g.locator_to_label(kvp.first)
            << " " << kvp.second.edge_triangle_count << std::endl;
      }
    }
  }

  //
  // Count total triangle edge count to return
  uint64_t count_to_return(0);
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    for (auto& kvp : dod_graph_truss[*vitr]) {
      count_to_return += kvp.second.edge_triangle_count;
    }
  }
  for (auto vitr = g.controller_begin(); vitr != g.controller_end(); ++vitr) {
    for (auto& kvp : dod_graph_truss[*vitr]) {
      count_to_return += kvp.second.edge_triangle_count;
    }
  }

  return comm_world().all_reduce(count_to_return, MPI_SUM);
}

}  // end namespace havoqgt
