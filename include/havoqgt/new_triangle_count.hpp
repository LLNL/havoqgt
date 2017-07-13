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

#ifndef HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
#define HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>
#include <fstream>
#include <boost/container/flat_map.hpp>


namespace havoqgt {


template <typename Visitor>
class triangle_priority_queue
{

protected:
  std::priority_queue< Visitor, std::deque<Visitor>,
                               std::greater<Visitor> > m_data;
public:
  triangle_priority_queue() { }

  bool push(Visitor const & task)
  {
    m_data.push(task);
    return true;
  }

  void pop()
  {
    m_data.pop();
  }

  Visitor const & top() //const
  {
    return m_data.top();
  }

  size_t size() const
  {
    return m_data.size();;
  }

  bool empty() const
  {
    return m_data.empty();
  }

  void clear()
  {
    m_data.clear();
  }
};


template<typename Graph>
class triangle_count_visitor {
public:
  typedef triangle_count_visitor<Graph>    my_type;
  typedef typename Graph::vertex_locator                 vertex_locator;

  triangle_count_visitor(): vertex(), first(), second() { }

  triangle_count_visitor(vertex_locator v): vertex(v), first(), second() { }

  triangle_count_visitor(vertex_locator v, vertex_locator f): vertex(v), first(f), second() { }

  triangle_count_visitor(vertex_locator v, vertex_locator f, vertex_locator s) : vertex(v), first(f), second(s) { }

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    assert(!first.is_valid() && !second.is_valid());
    //std::cout << "Init vertex_count_lower_degree = " << vertex_count_lower_degree(alg_data, vertex) << std::endl;
    if(vertex_count_lower_degree(alg_data, vertex) > 0) {
      for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        if(do_forward(alg_data, eitr, g)) {
          auto neighbor = eitr.target();
          my_type new_visitor( neighbor, vertex );
          //std::cout << "First clone" << std::endl;
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return true;
    }
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    assert(first.is_valid());

    if(!second.is_valid()) {
      //std::cout << "Second" << std::endl;
      if(vertex_count_lower_degree(alg_data, vertex) > 0) {
        for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
          if(do_forward(alg_data, eitr, g)) {
            auto neighbor = eitr.target();
            my_type new_visitor( neighbor, vertex, first );
            vis_queue->queue_visitor(new_visitor);
          }
        }
        return true;
      }
      return false;
    } else {
      //std::cout << "Looking for closing edge!" << std::endl;
      for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        if(g.locator_to_label(neighbor) == g.locator_to_label(second)) {
          std::get<4>(alg_data)[vertex] = std::get<4>(alg_data)[vertex] + 1;
        }
      }
      return true;
    }
    assert(false);
  }

  template< typename AlgData>
  bool do_forward(AlgData& ad, typename Graph::edge_iterator e, const Graph& g) const {
    if(vertex_degree(ad,vertex) > edge_target_degree(ad,e)) {
      //std::cout << "vertex_degree(ad,vertex) > edge_target_degree(ad,e) : " << edge_target_count_lower_degree(ad,e) << std::endl;
      return edge_target_count_lower_degree(ad,e) > 0;
    } else if(vertex_degree(ad,vertex) == edge_target_degree(ad,e) && g.locator_to_label(vertex) > g.locator_to_label(e.target())) {
      return edge_target_count_lower_degree(ad,e) > 0;
    }
    return false;
  }

  //std::forward_as_tuple(vertex_count_lower_degree, vertex_degree, edge_target_degree, edge_target_count_lower_degree);
  template<typename AlgData>
  uint32_t vertex_count_lower_degree(const AlgData& ad, vertex_locator v) const { return std::get<0>(ad)[v]; }

  template<typename AlgData>
  uint32_t vertex_degree(const AlgData& ad, vertex_locator v) const { return std::get<1>(ad)[v]; }

  template<typename AlgData>
  uint32_t edge_target_degree(const AlgData& ad, typename Graph::edge_iterator e) const { return std::get<2>(ad)[e]; }

  template<typename AlgData>
  uint32_t edge_target_count_lower_degree(const AlgData& ad, typename Graph::edge_iterator e) const { return std::get<3>(ad)[e]; }

  int get_state() const {
    if(!first.is_valid()) {
      return 0;
    }
    if(!second.is_valid()) {
      return 1;
    }
    return 2;
  }

  friend inline bool operator>(const triangle_count_visitor& v1, const triangle_count_visitor& v2) {
    return v1.get_state() < v2.get_state();
  }

  friend inline bool operator<(const triangle_count_visitor& v1, const triangle_count_visitor& v2) {
    return v1.get_state() > v2.get_state();
  }

  vertex_locator vertex, first, second;
};


template<typename Graph>
class core2_visitor {
public:
  typedef core2_visitor<Graph>    my_type;
  typedef typename Graph::vertex_locator                 vertex_locator;

  core2_visitor() : vertex(), init(true) { }

  core2_visitor(vertex_locator v) : vertex(v), init(true) { }

    core2_visitor(vertex_locator v, bool _init) : vertex(v), init(_init) { }

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    //Retrun true if alive
    return std::get<1>(alg_data)[vertex];
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
      if(std::get<1>(alg_data)[vertex]) {
        //   --(std::get<0>(alg_data)[vertex]);
        if(std::get<0>(alg_data)[vertex] < 2) {
          //remove from 2 core
          std::get<1>(alg_data)[vertex] = false;
          std::get<0>(alg_data)[vertex] = 0;
          for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
            vertex_locator neighbor = eitr.target();
            my_type new_visitor( neighbor, false);
            vis_queue->queue_visitor(new_visitor);
          }
        }
        return true;
      }
      return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(init) {
      if(std::get<0>(alg_data)[vertex] < 2) {
        //remove from 2 core
        std::get<1>(alg_data)[vertex] = false;
        std::get<0>(alg_data)[vertex] = 0;
        for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type new_visitor( neighbor, false);
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return true;
    }


    if(std::get<1>(alg_data)[vertex]) {
      --(std::get<0>(alg_data)[vertex]);
      if(std::get<0>(alg_data)[vertex] < 2) {
        //remove from 2 core
        std::get<1>(alg_data)[vertex] = false;
        std::get<0>(alg_data)[vertex] = 0;
        for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type new_visitor( neighbor, false );
          vis_queue->queue_visitor(new_visitor);

        }
      }
      return true;
    }
    return false;
  }

  friend inline bool operator>(const core2_visitor& v1, const core2_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const core2_visitor& v1, const core2_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  bool init;
};





template<typename Graph>
class directed_core2 {
public:
  typedef directed_core2<Graph>    my_type;
  typedef typename Graph::vertex_locator                 vertex_locator;

  directed_core2() : vertex(), init(true) { }

  directed_core2(vertex_locator v) : vertex(v), init(true) { }

  directed_core2(vertex_locator v, uint64_t _from, uint32_t _from_degree) : vertex(v), from_label(_from), from_degree(_from_degree), init(false) { }

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(std::get<0>(alg_data)[vertex] >= 2) {
      if(from_degree >= std::get<2>(alg_data).degree(vertex)) {
        return true;
      }
    }
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
      if(std::get<0>(alg_data)[vertex] >= 2 ) {
        //if in 2core, send degree to neighbors
        uint64_t my_label  = g.locator_to_label(vertex);
        uint32_t my_degree = g.degree(vertex);
        for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type new_visitor( neighbor, my_label, my_degree);
          vis_queue->queue_visitor(new_visitor);
        }
        return true;
      }
      return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(init) {
      if(std::get<0>(alg_data)[vertex] >= 2 ) {
        //if in 2core, send degree to neighbors
        uint64_t my_label  = g.locator_to_label(vertex);
        uint32_t my_degree = g.degree(vertex);
        for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
          vertex_locator neighbor = eitr.target();
          my_type new_visitor( neighbor, my_label, my_degree);
          vis_queue->queue_visitor(new_visitor);
        }
        return true;
      }
    }

    if(from_degree > g.degree(vertex) ||
      (from_degree == g.degree(vertex) && from_label > g.locator_to_label(vertex))) {

      std::get<1>(alg_data)[vertex][from_label] = from_degree;
    }
    return false;
  }

  friend inline bool operator>(const directed_core2& v1, const directed_core2& v2) {
    return false;
  }

  friend inline bool operator<(const directed_core2& v1, const directed_core2& v2) {
    return false;
  }

  vertex_locator vertex;
  uint64_t from_label;
  uint32_t from_degree;
  bool init;
};



template<typename Graph>
class core2_wedges {
public:
  typedef core2_wedges<Graph>    my_type;
  typedef typename Graph::vertex_locator                 vertex_locator;

  core2_wedges() : vertex() { }

  core2_wedges(vertex_locator v) : vertex(v) { }

  core2_wedges(vertex_locator v, uint64_t _check_close) : vertex(v), check_close(_check_close) { }

  template <typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(vertex.is_delegate()) {
      if(vertex.is_delegate_master()) {
        ++std::get<2>(alg_data);
        if(std::get<0>(alg_data)[vertex].count(check_close) > 0) {
          ++std::get<1>(alg_data);
          //std::cout << "FOUND" << std::endl;
        }
        return false;
      }
      return true;
    }
    ++std::get<2>(alg_data);
    if(std::get<0>(alg_data)[vertex].count(check_close) > 0) {
      ++std::get<1>(alg_data);
      //std::cout << "FOUND" << std::endl;
    }
    return false;

  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(std::get<0>(alg_data)[vertex].size() > 1) {
      for(auto pair_a : std::get<0>(alg_data)[vertex]) {
        for(auto pair_b : std::get<0>(alg_data)[vertex]) {

          if(pair_a.second < pair_b.second ||
            (pair_a.second == pair_b.second && pair_a.first < pair_b.first)) {
              my_type new_visitor(g.label_to_locator(pair_a.first), pair_b.first);
              vis_queue->queue_visitor(new_visitor);
            }
        }
      }
    }
      return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    ++std::get<2>(alg_data);
    if(std::get<0>(alg_data)[vertex].count(check_close) > 0) {
      ++std::get<1>(alg_data);
      //std::cout << "FOUND" << std::endl;
    }
    return false;
  }

  friend inline bool operator>(const core2_wedges& v1, const core2_wedges& v2) {
    return false;
  }

  friend inline bool operator<(const core2_wedges& v1, const core2_wedges& v2) {
    return false;
  }

  vertex_locator vertex;
  uint64_t check_close;
};


void output_degree_distribution(const std::map<uint64_t, uint64_t>& local_degree_count, const char* fname) {
  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  std::vector< std::vector<std::pair<uint64_t, uint64_t>> > send_p_vec(mpi_size);
  std::vector< std::vector<std::pair<uint64_t, uint64_t>> > recv_p_vec(mpi_size);
  for(auto deg_cnt : local_degree_count) {
    send_p_vec[deg_cnt.first % uint64_t(mpi_size)].push_back(deg_cnt);
  }

  mpi_all_to_all(send_p_vec, recv_p_vec, MPI_COMM_WORLD);

  std::map<uint64_t, uint64_t> partitioned_deg_count;
  for(auto& vec_deg_cnt : recv_p_vec) {
    for(auto deg_cnt : vec_deg_cnt) {
      partitioned_deg_count[deg_cnt.first] += deg_cnt.second;
    }
  }


  //
  // Send all to Rank 0 -- this is not efficient way
  send_p_vec.clear();
  send_p_vec.resize(mpi_size);
  recv_p_vec.clear();
  recv_p_vec.resize(mpi_size);
  for(auto vec_deg_cnt : partitioned_deg_count) {
    send_p_vec[0].push_back(vec_deg_cnt);
  }
  mpi_all_to_all(send_p_vec, recv_p_vec, MPI_COMM_WORLD);

  if(mpi_rank == 0) {
    std::vector<std::pair<uint64_t, uint64_t>> all_sorted;
    for(auto & vec_deg_cnt : recv_p_vec) {
      for(auto deg_cnt : vec_deg_cnt) {
        all_sorted.push_back(deg_cnt);
      }
    }
    std::sort(all_sorted.begin(), all_sorted.end());
    std::ofstream outfile(fname);
    for(auto deg_cnt : all_sorted) {
      outfile << deg_cnt.first << "\t" << deg_cnt.second << std::endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

}


template <typename TGraph>
uint64_t new_triangle_count(TGraph& g, const char* deg_output_fname) {
  typedef TGraph                                             graph_type;
  typedef triangle_count_visitor<TGraph>             visitor_type;
  typedef typename TGraph::vertex_locator                 vertex_locator;
  typename graph_type::template vertex_data<uint64_t, std::allocator<uint64_t> >   tc_data(g);


  int mpi_size(0), mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

  //
  // 1)  Calculate 2core degree.   0 degree means that vertex is not in 2core
  typename graph_type::template vertex_data<uint32_t, std::allocator<uint32_t> >   core2_degree(g);
  {

    typename graph_type::template vertex_data<bool, std::allocator<bool> >           core2_alive(g);
    core2_alive.reset(true);
    for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
      core2_degree[*vitr] = g.degree(*vitr);
    }
    for (auto ditr = g.delegate_vertices_begin(); ditr != g.delegate_vertices_end(); ++ditr) {
      core2_degree[*ditr] = g.degree(*ditr);
    }


    double start_time = MPI_Wtime();
    {
      auto alg_data = std::forward_as_tuple(core2_degree, core2_alive);
      auto vq = create_visitor_queue<core2_visitor<graph_type>, detail::visitor_priority_queue>(&g, alg_data);
      vq.init_visitor_traversal_new();
    }
    double end_time = MPI_Wtime();
    if(mpi_rank == 0) {
      std::cout << "2Core time = " << end_time - start_time << std::endl;
    }



    // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    //   std::cout << "Vertex " << g.locator_to_label(*vitr) << " has core2_degree = " << core2_degree[*vitr] << std::endl;
    // }
    // for (auto ditr = g.delegate_vertices_begin(); ditr != g.delegate_vertices_end(); ++ditr) {
    //   std::cout << "Delegate " << g.locator_to_label(*ditr) << " has core2_degree = " <<  core2_degree[*ditr] << std::endl;
    // }
  }

  //
  // 2)  Calculate directed 2core edges
  typename graph_type::template vertex_data<boost::container::flat_map<uint64_t,uint32_t>, std::allocator<boost::container::flat_map<uint64_t,uint32_t>> >   core2_directed(g);
  double start_time = MPI_Wtime();
  {
    auto alg_data = std::forward_as_tuple(core2_degree, core2_directed, g);
    auto vq = create_visitor_queue<directed_core2<graph_type>, detail::visitor_priority_queue>(&g, alg_data);
    vq.init_visitor_traversal_new();
  }
  double end_time = MPI_Wtime();
  if(mpi_rank == 0) {
    std::cout << "Directed 2Core time = " << end_time - start_time << std::endl;
  }

  // for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
  //   std::cout << "Vertex " << g.locator_to_label(*vitr) << " :: ";
  //     for(auto p : core2_directed[*vitr]) {
  //       std::cout << p.first << " ";
  //     }
  //     std::cout << std::endl;
  // }

  {
  //
  // 4)  Compute distributions
  std::map<uint64_t, uint64_t> local_orig_degree;
  std::map<uint64_t, uint64_t> local_2core_degree;
  std::map<uint64_t, uint64_t> local_2core_out_degree;

  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    local_orig_degree[g.degree(*vitr)]++;
    local_2core_degree[core2_degree[*vitr]]++;
    local_2core_out_degree[core2_directed[*vitr].size()]++;

  }
  for (auto citr = g.controller_begin(); citr != g.controller_end(); ++citr) {
    local_orig_degree[g.degree(*citr)]++;
    local_2core_degree[core2_degree[*citr]]++;
    local_2core_out_degree[core2_directed[*citr].size()]++;
  }

  std::stringstream orig_degree, tcore_degree, tcore_out_degree;
  orig_degree << deg_output_fname << "_orig_degree.txt";
  tcore_degree << deg_output_fname << "_2core_degree.txt";
  tcore_out_degree << deg_output_fname << "_2core_out_degree.txt";

  output_degree_distribution(local_orig_degree, orig_degree.str().c_str());
  output_degree_distribution(local_2core_degree, tcore_degree.str().c_str());
  output_degree_distribution(local_2core_out_degree, tcore_out_degree.str().c_str());


  }

  //
  // 3)  Build wedges & count
  uint64_t local_triangle_count(0), local_wedge_count(0);
  start_time = MPI_Wtime();
  {
    auto alg_data = std::forward_as_tuple(core2_directed, local_triangle_count, local_wedge_count);
    auto vq = create_visitor_queue<core2_wedges<graph_type>, detail::visitor_priority_queue>(&g, alg_data);
    vq.init_visitor_traversal_new();
  }
  end_time = MPI_Wtime();
  uint64_t global_wedge_count = mpi_all_reduce(local_wedge_count,std::plus<uint64_t>(), MPI_COMM_WORLD);
  if(mpi_rank == 0) {
    std::cout << "TC on directed 2core time = " << end_time - start_time << std::endl;
    std::cout << "Total wedges checked = " << global_wedge_count << std::endl;
  }


  uint64_t global_triangle_count = mpi_all_reduce(local_triangle_count,std::plus<uint64_t>(), MPI_COMM_WORLD);
  return global_triangle_count;

//  for (auto ditr = g.delegate_vertices_begin(); ditr != g.delegate_vertices_end(); ++ditr) {
//    std::cout << "Delegate " << g.locator_to_label(*ditr) << " has core2_degree = " <<  core2_degree[*ditr] << std::endl;
//  }



  /*
  //
  // 2)  count neighbors with lower degree
  typename graph_type::template vertex_data<uint32_t, std::allocator<uint32_t> >   vertex_count_lower_degree(g);
  for (auto vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    for(auto eitr = g.edges_begin(*vitr); eitr != g.edges_end(*vitr); ++eitr) {
      auto neighbor = eitr.target();
      uint32_t neighbor_degree = edge_target_degree[eitr];
      //std::cout << "neighbor_degree = " << neighbor_degree << std::endl;
      if(neighbor_degree < vertex_degree[*vitr]) {
        ++vertex_count_lower_degree[*vitr];
      } else if(neighbor_degree == vertex_degree[*vitr] && g.locator_to_label(neighbor) < g.locator_to_label(*vitr)) {
        ++vertex_count_lower_degree[*vitr];
      }
    }
  }
  for (auto ditr = g.delegate_vertices_begin(); ditr != g.delegate_vertices_end(); ++ditr) {
    for(auto eitr = g.edges_begin(*ditr); eitr != g.edges_end(*ditr); ++eitr) {
      auto neighbor = eitr.target();
      uint32_t neighbor_degree = edge_target_degree[eitr];
      if(neighbor_degree < vertex_degree[*ditr]) {
        ++vertex_count_lower_degree[*ditr];
      } else if(neighbor_degree == vertex_degree[*ditr] && g.locator_to_label(neighbor) < g.locator_to_label(*ditr)) {
        ++vertex_count_lower_degree[*ditr];
      }
    }
  }
  vertex_count_lower_degree.all_reduce();
  //
  // 3)  store lower degree count on edges
  typename graph_type::template edge_data<uint32_t, std::allocator<uint32_t> >   edge_target_count_lower_degree(g);
  {
    auto alg_data = std::forward_as_tuple(vertex_count_lower_degree, edge_target_count_lower_degree);
    auto vq = create_visitor_queue<clone_vertex_data_visitor<graph_type, uint32_t>, detail::visitor_priority_queue>(&g, alg_data);
    vq.init_visitor_traversal_new();
  }

  //
  // 4)
  {
    auto alg_data = std::forward_as_tuple(vertex_count_lower_degree, vertex_degree, edge_target_degree, edge_target_count_lower_degree, tc_data);
    auto vq = create_visitor_queue<triangle_count_visitor<graph_type>, detail::visitor_priority_queue>(&g, alg_data);
    vq.init_visitor_traversal_new();
  }

  return tc_data.tc_dataumulate();  */
  //return 0;
}

} //end namespace havoqgt

#endif //HAVOQGT_MPI_TRIANGLE_COUNT_HPP_INCLUDED
