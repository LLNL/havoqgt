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
 

#ifndef HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED
#define HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED



#include <havoqgt/visitor_queue.hpp>
#include <boost/container/deque.hpp>
#include <vector>

namespace havoqgt { namespace mpi {

typedef double edge_data_type;
typedef std::tuple<std::pair<uint64_t, uint64_t>, edge_data_type> edge_type;

template <typename Visitor>
class pr_queue
{

protected:
  std::vector< Visitor > m_data;
public:
  pr_queue() { }

  bool push(Visitor const & task)
  {
    m_data.push_back(task);
    return true;
  }

  void pop()
  {
    m_data.pop_back();
  }

  Visitor const & top() //const
  {
    return m_data.back();
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
class pr_visitor {
public:
  typedef typename Graph::vertex_locator                 vertex_locator;
  pr_visitor(): rank(std::numeric_limits<double>::min())  { }

  pr_visitor(vertex_locator _vertex, double _rank)
    : vertex(_vertex)
    , rank(_rank) { }

  pr_visitor(vertex_locator _vertex)
    : vertex(_vertex)
    , rank(std::numeric_limits<double>::min()) { }      

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if(rank == std::numeric_limits<double>::min()) {
      HAVOQGT_ERROR_MSG("This is a damn logic error!");
      return true;
    }
    std::get<1>(alg_data)[vertex] += rank; //change to next_rank 
    return false;
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //std::cout << "PageRank visit: vertex: " << vertex.local_id() << " degree: " << g.degree(vertex) 
    //          << " local degree: " << g.local_degree(vertex) << std::endl;
    typedef typename Graph::edge_iterator eitr_type;
    //eitr_type eitrA = g.edges_begin(vertex);
    //std::cout << vertex.local_id() << " Neighbour :" << eitrA.source().local_id() << std::endl;
    //eitr_type eitrB = g.edges_end(vertex);
    //std::cout << vertex.local_id() << " Neighbour :" << eitrB.source().local_id() << std::endl;
    
    //for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
    //  vertex_locator neighbor = eitr.target();
    //  std::cout << vertex.local_id() << "Neighbour :" << neighbor.local_id() << std::endl; 
    //}
    
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    //std::cout << "MPI Rank -> " << mpi_rank << " initiating visit ...... "<< std::endl;

    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //change to cur_rank
    double old_rank = std::get<0>(alg_data)[vertex];    
    uint64_t degree = g.degree(vertex);
    double send_rank = old_rank / double(degree);

    //std::cout << "PageRank visit: vertex: " << vertex.local_id() << " degree: " << degree << std::endl;

    typedef typename Graph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      int mpi_rank(0);
      CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

      double edge_data = (double)eitr.edge_data(); 
      //if (vertex.is_delegate()) {
        //std::cout << "# MPI Rank -> " << mpi_rank << " inside loop body " << std::endl;
        std::cout << "MPI Rank -> " << mpi_rank << " Source: " << g.locator_to_label(vertex) << " is delegate " << vertex.is_delegate() << " Neighbour : " << g.locator_to_label(neighbor) << " Edge data: " << edge_data << std::endl;
      //}
      pr_visitor new_visitor( neighbor, send_rank);
      vis_queue->queue_visitor(new_visitor);       
    }
    return true;
  }


  friend inline bool operator>(const pr_visitor& v1, const pr_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const pr_visitor& v1, const pr_visitor& v2) {
    return false;
  }

  vertex_locator   vertex;
  double           rank;
};

template <typename TGraph, typename PRData>
void page_rank(TGraph& g, PRData& cur_rank, PRData& next_rank, bool initial) {
  typedef  pr_visitor<TGraph>    visitor_type;

  int mpi_rank(0);
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  std::cout << "MPI Rank -> " << mpi_rank << " starting process ...... "<< std::endl;

  // gather local edges in the rank 

  //std::vector<edge_type> local_edges(g.num_local_edges());
  //std::vector<edge_type> global_edges(44);

  std::vector<uint64_t> vec_local_mpi_ranks;
  std::vector<uint64_t> vec_global_mpi_ranks; 
  std::vector<uint64_t> vec_local_edge_source;//(g.num_local_edges());
  std::vector<uint64_t> vec_global_edge_source;//(44);
  std::vector<uint64_t> vec_local_edge_target;//(g.num_local_edges());
  std::vector<uint64_t> vec_global_edge_target;//(44);
  std::vector<uint64_t> vec_local_edge_data;//(g.num_local_edges());
  std::vector<uint64_t> vec_global_edge_data;//(44);

  //std::cout << "MPI Rank -> " << mpi_rank << " #local edges " << local_edges.size() << std::endl; 

  typedef typename TGraph::vertex_iterator vitr_type;   
  typedef typename TGraph::vertex_locator vloc_type;
  for (vitr_type vitr = g.vertices_begin(); vitr != g.vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    //std::cout << "MPI Rank -> " << mpi_rank << " local vertex " << g.locator_to_label(vertex) << std::endl; 
    //create_local_edge_list(g, vertex, local_edges); 
    create_local_edge_data_list(g, vertex, vec_local_edge_source, vec_local_edge_target, vec_local_edge_data, vec_local_mpi_ranks); 
  }

  //typedef typename TGraph::controller_iterator citr_type; 
  //for (citr_type citr = g.controller_begin(); citr != g.controller_end(); ++citr) {
  for (vitr_type vitr = g.delegate_vertices_begin(); vitr != g.delegate_vertices_end(); ++vitr) { 
    vloc_type vertex = *vitr;
    //std::cout << "MPI Rank -> " << mpi_rank << " controller vertex " << g.locator_to_label(vertex) << std::endl;
    //create_local_edge_list(g, vertex, local_edges); 
    create_local_edge_data_list(g, vertex, vec_local_edge_source, vec_local_edge_target, vec_local_edge_data, vec_local_mpi_ranks);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "MPI Rank -> " << mpi_rank << " vec_local_edge_data.size " << vec_local_edge_data.size() << std::endl;
    for (size_t i; i < vec_local_edge_data.size(); ++i) {
      //std::pair<uint64_t, uint64_t> edge = std::get<0>(local_edges[i]);
      uint64_t ed = vec_local_edge_data[i];//std::get<1>(local_edges[i]);
      uint64_t s = vec_local_edge_source[i];//edge.first;
      uint64_t t = vec_local_edge_target[i];//edge.second;  
      //std::cout << "MPI Rank -> " << mpi_rank << " local " << s << " " << t << " " << ed << std::endl;
    }
  }

  // gather local edges in the rank 
  
  // gather global edges
  mpi_all_gather(vec_local_edge_source, vec_global_edge_source, MPI_COMM_WORLD); 
  mpi_all_gather(vec_local_edge_target, vec_global_edge_target, MPI_COMM_WORLD); 
  mpi_all_gather(vec_local_edge_data, vec_global_edge_data, MPI_COMM_WORLD);    
  mpi_all_gather(vec_local_mpi_ranks, vec_global_mpi_ranks, MPI_COMM_WORLD);
  // gather global edges
  
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> vec_input_graph;
  
  vec_input_graph.push_back(std::make_tuple(0, 1, 10));
  vec_input_graph.push_back(std::make_tuple(0, 5, 5));
/*0 1 10
0 5 5
1 0 21
1 2 31
1 6 51
2 1 16
2 3 19
2 7 53
3 2 12
3 4 10
3 8 57
4 3 13
4 9 55
5 0 57
5 6 18
5 10 59
6 1 51
6 5 12
6 7 14
6 11 50
7 2 55
7 6 14
7 8 13
7 12 51
8 3 52
8 7 13
8 9 14
8 13 55
9 4 55
9 8 16
9 14 57
10 5 58
10 11 19
11 6 54
11 10 13
11 12 12
12 7 51
12 11 12
12 13 13
13 8 54
13 12 15
13 14 18
14 9 53
14 13 19*/

  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> vec_global_edges;  
  
  if (mpi_rank == 0) {

    
  

    std::cout << "MPI Rank -> " << mpi_rank << " vec_global_edge_data.size " << vec_global_edge_data.size() << std::endl;
    for (size_t i; i < vec_global_edge_data.size(); ++i) {
      //std::pair<uint64_t, uint64_t> edge = std::get<0>(global_edges[i]);
      uint64_t ed = vec_global_edge_data[i]; //std::get<1>(global_edges[i]);
      uint64_t s = vec_global_edge_source[i]; //edge.first;
      uint64_t t = vec_global_edge_target[i]; //edge.second;
      uint64_t mr = vec_global_mpi_ranks[i]; 
      std::cout << "MPI Rank -> " << mr << " global " << s << " " << t << " " << ed << std::endl;
      vec_global_edges.push_back(std::make_tuple(s, t, ed));
    }
  

  std::stable_sort(vec_global_edges.begin(),vec_global_edges.end(),
       [](const std::tuple<uint64_t, uint64_t, uint64_t>& a,
       const std::tuple<uint64_t, uint64_t, uint64_t>& b) -> bool
       {
         return std::get<0>(a) < std::get<0>(b);
       });

  for (size_t i = 0; i < vec_global_edges.size(); i++) {
    uint64_t s = std::get<0>(vec_global_edges[i]);
    uint64_t t = std::get<1>(vec_global_edges[i]);
    uint64_t ed = std::get<2>(vec_global_edges[i]);
    std::cout << s << " " << t << " " << ed << std::endl;
  }

  // compare container sizes
  // compare element by element 
  for (size_t i = 0; i < vec_input_graph.size(); ++i ) {
    if (vec_input_graph[i] == vec_global_edges[i]) {
      std::cout << i << " Pass!" << std::endl;
    } else {
      std::cout << i << " Fail!" << std::endl;  
    }
  }
  


  } // mpi_rank = 0

  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "MPI Rank -> " << mpi_rank << " Barrier ...... Barrier ...... Barrier ...... Barrier ...... "<< std::endl;

  MPI_Barrier(MPI_COMM_WORLD); 
 
  auto alg_data = std::forward_as_tuple(cur_rank, next_rank);
  
  if(initial) {
    cur_rank.reset(double(1)/double(g.max_global_vertex_id()));
  }
   
  auto vq = create_visitor_queue<visitor_type, detail::visitor_priority_queue>(&g, alg_data);
  vq.init_visitor_traversal_new();
  next_rank.all_reduce();
}

template <typename TGraph, typename vertex_locator>
void create_local_edge_data_list(TGraph& g, vertex_locator vloc, std::vector<uint64_t>& edge_source, 
                                 std::vector<uint64_t>& edge_target, std::vector<uint64_t>& edge_data, std::vector<uint64_t>& mpi_ranks) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
//    std::cout << "## MPI Rank -> " << mpi_rank << " Source: " << g.locator_to_label(vloc) << " is delegate " << vloc.is_delegate() << std::endl;
    typedef typename TGraph::edge_iterator eitr_type;
    vertex_locator vertex = g.label_to_locator(g.locator_to_label(vloc));  
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
//      std::cout << "# MPI Rank -> " << mpi_rank << " inside loop body " << std::endl;
      vertex_locator neighbor = eitr.target();
      uint64_t ed = (uint64_t)eitr.edge_data();
//      std::cout << "# MPI Rank -> " << mpi_rank << " Source: " << g.locator_to_label(vertex) << " is target delegate " << neighbor.is_delegate() << " Neighbour : " << g.locator_to_label(neighbor) << " Edge data: " << ed << std::endl;
      edge_source.push_back(g.locator_to_label(vertex));
      edge_target.push_back(g.locator_to_label(neighbor));
      edge_data.push_back(ed);
      mpi_ranks.push_back(mpi_rank);
    }
}

template <typename TGraph, typename vertex_locator>
void create_local_edge_list(TGraph& g, vertex_locator vertex, std::vector<edge_type>& local_edges) {
    int mpi_rank(0);
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    //std::cout << "## MPI Rank -> " << mpi_rank << " Source: " << g.locator_to_label(vertex) << " is delegate " << vertex.is_delegate() << std::endl;
    typedef typename TGraph::edge_iterator eitr_type;
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      double edge_data = (double)eitr.edge_data();
      std::cout << "# MPI Rank -> " << mpi_rank << " Source: " << g.locator_to_label(vertex) << " is delegate " << vertex.is_delegate() << " Neighbour : " << g.locator_to_label(neighbor) << " Edge data: " << edge_data << std::endl;	

      local_edges.push_back(std::forward_as_tuple(
                            std::make_pair(g.locator_to_label(vertex), 
                                           g.locator_to_label(neighbor)),
                            edge_data) );	
    }

}

}} //end namespace havoqgt::mpi




#endif //HAVOQGT_MPI_PAGE_RANK_HPP_INCLUDED
