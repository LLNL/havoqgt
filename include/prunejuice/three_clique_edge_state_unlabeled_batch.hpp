#pragma once

//#include <array>
//#include <chrono>
#include <deque>
//#include <limits>
//#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace prunejuice { namespace indexing {

const size_t max_msg_items = 16;
static size_t max_msg_count = 16;
static uint64_t pattern_count = 0;

template<typename Visitor>
class three_clique_queue {

public:

  three_clique_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    //data.pop_back();
    data.pop_front();
  }

  Visitor const& top() {
    //return data.back();
    return data.front();
  }

  size_t size() const {
    return data.size();;
  }

  bool empty() const {
    return data.empty();
  }

  void clear() {
    data.clear();
  }

protected:

  //std::vector<Visitor> data;
  std::deque<Visitor> data;
};

// visitor class
template<typename Graph, typename Vertex, typename VertexUintMap>
class three_clique_visitor {

public:

  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  three_clique_visitor() {}

  three_clique_visitor(vertex_locator _vertex) :
    vertex(_vertex) {}

  template <typename NeighborArrayStatic>
  three_clique_visitor(vertex_locator _vertex,
    vertex_locator _parent, NeighborArrayStatic& _parent_neighbors) : //, Vertex _third_vertex_label) : 
    vertex(_vertex), parent(_parent) //, 
    //third_vertex_label(_third_vertex_label) 
    {      
      std::copy(std::begin(_parent_neighbors), std::end(_parent_neighbors),  
        std::begin(parent_neighbors));    
    }

  ~three_clique_visitor() {}

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {

    int mpi_rank = havoqgt::comm_world().rank();

    auto g = std::get<0>(alg_data);

    if (vertex.is_delegate() && g->master(vertex) != mpi_rank) {
      return true;
    }

    //typedef std::unordered_map<Vertex, uint8_t> VertexUintMap;

    for (auto& q : parent_neighbors) {

      auto find_q = std::get<1>(alg_data)[vertex].find(q);
      if (find_q != std::get<1>(alg_data)[vertex].end()) {
        pattern_count++;
      
        Vertex third_vertex_label = q; 

    // u --w-- v   
    {
    auto find_edge = std::get<2>(alg_data)[vertex].find(g->locator_to_label(parent));
    if (find_edge == std::get<2>(alg_data)[vertex].end()) {
      auto insert_status = std::get<2>(alg_data)[vertex].insert({g->locator_to_label(parent), VertexUintMap()});
      if (!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return false;
      }

      find_edge = insert_status.first;
      find_edge->second.insert({third_vertex_label, 1});
    } else {
      auto find_w = find_edge->second.find(third_vertex_label);
      if (find_w == find_edge->second.end()) {
        find_edge->second.insert({third_vertex_label, 1});
      }
    }  
    }
           
    // u --v-- w
    {  
    auto find_edge = std::get<2>(alg_data)[vertex].find(third_vertex_label);
    if (find_edge == std::get<2>(alg_data)[vertex].end()) {
      auto insert_status = std::get<2>(alg_data)[vertex].insert({third_vertex_label, VertexUintMap()});
      if (!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return false;
      }

      find_edge = insert_status.first;
      find_edge->second.insert({g->locator_to_label(parent), 1});
    } else {
      auto find_w = find_edge->second.find(g->locator_to_label(parent));
      if (find_w == find_edge->second.end()) {
        find_edge->second.insert({g->locator_to_label(parent), 1});
      }
    }  
    }

      } // if
    } // for

    // Test
    //if (mpi_rank == 0) {
    //  std::cout << std::endl;  
    //  for (size_t i = 0; i < parent_neighbors.size(); i++) {
    //    std::cout << parent_neighbors[i] << " "; 
    //  }     
    //}
    // Test

    //for (auto& q : parent_neighbors) { 
    //  auto find_q = std::get<1>(alg_data)[vertex].find(q);
    //  if (find_q != std::get<1>(alg_data)[vertex].end()) {
    //    pattern_count++;
    //  }   
    //}

    //pattern_count++;

    return false; 
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    //if(!std::get<13>(alg_data)[vertex]) {
    //  return false;
    //} else {
      return visit(g, vis_queue, alg_data);
    //}
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    
    int mpi_rank = havoqgt::comm_world().rank();

    if (vertex.is_delegate() && g.master(vertex) != mpi_rank) {
      return false;
    }  

    size_t vertex_degree = std::get<1>(alg_data)[vertex].size(); 

    for (auto& p : std::get<1>(alg_data)[vertex]) {

      auto find_p = std::get<3>(alg_data).find(p.first);
      if (find_p == std::get<3>(alg_data).end()) {
        std::cerr << "Error: did not find the expected item in the map." 
        << std::endl;
        return false;
      } else if (vertex_degree > find_p->second) {
        continue;
      } else if ((vertex_degree == find_p->second) && (g.locator_to_label(vertex) > p.first)) {
        continue;
      } 

      vertex_locator neighbor = g.label_to_locator(p.first);

      size_t msg_count = 0;
      size_t item_count = 0;

      std::array<Vertex, max_msg_items> q_neighbors;
      q_neighbors.fill(std::numeric_limits<std::uint64_t>::max()); 

      for (auto& q : std::get<1>(alg_data)[vertex]) {
        if (p.first == q.first) {
          item_count++;
          continue;
        }

        auto find_q = std::get<3>(alg_data).find(q.first);
        if (find_q == std::get<3>(alg_data).end()) {
           std::cerr << "Error: did not find the expected item in the map."
          << std::endl;
          return false;
        } else if (find_p->second > find_q->second) {
          item_count++;
          continue;
        } else if ((find_p->second == find_q->second) && (p.first > q.first)) {
          item_count++;
          continue;
        }

        // Test
        //q_neighbors[msg_count] = q.first;
        //three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
        //vis_queue->queue_visitor(new_visitor);
        //q_neighbors.fill(std::numeric_limits<std::uint64_t>::max());
        // Test 

        // send in chunks  
        if (msg_count >= (max_msg_items - 1)) {
          // add
          q_neighbors[msg_count] = q.first; 
          // send
          three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
          vis_queue->queue_visitor(new_visitor);
          // reset
          q_neighbors.fill(std::numeric_limits<std::uint64_t>::max()); 
          msg_count = 0;  
        } //else if (remaining < 1) {
          // add
          //q_neighbors[msg_count] = q.first;  
          // send
          //three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
          //vis_queue->queue_visitor(new_visitor); 
          //q_neighbors.fill(std::numeric_limits<std::uint64_t>::max()); 
          //break; 
        //} 
        else if (msg_count < (max_msg_items - 1)) {
          // add
          q_neighbors[msg_count] = q.first;
          //if (item_count >= (vertex_degree - 2)) {
          //  three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
          //  vis_queue->queue_visitor(new_visitor);  
          //} 
          msg_count++;      
        }
         
        item_count++;   

        /*if ((item_count == vertex_degree - 1) && (msg_count < max_msg_count - 1)) {
          //q_neighbors[msg_count] = q.first;  
          three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
          vis_queue->queue_visitor(new_visitor);
        }*/  
         
      } // for

      if ((msg_count > 0) && (msg_count <= max_msg_items - 1)) {
        three_clique_visitor new_visitor(neighbor, vertex, q_neighbors);
        vis_queue->queue_visitor(new_visitor);
      } 

      assert(item_count == (vertex_degree - 1)); 

    } // for   

    return false;
  }  

  friend inline bool operator>(const three_clique_visitor& v1, 
    const three_clique_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator parent;
  //Vertex third_vertex_label;
  std::array<Vertex, max_msg_items> parent_neighbors;
  //uint8_t msg_type; // 0 - init, 1 - forward

}; // class three_clique_visitor

template <typename TGraph, typename Vertex, typename Edge, typename VertexData,
  typename EdgeData, typename VertexMetadata, typename EdgeMetadata,
  typename VertexActive, typename VertexUint8MapCollection,
  typename TemplateVertex, typename VertexStateMap, typename PatternGraph,
  typename PatternUtilities, typename VertexUint8Map,
  typename VertexSetCollection,
  template<typename> class DelegateGraphVertexDataSTDAllocator,
  typename Boolean, typename BitSet, typename VertexNonlocalConstraintMatches,
  typename EdgeStateCollection, typename VertexSizeMap>

void three_clique_edge_state_batch(TGraph* g, VertexMetadata& vertex_metadata,
  VertexActive& vertex_active,
  VertexUint8MapCollection& vertex_active_edges_map,
  TemplateVertex& template_vertices, VertexStateMap& vertex_state_map,   
  VertexSetCollection& vertex_token_source_set,
  std::uint64_t batch_size,
  uint64_t& message_count,  
  VertexNonlocalConstraintMatches& vertex_nlc_matches,
  EdgeStateCollection& vertex_active_edge_state_map,
  VertexSizeMap& active_edge_target_vertex_degree_map) {


  int mpi_rank = havoqgt::comm_world().rank();
  int mpi_size = havoqgt::comm_world().size();

  if (mpi_rank == 0) {
    std::cout << "Three Clique Edge State ... " << std::endl;
  }

  typedef typename TGraph::vertex_locator vertex_locator;
  typedef typename TGraph::vertex_iterator vertex_iterator;
 
  typedef three_clique_visitor<TGraph, Vertex, VertexUint8Map> visitor_type;

  //////////////////////////////////////////////////////////////////////////////

  // batching

  uint64_t max_batch_size = mpi_size;
  uint64_t max_ranks_per_itr = batch_size;
  //uint64_t max_ranks_per_itr = 10;
  //uint64_t max_ranks_per_itr = mpi_size / 36; // quartz

  uint64_t batch_count = 0;

  if (mpi_rank == 0) {
    std::cout << " - | Batch Size : "
      << max_ranks_per_itr  
      << std::endl;
  }
   
  // batch processing

  for (auto batch_mpi_rank = 0;  batch_mpi_rank < max_batch_size;
    batch_mpi_rank+=max_ranks_per_itr) {

    auto batch_max_mpi_rank = (batch_mpi_rank + max_ranks_per_itr - 1) <= (mpi_size - 1) ?
      (batch_mpi_rank + max_ranks_per_itr - 1) : (mpi_size - 1);

    double time_start = MPI_Wtime();
    
    assert(batch_mpi_rank + max_ranks_per_itr <= mpi_size);
   
    if (mpi_rank == 0) {
      std::cout << " - | Batch #" << batch_count
      << " | MPI Ranks : " << batch_mpi_rank << " - " << batch_max_mpi_rank
      <<  " | Asynchronous Traversal ..." << std::endl;
    }

    ////////////////////////////////////////////////////////////////////////////

    // asynchronous traversal
 
    auto alg_data = std::forward_as_tuple(g, vertex_active_edges_map, 
      vertex_active_edge_state_map, active_edge_target_vertex_degree_map);

    auto vq = havoqgt::create_visitor_queue<visitor_type, 
      three_clique_queue>(g, alg_data); 

    ///vq.init_visitor_traversal_new();
    //vq.init_visitor_traversal_new_batch();
    //vq.init_visitor_traversal_new_alt();
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

    // asynchronous traversal

    ////////////////////////////////////////////////////////////////////////////
  
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << " - | Batch #" << batch_count
      << " | MPI Ranks : " << batch_mpi_rank << " - " << batch_max_mpi_rank      
      <<  " | Time : " << time_end - time_start << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

    batch_count++;

  } // for

  // batch processing

  // batching

  //////////////////////////////////////////////////////////////////////////////
  
  uint64_t global_pattern_count = havoqgt::mpi_all_reduce(pattern_count,
    std::plus<uint64_t>(), MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "Three Clique | Global Pattern Count : "
      << global_pattern_count << std::endl;
  }
  
  //return global_pattern_count;

  //MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

}

}} // end namespace prunejuice::indexing
