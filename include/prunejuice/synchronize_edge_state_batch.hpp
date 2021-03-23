#pragma once

//#include <array>
//#include <chrono>
#include <deque>
//#include <limits>
//#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace prunejuice { namespace indexing {

template<typename Visitor>
class sync_edge_state_queue {

public:
  sync_edge_state_queue() {}

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
class sync_edge_state_visitor {

public:

  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  sync_edge_state_visitor() {}

  sync_edge_state_visitor(vertex_locator _vertex) :
    vertex(_vertex) {}

  sync_edge_state_visitor(vertex_locator _vertex,
    vertex_locator _parent, Vertex _third_vertex_data) : 
    vertex(_vertex), parent(_parent), 
    third_vertex_data(_third_vertex_data) {}

  ~sync_edge_state_visitor() {}

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {

    int mpi_rank = havoqgt::comm_world().rank();

    auto g = std::get<0>(alg_data);

    if (vertex.is_delegate() && g->master(vertex) != mpi_rank) {
      return true;
    }

    //typedef std::unordered_map<Vertex, uint8_t> VertexUintMap;

    // u --w-- v   
    {
    auto find_edge = std::get<1>(alg_data)[vertex].find(g->locator_to_label(parent));
    if (find_edge == std::get<1>(alg_data)[vertex].end()) {
      auto insert_status = std::get<1>(alg_data)[vertex].insert({g->locator_to_label(parent), VertexUintMap()});
      if (!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return false;
      }

      find_edge = insert_status.first;
      find_edge->second.insert({third_vertex_data, 1});
    } else {
      auto find_w = find_edge->second.find(third_vertex_data);
      if (find_w == find_edge->second.end()) {
        find_edge->second.insert({third_vertex_data, 1});
      }
    }  
    }
    
    // not required when veretx data is stored on the edge        
    // u --v-- w
    /*{  
    auto find_edge = std::get<1>(alg_data)[vertex].find(third_vertex_label);
    if (find_edge == std::get<1>(alg_data)[vertex].end()) {
      auto insert_status = std::get<1>(alg_data)[vertex].insert({third_vertex_label, VertexUintMap()});
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
    }*/

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

    // u = vertex, v = neighbor, u --w--> v 
    for (auto& v : std::get<1>(alg_data)[vertex]) {
      vertex_locator neighbor = g.label_to_locator(v.first);
      for (auto&  w : v.second) {
        //vertex_locator neighbor = g.label_to_locator(w.first);
        sync_edge_state_visitor new_visitor(neighbor, vertex, w.first); // TODO: send in chunks
        vis_queue->queue_visitor(new_visitor);
      }  
    }  
 
    /*for (auto& v : std::get<1>(alg_data)[vertex]) {
      //vertex_locator neighbor = g.label_to_locator(v.first);
      for (auto&  w : v.second) {
        vertex_locator neighbor = g.label_to_locator(w.first);
        sync_edge_state_visitor new_visitor(neighbor, vertex, v.first);
        vis_queue->queue_visitor(new_visitor); 
      }  
    }*/  
 
    return false;
  }  

  friend inline bool operator>(const sync_edge_state_visitor& v1, 
    const sync_edge_state_visitor& v2) {
    return false;
  }

  vertex_locator vertex;
  vertex_locator parent;
  Vertex third_vertex_data;
  //uint8_t msg_type; // 0 - init, 1 - forward

}; // class sync_edge_state_visitor

template <typename TGraph, typename Vertex, typename Edge, typename VertexData,
  typename EdgeData, typename VertexMetadata, typename EdgeMetadata,
  typename VertexActive, typename VertexUint8MapCollection,
  typename TemplateVertex, typename VertexStateMap, typename PatternGraph,
  typename PatternUtilities, typename VertexUint8Map,
  typename VertexSetCollection,
  template<typename> class DelegateGraphVertexDataSTDAllocator,
  typename Boolean, typename BitSet, typename VertexNonlocalConstraintMatches,
  typename EdgeStateCollection>

void synchronize_edge_state_batch(TGraph* g, VertexMetadata& vertex_metadata,
  VertexActive& vertex_active,
  VertexUint8MapCollection& vertex_active_edges_map,
  TemplateVertex& template_vertices, VertexStateMap& vertex_state_map,  
  VertexSetCollection& vertex_token_source_set,  
  std::uint64_t batch_size,
  uint64_t& message_count,  
  VertexNonlocalConstraintMatches& vertex_nlc_matches,
  EdgeStateCollection& vertex_active_edge_state_map) {


  int mpi_rank = havoqgt::comm_world().rank();
  int mpi_size = havoqgt::comm_world().size();

  if (mpi_rank == 0) {
    std::cout << "Synchrnizing Edge State ... " << std::endl;
  }

  typedef typename TGraph::vertex_locator vertex_locator;
  typedef typename TGraph::vertex_iterator vertex_iterator;
 
  typedef sync_edge_state_visitor<TGraph, Vertex, VertexUint8Map> visitor_type;

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
 
    auto alg_data = std::forward_as_tuple(g, vertex_active_edge_state_map);

    auto vq = havoqgt::create_visitor_queue<visitor_type, 
      sync_edge_state_queue>(g, alg_data); 

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
  
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

}

////////////////////////////////////////////////////////////////////////////////

template <typename TGraph, typename Vertex, typename Edge, typename VertexData,
  typename EdgeData, typename VertexMetadata, typename EdgeMetadata,
  typename VertexActive, typename VertexUint8MapCollection,
  typename TemplateVertex, typename VertexStateMap, typename PatternGraph,
  typename PatternUtilities, typename VertexUint8Map,
  typename VertexSetCollection,
  template<typename> class DelegateGraphVertexDataSTDAllocator,
  typename Boolean, typename BitSet, typename VertexNonlocalConstraintMatches,
  typename EdgeStateCollection>

void synchronize_edge_state_batch(TGraph* g, VertexMetadata& vertex_metadata,
  VertexActive& vertex_active,
  VertexUint8MapCollection& vertex_active_edges_map,
  TemplateVertex& template_vertices, VertexStateMap& vertex_state_map,
  PatternGraph& pattern_graph, PatternUtilities& pattern_utilities, size_t pl,
  VertexUint8Map& token_source_map,
  VertexSetCollection& vertex_token_source_set,
  std::vector<uint8_t>::reference pattern_found,
  std::uint64_t batch_size,
  std::ofstream& paths_result_file, uint64_t& message_count,
  uint64_t pattern_constraint,
  VertexNonlocalConstraintMatches& vertex_nlc_matches,
  EdgeStateCollection& vertex_active_edge_state_map) {


  int mpi_rank = havoqgt::comm_world().rank();
  int mpi_size = havoqgt::comm_world().size();

  if (mpi_rank == 0) {
    std::cout << "Synchrnizing Edge State ... " << std::endl;
  }

  typedef typename TGraph::vertex_locator vertex_locator;
  typedef typename TGraph::vertex_iterator vertex_iterator;
 
  typedef sync_edge_state_visitor<TGraph, Vertex, VertexUint8Map> visitor_type;

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
 
    auto alg_data = std::forward_as_tuple(g, vertex_active_edge_state_map);

    auto vq = havoqgt::create_visitor_queue<visitor_type, 
      sync_edge_state_queue>(g, alg_data); 

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
  
  //MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

}

}} // end namespace prunejuice::indexing
