#pragma once

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <independent_set/utility.hpp>

namespace maximal_independent_set {

static uint8_t not_finished = 0; // true 1, false 0 

template<typename Visitor>
class mis_queue {

public:

  mis_queue() {}

  bool push(Visitor const& element) {
    data.push_back(element);
    return true;
  }

  void pop() {
    data.pop_back();	   
  } 	   

  Visitor const& top() {
    return data.back();
  }

  size_t size() const {
    return data.size();	   
  }	  

  bool empty() const {
    return data.empty();	  
  }	  

  void clear() {
    data.clear();	  
  }	  

protected:

  std::vector<Visitor> data;      
};

// visitor 
template<typename Graph, typename VertexID, typename VertexPriority>
class mis_visitor {

public:

  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  mis_visitor()
    {}

  mis_visitor(vertex_locator _vertex) : 
    vertex(_vertex) {}

  mis_visitor(vertex_locator _vertex, VertexID _neighbor_ID, 
    VertexPriority _neighbor_priority) :
    vertex(_vertex), neighbor_ID(_neighbor_ID), 
    neighbor_priority(_neighbor_priority) {}

  ~mis_visitor() {}

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    	  
    //int mpi_rank = havoqgt::comm_world().rank();
    //if (mpi_rank == 0) { 
    //  std::cout << neighbor_priority << std::endl;	    
    //}

    // temporary hack
    // delegates will not get here, they are not in the neighbor maps

    //auto graph = std::get<0>(alg_data);
    auto& higher_priority_neighbor_map = std::get<2>(alg_data)[vertex];
    auto& lower_priority_neighbor_map = std::get<3>(alg_data)[vertex];

    if ( (std::get<1>(alg_data)[vertex] == maximal_independent_set::in) || 
      (std::get<1>(alg_data)[vertex] == maximal_independent_set::out) ) {
      return false; 	    
    } else {
      // vertex is undecided
     
      if (neighbor_priority == maximal_independent_set::in) {
        std::get<1>(alg_data)[vertex] = maximal_independent_set::out;
        // TODO: notify neighbors in the next round before empting maps 
      }
      
      // remove neighbor from neighbor maps, neighbor is in or out 

      if (higher_priority_neighbor_map.size() > 0) {
        auto erase_status_h = higher_priority_neighbor_map.erase(neighbor_ID);
	if (erase_status_h < 1) {
	  //std::cerr << "Error: failed to remove an element from the map." << 
          //  std::endl;  	
	}	
      }
      
      if (lower_priority_neighbor_map.size() > 0) {
        auto erase_status_l = lower_priority_neighbor_map.erase(neighbor_ID);	      
        if (erase_status_l < 1) {
          //std::cerr << "Error: failed to remove an element from the map." << 
          //  std::endl; 
        }
      }

      not_finished = 1;
    } // else  

    return false;	  
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }    

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
	  
    //int mpi_rank = havoqgt::comm_world().rank();
    //if (mpi_rank == 0) {
    //  std::cout << g.locator_to_label(vertex) << " " <<
    //    std::get<1>(alg_data) << std::endl;   	    
    //}
   
    // std::get<1>(alg_data) vertex priority list
    auto& higher_priority_neighbor_map = std::get<2>(alg_data)[vertex];
    auto& lower_priority_neighbor_map = std::get<3>(alg_data)[vertex];

    if (std::get<1>(alg_data)[vertex] == maximal_independent_set::in) {
      return false;  
    } else if (std::get<1>(alg_data)[vertex] == maximal_independent_set::out) {
      
      // temporary hack
      if (vertex.is_delegate()) { 
        return false;	      
      }  	      
      // temporary hack 

      // notify neighbors 

      for (auto& item : higher_priority_neighbor_map) {
        vertex_locator neighbor = g.label_to_locator(item.first);
        mis_visitor new_visitor(neighbor, g.locator_to_label(vertex),
          std::get<1>(alg_data)[vertex]);
	vis_queue->queue_visitor(new_visitor);
      } // for	      

      for (auto& item : lower_priority_neighbor_map) {
        vertex_locator neighbor = g.label_to_locator(item.first);
        mis_visitor new_visitor(neighbor, g.locator_to_label(vertex),
          std::get<1>(alg_data)[vertex]);
	vis_queue->queue_visitor(new_visitor);
      } // for

      higher_priority_neighbor_map.clear();
      lower_priority_neighbor_map.clear();  
      return false; 	   

    } else { // undecided
      not_finished = 1; // TODO: may not need this	    
      bool vertex_has_max_prioriy = true;

      for (auto& item : higher_priority_neighbor_map) {
        //vertex_locator neighbor = g.label_to_locator(item.first);
        if (item.second > std::get<1>(alg_data)[vertex]) { // item.second neighbor's priority
	  vertex_has_max_prioriy = false;
          break;   	  
	}           	
      } // for

      if (vertex_has_max_prioriy) {
        std::get<1>(alg_data)[vertex] = maximal_independent_set::in;

        // notify neighbors
	 
        for (auto& item : higher_priority_neighbor_map) {
          vertex_locator neighbor = g.label_to_locator(item.first);
          mis_visitor new_visitor(neighbor, g.locator_to_label(vertex),
            std::get<1>(alg_data)[vertex]); 
    	  vis_queue->queue_visitor(new_visitor);
        } // for	      

        for (auto& item : lower_priority_neighbor_map) {
          vertex_locator neighbor = g.label_to_locator(item.first);
          mis_visitor new_visitor(neighbor, g.locator_to_label(vertex),
            std::get<1>(alg_data)[vertex]);
	  vis_queue->queue_visitor(new_visitor);
        } // for

        higher_priority_neighbor_map.clear();
        lower_priority_neighbor_map.clear();  
        return false; 	   

      } // if

    } // else	    
 
    return true; // controller invoking delegate vertices 
    //return false;	  
  }	  

  friend inline bool operator>(const mis_visitor& v1, const mis_visitor& v2) {
    return false;	  
  }

  friend inline bool operator<(const mis_visitor& v1, const mis_visitor& v2) {
    return false;	  
  }	  

  vertex_locator vertex;
  VertexID neighbor_ID;
  VertexPriority neighbor_priority;
};

template <typename TGraph, typename VertexID, typename VertexPriority, 
  typename VertexPriorityCollection, typename VertexIDPriorityMapCollection>
void maximal_independent_set(TGraph* graph, 
  VertexPriorityCollection& vertex_priority_list, 
  VertexIDPriorityMapCollection& higher_priority_neighbor_map,
  VertexIDPriorityMapCollection& lower_priority_neighbor_map) {

  typedef typename TGraph::vertex_locator vloc_type;
  typedef typename TGraph::vertex_iterator vitr_type;
  typedef typename TGraph::controller_iterator citr_type;
  typedef typename TGraph::edge_iterator eitr_type;

  int mpi_rank = havoqgt::comm_world().rank();
  if (mpi_rank == 0) {
    std::cout << "Maximal independent set" << std::endl;   
  }

  // visitor

  typedef mis_visitor<TGraph, VertexID, VertexPriority> visitor_type;
  auto alg_data = std::forward_as_tuple(
    graph, // 0
    vertex_priority_list, // 1, initialized 
    higher_priority_neighbor_map, // 2 
    lower_priority_neighbor_map // 3
    );
  auto vq = havoqgt::create_visitor_queue<visitor_type,
    havoqgt::detail::visitor_priority_queue>(graph, alg_data);

  do {
    not_finished = 0; // false	  
    vq.init_visitor_traversal();
    MPI_Barrier(MPI_COMM_WORLD); // TODO: remove?

    not_finished = havoqgt::mpi_all_reduce(not_finished, 
      std::greater<uint8_t>(), MPI_COMM_WORLD);

    if (mpi_rank == 0) { 
      std::cout << "Undecided vertices exist: " << 
       (not_finished > 0 ? "true" : "false") << std::endl;
    }

    //MPI_Barrier(MPI_COMM_WORLD);
  } while(not_finished); 

  //MPI_Barrier(MPI_COMM_WORLD);
}	

} // end namespace maximal_independent_set	
