#pragma once

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <independent_set/utility.hpp>

namespace maximal_independent_set {

template<typename Visitor>
class init_prio_queue {

public:

  init_prio_queue() {}

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
class init_prio_visitor {

public:

  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  typedef VertexID EdgeID; 

  init_prio_visitor()
    {}

  init_prio_visitor(vertex_locator _vertex) : 
    vertex(_vertex) {}

  init_prio_visitor(vertex_locator _vertex, VertexID _neighbor_ID, 
    VertexPriority _neighbor_priority) :
    vertex(_vertex), neighbor_ID(_neighbor_ID), 
    neighbor_priority(_neighbor_priority) {}

  ~init_prio_visitor() {}

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    	  
    //int mpi_rank = havoqgt::comm_world().rank();
    //if (mpi_rank == 0) { 
    //  std::cout << neighbor_priority << std::endl;	    
    //}

    // temporary hack	  
    if (vertex.is_delegate()) {
      std::get<3>(alg_data)[vertex] = maximal_independent_set::out;	    
      return false;
    }
    // temporary hack

    auto graph = std::get<2>(alg_data);
    auto& higher_priority_neighbor_map = std::get<4>(alg_data)[vertex];
    auto& lower_priority_neighbor_map = std::get<5>(alg_data)[vertex];

    auto vertex_ID = graph->locator_to_label(vertex);

    if (std::get<3>(alg_data)[vertex] < 1) { // vertex_priority_list
      // vertex priority is 0, i.e., OUT     
    
      auto vertex_degree = graph->degree(vertex);

      std::get<3>(alg_data)[vertex] = 
        maximal_independent_set::vertex_priority
 	  <VertexID, EdgeID, VertexPriority>(vertex_ID, vertex_degree);
    }

    if ((std::get<3>(alg_data)[vertex] > neighbor_priority) || 
      ((std::get<3>(alg_data)[vertex] == neighbor_priority) && 
      (vertex_ID > neighbor_ID)) ) {
      // add to lower priority neighbor map
      auto find_item = lower_priority_neighbor_map.find(neighbor_ID);
      if (find_item == lower_priority_neighbor_map.end()) {
        auto insert_status = lower_priority_neighbor_map.insert({neighbor_ID, 
          neighbor_priority});
	if (!insert_status.second) {
	  std::cerr << "Error: failed to add item to the map." << 
            std::endl;
	  return false;
	}       	
      } else {
//        std::cerr << "Error: unexpected item in the map." << std::endl;
        return false;     	
      }   	      
    } else {
      // add to higher priority neighbor map          
      auto find_item = higher_priority_neighbor_map.find(neighbor_ID);
      if (find_item == higher_priority_neighbor_map.end()) {
        auto insert_status = higher_priority_neighbor_map.insert({neighbor_ID, 
          neighbor_priority});
	if (!insert_status.second) {
	  std::cerr << "Error: failed to add item to the map." << 
            std::endl;
	  return false;
	}       	
      } else {
//        std::cerr << "Error: unexpected item in the map." << std::endl;
        return false;     	
      }   	      
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
    //    std::get<0>(alg_data) << std::endl;   	    
    //}
    
    // temporary hack	  
    if (vertex.is_delegate()) {
      std::get<3>(alg_data)[vertex] = maximal_independent_set::out;	    
      return true;	    
    }
    // temporary hack

    if (std::get<3>(alg_data)[vertex] < 1) { // vertex_priority_list
      // vertex priority is 0, i.e., OUT	    

      auto vertex_ID = g.locator_to_label(vertex);
      auto vertex_degree = g.degree(vertex);

      std::get<3>(alg_data)[vertex] = 
        maximal_independent_set::vertex_priority
 	  <VertexID, EdgeID, VertexPriority>(vertex_ID, vertex_degree); 
    }  
    
    for(eitr_type eitr = g.edges_begin(vertex);
      eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      init_prio_visitor new_visitor(neighbor, g.locator_to_label(vertex), 
        std::get<3>(alg_data)[vertex]); // vertex_priority_list 
      vis_queue->queue_visitor(new_visitor);    
    } // for	    

    return true;	  
  }	  

  friend inline bool operator>(const init_prio_visitor& v1, 
    const init_prio_visitor& v2) {
    return false;	  
  }

  friend inline bool operator<(const init_prio_visitor& v1, 
    const init_prio_visitor& v2) {
    return false;	  
  }	  

  vertex_locator vertex;
  VertexID neighbor_ID;
  VertexPriority neighbor_priority;
};

template <typename TGraph, typename VertexID, typename VertexPriority, 
  typename VertexPriorityCollection, typename VertexIDPriorityMapCollection>
void initialize_priority(TGraph* graph, 
  VertexPriorityCollection& vertex_priority_list, 
  VertexIDPriorityMapCollection& higher_priority_neighbor_map,
  VertexIDPriorityMapCollection& lower_priority_neighbor_map) {

  typedef typename TGraph::vertex_locator vloc_type;
  typedef typename TGraph::vertex_iterator vitr_type;
  typedef typename TGraph::controller_iterator citr_type;
  typedef typename TGraph::edge_iterator eitr_type;

  int mpi_rank = havoqgt::comm_world().rank();
  if (mpi_rank == 0) {
    std::cout << "Initialize vertex priority" << std::endl;   
  }

  const size_t average_degree = maximal_independent_set::average_degree;
  const double scaled_average_degree = 
    maximal_independent_set::scaled_average_degree; 

  // visitor

  typedef init_prio_visitor<TGraph, VertexID, VertexPriority> visitor_type;
  auto alg_data = std::forward_as_tuple(
    average_degree, // 0
    scaled_average_degree, // 1
    graph, // 2
    vertex_priority_list, // 3, initialized to 0 
    higher_priority_neighbor_map, // 4 
    lower_priority_neighbor_map // 5
    );
  auto vq = havoqgt::create_visitor_queue<visitor_type,
    havoqgt::detail::visitor_priority_queue>(graph, alg_data);
  vq.init_visitor_traversal();
  //MPI_Barrier(MPI_COMM_WORLD);

  // temporary hack
  //for(vitr_type vitr = graph->delegate_vertices_begin();
  //  vitr != graph->delegate_vertices_end(); ++vitr) {
  //  vloc_type vertex = *vitr;
  //  vertex_priority_list[vertex] = maximal_independent_set::out; 
  //} // for
  // temporary hack

  // verification 
  
  /*for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
    ++vitr) {
    vloc_type vertex = *vitr;
    auto s = higher_priority_neighbor_map[vertex].size() + 
      lower_priority_neighbor_map[vertex].size();    
    assert(s == graph->degree(vertex));

    //if (mpi_rank == 0) {
    if (static_cast<VertexID>(graph->locator_to_label(vertex)) == 
      static_cast<VertexID>(218128)) {
      std::cout << graph->locator_to_label(vertex) << " " << 
        s << " " << graph->degree(vertex) << std::endl;
    }
  } 

  for(vitr_type vitr = graph->delegate_vertices_begin();
    vitr != graph->delegate_vertices_end(); ++vitr) {
    vloc_type vertex = *vitr;
    auto s = higher_priority_neighbor_map[vertex].size() + 
      lower_priority_neighbor_map[vertex].size(); 
    assert(s == graph->degree(vertex));

    //if (mpi_rank == 0) {
    if (static_cast<VertexID>(graph->locator_to_label(vertex)) == 
      static_cast<VertexID>(218128)) {
      std::cout << graph->locator_to_label(vertex) << " " << 
        s << " " << graph->degree(vertex) << std::endl;
    }
  }*/

  //for (citr_type citr = graph->controller_begin(); 
  //  citr != graph->controller_end(); ++citr) {
  //  vloc_type vertex = *citr;
  //  ++local_vertex_count;
  //  local_edge_count+=graph->degree(vertex); 
  //}   

  // verification
 
}	

} // end namespace maximal_independent_set	
