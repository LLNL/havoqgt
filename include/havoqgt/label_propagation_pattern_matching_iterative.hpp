#ifndef HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED
#define HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename IntegralType>
class vertex_state {
public:
  vertex_state() :
  is_active(false),
  vertex_pattern_index(0) {}

  bool is_active;
  size_t vertex_pattern_index; // TODO: change type
  std::unordered_map<size_t, IntegralType> pattern_vertex_itr_count_map; // TODO: not itr_count anymore, more like true / false 
};

template<typename Visitor>
class lppm_queue {

public:
  lppm_queue() {}

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
    return data.size();;
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

// label propagation pattern matching visitor class
template<typename Graph, typename VertexData>
class lppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  lppm_visitor() : 
    msg_type(0) {}

  lppm_visitor(vertex_locator _vertex, uint8_t _msg_type = 0) : 
    vertex(_vertex), 
    msg_type(_msg_type) {}

  lppm_visitor(vertex_locator _vertex, size_t _parent_pattern_index,
    uint8_t _msg_type) :
    vertex(_vertex),
    parent_pattern_index(_parent_pattern_index),
    msg_type(_msg_type) {}

  ~lppm_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep 
    auto g = std::get<10>(alg_data); 
 
    bool match_found = false;
    bool valid_parent_found = false;

    size_t vertex_pattern_index = 0;

    int mpi_rank = havoqgt_env()->world_comm().rank();     

      if (vertex.is_delegate() && g->master(vertex) != mpi_rank && msg_type == 1) { // a delegate but not the controller
       // the vertex_state is only maintained on the controller

        // TODO: avoid figuring out vertex_pattern_index everytime
        // does vertex_data match any entry in the query pattern
        for (vertex_pattern_index = 0;
          vertex_pattern_index < pattern_graph.vertex_data.size();
          vertex_pattern_index++) {
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
            match_found = true;   

            // verify if heard from a valid parent
            if (msg_type == 1 && match_found) {
              //match_found = false;
              for (auto e = pattern_graph.vertices[vertex_pattern_index];
                e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
                if (pattern_graph.edges[e] == parent_pattern_index) {
                  //match_found = true;
                  valid_parent_found = true; 
                  break; 
                }
              } // for
 
              if (!valid_parent_found) {
                return false; 
              }    

            } // if            
 
          } // if

          if (valid_parent_found) {
            break;
          }
        } // for 

        // initial case - return true to handle delegates // TODO: try this one too
        //if (std::get<4>(alg_data)[vertex] && msg_type == 0 && !match_found) {
        //  return true;  
        //}

        if (!match_found) {
          std::get<4>(alg_data)[vertex] = false; 
          //return true; // send to the controller
          return false;
        } else {
          return true; // send to the controller
        } 
 
      } 

      // for local vertex and controller only
 
      // first LP superstep of the first iteration 
      if (std::get<8>(alg_data) == 0  && std::get<9>(alg_data)) {
        for (vertex_pattern_index = 0; 
          vertex_pattern_index < pattern_graph.vertex_data.size(); 
          vertex_pattern_index++) { 
          if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
             match_found = true;
             break; // Important 
          }
        }

        if (!match_found) { // TODO: controller return true?
          std::get<4>(alg_data)[vertex] = false;
          if (vertex.is_delegate() && g->master(vertex) == mpi_rank) { // controller
            //return true; // to invalidate the delegates
            return false;
          }           
          return false; 
        }
      }
     
      // if the vertex is not in the global map after the first LP 
      // superstep of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) { 
        auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false;
        } else {
          vertex_pattern_index = find_vertex->second.vertex_pattern_index;
        }  
      }
 
      if (msg_type == 1) {
        verify_and_update_vertex_state_map(alg_data, vertex_pattern_index);
      }
   
      return false;     
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {    
    if (!std::get<4>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank(); 

    // Important : skip this verification for the delegates as the 
    // vertex_state is only maintained on the contrller
    if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
      // if the vertex is not in the global map after the first LP superstep 
      // of the first iteration, ignore it
      if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) {
        auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
        if (find_vertex == std::get<6>(alg_data).end()) {
          return false; 
        }
      } 
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    // std::get<4>(alg_data) - vertex_active
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep
    // std::get<10>(alg_data) - g

    // does vertex_data match an entry in the query pattern
    bool match_found = false;

    // TODO: do you want to compute this every time or store in the memory? Overhead is not noticable though.
    std::vector<size_t> vertex_pattern_indices(0); // a vertex label could be a match for multiple pattern labels
    for (size_t vertex_pattern_index = 0; 
      vertex_pattern_index < pattern_graph.vertex_data.size(); 
      vertex_pattern_index++) { 
      if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
        vertex_pattern_indices.push_back(vertex_pattern_index);
        // TODO: compare with the entry in pattern_indices to detect loop or 
        // use token passing
        match_found = true;
        //break; 
      }       
    }

    if (!match_found) {
      std::get<4>(alg_data)[vertex] = false;
      return false;
      //return true; // TODO: ask Roger?
    } 

    if (msg_type == 0 && match_found) {
      // send to all the neighbours
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();

        // TODO: only handling undirected grpahs

        for (auto vertex_pattern_index : vertex_pattern_indices) {
          // do this for all the pattern indices for this vertex
          lppm_visitor new_visitor(neighbor, vertex_pattern_index, 1);
          vis_queue->queue_visitor(new_visitor);
        } // for
      } // for
      return true;
    } else if (msg_type == 1 && match_found) {        
      // must go all the way to the controller 
      //return true; // false?
      return false; // TODO: ask Roger?
    } else {
      return false; 
    } 
    return true;
  }

  friend inline bool operator>(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  template<typename AlgData>
  uint64_t verify_and_update_vertex_state_map(AlgData& alg_data, 
    size_t vertex_pattern_index) const {

    typedef vertex_state<uint8_t> VertexState; // TODO: use Vertex type

    //auto& pattern = std::get<1>(alg_data);  
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data); 
    auto g = std::get<10>(alg_data);
 
    bool match_found = false;

    // verify if parent_pattern_index is valid   
    for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
      e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
      if (pattern_graph.edges[e] == parent_pattern_index) {
        match_found = true;
        break; 
      }  
    }

    if (!match_found) {
      return 0;
    }

    // vertex heard from a valid neighbour 
    // create an entry for this vertex in the vertex_state_map or 
    // update, if exists already 
    auto find_vertex = std::get<6>(alg_data).find(g->locator_to_label(vertex));
    if (find_vertex == std::get<6>(alg_data).end()) {
      auto insert_status = std::get<6>(alg_data).insert({g->locator_to_label(vertex), VertexState()});
      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map." << std::endl;
        return 0;
      }     	
      find_vertex = insert_status.first;
      find_vertex->second.vertex_pattern_index = vertex_pattern_index; // ID of the vertex in the pattern_graph 
    }

    if (std::get<6>(alg_data).size() < 1) {
      return 0;
    }

    // figure out what pattern indices are expected and add them to pattern_vertex_itr_count_map
    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      for (auto e = pattern_graph.vertices[vertex_pattern_index];
        e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
              
        auto pattern_index = pattern_graph.edges[e];   
        
          auto find_pattern_vertex =  find_vertex->second.pattern_vertex_itr_count_map.find(pattern_index);
          if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
            auto insert_status = find_vertex->second.pattern_vertex_itr_count_map.insert({pattern_index, 0});
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl;
              return 0;
            }
          } 

      } // for
      
    } // if

    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      return 0;
    }

    // set status of parent_pattern_index to 1
    auto find_pattern_vertex = find_vertex->second.pattern_vertex_itr_count_map.find(parent_pattern_index);  
    if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
      std::cerr << "Error: did not find the expected item in the map." << std::endl;
      return 0;
    }      
   
    // update status of the pattern vertex 
    if (find_pattern_vertex->second < 1) {
      find_pattern_vertex->second = 1;
    }  
 
    return 1; 
  } 

  vertex_locator vertex;
  size_t parent_pattern_index; // TODO: pass type as template argument
  uint8_t msg_type; // 0 - init, 1 - alive
};

template <typename TGraph, typename AlgData, typename VertexStateMap, 
  typename PatternGraph, typename VertexActive, typename VertexIteration>
void verify_and_update_vertex_state_map(TGraph* g, AlgData& alg_data, 
  VertexStateMap& vertex_state_map, PatternGraph& pattern_graph, 
  VertexActive& vertex_active, 
  VertexIteration& vertex_iteration, uint64_t superstep, bool& global_not_finished) {

  typedef typename TGraph::vertex_locator vertex_locator;

  //auto vertex_temp = vertex_state_map.begin()->first;
  //std::vector<decltype(vertex_temp)> vertex_remove_from_map_list(0); // hack
  std::vector<uint64_t> vertex_remove_from_map_list; // TODO: use Vertex type

  for (auto& v : vertex_state_map) { // TODO: use C++11 approach to remove item from map
    auto v_locator = g->label_to_locator(v.first);
    
    for (auto& p : v.second.pattern_vertex_itr_count_map) {
      if (p.second < 1) {
        vertex_remove_from_map_list.push_back(v.first);
        vertex_active[v_locator] = false; 
        break;
      } else {
        p.second = 0; // reset for next iteration   
      }   
    }   
  } // for

  if (vertex_remove_from_map_list.size() > 0) {
    global_not_finished = true;
  }

  for (auto v : vertex_remove_from_map_list) {
    if (vertex_state_map.erase(v) < 1) {
      std::cerr << "Error: failed to remove an element from the map." 
        << std::endl;  
    }    
  }

  vertex_active.all_min_reduce(); 
  MPI_Barrier(MPI_COMM_WORLD); 
}   

template <typename TGraph, typename VertexMetaData, typename VertexData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename VertexActive, 
  typename VertexIteration, typename VertexStateMap, typename PatternGraph>
void label_propagation_pattern_matching_bsp(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, VertexRank& vertex_rank,
  VertexActive& vertex_active, VertexIteration& vertex_iteration, VertexStateMap& vertex_state_map, 
  PatternGraph& pattern_graph, bool initstep, bool& global_not_finished, size_t global_itr_count, 
  std::ofstream& superstep_result_file, std::ofstream& active_vertices_count_result_file) {

  int mpi_rank = havoqgt_env()->world_comm().rank();
  uint64_t superstep_var = 0;
  uint64_t& superstep_ref = superstep_var; 

  typedef lppm_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank,
    vertex_active, vertex_iteration, vertex_state_map, pattern_graph, superstep_var, initstep, g);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // beiginning of BSP execution
  // TODO: change for loop to use a local termination detection at the end of a seperstep
  //bool not_finished = false;
 
  for (uint64_t superstep = 0; superstep < pattern_graph.diameter; superstep++) {
    superstep_ref = superstep;
    if (mpi_rank == 0) { 
      //std::cout << "Superstep #" << superstep << std::endl;
      std::cout << "Label Propagation | Superstep #" << superstep;
    }

    //MPI_Barrier(MPI_COMM_WORLD); 
    double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); 
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {    
      //std::cout << "Superstep #" << superstep <<  " Synchronizing ... " << std::endl;
      std::cout <<  " | Synchronizing ...";
    }

    //vertex_active.all_min_reduce(); // do not need this here
    ///MPI_Barrier(MPI_COMM_WORLD);
 
    verify_and_update_vertex_state_map(g, alg_data, vertex_state_map, pattern_graph, 
      vertex_active, vertex_iteration, superstep, global_not_finished);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      //std::cout << "Superstep #" << superstep <<  " Time " << time_end - time_start << std::endl;
      std::cout << " | Time : " << time_end - time_start << std::endl;
    }

    // result
    if (mpi_rank == 0) { 
      superstep_result_file << global_itr_count << ", LP, "
        << superstep << ", "
        << time_end - time_start << "\n"; 
    }

    // Important : This may slow things down -only for presenting results
    uint64_t active_vertices_count = 0; 
    for (auto& v : vertex_state_map) {
      auto v_locator = g->label_to_locator(v.first);
      if (v_locator.is_delegate() && (g->master(v_locator) == mpi_rank)) {
        active_vertices_count++;  
      } else if (!v_locator.is_delegate()) {
        active_vertices_count++;
      }
    }
 
    active_vertices_count_result_file << global_itr_count << ", LP, "
      << superstep << ", " 
      << active_vertices_count << "\n";

    // Test
/*    if (superstep == 1) {
      std::string active_vertices_result_filename = "/p/lscratchf/havoqgtu/reza2_tmp/rmat_tmp_result/0/all_ranks_active_vertices_lp/active_vertices_" + std::to_string(mpi_rank);
      std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);
      for (auto& v : vertex_state_map) {
        auto v_locator = g->label_to_locator(v.first);
        if (v_locator.is_delegate() && (g->master(v_locator) == mpi_rank)) {
          active_vertices_result_file << mpi_rank << ", c, " << vertex_metadata[v_locator] << ", "
            << v.first << ", " 
            << v.second.vertex_pattern_index << "\n";
        } else if (v_locator.is_delegate() && (g->master(v_locator) != mpi_rank)) {
          //active_vertices_result_file << mpi_rank << ", d, "
          // << v.first << ", "
          // << v.second.vertex_pattern_index << "\n"; 
        } else if (!v_locator.is_delegate()) {
          active_vertices_result_file << mpi_rank << ", l, " << vertex_metadata[v_locator] << ", "
            << v.first << ", "
            << v.second.vertex_pattern_index << "\n";
        }
      }		   
    }*/
    // Test

    // TODO: global reduction on global_not_finished before next iteration

  } // for 
  // end of BSP execution  
}  

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_ITERATIVE_HPP_INCLUDED 
