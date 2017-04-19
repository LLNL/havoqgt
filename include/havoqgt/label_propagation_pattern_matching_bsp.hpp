#ifndef HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_BSP_HPP_INCLUDED
#define HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_BSP_HPP_INCLUDED

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

template<typename IntegralType>
class vertex_state {
public:
  vertex_state() :
  is_active(false),
  vertex_pattern_index(0),
  global_itr_count(0),
  local_itr_count(0) {
  }

  bool is_active;
  size_t vertex_pattern_index; // TODO: change type
  IntegralType global_itr_count; // TODO: remove for BSP
  IntegralType local_itr_count; // TODO: remove for BSP
  // map<pattern_vertex, itr_count>
  // TODO: cahnge size_t to pattern_vertex type  
  std::unordered_map<size_t, IntegralType> pattern_vertex_itr_count_map; 
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
    itr_count(0),
    do_update_vertex_pattern_id(false) {}

  lppm_visitor(vertex_locator _vertex, uint64_t _itr_count = 0, 
    bool _do_update_vertex_pattern_id = false) : 
    vertex(_vertex), 
    itr_count(_itr_count), 
    do_update_vertex_pattern_id(_do_update_vertex_pattern_id) {}

  lppm_visitor(vertex_locator _vertex, vertex_locator _parent, 
    VertexData _parent_vertex_data, size_t _parent_pattern_index, 
    uint64_t _itr_count, bool _do_update_vertex_pattern_id = false) :
    vertex(_vertex), 
    parent(_parent),
    parent_vertex_data(_parent_vertex_data),
    parent_pattern_index(_parent_pattern_index), 
    itr_count(_itr_count), 
    do_update_vertex_pattern_id(_do_update_vertex_pattern_id) {}

  ~lppm_visitor() {}

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    if (!std::get<4>(alg_data)[vertex]) { // vertex_active
      return false;
    }
 
    // TODO: pass pointer to the delegate graph to pre_visit
    // if the vertex is not in the global map after the first superstep, ignore it
    //if (std::get<8>(alg_data) > 0) { // supestep #
    //  auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
    //  if (find_vertex == std::get<6>(alg_data).end()) {
    //    return false; 
    //  }
    //}    

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data);
    // std::get<10>(alg_data) - g

    bool match_found = false;
    bool valid_parent_found = false;

    // does vertex_data match any entry in the query pattern
    for (size_t vertex_pattern_index = 0;
      vertex_pattern_index < pattern_graph.vertex_data.size();
      vertex_pattern_index++) {
      if (pattern_graph.vertex_data[vertex_pattern_index] == vertex_data) {
        match_found = true;   

        // verify if heard from a valid parent
        if ((itr_count == 1 || itr_count == 2) && match_found) {
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

    // initial case - return true to handle delegates
    if (std::get<4>(alg_data)[vertex] && itr_count == 0 && !match_found) {
      return true;  
    }
 
    // TODO: update veretx_pattern_id
    // need to write a new constructor to do update only

    return match_found; //true;
  }
  
  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    
    if (!std::get<4>(alg_data)[vertex]) { // TODO: we might not need to maintain this state at all
      return false;
    }

    // if the vertex is not in the global map after the first superstep 
    // (first LP iteration of first global step), ignore it
    // std::get<8>(alg_data) - superstep
    // std::get<9>(alg_data) - initstep
    // std::get<10>(alg_data) - g
    if (std::get<8>(alg_data) > 0 || !std::get<9>(alg_data)) {
      auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
      if (find_vertex == std::get<6>(alg_data).end()) {
        return false; 
      }
    } 

    auto vertex_data = std::get<0>(alg_data)[vertex];
    //auto& pattern = std::get<1>(alg_data);
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data);

    int mpi_rank = havoqgt_env()->world_comm().rank(); 

    // does vertex_data match any entry in the query pattern
    bool match_found = false;

    // TODO: do you want to compute this every time or store in the memory?
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
      //return false;
      return true; // TODO: ask Roger?
    } 

    //std::cout << g.locator_to_label(vertex) << std::endl; // Test

    if (itr_count == 0 && match_found) { // int the BSP implementaion, it_count is the message ID 
      // send to all the neighbours
      for(eitr_type eitr = g.edges_begin(vertex); 
        eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();

        // TODO: if the grpah is undirected, for BSP, you can have only one exchange, otherwise, two exchanges 
        // TODO: avoid creating visitors with the "default" constructor ?  
        for (auto vertex_pattern_index : vertex_pattern_indices) {
          // do this for all the pattern indices for this vertex
          //if (g.locator_to_label(vertex) == 28) // Test
	  //  std::cout << g.locator_to_label(vertex) << " sending to " << g.locator_to_label(neighbor) << std::endl; // Test	 
          lppm_visitor new_visitor(neighbor, vertex, vertex_data, vertex_pattern_index, 1); 
          vis_queue->queue_visitor(new_visitor); 
        } // for
      } // for
      return true;
    }

    // itr_count == 1 as a result of BCAST to out-edges
    // itr_count == 2 as a result of P2P to in-edges 

    else if ((itr_count == 1 || itr_count == 2) && match_found) {
    //if (itr_count == 1 && match_found) {  
      ///if (g.locator_to_label(vertex) == 28 || g.locator_to_label(vertex) == 89) // Test
      ///  std::cout << g.locator_to_label(vertex) << " receiving from " << g.locator_to_label(parent) << " " << itr_count << std::endl; // Test 
      // heard from a neighbour (parent) whose vertex_data matches an entry in the query pattern
      
      // verify if parent_vertex_data meet constrains for current vertex
      // what pattern label the vartex_data of this vertex corresponds to  
      // what are the valid parent labels for this vertex
      match_found = false; // TODO: Ask Roger about how to set the visitor property

      for (auto vertex_pattern_index : vertex_pattern_indices) { 

        // verify and decide if vertex should reply to the parent and 
        // update vertex_state_map accordingly  
        uint64_t next_itr_count = verify_and_update_vertex_state_map(g, vis_queue, alg_data, vertex_pattern_index);

        if (next_itr_count == 0) {
          //std ::cout << "."; // Test 
          continue;
        }

        // TODO: not very correct
        //else if (/*itr_count*/next_itr_count >= pattern.size()) { // reached max iterations
          // TODO: Output vertex_ID and pattern_index at the end of nth iteration, similar to updating vertex rank
          //std ::cout << itr_count; // Test          
          //match_found = true;
          //continue;

        //}

        else {
          //std ::cout << itr_count; // Test
          ///if (g.locator_to_label(vertex) == 28 || g.locator_to_label(vertex) == 89) // Test
          ///  std::cout << g.locator_to_label(vertex) << " | " << vertex_pattern_index << " sending to " << g.locator_to_label(parent) << " " << (itr_count + 1) << std::endl; // Test

          match_found = true;
          if (itr_count == 1) { // TODO: an undirected graph should not reply on the back edge
            //match_found = true;
            //if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
            //  return true;
            //} 
            //else {
            //  lppm_visitor new_visitor(parent, vertex, vertex_data, vertex_pattern_index, itr_count + 1/*next_itr_count*/);
            //  vis_queue->queue_visitor(new_visitor);
            //}
            //std::cout << " >> Superstep # " << std::get<8>(alg_data) << std::endl; // Test 
          }            

        } // else  
 
      } // for all vertex patterns 

      if (!match_found) { 
        //std::get<4>(alg_data)[vertex] = false;
        return false;
      } else {
        return true;
      }

    } else {
      // exceeded max iterations  
      return false; 
    }
 
    return true; //true;
  }

  friend inline bool operator>(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  friend inline bool operator<(const lppm_visitor& v1, const lppm_visitor& v2) {
    return false;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  uint64_t verify_and_update_vertex_state_map(Graph& g, 
    VisitorQueueHandle vis_queue, AlgData& alg_data, 
    size_t vertex_pattern_index) const {

    typedef vertex_state<uint8_t> VertexState;

    //auto& pattern = std::get<1>(alg_data);  
    //auto& pattern_indices = std::get<2>(alg_data);
    auto& pattern_graph = std::get<7>(alg_data); 
 
    bool match_found = false;

    // verify if parent_pattern_index is valid   
    for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
      e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
      if (pattern_graph.edges[e] == parent_pattern_index) {
        match_found = true; 
      }  
    }

    if (!match_found) {
      return 0;
    }

    //match_found = false; 
    
    // TODO: this only works for chains, better move to adjacency list representation 
/*    if (vertex_pattern_index == 0 && 
      !(parent_pattern_index == vertex_pattern_index + 1)) {
      //continue;
      return 0;
    } else if (vertex_pattern_index == pattern.size() - 1 && 
      !(parent_pattern_index == vertex_pattern_index - 1)) {
      //continue;
      return 0;
    } else if (!(parent_pattern_index == vertex_pattern_index - 1) &&
      !(parent_pattern_index == vertex_pattern_index + 1)) {
      //continue;
      return 0;
    }*/

    // vertex heard from a valid neighbour (possibly) 
    // create an entry for this vertex in the vertex_state_map or 
    // update, if exists already 
    auto find_vertex = std::get<6>(alg_data).find(g.locator_to_label(vertex));
    if (find_vertex == std::get<6>(alg_data).end()) {
      auto insert_status = std::get<6>(alg_data).insert({g.locator_to_label(vertex), VertexState()});
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

    // Important: a vertex could receive messages after it has met all the 
    // constraints. In that case, do not update vertex state but return 1, 
    // since the parent is a valid one. 
    if (find_vertex->second.is_active) {
      return 1;
    }

    // figure out what pattern indices are expected and add them to pattern_vertex_itr_count_map
    //if (itr_count == 1) {
    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      //match_found = false; 
      //for (size_t pattern_index = 0;  pattern_index < pattern_indices.size(); pattern_index++) {
      for (auto e = pattern_graph.vertices[vertex_pattern_index];
        e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {
              
        //match_found = false;

        auto pattern_index = pattern_graph.edges[e];   

        //if (vertex_pattern_index - 1 == pattern_index) {
        //  match_found = true;
        //} else if (vertex_pattern_index + 1 == pattern_index) {
        //  match_found = true;
        //}
        
        //if (match_found) {
          auto find_pattern_vertex =  find_vertex->second.pattern_vertex_itr_count_map.find(pattern_index);
          if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
            auto insert_status = find_vertex->second.pattern_vertex_itr_count_map.insert({pattern_index, 0});
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl;
              return 0;
            }
            //find_pattern_vertex = insert_status.first;
          } 
        //} 

      } // for
      
    } // if

    if (find_vertex->second.pattern_vertex_itr_count_map.size() < 1) {
      return 0;
    }

    auto find_pattern_vertex = find_vertex->second.pattern_vertex_itr_count_map.find(parent_pattern_index);  
    if (find_pattern_vertex == find_vertex->second.pattern_vertex_itr_count_map.end()) {
      //auto insert_status = find_vertex->second.pattern_vertex_itr_count_map.insert({parent_pattern_index, 0});
      //if(!insert_status.second) {
      //  std::cerr << "Error: failed to add an element to the map." << std::endl;
      //  return 0;
      //}
      //find_pattern_vertex = insert_status.first; 
      std::cerr << "Error: did not find the expected item in the map." << std::endl;
      return 0;
    }      
   
    // update itr_count of the pattern vertex 
    if (find_pattern_vertex->second < itr_count) {
      find_pattern_vertex->second = itr_count;
    }   

    // verify if heard from all the neighbours
    bool did_heard_from_all = true; 
    // figure out current iteration count for this vertex      
    uint64_t min_itr_count = find_vertex->second.pattern_vertex_itr_count_map.begin()->second;

    // TODO: verify if v.first's are valid
    ///if (find_vertex->first == g.locator_to_label(vertex) && g.locator_to_label(vertex) == 89 
    ///    || g.locator_to_label(vertex) == 28
    ///    || g.locator_to_label(vertex) == 55681294
    ///    || g.locator_to_label(vertex) == 55681295
    ///    || g.locator_to_label(vertex) == 55668937
    ///    || g.locator_to_label(vertex) == 55668938) // Test    
    ///  std::cout << " > " << g.locator_to_label(vertex) << " : "; // Test
 
    for (auto& v : find_vertex->second.pattern_vertex_itr_count_map) {

      ///if (find_vertex->first == g.locator_to_label(vertex) && g.locator_to_label(vertex) == 89 
      ///    || g.locator_to_label(vertex) == 28
      ///    || g.locator_to_label(vertex) == 55681294
      ///    || g.locator_to_label(vertex) == 55681295
      ///    || g.locator_to_label(vertex) == 55668937
      ///    || g.locator_to_label(vertex) == 55668938) // Test   
      ///  std::cout << "(" << v.first << ", " << v.second << ") "; // Test

      if (v.second == 0) {
        did_heard_from_all = false;
        break; // TODO: This will not work with min_itr_count !!!!!!!!!!!!!!!!!!!!!!!
      } 

      if (v.second < min_itr_count) { // TODO: should be min
        min_itr_count = v.second;  
      }
 
    } // for
    find_vertex->second.is_active = did_heard_from_all; // update is_active status for this vertex

    // TODO: probably do not need the itr_count fields  
        
    ///if (find_vertex->first == g.locator_to_label(vertex) && g.locator_to_label(vertex) == 89 
    ///    || g.locator_to_label(vertex) == 28
    ///    || g.locator_to_label(vertex) == 55681294 
    ///    || g.locator_to_label(vertex) == 55681295
    ///    || g.locator_to_label(vertex) == 55668937 
    ///    || g.locator_to_label(vertex) == 55668938) // Test
    ///  std::cout << std::endl; // Test

    // update current iteration count for this vertex
//    if (find_vertex->second.local_itr_count <= min_itr_count) {  
//      find_vertex->second.local_itr_count = min_itr_count + 1;    
//    } 

    // update vertex_iteration
//    if (find_vertex->second.local_itr_count > std::get<5>(alg_data)[vertex]) {
//      std::get<5>(alg_data)[vertex] = find_vertex->second.local_itr_count;    
//    }

    // Test
    ///if (find_vertex->first == g.locator_to_label(vertex) && g.locator_to_label(vertex) == 89 || g.locator_to_label(vertex) == 28) {// Test
    ///  std::cout << " > " << g.locator_to_label(vertex) << " : current local_itr_count " << find_vertex->second.local_itr_count << std::endl; // Test   

      ///for (auto e = pattern_graph.vertices[vertex_pattern_index]; e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
      ///  auto v_nbr = pattern_graph.edges[e];
      ///  std::cout << v_nbr << ", " ; 
      ///}
      ///std::cout << std::endl;

    ///} // Test

    // Importnt: reset vertex state for the next iteration
    if (find_vertex->second.is_active) {
      for (auto& v : find_vertex->second.pattern_vertex_itr_count_map) {
        v.second = 0;
      }   
    }  

    return 1; //find_vertex->second.local_itr_count; 
 
  } 

  vertex_locator vertex;
  vertex_locator parent;
  VertexData parent_vertex_data; 
  size_t parent_pattern_index;
  uint64_t itr_count;
  bool do_update_vertex_pattern_id;
};

template <typename TGraph, typename AlgData, typename VertexStateMap, 
  typename PatternGraph, typename VertexActive, typename VertexIteration>
void verify_and_update_vertex_state_map(TGraph* g, AlgData& alg_data, 
  VertexStateMap& vertex_state_map, PatternGraph& pattern_graph, 
  VertexActive& vertex_active, 
  VertexIteration& vertex_iteration, uint64_t superstep, bool& global_not_finished) {

  int mpi_rank = havoqgt_env()->world_comm().rank();
  
  uint64_t max_superstep = pattern_graph.vertex_data.size() - 1; // Test
  
  if (vertex_state_map.size() > 0 && superstep == max_superstep) {
  // std::cout << mpi_rank << " Superstep #" << superstep << " vertex state map size " 
  //   << vertex_state_map.size() << std::endl; // Test  
  }
 
  auto vertex_temp = vertex_state_map.begin()->first;
  std::vector<decltype(vertex_temp)> vertex_remove_from_map_list(0); // TODO: use Vetex type instead of uint64_t

  for (auto& v : vertex_state_map) {
    auto v_locator = g->label_to_locator(v.first);
    auto vertex_data = std::get<0>(alg_data)[v_locator]; // Test 
    
    //for (auto& p : v.second.pattern_vertex_itr_count_map) {
    //  if (p.second < 1) {
    //    vertex_active[v_locator] = false; 
    //    vertex_remove_from_map_list.push_back(v.first);
    //    break;
    //  } else {
    //    p.second = 0; // reset for next iteration   
    //  }   
    //}
   
    if (vertex_state_map.size() > 0 && superstep == max_superstep) { 
    //  std::cout << mpi_rank << " Superstep #" << superstep << " | vertex " << v.first <<  " data " << vertex_data  
    //    << " is in the map" << std::endl; // Test   
    }
 
    if (!v.second.is_active) {
      vertex_remove_from_map_list.push_back(v.first);
      vertex_active[v_locator] = false;  
    } else {
      v.second.is_active = false; // reset for next iteration
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

  if (vertex_state_map.size() > 0 && superstep == max_superstep) {
  // std::cout << mpi_rank << " Superstep #" << superstep << " | " 
  //    << vertex_remove_from_map_list.size() 
  //    << " vertices were removed from the vertex state map, new map size " 
  //    << vertex_state_map.size()
  //    << std::endl; // Test     
  }

  vertex_active.all_max_reduce(); 
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
  //std::cout << "label_propagation_pattern_matching_bsp.hpp" << std::endl;

  int mpi_rank = havoqgt_env()->world_comm().rank();
  uint64_t superstep_var = 0;
  uint64_t& superstep_ref = superstep_var; 

  typedef lppm_visitor<TGraph, VertexData> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank,
    vertex_active, vertex_iteration, vertex_state_map, pattern_graph, superstep_var, initstep, g);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

  // beiginning of BSP execution
  // TODO: change this to use a local termination detection at the end of a seperstep
  bool not_finished = false;
 
  for (uint64_t superstep = 0; superstep < pattern_graph.diameter; superstep++) {
    superstep_ref = superstep;
    if (mpi_rank == 0) { 
      std::cout << "Superstep #" << superstep << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    double time_start = MPI_Wtime();
    vq.init_visitor_traversal_new(); 
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {    
      std::cout << "Superstep #" << superstep <<  " Synchronizing ... " << std::endl;
    }

    vertex_active.all_max_reduce();
    MPI_Barrier(MPI_COMM_WORLD);
 
    verify_and_update_vertex_state_map(g, alg_data, vertex_state_map, pattern_graph, 
      vertex_active, vertex_iteration, superstep, global_not_finished);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    double time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Superstep #" << superstep <<  " Time : " << time_end - time_start << std::endl;
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
    /*if (superstep == 1) {
      std::string active_vertices_result_filename = "/p/lscratchf/havoqgtu/reza2_tmp/rmat_tmp_result/0/all_ranks_active_vertices_lp/active_vertices_" + std::to_string(mpi_rank);
      std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);
      for (auto& v : vertex_state_map) {
        auto v_locator = g->label_to_locator(v.first);
        if (v_locator.is_delegate() && (g->master(v_locator) == mpi_rank)) {
          active_vertices_result_file << mpi_rank << ", c, " << vertex_metadata[v_locator] << ", " 
            << v.first << ", "
            << v.second.vertex_pattern_index << "\n";
        } else if (v_locator.is_delegate() && (g->master(v_locator) != mpi_rank)) {
          active_vertices_result_file << mpi_rank << ", d, " << vertex_metadata[v_locator] << ", "
            << v.first << ", "
            << v.second.vertex_pattern_index << "\n";
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
  
  //vertex_rank.all_reduce();
  //vertex_iteration.all_max_reduce();
  //MPI_Barrier(MPI_COMM_WORLD);
}  

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_LABEL_PROPAGATION_PATTERN_MATCHING_BSP_HPP_INCLUDED 
