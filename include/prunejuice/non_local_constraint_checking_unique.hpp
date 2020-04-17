#pragma once

#include <deque>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

//using namespace havoqgt;

///namespace havoqgt { ///namespace mpi {

namespace prunejuice {

static bool enable_vertex_token_source = true;
static uint64_t tp_visitor_count = 0;

template<typename Visitor>
class tppm_queue {

public:
  tppm_queue() {}

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

// token passing pattern matching visitor class
template<typename Graph, typename BitSet>
class tppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  tppm_visitor() : 
  itr_count(0), 
  do_pass_token(false),
  is_init_step(true),
  ack_success(false),
  source_index_pattern_indices(0), 
  parent_pattern_index(0) {}

  tppm_visitor(vertex_locator _vertex) :  
    vertex(_vertex), 
    itr_count(0),
    do_pass_token(false), 
    is_init_step(true),
    ack_success(false),
    source_index_pattern_indices(0), 
    parent_pattern_index(0) {}
   
  tppm_visitor(vertex_locator _vertex, 
    vertex_locator _parent, 
    vertex_locator _target_vertex, 
    size_t _itr_count, 
    size_t _max_itr_count, 
    size_t _source_index_pattern_indices, 
    size_t _parent_pattern_index, 
    bool _expect_target_vertex = true, 
    bool _do_pass_token = true, 
    bool _is_init_step = false, 
    bool _ack_success = false) : 
    vertex(_vertex),
    parent(_parent),
    target_vertex(_target_vertex), 
    itr_count(_itr_count), 
    max_itr_count(_max_itr_count), 
    expect_target_vertex(_expect_target_vertex), 
    do_pass_token(_do_pass_token), 
    is_init_step(_is_init_step),
    ack_success(_ack_success), 
    source_index_pattern_indices(_source_index_pattern_indices), 
    parent_pattern_index(_parent_pattern_index) {}  

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    tp_visitor_count++;
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);
    //auto& pattern_graph = std::get<4>(alg_data);
    // std::get<5>(alg_data); // vertex_state_map
    auto g = std::get<11>(alg_data); // graph
    // std::get<12>(alg_data); // vertex_token_source_set   
    // std::get<13>(alg_data); // vertex_active
    // std::get<14>(alg_data); // template_vertices
    // std::get<15>(alg_data); // vertex_active_edges_map
    // std::get<16>(alg_data); // pattern_selected_vertices
    // std::get<17>(alg_data); // pattern_selected_edges
    // std::get<18>(alg_data); // pattern_mark_join_vertex
    // std::get<19>(alg_data); // pattern_ignore_join_vertex
    // std::get<20>(alg_data); // pattern_join_vertex      

    ///int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_rank = havoqgt::comm_world().rank();

    if (ack_success) {
      // TODO: if vertex == target_vertex
      // it is in the token_source map
      // !bitset.none
      // if delegate, return true otherwise false
      return true; // Important : must return true to handle delegates 
    } 
   
    if(enable_vertex_token_source) {    
    // verify if this vertex have already forwarded a token from the originating vertex
    if (!is_init_step && max_itr_count > itr_count) {
      auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g->locator_to_label(target_vertex));
      if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
        return false;
      }
    }
    } //enable_vertex_token_source
/* 
    // TODO: use vertex_active
    // TODO: only start token passing from the controller 
    auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
    if (find_vertex == std::get<5>(alg_data).end()) {    
      return false;
    }
*/    
    if (!do_pass_token && is_init_step && itr_count == 0) { // probably it never gets here
      // TODO: this is probbaly wrong
//      if (vertex.is_delegate() && (g->master(vertex) != mpi_rank)) {
//        return false;
//      }

//      if (!(find_vertex->second.vertex_pattern_index == pattern_indices[0] && vertex_data == pattern[0])) {
      if (std::get<16>(alg_data)) {
         if ((vertex_data != pattern[0]) && (vertex_data != pattern[pattern.size() - 1])) {   
           return false;  
         } else if (vertex_data == pattern[pattern.size() - 1]) { 
           // TODO: add vertex to token_source_map and return false
           return true;  
         } else if (vertex_data == pattern[0]) {
           return true; 
         }   
      } else {
        if(vertex_data == pattern[0]) {
          return true;
        } else {
          return false;
        }
      }     
         
    } else if (!is_init_step) { // relay the token
       // target_vertex cannot relay the token      
       if (do_pass_token && (max_itr_count > itr_count) && 
         (g->locator_to_label(vertex) == g->locator_to_label(target_vertex))) {
         return false;
       }

       auto new_itr_count = itr_count + 1;
       auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index       
       auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;

       // TODO: read vertex_pattern_index from template_vertices
       // vertex_data == pattern[next_pattern_index] and vertex_pattern_index == pattern_indices[next_pattern_index] must hold  
       // otherwise return false
       if (vertex_data != pattern[next_pattern_index]) {
         return false;
       }

       bool match_found = false; 
       
       BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices
  
       if (vertex_template_vertices.none()) {
         return false;  
       }

       for (size_t i = 0; i < vertex_template_vertices.size(); i++) { // TODO: vertex_template_vertices.test(pattern_indices[next_pattern_index]) 
         if (vertex_template_vertices.test(i)) {
           if (i == pattern_indices[next_pattern_index]) {
             match_found = true;
             vertex_pattern_index = i; 
             break;  
           }    
         }   
       }
 
       if (!match_found) {
         return false;
       }    
   
/*//--       if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
         if (vertex_data == pattern[next_pattern_index]) {
           vertex_pattern_index = pattern_indices[next_pattern_index];
         } else {
           return false;
         } 
//--       } else {
//--         auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
//--         if (find_vertex == std::get<5>(alg_data).end()) {
//--           return false;
//--         }
//--         vertex_pattern_index = find_vertex->second.vertex_pattern_index;
         
//--       }*/

       // verify if received from a valid parent
       if (vertex_data == pattern[next_pattern_index] &&
         vertex_pattern_index == pattern_indices[next_pattern_index]) {
         if (vertex_data == pattern[next_pattern_index] && 
           parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 

//+           if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
//+             return true;
//+           }  

           if (do_pass_token && (max_itr_count > itr_count)) {

             // join vertex
             // std::get<18>(alg_data); // pattern_mark_join_vertex
             // std::get<19>(alg_data); // pattern_ignore_join_vertex
             // std::get<20>(alg_data); // pattern_join_vertex 
             // if vertex is not a delegate, delegate vertices are not added to the vertex_state_map
             if (!(vertex.is_delegate() && (g->master(vertex) != mpi_rank))) { 
               if (std::get<18>(alg_data) && (vertex_pattern_index == std::get<20>(alg_data))) { // mark join vertex
                 // find vertex in the vertex_state_map and mark it
                 auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
                 if (find_vertex == std::get<5>(alg_data).end()) {
                   std::cerr << "Error: did not find the expected item in the map." << std::endl;
                   return false;
                 } else {
                   if (!find_vertex->second.is_active) {
                     find_vertex->second.is_active = true;
                   }
                 }      
               } else if (std::get<19>(alg_data) && (vertex_pattern_index == std::get<20>(alg_data))) { // ignore join vertex
                 // find vertex in the vertex_state_map and ignore it 
                 auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
                 if (find_vertex == std::get<5>(alg_data).end()) {
                   std::cerr << "Error: did not find the expected item in the map." << std::endl;
                   return false;
                 } else {
                   if (find_vertex->second.is_active) {
                     return false;
                   }   
                 } 
               }
             }  

             if (enable_vertex_token_source) {  
             // OK to forwarded a token from a source, now update vertex_token_source_set
             auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g->locator_to_label(target_vertex));
             if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
               auto insert_status = std::get<12>(alg_data)[vertex].insert(g->locator_to_label(target_vertex));
               if(!insert_status.second) {
                 std::cerr << "Error: failed to add an element to the set." << std::endl;
                 return false;
               }
               // std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(target_vertex)
	       //  << " to the vertex set" << std::endl; // Test     
	     } else {
               std::cerr << "Error: unexpected item in the set." << std::endl;
               return false;
             }
             } // enable_vertex_token_source 
 
           } 

           //TODO: max_itr_count == itr_count // handle in pre_visit
           
           return true; // delegate forwarding to the controller 
         } else {
           return false;
         }    
       } else {
         return false;
       }  
    } else {
      return false;  
    }
//++    return true;
    return false; // Test
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    //tp_visitor_count++;
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    ///int mpi_rank = havoqgt_env()->world_comm().rank();
    int mpi_rank = havoqgt::comm_world().rank();

    if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) { // Test ++
      std::cerr << "Error: controller forwarding visitor to delegates." << std::endl;
      return false;
    } // Test ++


    if (ack_success) {
      auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
      if (find_token_source == std::get<6>(alg_data).end()) {
        std::cerr << "Error: did not find the expected item in the map." << std::endl;
        //return true; // false ?           
        return false;
      }
      find_token_source->second = 1; //true;   
//    }
 
      std::get<9>(alg_data) = 1; // true; // pattern_found

      //std::cerr << "ACK_SUCCESS " << g.locator_to_label(vertex) << " " << itr_count << std::endl;	// Test
         
//++      return true; // Important : must return true to handle delegates 
      return false;
    } 

    // if vertex is a delegate  
    if (!is_init_step && vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
      if (enable_vertex_token_source) { 
      // verify if this vertex have already forwarded a token from the originating vertex
      if (!is_init_step && max_itr_count > itr_count) {
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
          return false; 
        } 		
      }
      } // enable_vertex_token_source
    }
/*
    // TODO: verify if this vertex is alive
    auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
    if (find_vertex == std::get<5>(alg_data).end()) {
      return false;
    } //else { // Test
      //std::cout << find_vertex->first << " " << find_vertex->second.vertex_pattern_index 
      //<< " " << std::get<0>(alg_data)[vertex] << std::endl;
    //} // Test
*/    
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);   
    //auto& pattern_graph = std::get<4>(alg_data);
    // std::get<5>(alg_data); // vertex_state_map
    auto pattern_cycle_length = std::get<7>(alg_data);   
    auto pattern_valid_cycle = std::get<8>(alg_data);
    //auto& pattern_found = std::get<9>(alg_data);
    //auto& edge_metadata = std::get<10>(alg_data); 
    // <11> // graph
    // std::get<12>(alg_data); // vertex_token_source_set
    // std::get<13>(alg_data); // vertex_active 
    // std::get<14>(alg_data); // template_vertices
    // std::get<15>(alg_data); // vertex_active_edges_map
    // std::get<16>(alg_data); // pattern_selected_vertices
    // std::get<17>(alg_data); // pattern_selected_edges
    
    // std::get<18>(alg_data); // pattern_mark_join_vertex
    // std::get<19>(alg_data); // pattern_ignore_join_vertex
    // std::get<20>(alg_data); // pattern_join_vertex 

    if (!do_pass_token && is_init_step && itr_count == 0) {
      // create visitors only for the source vertices
       
      // Test 
//      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
//        vertex_locator neighbor = eitr.target();
//        if (g.locator_to_label(neighbor) != 46) { // TODO: put all the valid vertices form label propagation in a queue and pop here and create visitors like this 
//          continue;    
//        }  
        // 65 66 67 66 67 69 68 66 67 66 67 : 0 1 2 3 4 5 6 3 4 1 2 
//        tppm_visitor new_visitor(neighbor, g.label_to_locator(46), 0, 7, 1, pattern_indices[1], true, true, true); 
//        vis_queue->queue_visitor(new_visitor);              
//      }
      // Test
       
      // TODO: this is probbaly wrong
      // for delegate vertices only start token passing from the controller   
//      if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
//        return false;
//      }

      //if (!(find_vertex->second.vertex_pattern_index == pattern_indices[0] && vertex_data == pattern[0])) {
      if (std::get<16>(alg_data)) { // pattern_selected_vertices // unique vertex metadata
        if ((vertex_data != pattern[0]) && (vertex_data != pattern[pattern.size() - 1])) { 
          return false;  
        } else if (vertex_data == pattern[pattern.size() - 1]) {
          // add (destination) vertex to token source map
          // vertex will be deleted if it does not hear from at least one of
          // the vertices in its token_source_set 
          auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
          if (find_token_source == std::get<6>(alg_data).end()) {
            auto insert_status = std::get<6>(alg_data).insert({g.locator_to_label(vertex), false});
            if(!insert_status.second) {
              std::cerr << "Error: failed to add an element to the map." << std::endl;
              return false;
            }
            //std::cout << "Instrting " << vertex_data << " to token_source_map" << std::endl; // Test
          } 
          return false; // Important 
        }  
      } else {
        if (vertex_data != pattern[0]) {
           return false;     
        } 
      }   

      // pattern_selected_vertices and  vertex_token_source_set 
      //if (std::get<16>(alg_data) && std::get<12>(alg_data)[vertex].empty()) { 
      //  return false;          
      //}  

      BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices
 
      if (vertex_template_vertices.none() || !vertex_template_vertices.test(pattern_indices[0])) {
        return false;  
      }

      // TODO: disabling 13 Dec 2017, do we want to keep this
      // nonunique vertex metadata
      // initiate tokens only from vertices with multiple template vertex matches 
      if (!pattern_valid_cycle && !std::get<16>(alg_data))  { // path checking // not pattern_selected_vertices
        if (!(vertex_template_vertices.test(pattern_indices[0]) 
          && vertex_template_vertices.test(pattern_indices[pattern_indices.size() -1]) ) ) {
          return false;
        } 
        
        //std::cout << vertex_template_vertices << std::endl;      
      } 
      // nonunique vertex metadata    
 
      // token passing constraints go here
//      tppm_visitor new_visitor(vertex, g.label_to_locator(28), 0, 2, 0, pattern_indices[0], true, true, true); 
//      vis_queue->queue_visitor(new_visitor);
//      return true;
      
      //std::cout << "Found source vertex " << g.locator_to_label(vertex) << " vertex_data " << vertex_data << std::endl; // Test
      
      // Important: only the token_source_map on the controller contains the source vertex (wrong I think)
//      if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
      // local, controller and delegates are added to the token_source_map 
      if (!std::get<16>(alg_data)) { // not pattern_selected_vertices     
 
        auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
        if (find_token_source == std::get<6>(alg_data).end()) {
          auto insert_status = std::get<6>(alg_data).insert({g.locator_to_label(vertex), false});
          if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the map." << std::endl;
            return false;   
          } 
          //std::cout << "Instrting " << vertex_data << " to token_source_map" << std::endl; // Test    
        }

      }  	
//      }

      // initiate token passing from the source vertex
//--      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
//--        vertex_locator neighbour = eitr.target();

      for (auto& item : std::get<15>(alg_data)[vertex]) { // vertex_active_edges_map
        vertex_locator neighbour = g.label_to_locator(item.first);
 
        //std::get<10>(alg_data)[eitr] = 240; // Test
        //int neighbour_edge_data = std::get<10>(alg_data)[eitr]; // Test

        //std::cout << "found source vertex " << g.locator_to_label(vertex) 
        //<< " vertex_data " << vertex_data 
        //<< " neighbour " << g.locator_to_label(neighbour) 
        //<< " edge data " << neighbour_edge_data << std::endl; // Test
        
        // token passing constraints go here
//        tppm_visitor new_visitor(neighbour, g.label_to_locator(4902294), 0, 2, 0, pattern_indices[0], true, true, false);
//        tppm_visitor new_visitor(neighbour, g.label_to_locator(932746), 0, (pattern_indices.size() - 2), 0, pattern_indices[0], true, true, false); // Test              

        // loop detection - path back to the source vertex is valid
        //tppm_visitor new_visitor(neighbour, vertex, 0, 2, 0, pattern_indices[0], true, true, false);

//        if (std::get<16>(alg_data)) { // pattern_selected_vertices
          // vertex_token_source_set // the set was populate in the previous iteration 
          // std::cout << "Token source: " << g.locator_to_label(vertex) << " " << std::get<12>(alg_data)[vertex].size() << std::endl; // Test
//          for (auto v : std::get<12>(alg_data)[vertex]) {
//            tppm_visitor new_visitor(neighbour, g.label_to_locator(v), 0, pattern_cycle_length, 0, pattern_indices[0], pattern_valid_cycle, true, false);
//            vis_queue->queue_visitor(new_visitor);
//          }        
//        } else {  
//          tppm_visitor new_visitor(neighbour, vertex, 0, (pattern_indices.size() - 2), 0, pattern_indices[0], true, true, false);
          
          // pattern_selected_edges          
          if (std::get<17>(alg_data) && !item.second) {
            continue;
          }
          
          tppm_visitor new_visitor(neighbour, vertex, vertex, 0, pattern_cycle_length, 0, pattern_indices[0], pattern_valid_cycle, true, false);
        
          // loop detection - path back to the source vertex is invalid
          //tppm_visitor new_visitor(neighbour, vertex, 0, 2, 0, pattern_indices[0], false, true, false);
          vis_queue->queue_visitor(new_visitor);
//        } 
      }
//++      return true;
      return false; // controller not forwarding visitor to delegates
 
    } /*else if ((find_vertex->second.vertex_pattern_index == pattern_indices[0]) && itr_count == 0 && is_init_step) {
      // initiate token passing from the source vertex
      //std::cout << "found source vertex " << g.locator_to_label(vertex) << " vertex_data " << vertex_data << std::endl; // Test
      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        tppm_visitor new_visitor(neighbor, target_vertex, itr_count, max_itr_count, 
          source_index_pattern_indices, parent_pattern_index, true, true, false);
        vis_queue->queue_visitor(new_visitor);
      }  
      return true;
 
    }*/ else if (!is_init_step) { // else if      

      // target_vertex cannot relay the token
      if (do_pass_token && (max_itr_count > itr_count) &&
        (g.locator_to_label(vertex) == g.locator_to_label(target_vertex))) {
        return false;
      }     

      bool do_forward_token = false;
      auto new_itr_count = itr_count + 1;
      auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index  
      auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      //auto vertex_pattern_index = pattern_indices[source_index_pattern_indices + new_itr_count]; // TODO: read from the map
      
      // TODO: read vertex_pattern_index from template_vertices
      // vertex_data == pattern[next_pattern_index] and vertex_pattern_index == pattern_indices[next_pattern_index] must hold  
      // otherwise return false
      if (vertex_data != pattern[next_pattern_index]) {
        return false;
      } 

      bool match_found = false;

      BitSet vertex_template_vertices(std::get<14>(alg_data)[vertex]); // template_vertices

      if (vertex_template_vertices.none()) {
        return false;
      }

      for (size_t i = 0; i < vertex_template_vertices.size(); i++) { // TODO: vertex_template_vertices.test(pattern_indices[next_pattern_index])
        if (vertex_template_vertices.test(i)) {
          if (i == pattern_indices[next_pattern_index]) {
            match_found = true;
            vertex_pattern_index = i;
            break; 
          }
        }
      }

      if (!match_found) {
        return false;
      }

/*//--      if (vertex.is_delegate() && g.master(vertex) != mpi_rank) { // delegate but not the controller
        if (vertex_data == pattern[next_pattern_index]) {
          vertex_pattern_index = pattern_indices[next_pattern_index];  
        } else {
          return false;  
        }  
//--      } else {
//--        auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
//--        if (find_vertex == std::get<5>(alg_data).end()) {
//--          return false;
//--        }
//--        vertex_pattern_index = find_vertex->second.vertex_pattern_index;
//--      }*/      
      
      // TODO: verify next_pattern_index < pattern_indices.size() before anythin else

      //int mpi_rank = havoqgt_env()->world_comm().rank(); // Test
      //if (mpi_rank == 51) { 
      //  std::cout << g.locator_to_label(vertex) << " vertex_pattern_index "
      //  << vertex_pattern_index << " received from parent_pattern_index " 
      //  << parent_pattern_index << " itr_count " << itr_count 
      //  << " max_itr_count " << max_itr_count << " next_pattern_index " << next_pattern_index
      //  << std::endl; // Test	 
      //} // Test  

      if (max_itr_count > itr_count) {

        // are vertex_data and vertex_pattern_index valid
//        if (vertex_data == pattern[pattern_indices[next_pattern_index]] && 
	// TODO: Important: for token passing only        
        if (vertex_data == pattern[next_pattern_index] &&
          vertex_pattern_index == pattern_indices[next_pattern_index]) {
          // verify if received from a valid parent
//          for (auto e = pattern_graph.vertices[vertex_pattern_index]; 
//            e < pattern_graph.vertices[vertex_pattern_index + 1]; e++) {  
//            if (pattern_graph.edges[e] == parent_pattern_index) {
//              do_forward_token = true; 
//              break; 
//            }  
//          } // for      
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 
            do_forward_token = true;

            // join vertex
            // std::get<18>(alg_data); // pattern_mark_join_vertex
            // std::get<19>(alg_data); // pattern_ignore_join_vertex
            // std::get<20>(alg_data); // pattern_join_vertex 
            // if vertex is not a delegate, delegate vertices are not added to the vertex_state_map
            if (!(vertex.is_delegate() && (g.master(vertex) != mpi_rank))) { 
              if (std::get<18>(alg_data) && (vertex_pattern_index == std::get<20>(alg_data))) { // mark join vertex
                // find vertex in the vertex_state_map and mark it
                auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
                if (find_vertex == std::get<5>(alg_data).end()) {
                  std::cerr << "Error: did not find the expected item in the map." << std::endl;
                  return false;
                } else {
                  if (!find_vertex->second.is_active) {
                    find_vertex->second.is_active = true;
                  }
                }      
              } else if (std::get<19>(alg_data) && (vertex_pattern_index == std::get<20>(alg_data))) { // ignore join vertex
                // find vertex in the vertex_state_map and ignore it 
                auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
                if (find_vertex == std::get<5>(alg_data).end()) {
                  std::cerr << "Error: did not find the expected item in the map." << std::endl;
                  return false;
                } else {
                  if (find_vertex->second.is_active) {
                    return false;
                  }   
                } 
              }
            } 
            // join vertex
              
          }   
        } // if    
 
      } else if (max_itr_count == itr_count) {
        // are vertex_data and vertex_pattern_index valid
        //bool match_found = false;
        match_found = false;
//        if (vertex_data == pattern[pattern_indices[next_pattern_index]] &&
        // TODO: Important: for token passing only        
        if (vertex_data == pattern[next_pattern_index] &&        
          vertex_pattern_index == pattern_indices[next_pattern_index]) { 
          // verify if received from a valid parent
//          for (auto e = pattern_graph.vertices[vertex_pattern_index];
//            e < pattern_graph.vertices[vertex_pattern_index + 1]; e      return false;) {
//            if (pattern_graph.edges[e] == parent_pattern_index) {
//              match_found = true;
//              break;
//            }
//          } // for
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) {
            match_found = true; 
          }      
        } // if

        // path checking
        if (match_found && !expect_target_vertex) {
          if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex)) {
            //std::cout << "found invalid cycle - vertex " <<  " | parent_pattern_index "
            //  <<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex)
            //  <<  " vertex_pattern_index " << vertex_pattern_index << " itr "
            //  << itr_count << std::endl; // Test  
            return false; //true; 
                      
          } else {
            //std::cout << "found valid path - vertex " <<  " | parent_pattern_index "
            //  <<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex)
            //  <<  " vertex_pattern_index " << vertex_pattern_index << " itr "
            //  << itr_count << std::endl; // Test              
               
            if (std::get<16>(alg_data)) { // pattern_selected_vertices
              //std::cout << "found valid path " << vertex_data << " " 
              //  << std::get<12>(alg_data)[vertex].size() << std::endl; // Test
                
              // find the target_vertex in the vertex_token_source_set
              auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
              if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
                return false;    
              } else {     

                auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
                if (find_token_source == std::get<6>(alg_data).end()) {
                  std::cerr << "Error: did not find the expected item in the map." << std::endl;
                  //return true; // false ?           
                  return false;
                }

                find_token_source->second = 1; //true;   
                std::get<9>(alg_data) = 1; // true; // pattern_found
	
              }

              return false;  
            } else {  
              // valid path found, send ack to the token source so it could update its state
              tppm_visitor new_visitor(target_vertex, vertex, target_vertex, itr_count, max_itr_count,
                source_index_pattern_indices, 0, expect_target_vertex, false, false, true);
                vis_queue->queue_visitor(new_visitor);  
              return false; //true;            
            }   
          }  
        } 
 
        // cycle checking
        // is this the target vertex
        else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
          && match_found && expect_target_vertex) {
          // found cycle 
          //std::cout << "found valid cycle - vertex " << g.locator_to_label(vertex) <<  " | parent_pattern_index " 
          //<<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex) 
          //<<  " vertex_pattern_index " << vertex_pattern_index << " itr " 
          //<< itr_count << std::endl; // Test

          //return false; // TODO: true ?	

          // TODO: this is probbaly wrong
          // Important: only the token_source_map on the controller contains the source vertex (wrong I think)
//          if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
//            std::get<9>(alg_data) = 1; // true; // pattern_found
//            return true;		        
//          }  
//          if (!(vertex.is_delegate() && g.master(vertex) != mpi_rank)) {
          auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
          if (find_token_source == std::get<6>(alg_data).end()) {
            std::cerr << "Error: did not find the expected item in the map." << std::endl;
            //return true; // false ?           
            return false;
          }
          find_token_source->second = 1; //true;   
//          }
 
          std::get<9>(alg_data) = 1; // true; // pattern_found	

//++          return true; // Important : must return true to handle delegates
          
          // mark the edge received the successful token on     	  
	  //std::get<15>(alg_data)[vertex])  // vertex_active_edges_map  	 	
	  auto find_edge = std::get<15>(alg_data)[vertex].find(g.locator_to_label(parent));
	  if (find_edge == std::get<15>(alg_data)[vertex].end()) {
            std::cerr << "Error: did not find the expected item in the map." << std::endl;
            return false;
          } else {
            find_edge->second = 1; 
          }  		

	  return false;
        } 

        //else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
        //  && match_found && !expect_target_vertex) {
          // TODO: TBA 
          //std::cout << "found invalid cycle - vertex " <<  " | parent_pattern_index "
          //<<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex)
          //<<  " vertex_pattern_index " << vertex_pattern_index << " itr "
          //<< itr_count << std::endl; // Test  
        //  return true; 
        //} 

        else {
          // reached max iteration but did not find the target vertex or a loop
          //std::cout << "At " << g.locator_to_label(vertex) <<  ", did not find target " 
          //<< g.locator_to_label(target_vertex) <<  " after " << itr_count 
          //<< " iterations" <<std::endl; // Test
          return false;//true;//false;  
        }   
      } else {
        std::cerr << "Error: wrong code branch." << std::endl;    
        return false;
      }   

      if (!do_forward_token) {
        return false;
      } //else {
        //return false; // Test 
      //} 

      // all good, forward along the token
      
      // if vertex is a delegate
/*++      if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
        // forwarded a token from a source, now update vertex_token_source_set 
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded == std::get<12>(alg_data)[vertex].end()) {
          auto insert_status = std::get<12>(alg_data)[vertex].insert(g.locator_to_label(target_vertex));		
	  if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the set." << std::endl;
            return false;
          }
          //std::cout << g.locator_to_label(vertex) << " adding " << g.locator_to_label(target_vertex) 
          //  << " to the vertex set" << std::endl; // Test 
        } else {
          std::cerr << "Error: unexpected item in the set." << std::endl;
	  return false;
        }        
      } // if vertex is a delegate 
*/

      //if (vertex_pattern_index == 2) // Test 
        //std::cout << g.locator_to_label(vertex) << " vertex_pattern_index " 
        //<< vertex_pattern_index << " " << new_itr_count << " forwarding ... " 
        //<< g.locator_to_label(target_vertex) << std::endl; // Test
      //size_t max_nbr_count = 0; // Test 
//      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
//        vertex_locator neighbor = eitr.target();

      for (auto& item : std::get<15>(alg_data)[vertex]) { // vertex_active_edges_map
        vertex_locator neighbour = g.label_to_locator(item.first);

        // do not forward the token to the parent the vertex received it from
        if (g.locator_to_label(neighbour) == g.locator_to_label(parent)) {
          continue;
        }
      
        tppm_visitor new_visitor(neighbour, vertex, target_vertex, new_itr_count, max_itr_count, 
          source_index_pattern_indices, vertex_pattern_index, expect_target_vertex); 
        // vertex_pattern_index = parent_pattern_index for the neighbours 
        vis_queue->queue_visitor(new_visitor);
        // Test
        //if (max_nbr_count == 20) {
        //  break;     
        //} else {
        //  max_nbr_count++;
        //}  
        // Test    
      }
	     
//++      return true;
      return false; // controller not forwarding visitor to delegates

      // else if
    } else {
      return false;
    }
    return false;		 
  }

  friend inline bool operator>(const tppm_visitor& v1, const tppm_visitor& v2) {
    return false;
    /*if (v1.itr_count > v2.itr_count) {
      return true;
    } else if (v1.itr_count < v2.itr_count) {
      return false; 
    }
    if (v1.vertex == v2.vertex) {
      return false;
    }
    return !(v1.vertex < v2.vertex);*/
    
    /*if (v1.itr_count <= v2.itr_count) {
      return true;
    } else {
      return false;
    }*/ 
  }

  //friend inline bool operator<(const tppm_visitor& v1, const tppm_visitor& v2) {
    //return false;
    //if (v1.itr_count < v2.itr_count) {
    //  return true;
    //} else if (v1.itr_count > v2.itr_count) {
    //  return false;
    //}
    //if (v1.vertex == v2.vertex) {
    //  return false;
    //}
    //return !(v1.vertex < v2.vertex);    
  //}

  vertex_locator vertex;
  vertex_locator parent;  
  vertex_locator target_vertex; // for a cycle, this is also the originating vertex
  size_t itr_count; // TODO: change type
  size_t max_itr_count; // equal to diameter - 1 of the pattern as itr_count is initialized to 0 // TODO: change type
  bool expect_target_vertex;
  bool do_pass_token;
  bool is_init_step;
  bool ack_success; 
  size_t source_index_pattern_indices; // index of the token source in the pattern_indices container 
  size_t parent_pattern_index; // TODO: change to the same type as in the pattern_graph
};

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename PatternGraph, 
  typename VertexStateMap, typename TokenSourceMap, typename EdgeMetaData, 
  typename VertexSetCollection, typename VertexActive, typename TemplateVertex, 
  typename VertexUint8MapCollection, typename BitSet>
void token_passing_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, 
  VertexRank& vertex_rank, PatternGraph& pattern_graph, VertexStateMap& vertex_state_map,
  TokenSourceMap& token_source_map, size_t pattern_cycle_length, bool pattern_valid_cycle, 
  std::vector<uint8_t>::reference pattern_found, EdgeMetaData& edge_metadata, 
  VertexSetCollection& vertex_token_source_set, VertexActive& vertex_active, 
  TemplateVertex& template_vertices, VertexUint8MapCollection& vertex_active_edges_map, 
  bool pattern_selected_vertices, bool pattern_selected_edges,
  bool pattern_mark_join_vertex, bool pattern_ignore_join_vertex, size_t pattern_join_vertex, 
  uint64_t& message_count) { // TODO: bool& pattern_found does not work, why?
  //std::cout << "token_passing_pattern_matching_new.hpp" << std::endl;
  
  tp_visitor_count = 0;

  int mpi_rank = havoqgt::comm_world().rank();
  int mpi_size = havoqgt::comm_world().size();

  if (mpi_rank == 0) {
    std::cout << "Token Passing ... " << std::endl;
    //std::cout << "Token Passing [" << pl << "] ... " << std::endl;
  } 

  typedef tppm_visitor<TGraph, BitSet> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank, 
    pattern_graph, vertex_state_map, token_source_map, pattern_cycle_length, pattern_valid_cycle, pattern_found, 
    edge_metadata, g, vertex_token_source_set, vertex_active, template_vertices, vertex_active_edges_map, 
    pattern_selected_vertices, pattern_selected_edges, 
    pattern_mark_join_vertex, pattern_ignore_join_vertex, pattern_join_vertex);
  auto vq = havoqgt::create_visitor_queue<visitor_type, /*havoqgt::detail::visitor_priority_queue*/tppm_queue>(g, alg_data);
  ///vq.init_visitor_traversal_new();
  //vq.init_visitor_traversal_new_alt();
  vq.init_visitor_traversal();
  MPI_Barrier(MPI_COMM_WORLD);

  message_count = tp_visitor_count;
}

} ///} //end namespace havoqgt::mpi
