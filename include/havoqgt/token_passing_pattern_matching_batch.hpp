#ifndef HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_BATCH_HPP_INCLUDED
#define HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_BATCH_HPP_INCLUDED

#include <deque>
#include <unordered_set>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

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
template<typename Graph>
class tppm_visitor {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;
  tppm_visitor() : 
  itr_count(0), 
  do_pass_token(false), 
  is_init_step(true), 
  source_index_pattern_indices(0), 
  parent_pattern_index(0) {}

  tppm_visitor(vertex_locator _vertex) :  
    vertex(_vertex), 
    itr_count(0),
    do_pass_token(false), 
    is_init_step(true), 
    source_index_pattern_indices(0), 
    parent_pattern_index(0) {}
   
  tppm_visitor(vertex_locator _vertex, vertex_locator _target_vertex, 
    size_t _itr_count, 
    size_t _max_itr_count, 
    size_t _source_index_pattern_indices, 
    size_t _parent_pattern_index, 
    bool _expect_target_vertex = true, 
    bool _do_pass_token = true, 
    bool _is_init_step = false) : 
    vertex(_vertex),
    target_vertex(_target_vertex), 
    itr_count(_itr_count), 
    max_itr_count(_max_itr_count), 
    expect_target_vertex(_expect_target_vertex), 
    do_pass_token(_do_pass_token), 
    is_init_step(_is_init_step), 
    source_index_pattern_indices(_source_index_pattern_indices), 
    parent_pattern_index(_parent_pattern_index) {}  

  template<typename AlgData>
  bool pre_visit(AlgData& alg_data) const {
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);
    //auto& pattern_graph = std::get<4>(alg_data);
    auto g = std::get<11>(alg_data); // graph
    // std::get<12>(alg_data); // vertex_token_source_set   
    // std::get<13>(alg_data); // vertex_active

    int mpi_rank = havoqgt_env()->world_comm().rank();
    
    // verify if this vertex has already forwarded a copy of the token
    if (!is_init_step && max_itr_count > itr_count) {
      auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g->locator_to_label(target_vertex));
      if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
        return false;
      }
    }

    if (!do_pass_token && is_init_step && itr_count == 0) { // probably it never gets here
      if(vertex_data == pattern[0]) { 
        return true;
      } else {
        return false; 
      } 
    } else if (!is_init_step) { // relay token
       auto new_itr_count = itr_count + 1;
       auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index       
       auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      
       if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
         if (vertex_data == pattern[next_pattern_index]) {
           vertex_pattern_index = pattern_indices[next_pattern_index];
         } else {
           return false;
         } 
       } else {
         auto find_vertex = std::get<5>(alg_data).find(g->locator_to_label(vertex));
         if (find_vertex == std::get<5>(alg_data).end()) {
           return false;
         }
         vertex_pattern_index = find_vertex->second.vertex_pattern_index;
       }

       if (vertex_data == pattern[next_pattern_index] &&
         vertex_pattern_index == pattern_indices[next_pattern_index]) {
         if (vertex_data == pattern[next_pattern_index] && 
           parent_pattern_index == pattern_indices[next_pattern_index - 1]) {

           if (vertex.is_delegate() && g->master(vertex) != mpi_rank) { // delegate but not the controller
             // TODO: try to move code from visit to here
             return true;
           }  

           if (max_itr_count > itr_count) {
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
           } 

           //TODO: if (max_itr_count == itr_count) 
  
           return true; 
         } else {
           return false;
         }    
       } else {
         return false;
       }  
    } else {
      return false;  
    }
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue,
    AlgData& alg_data) const {
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if(!std::get<13>(alg_data)[vertex]) {
      return false;
    }

    int mpi_rank = havoqgt_env()->world_comm().rank();

    // if vertex is a delegate  
    if (!is_init_step && vertex.is_delegate() && (g.master(vertex) != mpi_rank)) { 
      // verify if this vertex has already forwarded a copy of the token
      if (!is_init_step && max_itr_count > itr_count) {
        auto find_token_source_forwarded = std::get<12>(alg_data)[vertex].find(g.locator_to_label(target_vertex));
        if (find_token_source_forwarded != std::get<12>(alg_data)[vertex].end()) {
          return false; 
        } 		
      }
    }
    
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto& pattern = std::get<1>(alg_data);
    auto& pattern_indices = std::get<2>(alg_data);   
    //auto& pattern_graph = std::get<4>(alg_data);

    auto pattern_cycle_length = std::get<7>(alg_data);   
    auto pattern_valid_cycle = std::get<8>(alg_data);
    //auto& pattern_found = std::get<9>(alg_data);
    //auto& edge_metadata = std::get<10>(alg_data); 
    // std::get<12>(alg_data); // graph
    // std::get<12>(alg_data); // vertex_token_source_set
    // std::get<13>(alg_data); // vertex_active 

    if (!do_pass_token && is_init_step && itr_count == 0) {

      // if vertex is in the token_source_map, initiate a token from this vertex
      auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
      if (find_token_source == std::get<6>(alg_data).end()) {
        return false;
      }

      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbour = eitr.target();
        tppm_visitor new_visitor(neighbour, vertex, 0, pattern_cycle_length, 0, 
          pattern_indices[0], pattern_valid_cycle, true, false);
        vis_queue->queue_visitor(new_visitor);
      }
      return true;
 
    } else if (!is_init_step) { 

      bool do_forward_token = false;
      auto new_itr_count = itr_count + 1;
      auto next_pattern_index = source_index_pattern_indices + new_itr_count; // expected next pattern_index  
      auto vertex_pattern_index = 0; //find_vertex->second.vertex_pattern_index;
      // TODO: verify next_pattern_index < pattern_indices.size() before anything else

      if (vertex.is_delegate() && g.master(vertex) != mpi_rank) { // delegate but not the controller
        if (vertex_data == pattern[next_pattern_index]) {
          vertex_pattern_index = pattern_indices[next_pattern_index];  
        } else {
          return false;  
        }  
      } else {
        auto find_vertex = std::get<5>(alg_data).find(g.locator_to_label(vertex));
        if (find_vertex == std::get<5>(alg_data).end()) {
          return false;
        }
        vertex_pattern_index = find_vertex->second.vertex_pattern_index;
      }           

      if (max_itr_count > itr_count) {
        // are vertex_data and vertex_pattern_index valid
        if (vertex_data == pattern[next_pattern_index] &&
          vertex_pattern_index == pattern_indices[next_pattern_index]) {
          // verify if received from a valid parent
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) { 
            do_forward_token = true; 
          }   
        } 
 
      } else if (max_itr_count == itr_count) {
        // are vertex_data and vertex_pattern_index valid
        bool match_found = false;
        if (vertex_data == pattern[next_pattern_index] &&        
          vertex_pattern_index == pattern_indices[next_pattern_index]) { 
          // verify if received from a valid parent
          if (parent_pattern_index == pattern_indices[next_pattern_index - 1]) {
            match_found = true; 
          }      
        }

        // is this the token source vertex
        if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
          && match_found && expect_target_vertex) {
          // found valid cycle 
          auto find_token_source = std::get<6>(alg_data).find(g.locator_to_label(vertex)); // token_source_map
          if (find_token_source == std::get<6>(alg_data).end()) {
            std::cerr << "Error: did not find the expected item in the map." << std::endl;    
            return false;
          }
          find_token_source->second = 1; //true;   
          std::get<9>(alg_data) = 1; // true; // pattern_found	
          return true; // Important : must return true to handle delegates
        } else if (g.locator_to_label(vertex) == g.locator_to_label(target_vertex) 
          && match_found && !expect_target_vertex) {
          // TODO: TBA 
          //std::cout << "found invalid cycle - vertex " <<  " | parent_pattern_index "
          //<<  parent_pattern_index <<  " | " << g.locator_to_label(target_vertex)
          //<<  " vertex_pattern_index " << vertex_pattern_index << " itr "
          //<< itr_count << std::endl; // Test  
          return true; 
        } else {
          // reached max iteration but did not find the target vertex or a cycle
          //std::cout << "At " << g.locator_to_label(vertex) <<  ", did not find target " 
          //<< g.locator_to_label(target_vertex) <<  " after " << itr_count 
          //<< " iterations" <<std::endl; // Test
          return true; //false;  
        }   
      } else {
        std::cerr << "Error: wrong code branch." << std::endl;    
        return false;
      }   

      if (!do_forward_token) {
        return false;
      } 

      // all good, forward along the token
      
      // if vertex is a delegate
      if (vertex.is_delegate() && (g.master(vertex) != mpi_rank)) {
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

      for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
        vertex_locator neighbor = eitr.target();
        tppm_visitor new_visitor(neighbor, target_vertex, new_itr_count, max_itr_count, 
          source_index_pattern_indices, vertex_pattern_index, expect_target_vertex); 
        // vertex_pattern_index = parent_pattern_index for the neighbours 
        vis_queue->queue_visitor(new_visitor);
      }
	     
      return true;

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
  vertex_locator target_vertex; // for a cycle, this is also the token source vertex
  size_t itr_count; // TODO: change type // initialized to 0
  size_t max_itr_count; // equal to diameter - 1 of the pattern as itr_count is initialized to 0 // TODO: change type
  bool expect_target_vertex;
  bool do_pass_token;
  bool is_init_step;
  size_t source_index_pattern_indices; // index of the token source in the pattern_indices container 
  size_t parent_pattern_index; // TODO: change to the same type as in the pattern_graph
};

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices, typename VertexRank, typename PatternGraph, 
  typename VertexStateMap, typename TokenSourceMap, typename EdgeMetaData, 
  typename VertexSetCollection, typename VertexActive>
void token_passing_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, 
  VertexRank& vertex_rank, PatternGraph& pattern_graph, VertexStateMap& vertex_state_map,
  TokenSourceMap& token_source_map, size_t pattern_cycle_length, bool pattern_valid_cycle, 
  std::vector<uint8_t>::reference pattern_found, EdgeMetaData& edge_metadata, 
  VertexSetCollection& vertex_token_source_set, VertexActive& vertex_active, 
  bool& global_not_finished) { 
  // TODO: bool& pattern_found does not work, why?

  typedef typename TGraph::vertex_iterator vertex_iterator;
  typedef typename TGraph::vertex_locator vertex_locator;

  int mpi_rank = havoqgt_env()->world_comm().rank();
  int mpi_size = havoqgt_env()->world_comm().size();

  // setup token sources
  std::unordered_set<uint64_t> token_source_set; // TODO: use Vertex type
  
  for (vertex_iterator vitr = g->vertices_begin(); 
    vitr != g->vertices_end(); ++vitr) {  
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex] && vertex_metadata[vertex] == pattern[0] ) {
      auto find_token_source = token_source_set.find(g->locator_to_label(vertex));
      if (find_token_source == token_source_set.end()) {
        auto insert_status = token_source_set.insert(g->locator_to_label(vertex));
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the set." << std::endl;
          //return false;
        }
      }      
    }  
  }

  for(vertex_iterator vitr = g->delegate_vertices_begin();
    vitr != g->delegate_vertices_end(); ++vitr) {
    vertex_locator vertex = *vitr;
    if (vertex_active[vertex] && vertex_metadata[vertex] == pattern[0] ) {
      auto find_token_source = token_source_set.find(g->locator_to_label(vertex));
      if (find_token_source == token_source_set.end()) {
        auto insert_status = token_source_set.insert(g->locator_to_label(vertex));
        if(!insert_status.second) {
          std::cerr << "Error: failed to add an element to the set." << std::endl;
          //return false;
        }
      } 
    }
  } 

  uint64_t global_token_source_set_size =
  havoqgt::mpi::mpi_all_reduce(token_source_set.size(), std::plus<size_t>(),
    MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    std::cout << "Token Passing | Global Token Source Count : "
      << global_token_source_set_size << std::endl;
  }

  // batch parameters
  uint64_t max_batch_size = 1;
  uint64_t batch_size = 0;
  uint64_t batch_count = 0;
  uint64_t max_ranks_per_itr = mpi_size / 4; //(mpi_size / 24) * 2;
  uint64_t max_itr_count = mpi_size;

  // pass around tokens in batches
  for (uint64_t itr_count = 0 ; itr_count < max_itr_count;
    itr_count+=max_ranks_per_itr) {

    double time_start = MPI_Wtime();
    batch_size = 0;
    token_source_map.clear(); // Important

    assert(iter_count + max_ranks_per_itr <= mpi_size);

    // setup token source map
    if (mpi_rank >= itr_count && mpi_rank < (itr_count + max_ranks_per_itr)) {
      for (auto v = token_source_set.cbegin(); v != token_source_set.cend();) {
        auto find_token_source = token_source_map.find(*v);
        if (find_token_source == token_source_map.end()) {
          auto insert_status = token_source_map.insert({*v, false});
          if(!insert_status.second) {
            std::cerr << "Error: failed to add an element to the map." << std::endl;
          }
          v = token_source_set.erase(v); // C++11
          //batch_size++;
        } else {
          std::cerr << "Error: unexpected item in the set." << std::endl;
        }
      } // for  
    } // if
 
    MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
  
    uint64_t global_batch_token_source_map_size = 
      havoqgt::mpi::mpi_all_reduce(token_source_map.size(), 
      std::plus<size_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // pass around tokens
    typedef tppm_visitor<TGraph> visitor_type;
    auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank, 
      pattern_graph, vertex_state_map, token_source_map, pattern_cycle_length, pattern_valid_cycle, pattern_found, 
      edge_metadata, g, vertex_token_source_set, vertex_active);
    auto vq = create_visitor_queue<visitor_type, /*havoqgt::detail::visitor_priority_quue*/tppm_queue>(g, alg_data);
    vq.init_visitor_traversal_new();
    //vq.init_visitor_traversal_new_alt();
    MPI_Barrier(MPI_COMM_WORLD);

    double time_end = MPI_Wtime();
    if(mpi_rank == 0) {
    std::cout << "Token Passing | Batch #" << batch_count
      << " | MPI Ranks " << itr_count << " - "
      << (itr_count + max_ranks_per_itr - 1)
      << " | Global Batch Token Source Count "
      << global_batch_token_source_map_size
      <<  " | Time : " << time_end - time_start << std::endl;
    }

    // remove the invalid (token source) vertices from the vertex_state_map
    uint64_t remove_count = 0;  

    for (auto& s : token_source_map) {
      if (!s.second) {
        vertex_active[g->label_to_locator(s.first)] = false;
        if (!global_not_finished) { 
          global_not_finished = true;
        }       
      }
    }

    vertex_active.all_min_reduce();
    MPI_Barrier(MPI_COMM_WORLD);
 
    // remove from vertex_state_map
    for (auto& s : token_source_map) {
      auto v_locator = g->label_to_locator(s.first);
      if (!vertex_active[v_locator]) {
        auto find_vertex = vertex_state_map.find(s.first);
 
        if (find_vertex != vertex_state_map.end()) { 
       
          if (vertex_state_map.erase(s.first) < 1) { // s.first is the vertex
            std::cerr << "Error: failed to remove an element from the map."
              << std::endl;
          } else {
            remove_count++;
          }  
 
        }
      }
    } // for 

    if(mpi_rank == 0) {
    std::cout << "Token Passing Batch [" << batch_count << "] | MPI Rank [" << mpi_rank 
      << "] | Removed " << remove_count << " vertices."<< std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    batch_count++;
  
  } // for 
}

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_TOKEN_PASSING_PATTERN_MATCHING_BATCH_HPP_INCLUDED

