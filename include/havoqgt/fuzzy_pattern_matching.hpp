#ifndef HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED
#define HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED

#include <array>

#include <havoqgt/visitor_queue.hpp>
#include <havoqgt/detail/visitor_priority_queue.hpp>

namespace havoqgt { namespace mpi {

static constexpr size_t max_walk_history_size = 15;

template<typename Visitor>
class fpm_queue {

public:
  fpm_queue() {}

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

template<typename Graph>
class fpm_walker {
public:
  typedef typename Graph::vertex_locator vertex_locator;
  typedef typename Graph::edge_iterator eitr_type;

  fpm_walker() : 
    current_pattern_index(0), do_update_vertex_rank(false) {}

  fpm_walker(vertex_locator _vertex, bool _do_update_vertex_rank = false) :
    vertex(_vertex),
    current_pattern_index(0), 
    do_update_vertex_rank(_do_update_vertex_rank) {
    // has no parent
    walk_history[current_pattern_index] = _vertex;
  }

  template<typename WalkHistory>
  fpm_walker(vertex_locator _vertex, size_t _current_pattern_index, 
    WalkHistory& _parent_walk_history) :
    vertex(_vertex), 
    current_pattern_index(_current_pattern_index), 
    do_update_vertex_rank(false) {

    if (current_pattern_index == 0) {
      // has no parent
      // Important: I think it should never get here
      walk_history[current_pattern_index] = _vertex;
    } else {
      // copy _parent_walk_history to walk_history
      std::copy(_parent_walk_history.begin(), 
        _parent_walk_history.begin() + current_pattern_index, 
        walk_history.begin());
      walk_history[current_pattern_index] = _vertex;
    } 
  }

  ~fpm_walker() {
  } 

  template<typename AlgData> 
  bool pre_visit(AlgData& alg_data) const {
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);
    
    if (do_update_vertex_rank) {
      std::get<3>(alg_data)[vertex]++;   
      std::cout << "Pre-visit Updating rank " << std::get<3>(alg_data)[vertex] << std::endl;
      return false;
    } 
    
    if (vertex_data != pattern[current_pattern_index]) {
      return false;
    } else { 
      if (current_pattern_index == pattern.size() - 1) {
        std::cout << "Pre-visit Found" << std::endl;
        //return false;
        // TODO: can not update the ranks here, do not have reference to the graph or visitor_queue
        // Requires a new version of the pre_visit function in the visitor_queue class 
        return true;
      }
      // TODO: degree > 0 
      return true;
    }
    return true;
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool init_visit(Graph& g, VisitorQueueHandle vis_queue, 
    AlgData& alg_data) const {
     if (do_update_vertex_rank) {
       std::cerr << "Wrong code branch." << std::endl;
       return false;
     }
    // TODO: verifications? degree > 0 
    return visit(g, vis_queue, alg_data);
  }

  template<typename VisitorQueueHandle, typename AlgData>
  bool visit(Graph& g, VisitorQueueHandle vis_queue, AlgData& alg_data) const {
    if (do_update_vertex_rank) {
       std::cerr << "Wrong code branch." << std::endl;
       return false;
    }
    auto vertex_data = std::get<0>(alg_data)[vertex];
    auto pattern = std::get<1>(alg_data);
    auto pattern_indices = std::get<2>(alg_data);

    int mpi_rank(0); // test
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)); // test

    // Important: metadata verification must be doen here
    if (vertex_data != pattern[current_pattern_index]) {
      return false;
    } else {
      if (current_pattern_index == pattern.size() - 1) {

        for (size_t i = 0; i <= current_pattern_index; i++) { // TODO: loop termination condition i < current_pattern_index this will creat one less visitor
          std::cout << g.locator_to_label(walk_history[i]) << " " << std::get<0>(alg_data)[walk_history[i]] << std::endl; // test 
          // add vertex to the visitor queue to update rank
          fpm_walker new_visitor_update_rank(walk_history[i], true);
          vis_queue->queue_visitor(new_visitor_update_rank);
        } 
 
        //std::get<3>(alg_data)[vertex]++; // TODO: enable this ?        

        std::cout << "Found " << std::endl;

        return false;
      }
    }

    size_t next_pattern_index = current_pattern_index + 1;
    if (next_pattern_index > pattern.size() - 1) {
      return false; 
    }
     
    for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
      vertex_locator neighbor = eitr.target();
      auto neighbor_data = std::get<0>(alg_data)[neighbor];      
      
      //if (g.locator_to_label(vertex) < 90) // test
      //  std::cout << mpi_rank << " | " << g.locator_to_label(vertex) << " > " << g.locator_to_label(neighbor) << " : " << neighbor_data << std::endl; // test

      // pre-cloning verifications
      
      // TODO: neighbour degree > 0 if next_pattern_index < pattern.size() - 1 
     
      bool invalid_vertex = false;
      if (pattern_indices[next_pattern_index] == next_pattern_index) { // verify duplicate vertex
        for (size_t i = 0; i < next_pattern_index; i++) {
          //if (walk_history[i] == neighbor) { // does neighbor exist in parent's walk history? // doe not work
          if (g.locator_to_label(walk_history[i]) == g.locator_to_label(neighbor)) {
            invalid_vertex = true; 
            break; // duplicate exists, do not clone neighbor 
          }
        } // for        

        if (invalid_vertex) {
          continue;
        }

      } else if (pattern_indices[next_pattern_index] < next_pattern_index) { // loop closing detection
        //if (walk_history[pattern_indices[next_pattern_index]] != neighbor) { // does not work 
        if (g.locator_to_label(walk_history[pattern_indices[next_pattern_index]]) != g.locator_to_label(neighbor)) {
          continue;
        } 
      } else {
        std::cerr << "Wrong code branch." << std::endl;
        continue; 
      }

      // Important: you cannot do this verification if the neighbour is not a local vertex
      //if (neighbor_data == pattern[next_pattern_index]) {
        //if (next_pattern_index == pattern.size() - 1) {

          // test print
          //for (size_t i = 0; i < next_pattern_index; i++) {
            //std::cout << g.locator_to_label(walk_history[i]) << " " << std::get<0>(alg_data)[walk_history[i]] << std::endl;  
          //}
          //std::cout << mpi_rank << " | " << g.locator_to_label(neighbor) << " " << std::get<0>(alg_data)[neighbor] << std::endl;
 
          //std::cout << "Found " << walk_history.size() << std::endl;
          //continue;
        //}

        // create clone
        //if (g.locator_to_label(vertex) < 90) // test
        //  std::cout << g.locator_to_label(vertex) << " >> clone >> " << g.locator_to_label(neighbor) << std::endl; // test

        fpm_walker new_visitor(neighbor, next_pattern_index, walk_history);           
        vis_queue->queue_visitor(new_visitor);   
      //}
  
    }
    return true;
  }

  friend inline bool operator>(const fpm_walker& v1, const fpm_walker& v2) {
    return false;
  }

  friend inline bool operator<(const fpm_walker& v1, const fpm_walker& v2) {
    return false;
  }

  vertex_locator vertex;
  size_t current_pattern_index;
  std::array<vertex_locator, max_walk_history_size> walk_history;
  bool do_update_vertex_rank;
}; //__attribute__ ((packed));

template <typename TGraph, typename VertexMetaData, typename PatternData, 
  typename PatternIndices, typename VertexRank>
void  fuzzy_pattern_matching(TGraph* g, VertexMetaData& vertex_metadata, 
  PatternData& pattern, PatternIndices& pattern_indices, VertexRank& vertex_rank) {

  typedef fpm_walker<TGraph> visitor_type;
  auto alg_data = std::forward_as_tuple(vertex_metadata, pattern, pattern_indices, vertex_rank);
  auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);
  vq.init_visitor_traversal_new(); 
  MPI_Barrier(MPI_COMM_WORLD);
  vertex_rank.all_reduce();
  MPI_Barrier(MPI_COMM_WORLD);  
} 

}} //end namespace havoqgt::mpi

#endif //HAVOQGT_MPI_FUZZY_PATTERN_MATCHING_HPP_INCLUDED
