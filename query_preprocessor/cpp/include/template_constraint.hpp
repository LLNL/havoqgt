#pragma once

namespace prunejuice { namespace pattern {  

template <typename Vertex, typename Edge, typename EdgeSet,
  typename EdgeSetVector, typename EdgeListTupleVector>
class template_constraint {

public:

  typedef size_t ConstraintType;
  static constexpr ConstraintType PATH = 0;
  static constexpr ConstraintType CYCLE = 1;
  static constexpr ConstraintType TDS = 2;
  static constexpr ConstraintType ENUMERATION = 3;
  
  template_constraint() {}
 
  template_constraint(EdgeSet _edgeset, Edge _edgeset_hash, 
    size_t _constraint_ID, ConstraintType _constraint_type) : 
    edgeset(_edgeset), edgeset_hash(_edgeset_hash), 
    constraint_ID(_constraint_ID), constraint_type(_constraint_type) {      
  }

  // CYCLE
  template_constraint(EdgeSet _edgeset, Edge _edgeset_hash,
    EdgeSetVector _edgeset_vector, size_t _constraint_ID, 
    ConstraintType _constraint_type) :
    edgeset(_edgeset), edgeset_hash(_edgeset_hash), 
    edgeset_vector(_edgeset_vector),
    constraint_ID(_constraint_ID), constraint_type(_constraint_type) {
  }

  // TDS
  template_constraint(EdgeSet _edgeset, Edge _edgeset_hash,
    EdgeListTupleVector _edgelist_vector, size_t _constraint_ID, 
    ConstraintType _constraint_type) :
    edgeset(_edgeset), edgeset_hash(_edgeset_hash), 
    edgelist_vector(_edgelist_vector),
    constraint_ID(_constraint_ID), constraint_type(_constraint_type) {
  }

  // ENUMERATION
  template_constraint(Edge _edgeset_hash,
    EdgeListTupleVector _edgelist_vector, size_t _constraint_ID, 
    ConstraintType _constraint_type) :
    edgeset_hash(_edgeset_hash), 
    edgelist_vector(_edgelist_vector),
    constraint_ID(_constraint_ID), constraint_type(_constraint_type) {
  }

  ~template_constraint() {}

  EdgeSet edgeset; // for a unique path - the type is same for all ConstraintTypes
  EdgeSetVector edgeset_vector; // for all paths - CYCLE 
  EdgeListTupleVector edgelist_vector; // for all paths - TDS, ENUMERATION
  Edge edgeset_hash; // not useful for TDS and ENUMERATION 
  ConstraintType constraint_type;
  ConstraintType algorithm_type; // TODO: batch, no aggregation ?    
  size_t constraint_ID;

};  

}} // end namespace prunejuice::pattern
