#include <gtest/gtest.h>
#include <include/create_delegate_graph.hpp>
#include <include/havoqgt_setup.hpp>
#include <include/input_graph.hpp>
#include <include/read_delegate_graph.hpp>
#include <include/util.hpp>

namespace havoqgt { namespace test {

typedef graph_type::edge_data<edge_data_type, 
  bip::allocator<edge_data_type, segment_manager_t>> edge_data_t;

const std::string graph_unique_instance_name = "graph_obj";
const std::string edge_data_unique_instance_name = "graph_edge_data_obj";
const std::string input_graph_file_name = "/dev/shm/test_havoqgt_graph_6";

std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> input_graph;
std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> vec_global_edges; 
std::set<uint64_t> delegate_vertices;

void create_local_edge_list(graph_type& g, edge_data_t& edge_data_ptr, 
                            vloc_type vertex,  
                            std::vector<uint64_t>& edge_source, 
                            std::vector<uint64_t>& edge_target, 
                            std::vector<uint64_t>& edge_data) {
  for(eitr_type eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); 
      ++eitr) {
    edge_source.push_back(g.locator_to_label(vertex));
    edge_target.push_back(g.locator_to_label(eitr.target()));
    //edge_data.push_back((uint64_t)eitr.edge_data());   
    edge_data.push_back((uint64_t)edge_data_ptr[eitr]); 
  }
}

void test_Delegate_Graph_Weighted_Edges() {
  int mpi_rank(0); 
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));

  input_graph = grid_graph_weighted_edges();

  //MPI_Barrier(MPI_COMM_WORLD);

  create_delegate_graph(input_graph, input_graph_file_name, 
                        graph_unique_instance_name, 
                        edge_data_unique_instance_name, mpi_rank);

  MPI_Barrier(MPI_COMM_WORLD);

  //graph_type* g = read_delegate_graph("testD", mpi_rank);
  
  havoqgt::distributed_db ddb(havoqgt::db_open(), 
                              input_graph_file_name.c_str());
  graph_type *g = ddb.get_segment_manager()->
    find<graph_type>(graph_unique_instance_name.c_str()).first;
  assert(g != nullptr);

  edge_data_t* edge_data_ptr = ddb.get_segment_manager()->
    find<edge_data_t>(edge_data_unique_instance_name.c_str()).first;
  assert(edge_data_ptr != nullptr);  

  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    std::cout << "MPI Rank: " << mpi_rank << " Graph Loaded Ready." 
    << std::endl;
  }
  //graph->print_graph_statistics(); 
  // causes MPI runtime exception on flash.llnl.gov 

  MPI_Barrier(MPI_COMM_WORLD);
  
  // build a local copy of the entire graph on each rank 
  
  std::vector<uint64_t> vec_local_edge_source;
  std::vector<uint64_t> vec_global_edge_source;
  std::vector<uint64_t> vec_local_edge_target;
  std::vector<uint64_t> vec_global_edge_target;
  std::vector<uint64_t> vec_local_edge_data;
  std::vector<uint64_t> vec_global_edge_data;  

  for (vitr_type vitr = g->vertices_begin(); vitr != g->vertices_end(); 
       ++vitr) {
    vloc_type vertex = *vitr;
    if (vertex.is_delegate()) {
      delegate_vertices.insert(g->locator_to_label(vertex));
    }
    create_local_edge_list(*g, *edge_data_ptr, vertex, vec_local_edge_source, 
                           vec_local_edge_target, vec_local_edge_data); 
  }

  for (vitr_type vitr = g->delegate_vertices_begin(); 
       vitr != g->delegate_vertices_end(); ++vitr) { 
    vloc_type vertex = *vitr;
    if (vertex.is_delegate()) {
      delegate_vertices.insert(g->locator_to_label(vertex));
    } 
    create_local_edge_list(*g, *edge_data_ptr, vertex, vec_local_edge_source,
                           vec_local_edge_target, vec_local_edge_data);
  }
 
  MPI_Barrier(MPI_COMM_WORLD);

  // gather global edges
  mpi_all_gather(vec_local_edge_source, vec_global_edge_source, MPI_COMM_WORLD); 
  mpi_all_gather(vec_local_edge_target, vec_global_edge_target, MPI_COMM_WORLD);
  mpi_all_gather(vec_local_edge_data, vec_global_edge_data, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  for (size_t i; i < vec_global_edge_source.size(); ++i) {
    vec_global_edges.push_back(std::make_tuple(vec_global_edge_source[i], 
                                               vec_global_edge_target[i], 
                                               vec_global_edge_data[i]));
  }

  // sort edges
  std::stable_sort(vec_global_edges.begin(),vec_global_edges.end(),
    [](const std::tuple<uint64_t, uint64_t, uint64_t>& a,
       const std::tuple<uint64_t, uint64_t, uint64_t>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });
  
  for (size_t i = 0; i < grid_graph_weighted_edges_offset.size() - 1; ++i) {
     size_t start = grid_graph_weighted_edges_offset[i];
     size_t end = grid_graph_weighted_edges_offset[i+1];
      
     std::stable_sort(vec_global_edges.begin() + start, 
                      vec_global_edges.begin() + end,
       [](const std::tuple<uint64_t, uint64_t, uint64_t>& a,
          const std::tuple<uint64_t, uint64_t, uint64_t>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });    
  } // sort edges  

  MPI_Barrier(MPI_COMM_WORLD);
}

TEST(my_test, test_Delegate_Graph_Weighted_Edges_Verify_Edges) {
  test_Delegate_Graph_Weighted_Edges(); 
  MPI_Barrier(MPI_COMM_WORLD);  
  EXPECT_EQ(input_graph.size(), vec_global_edges.size());
  EXPECT_TRUE(std::equal(input_graph.begin(), input_graph.end(),
              vec_global_edges.begin()));
}

TEST(my_test, test_Delegate_Graph_Weighted_Edges_Verify_Delegates) {
  MPI_Barrier(MPI_COMM_WORLD);
  std::set<uint64_t> expected = {6, 7, 8};
  EXPECT_EQ(expected.size(), delegate_vertices.size());
  EXPECT_TRUE(std::equal(expected.begin(), expected.end(),
              delegate_vertices.begin()));
}

}} //end namespace apgaf::test

//mpi main for gteset
GTEST_API_ int main(int argc, char **argv) {
  // set up environment
  int mpi_rank(0), mpi_size(0);
  havoqgt::havoqgt_init(&argc, &argv);
  {
  CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
  CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
  havoqgt::get_environment();

  if (mpi_rank == 0) {
    std::cout << "MPI initialized with " << mpi_size << " ranks." << std::endl;
    //havoqgt::get_environment().print();
    //print_system_info(false); 
  }
  MPI_Barrier(MPI_COMM_WORLD);
   
  // execute tests
  testing::InitGoogleTest(&argc, argv);
  int to_return RUN_ALL_TESTS();
  }

  // delete the generated files

  havoqgt::havoqgt_finalize();

  //return to_return;
}
