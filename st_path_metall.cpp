#include <iostream>
#include <string>
#include <metall/metall.hpp>
#include <metall/container/experimental/jgraph/jgraph.hpp>

using namespace metall::container::experimental;

using graph_type = jgraph::jgraph<metall::manager::allocator_type<std::byte>>;

int main() {

  // Assumes that the format of the input file is JSON Lines
  //const std::string input_json_file_name("/usr/workspace/reza2/metall/metalljson/spoke.json");//("../spoke.json");
  const std::string input_json_file_name("/p/lustre2/havoqgtu/Spoke/spoke-20201130.json");

  {
    std::cout << "--- Create ---" << std::endl;
    metall::manager manager(metall::create_only, "/usr/workspace/reza2/metall/jgraph_spoke_obj");

    std::ifstream ifs(input_json_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open file: " << input_json_file_name << std::endl;
      std::abort();
    }

    std::string json_string;

    auto *graph = manager.construct<graph_type>(metall::unique_instance)(manager.get_allocator());

    // Parse each line of the input file one by one
    while (std::getline(ifs, json_string)) {

      try {

      // Parse the JSON string and allocate a JSON value object.
      // Pass a Metall allocator so that the contents of the object is allocated in Metall space.
      auto json_value = json::parse(json_string, graph->get_allocator());

      if (json_value.as_object()["type"].as_string() == "node") {
        const auto &vertex_id = json_value.as_object()["id"].as_string();
        graph->register_vertex(vertex_id);
        graph->vertex_value(vertex_id) = std::move(json_value);
      } else if (json_value.as_object()["type"].as_string() == "relationship") {
        const auto &src_id = json_value.as_object()["start"].as_object()["id"].as_string();
        const auto &dst_id = json_value.as_object()["end"].as_object()["id"].as_string();
        const auto &edge_id = json_value.as_object()["id"].as_string();
        graph->register_edge(src_id, dst_id, edge_id);
        graph->edge_value(edge_id) = std::move(json_value);
      }

      } catch (const std::exception& e) {
        std::cout << " a standard exception was caught, with message '" << e.what() << std::endl;
      }

    } // while

    std::cout << "#of vertices: " << graph->num_vertices() << std::endl;
    std::cout << "#of edges: " << graph->num_edges() << std::endl;
  }

  {
    std::cout << "\n--- Open ---" << std::endl;
    metall::manager manager(metall::open_read_only, "/usr/workspace/reza2/metall/jgraph_spoke_obj");

    const auto *graph = manager.find<graph_type>(metall::unique_instance).first;

    std::cout << "#of vertices: " << graph->num_vertices() << std::endl;
    std::cout << "#of edges: " << graph->num_edges() << std::endl;

    /*std::cout << "<Vertices>" << std::endl;
    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); vitr != vend; ++vitr) {
      std::cout << "Vertex ID = " << vitr->key() << std::endl;
      std::cout << "Vertex value = " << vitr->value() << std::endl;
    }

    // Access vertex values and edge values using the iterators
    std::cout << "\n<Edges>" << std::endl;
    for (auto vitr = graph->vertices_begin(), vend = graph->vertices_end(); vitr != vend; ++vitr) {
      for (auto eitr = graph->edges_begin(vitr->key()), eend = graph->edges_end(vitr->key()); eitr != eend; ++eitr) {
        std::cout << "Edge ID = " << eitr->key() << std::endl;
        std::cout << "Edge value = " << eitr->value() << std::endl;
      }
    }*/
    
  }

  return 0;
}
