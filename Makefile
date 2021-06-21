CXX = g++-8
#BOOST_INCLUDE = ~/boost

.PHONY: all connection_subgraph_undirected

all: connection_subgraph_undirected

connection_subgraph_undirected:
	#$(CXX_COMPILER) -std=c++17 -O3 -g -fopenmp -lrt -lpthread -I$(BOOST_INCLUDE) -I include -lstdc++fs connection_subgraph_undirected.cpp -o connection_subgraph_undirected
	$(CXX) -std=c++17 -O3 -g -fopenmp -lrt -lpthread -I include -lstdc++fs connection_subgraph_undirected.cpp -o connection_subgraph_undirected
.PHONY: clean
clean:	
	rm connection_subgraph_undirected
