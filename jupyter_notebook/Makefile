CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++
BOOST_INCLUDE=/usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-4.9.3

.PHONY: all connection_subgraph_undirected

all: connection_subgraph_undirected

connection_subgraph_undirected:
	$(CXX_COMPILER) -std=c++17 -O3 -g -fopenmp -lrt -lpthread -I$(BOOST_INCLUDE) -I include -lstdc++fs connection_subgraph_undirected.cpp -o connection_subgraph_undirected

.PHONY: clean
clean:	
	rm connection_subgraph_undirected
