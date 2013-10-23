BOOST_INC = /usr/include
CC      = g++

CFLAGS  = -I../include -I$(BOOST_INC) -fopenmp -O3

default: generate_graph run_bfs

generate_graph:
	$(CC) generate_graph.cpp -o generate_graph $(CFLAGS)

run_bfs:
	$(CC) run_bfs.cpp -o run_bfs $(CFLAGS)	

clean:
	rm generate_graph run_bfs
