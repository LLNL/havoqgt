CXX_COMPILER=/usr/tce/packages/gcc/gcc-8.1.0/bin/g++
BOOST_INCLUDE=/usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-4.9.3

.PHONY: all setup_application_edit_distance

all: setup_application setup_application_edit_distance

# -lstdc++fs must be included after the source file

setup_application_edit_distance:
	$(CXX_COMPILER) -std=c++17 -O3 -g -fopenmp -lrt -lpthread -ltcmalloc_minimal -I$(BOOST_INCLUDE) -I include setup_application_edit_distance.cpp -lstdc++fs -o setup_application_edit_distance

.PHONY: clean
clean:	
	rm setup_application_edit_distance
