#include <iostream>
#include <fstream>
#include <random>
#include <unordered_set>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <havoqgt/environment.hpp>
#include <havoqgt/parallel_edge_list_reader.hpp>
#include <havoqgt/rmat_edge_generator.hpp>

#include <havoqgt/graphstore/rhh/rhh_defs.hpp>
#include <havoqgt/graphstore/rhh/rhh_utilities.hpp>
#include <havoqgt/graphstore/rhh/rhh_allocator_holder.hpp>

#include <havoqgt/graphstore/graphstore_utilities.hpp>

using mapped_file_type     = boost::interprocess::managed_mapped_file;
using segment_manager_type = boost::interprocess::managed_mapped_file::segment_manager;

using key_type    = uint64_t;
using value_type  = uint64_t;

#if 0
#include <havoqgt/graphstore/rhh/rhh_container.hpp>
using rhh_type = graphstore::rhh_container<key_type, value_type, size_t, segment_manager_type>;
#else
#include <havoqgt/graphstore/rhh/blocked_rhh_container.hpp>
using rhh_type = graphstore::blocked_rhh_container<key_type, value_type, size_t, segment_manager_type>;
#endif

void usage()  {
  if(havoqgt::havoqgt_env()->world_comm().rank() == 0) {
    std::cerr << "Usage: -i <string> -s <int>\n"
              << " -s <string>   - output graph base filename (default is /dev/shm/segment_file)\n"
                 //         << " -e <string>   - filename that has a list of edgelist files\n"
              << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& segmentfile_name, std::vector<std::string>& edgelist_files) {
  std::cout << "CMD line:";
  for (int i=0; i<argc; ++i) {
    std::cout << " " << argv[i];
  }
  std::cout << std::endl;

  segmentfile_name = "/dev/shm/segment_file";

  bool found_segmentfile_name_ = true;
  bool found_edgelist_filename_ = true;

  char c;
  bool prn_help = false;
  while ((c = getopt(argc, argv, "s:e:h")) != -1) {
    switch (c) {
      case 'h':
        prn_help = true;
        break;
      case 's':
        found_segmentfile_name_ = true;
        segmentfile_name = optarg;
        break;
      case 'e':
      {
        found_edgelist_filename_ = true;
        std::string fname(optarg);
        std::ifstream fin(fname);
        std::string line;
        if (!fin.is_open()) {
          std::cerr << fname << std::endl;
          HAVOQGT_ERROR_MSG("Unable to open a file");
        }
        while (std::getline(fin, line)) {
          edgelist_files.push_back(line);
        }
        break;
      }
      default:
        std::cerr << "Unrecognized option: "<<c<<", ignore."<<std::endl;
        prn_help = true;
        break;
    }
  }
  if (prn_help || !found_segmentfile_name_ || !found_edgelist_filename_) {
    usage();
    exit(-1);
  }
}


/// basic test cases
void test1(rhh_type** _rhh, const size_t num_keys, const size_t num_duplicates)
{

  rhh_type* rhh = *_rhh;
  /// insertion and deletion with value iterator
  {
    /// init
    for (uint64_t k = 0; k <= num_duplicates; ++k) {
      for (uint64_t i = 0; i < num_keys; ++i) {
        graphstore::rhh::insert(&rhh, i, k);
      }
    }
    /// num edges
    assert(rhh->size() == num_keys * (num_duplicates + 1));

//    rhh->print_all_element();

    /// delete all edges with value iterator
    for (uint64_t i = 0; i < num_keys; ++i) {
      size_t cnt = 0;
      for (auto itr = rhh->find(i), end = rhh->find_end(); itr != end; ++itr) {
        rhh->erase(itr);
        ++cnt;
      }
      assert(cnt == (num_duplicates + 1));
    }
    /// num edges
    assert(rhh->size() == 0);
  }

  /// insertion and deletion with whole iterator
  {
    /// init
    for (uint64_t k = 0; k <= num_duplicates; ++k) {
      for (uint64_t i = 0; i < num_keys; ++i) {
        graphstore::rhh::insert(&rhh, i, k);
      }
    }
    /// num edges
    assert(rhh->size() == num_keys * (num_duplicates + 1));

    /// delete all edges
    size_t cnt = 0;
    for (auto itr = rhh->begin(), end = rhh->end(); itr != end; ++itr) {
      rhh->erase(itr);
      ++cnt;
    }
    assert(cnt == num_keys * (num_duplicates + 1));

    /// num edges
    assert(rhh->size() == 0);
  }

  *_rhh = rhh;
}


int main(int argc, char** argv) {

  /// --- parse argments ---- ///
  std::string segmentfile_name;
  std::vector<std::string> edgelist_files;
  parse_cmd_line(argc, argv, segmentfile_name, edgelist_files);


  /// --- create a segument file --- ///
  size_t graph_capacity = std::pow(2, 25);
  graphstore::utility::interprocess_mmap_manager::delete_file(segmentfile_name);
  graphstore::utility::interprocess_mmap_manager mmap_manager(segmentfile_name, graph_capacity);

  /// --- allocate a graphstore --- ///
  graphstore::rhh::init_allocator<typename rhh_type::allocator, segment_manager_type>(mmap_manager.get_segment_manager());
  rhh_type* rhh = rhh_type::allocate(2);

  std::cout << "num keys, num_duplicates" << std::endl;
  for (size_t num_duplicates = 0; num_duplicates <= 64; ++num_duplicates) {
    for (size_t num_keys = 1; num_keys <= (1ULL << 4); ++num_keys) {
      std::cout << "test: " << num_keys << ", " << num_duplicates << std::endl;
      test1(&rhh, num_keys, num_duplicates);
    }
  }

  std::cout << "All tests successed!!!" << std::endl;
  return 0;
}


