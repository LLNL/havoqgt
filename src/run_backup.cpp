#include <iostream>
#include <sstream>

#include <havoqgt/distributed_db.hpp>

using namespace havoqgt;

int main(int argc, char** argv) {
  if( argc < 3 ) {
    std::cout << "Usage: run_backup <Edgelist> <Metadata> <Backup EdgeList> <Backup Metadata>" << std::endl;
    exit(1);
  }
  std::string output_filename = std::string( argv[1]);
  std::string backup_filename = std::string( argv[2] );
  std::string metadata_output_filename;
  std::string backup_metadata_filename;
  if( argc > 3 ) {
    metadata_output_filename = std::string( argv[3] );
    backup_metadata_filename = std::string( argv[4] );
  }
 
  havoqgt_init(&argc, &argv);
  havoqgt::get_environment();

  havoqgt::distributed_db::transfer(output_filename.c_str(), backup_filename.c_str());
  if( argc > 3 ) {
    havoqgt::distributed_db::transfer(metadata_output_filename.c_str(), backup_metadata_filename.c_str());
  }
  
  havoqgt::havoqgt_finalize();
}
