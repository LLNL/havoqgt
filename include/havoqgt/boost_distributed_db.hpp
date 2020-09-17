// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_BOOST_DISTRIBUTED_DB_HPP_INCLUDED
#define HAVOQGT_BOOST_DISTRIBUTED_DB_HPP_INCLUDED

#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

namespace havoqgt {

/**
 *
 */
class db_create {};
class db_open {};

/**
 *
 */
class boost_distributed_db {
 private:
  struct header {
    int                comm_rank;
    int                comm_size;
    boost::uuids::uuid uuid;
    bool               clean_close;
  };

  typedef boost::interprocess::basic_managed_mapped_file<
      char, boost::interprocess::rbtree_best_fit<
                boost::interprocess::null_mutex_family>,
      boost::interprocess::iset_index>
      mapped_type;

 public:
  typedef mapped_type::segment_manager segment_manager_type;

  /**
   *
   */
  static void transfer(const char* src_base_fname,
                       const char* dest_base_fname) {
    std::string src_fname  = generate_filename(src_base_fname);
    std::string dest_fname = generate_filename(dest_base_fname);

    comm_world().barrier();
    double transfer_time = MPI_Wtime();
    {
      /*std::ifstream  src(src_fname.c_str(), std::ios::binary);
      std::ofstream  dst(dest_fname.c_str(),   std::ios::binary);
      dst << src.rdbuf();*/

      void*    buf;
      uint64_t src_total_size(0);
      bool     dest_direct_padded = false;
      size_t   buff_size          = 1024 * 1024 * 1;
      posix_memalign(&buf, 4096, buff_size);
      int source = -1;
#ifdef O_DIRECT
      source = open(src_fname.c_str(), O_RDONLY | O_DIRECT, 0);
      if (source == -1) {  // Attempt w/o O_DIRECT
        source = open(src_fname.c_str(), O_RDONLY /*| O_DIRECT*/, 0);
      }
#endif
      if (source == -1) {
        HAVOQGT_ERROR_MSG("Unable to open source file");
      }
      bool dest_direct = true;
      int  dest        = -1;
#ifdef O_DIRECT
      dest = open(dest_fname.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT,
                  0644);
      if (dest == -1) {  // Attempt w/o O_DIRECT
        dest_direct = false;
        dest        = open(dest_fname.c_str(),
                    O_WRONLY | O_CREAT | O_TRUNC /*| O_DIRECT*/, 0644);
      }
#endif
      if (dest == -1) {
        HAVOQGT_ERROR_MSG("Unable to open dest file");
      }
      size_t rsize;
      while ((rsize = read(source, buf, buff_size)) > 0) {
        src_total_size += rsize;
        if (dest_direct && (rsize % 4096 != 0)) {
          dest_direct_padded = true;
          char* bytes        = (char*)buf;
          for (size_t i = rsize; i < buff_size; ++i) {
            bytes[i] = 0;
          }
          rsize += 4096 - (rsize % 4096);
        }
        size_t wsize = write(dest, buf, rsize);
        if (wsize != rsize) {
          HAVOQGT_ERROR_MSG("Write Error");
        }
      }
      if (rsize == -1) {
        HAVOQGT_ERROR_MSG("Read Error");
      }
      //@todo need to truncate dest if used O_DIRECT
      close(source);
      close(dest);
      free(buf);
      if (dest_direct_padded) {
        // Truncate dest file if used O_DIRECT, because dest was padded;
        off_t off_t_src_total_size = src_total_size;
        if (off_t_src_total_size == src_total_size) {
          truncate(dest_fname.c_str(), src_total_size);
        }
      }
    }
    comm_nl().barrier();
    if (comm_nl().rank() == 0) {
      sync();
    }
    comm_world().barrier();
    if (comm_world().rank() == 0) {
      std::cout << "Transfer time = " << MPI_Wtime() - transfer_time
                << std::endl;
    }
    //    if(shrink) {
    //      mapped_type::shrink_to_fit(dest_fname.c_str());
    //    }
  }

  /**
   *
   */
  boost_distributed_db(db_create, const char* base_fname, double gbyte_per_rank) {
    int mpi_rank = comm_world().rank();
    int mpi_size = comm_world().size();

    uint64_t file_size = uint64_t(gbyte_per_rank * 1024 * 1024) * 1024ULL;

    init_rank_filename(base_fname);
    // std::cout << "file_size = file_size, " << file_size << ", fname = " <<
    // m_rank_filename << std::endl;
    if (rank_file_exists()) {
      // HAVOQGT_ERROR_MSG("File already exists.");
      ::remove(m_rank_filename.c_str());
    }

    m_pm = new mapped_type(boost::interprocess::create_only,
                           m_rank_filename.c_str(), file_size);

/*#ifdef HAVE_POSIX_FALLOCATE
    {
      int fd = open(m_rank_filename.c_str(), O_RDWR);
      if (fd == -1) {
        HAVOQGT_ERROR_MSG("Error opening file.");
      }
      int ret = posix_fallocate(fd, 0, file_size);
      if (ret != 0) {
        HAVOQGT_ERROR_MSG("posix_fallocate failed.");
      }
      close(fd);
    }
#else
//#warning posix_fallocate not found;  OSX?
#endif*/

    //
    // Create header
    header* phead =
        m_pm->construct<header>(boost::interprocess::unique_instance)();
    if (mpi_rank == 0) {
      phead->uuid = boost::uuids::random_generator()();
    }
    mpi_bcast(phead->uuid, 0, comm_world().mpi_comm());
    // std::cout << "Rank = " << mpi_rank << ", UUID = " << phead->uuid <<
    // std::endl;
    phead->comm_rank   = mpi_rank;
    phead->comm_size   = mpi_size;
    phead->clean_close = false;
  }

  /**
   *
   */
  boost_distributed_db(db_open, const char* base_fname) {
    int mpi_rank = comm_world().rank();
    int mpi_size = comm_world().size();

    init_rank_filename(base_fname);
    if (!rank_file_exists()) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": file not found.";
      throw std::runtime_error(error.str());
    }

    m_pm = new mapped_type(boost::interprocess::open_only,
                           m_rank_filename.c_str());

    //
    // Check for header
    std::pair<header*, std::size_t> ret =
        m_pm->find<header>(boost::interprocess::unique_instance);
    if (ret.second == 0) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__
            << ": header now found.";
      throw std::runtime_error(error.str());
    }

    //
    // Check UUID
    boost::uuids::uuid uuid;
    if (mpi_rank == 0) {
      uuid = ret.first->uuid;
    }
    mpi_bcast(uuid, 0, comm_world().mpi_comm());
    if (uuid != ret.first->uuid) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__
            << ": UUIDs don't match.";
      throw std::runtime_error(error.str());
    }
    // std::cout << "Rank = " << mpi_rank << ", UUID = " << ret.first->uuid <<
    // std::endl;
    if (ret.first->comm_rank != mpi_rank || ret.first->comm_size != mpi_size ||
        !ret.first->clean_close) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": DB corrupt.";
      throw std::runtime_error(error.str());
    }
  }

  /**
   *
   */
  ~boost_distributed_db() noexcept(false) {
    //
    // Mark clean close
    std::pair<header*, std::size_t> ret =
        m_pm->find<header>(boost::interprocess::unique_instance);
    if (ret.second == 0) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__
            << ": header now found.";
      throw std::runtime_error(error.str());
    }
    ret.first->clean_close = true;

    delete m_pm;
    m_pm = nullptr;
    comm_nl().barrier();
    bool shrink_ret = mapped_type::shrink_to_fit(m_rank_filename.c_str());
    comm_nl().barrier();
  }

  segment_manager_type* get_segment_manager() {
    return m_pm->get_segment_manager();
  }

 private:
  /**
   *
   */
  /*uint64_t get_file_size()
  {
    const char* fsize = getenv("HAVOQGT_DB_SIZE");
    if(fsize == NULL) {
      //return 700ULL * 1024ULL * 1024ULL * 1024ULL / 24;
      //return 799700000000ULL;
      //return 500ULL * 1024ULL * 1024ULL * 1024ULL;
      // causing probs return 799700000000ULL / 24ULL;
      // MAXES OUT CATALST!! return 799595142400 / 24ULL;
      //return 512*1024*1024 / comm_nl().size();
      return 64ULL * 1024 *1024 *1024 / 24;
    }
    return boost::lexical_cast<uint64_t>(fsize);
  }*/

  /**
   *
   */
  bool rank_file_exists() {
    std::ifstream fin(m_rank_filename.c_str());
    return fin.good();
  }

  /**
   *
   */
  void init_rank_filename(const char* base_fname) {
    m_rank_filename = generate_filename(base_fname);
  }

  static std::string generate_filename(const char* base_fname) {
    int               mpi_rank = comm_world().rank();
    int               mpi_size = comm_world().size();
    std::stringstream sstr;
    sstr << base_fname << "_" << mpi_rank << "_of_" << mpi_size;
    return sstr.str();
  }

 private:
  mapped_type* m_pm;
  std::string  m_rank_filename;
};

}  // havoqgt

#endif  // HAVOQGT_BOOST_DISTRIBUTED_DB_HPP_INCLUDED
