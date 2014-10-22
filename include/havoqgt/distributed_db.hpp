/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or process
 * disclosed, or represents that its use would not infringe privately-owned rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring by
 * the United States Government or Lawrence Livermore National Security, LLC. The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */



#ifndef HAVOQGT_DISTRIBUTED_DB_HPP_INCLUDED
#define HAVOQGT_DISTRIBUTED_DB_HPP_INCLUDED

#include <boost/interprocess/managed_mapped_file.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>            
#include <boost/uuid/uuid_generators.hpp> 
#include <boost/uuid/uuid_io.hpp>         

#include <string>
#include <sstream>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


namespace havoqgt {

/**
 *
 */
class db_create {};
class db_open {};

/**
 *
 */
class distributed_db {
private:
  struct header {
    int comm_rank;
    int comm_size;
    boost::uuids::uuid uuid;
    bool clean_close;
  };

  typedef boost::interprocess::managed_mapped_file mapped_type;

public:
  typedef mapped_type::segment_manager segment_manager_type;

  /**
   *
   */
  distributed_db(db_create, MPI_Comm comm, const char* base_fname)
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
    CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

    init_rank_filename(comm, base_fname);
    if(rank_file_exists())
    {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": file already exists.";
      throw std::runtime_error(error.str());
    }

    m_pm = new mapped_type(boost::interprocess::create_only, m_rank_filename.c_str(), get_file_size()); 

    {
      int fd  = open(m_rank_filename.c_str(), O_RDWR);
      if(fd == -1) {
        std::cerr << "ERROR opening file" << std::endl; exit(-1);
      }
      int ret = posix_fallocate(fd,0,get_file_size());
      if(ret != 0)
      {
        std::cerr << "ERROR:  posix_fallocate failed" << std::endl;
      }
      close(fd);
    }

    //
    // Create header
    header* phead = m_pm->construct<header>(boost::interprocess::unique_instance)();
    if(mpi_rank == 0) {
      phead->uuid = boost::uuids::random_generator()();
    }
    mpi::mpi_bcast(phead->uuid,0,comm);
    std::cout << "Rank = " << mpi_rank << ", UUID = " << phead->uuid << std::endl;
    phead->comm_rank = mpi_rank;
    phead->comm_size = mpi_size;
    phead->clean_close = false;
  }

  /**
   *
   */
  distributed_db(db_open, MPI_Comm comm, const char* base_fname)
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
    CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );

    init_rank_filename(comm, base_fname);
    if(!rank_file_exists())
    {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": file not found.";
      throw std::runtime_error(error.str());
    }
    
    m_pm = new mapped_type(boost::interprocess::open_only, m_rank_filename.c_str()); 

    //
    // Check for header
    std::pair<header*, std::size_t> ret = m_pm->find<header>(boost::interprocess::unique_instance);
    if(ret.second == 0) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": header now found.";
      throw std::runtime_error(error.str());
    }

    //
    // Check UUID
    boost::uuids::uuid uuid;
    if(mpi_rank == 0) {
      uuid = ret.first->uuid;
    }
    mpi::mpi_bcast(uuid,0,comm);
    if(uuid != ret.first->uuid)
    {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": UUIDs don't match.";
      throw std::runtime_error(error.str());
    }
    std::cout << "Rank = " << mpi_rank << ", UUID = " << ret.first->uuid << std::endl;
    if(ret.first->comm_rank != mpi_rank ||
       ret.first->comm_size != mpi_size ||
       !ret.first->clean_close) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": DB corrupt.";
      throw std::runtime_error(error.str());
    }
  }

  /**
   *
   */
  ~distributed_db()
  {
    //
    // Mark clean close
    std::pair<header*, std::size_t> ret = m_pm->find<header>(boost::interprocess::unique_instance);
    if(ret.second == 0) {
      std::stringstream error;
      error << "ERROR: " << __FILE__ << ":" << __LINE__ << ": header now found.";
      throw std::runtime_error(error.str());
    }
    ret.first->clean_close = true;

    delete m_pm;
    m_pm = nullptr;
    bool shrink_ret = mapped_type::shrink_to_fit(m_rank_filename.c_str());
  }

  segment_manager_type* get_segment_manager()
  {
    return m_pm->get_segment_manager();
  }

private:

  /**
   *
   */
  uint64_t get_file_size() 
  {
    const char* fsize = getenv("HAVOQGT_DB_SIZE");
    if(fsize == NULL) 
      //return 700ULL * 1024ULL * 1024ULL * 1024ULL / 24;
      //return 799700000000ULL;
      //return 500ULL * 1024ULL * 1024ULL * 1024ULL;
      // causing probs return 799700000000ULL / 24ULL;
      return 799595142400 / 24ULL;
    return boost::lexical_cast<uint64_t>(fsize);
  } 


  /**
   *
   */
  bool rank_file_exists()
  {
    std::ifstream fin(m_rank_filename.c_str());
    return fin.good();
  }

  /**
   *
   */
  void init_rank_filename(MPI_Comm comm, const char* base_fname)
  {
    int mpi_rank(0), mpi_size(0);
    CHK_MPI( MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank) );
    CHK_MPI( MPI_Comm_size( MPI_COMM_WORLD, &mpi_size) );
    std::stringstream sstr;
    sstr << base_fname << "_" << mpi_rank << "_of_" << mpi_size;
    m_rank_filename = sstr.str();
  }


private:
  mapped_type* m_pm;
  std::string  m_rank_filename;

};


} //havoqgt



#endif //HAVOQGT_DISTRIBUTED_DB_HPP_INCLUDED
