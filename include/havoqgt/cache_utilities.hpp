
/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Steven Feldman <feldman12@llnl.gov>.
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

#ifndef _HAVOQGT_CACHE_UTILI_HPP
#define _HAVOQGT_CACHE_UTILI_HPP

#include <sys/mman.h>
#include <unistd.h>

template<typename Vector>
char * get_address(Vector &vec) {
  uintptr_t temp = reinterpret_cast<uintptr_t>(&(vec[0]));
  temp -= temp % 4096;
  return reinterpret_cast<char *>(temp);
}


template<typename Vector>
size_t get_length(Vector &vec) {
  size_t length = vec.size() * sizeof(vec[0]);
  length += (4096 - length%4096);
  return length;
}


template<typename Vector>
void advise_vector_rand(Vector &vec) {
  char * addr = get_address(vec);
  size_t length = get_length(vec);

  int t = madvise(addr, length, MADV_RANDOM);
  assert(t == 0);
}


template<typename Vector>
void flush_advise_vector_dont_need(Vector &vec) {
  char * addr = get_address(vec);
  size_t length = get_length(vec);

  int t = msync(addr, length, MS_SYNC);
  if (t == 0) {
    int t2 = madvise(addr, length, MADV_DONTNEED);
    assert(t2 == 0);
  } else {
    assert(t == 0);
  }

}


template<typename Vector>
void flush_vector(Vector &vec) {
  char * addr = get_address(vec);
  size_t length = get_length(vec);

  int t = msync(addr, length, MS_SYNC);
  assert(t == 0);
}

template<typename Vector>
void flush_advise_vector(Vector &vec) {
  char * addr = get_address(vec);
  size_t length = get_length(vec);

  int t = msync(addr, length, MS_SYNC);
  if (t == 0) {
    int t2 = madvise(addr, length, MADV_DONTNEED);
    assert(t2 == 0);
    t2 = madvise(addr, length, MADV_RANDOM);
    assert(t2 == 0);
  } else {
    assert(t == 0);
  }
}

#ifndef DIRTY_THRESHOLD_GB
  #define DIRTY_THRESHOLD_GB 70
#endif

uint32_t get_disk_utilization() {
  uint32_t dirty_kb;

  //FILE *pipe;
  //pipe = popen("df -h /l/ssd | grep /dev/md0  | awk '{print $3}'", "r" );
  //fscanf(pipe, "%u", &dirty_kb);
  //pclose(pipe);

  return dirty_kb;
}

void print_system_info(bool print_dimmap) {
  printf("#################################################################\n");
  printf("System Information\n");
  printf("#################################################################\n");


  printf("\n-----------------------------------------------------------------\n");
  system("echo \"SLURM_NODELIST = $SLURM_NODELIST \"");
  printf("-----------------------------------------------------------------\n");

  printf("\n-----------------------------------------------------------------\n");
  printf("Tuned Info:\n");
  printf("-----------------------------------------------------------------\n");
  system("echo \"/proc/sys/vm/dirty_ratio = $(cat /proc/sys/vm/dirty_ratio)\"");
  system("echo \"/proc/sys/vm/dirty_background_ratio = $(cat /proc/sys/vm/dirty_background_ratio)\"");
  system("echo \"/proc/sys/vm/dirty_expire_centisecs = $(cat /proc/sys/vm/dirty_expire_centisecs)\"");

  printf("\n-----------------------------------------------------------------\n");
  printf("df -h /l/ssd\n");
  printf("-----------------------------------------------------------------\n");
  system("df -h /l/ssd");

  printf("\n-----------------------------------------------------------------\n");
  printf("ls /l/ssd\n");
  printf("-----------------------------------------------------------------\n");
  system("ls /l/ssd");


  printf("\n-----------------------------------------------------------------\n");
  printf("io-stat -m | grep md0 2>&1\n");
  printf("-----------------------------------------------------------------\n");
  system("iostat -m | grep Device 2>&1");
  system("iostat -m | grep md0 2>&1");


  if (print_dimmap) {
    printf("\n-----------------------------------------------------------------\n");
    system("echo \"/proc/di-mmap-runtimeA-stats = $(cat /proc/di-mmap-runtimeA-stats)\"");
    printf("-----------------------------------------------------------------\n");
  }



  printf("\n\n");

}

void print_dmesg() {
  printf("\n-----------------------------------------------------------------\n");
  printf("dmesg\n");
  printf("-----------------------------------------------------------------\n");

  system("dmesg");

  printf("\n\n");
}

uint32_t get_dirty_pages() {
  uint32_t dirty_kb;

  //FILE *pipe;
  //pipe = popen("grep Dirty /proc/meminfo | awk '{print $2}'", "r" );
  //fscanf(pipe, "%u", &dirty_kb);
  //pclose(pipe);

  return dirty_kb;
}

bool check_dirty_pages() {
  uint32_t dirty_kb = get_dirty_pages();
  const uint32_t dirty_threshold_kb = DIRTY_THRESHOLD_GB * 1000000;
  return (dirty_kb > dirty_threshold_kb);
}


void get_io_stat_info(int &r, int &w) {
  //FILE *pipe;
  //char str[250];
  //pipe = popen("iostat -m | grep md0 2>&1 | awk '{printf \"%d %d\\n\" , $5, $6}'", "r" );

  //fscanf(pipe, "%d %d", &r, &w);
  //pclose(pipe);
};



// template<typename mapped_t>
// void custom_flush(mapped_t * mapped) {

//   mapped->flush();

//   boost::interprocess::mapped_region::advice_types advise;
//   advise = boost::interprocess::mapped_region::advice_types::advice_dontneed;
//   assert(mapped->advise(advise));

//   advise = boost::interprocess::mapped_region::advice_types::advice_random;
//   assert(mapped->advise(advise));

// }



// template <typename SegmentManager>
// void
// delegate_partitioned_graph<SegmentManager>::
// try_flush(MPI_Comm comm) {
// #if 1
//   static uint64_t check_id = 0;

//   if ((m_mpi_rank  % 24) == (check_id++) % 24) {

//     uint32_t dirty_kb;
//     {
//       FILE *pipe;
//       pipe = popen("grep Dirty /proc/meminfo | awk '{print $2}'", "r" );
//       fscanf(pipe, "%u", &dirty_kb);
//       pclose(pipe);
//     }
//     const uint32_t dirty_threshold_kb = DIRTY_THRESHOLD_GB * 1000000;

//     if (dirty_kb > dirty_threshold_kb) {
//       m_flush_func();
//     }
//   }


// #else
//   bool do_flush;
//   if (m_mpi_rank == 0) {
//     uint32_t dirty_kb;
//     {
//       FILE *pipe;
//       pipe = popen("grep Dirty /proc/meminfo | awk '{print $2}'", "r" );
//       fscanf(pipe, "%u", &dirty_kb);
//       pclose(pipe);
//     }
//     const uint32_t dirty_threshold_kb = DIRTY_THRESHOLD_GB * 1000000;
//     do_flush = (dirty_kb > dirty_threshold_kb);
//   }

//   MPI_Bcast(&do_flush, 1, mpi_typeof(do_flush), 0, comm);


//   if (do_flush) {
//     m_flush_func();
//   }
// #endif
// }
#endif  // _HAVOQGT_CACHE_UTILI_HPP

