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

#ifndef GRAPHSTORE_UTILITIES_HPP
#define GRAPHSTORE_UTILITIES_HPP

#include <cstring>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <havoqgt/utilities.hpp>

namespace graphstore {

   class direct_file_reader {
    direct_file_reader(const char* const fname, const size_t buf_size, const size_t align_size = 4096)
    : m_fd(-1)
    , m_buffer(nulptr)
    , m_buf_size(buf_size)
    , m_align_size(align_size)
    {
      if (buf_size < align_size) {

      }
      init(fname);
    }

    ~direct_file_reader()
    {
        free(m_buffer);
    }

    unsigned char &operator[](const size_t i)
    {
        return m_buffer[i+m_offset];
    }

    ssize_t align_read(size_t count)
    {
        /// Alignment the lenght of read
        const size_t aligned_read_size = std::min(m_buf_size, count) & !(m_align_size - 1);

        /// Alignment the file offset of I/Os
        const off_t cur_pos = lseek(fd, 0, SEEK_CUR);
        const size_t off = cur_pos & !(m_align_size - 1);
        lseek(m_fd, -off, SEEK_CUR);

        const ssize_t actual_read_size = read(m_fd, m_buffer, aligned_read_size);
        if (actual_read_size == -1) {
           perror("Reading file");
        }
        m_offset = off;
        const ssize_t read_size = actual_read_size - off;

        return read_size;
    }

   private:
    void init(char* fname)
    {
      const int flags = O_RDONLY | O_DIRECT;
      const int fd = open(m_fd, flags);
      if (fd == -1) {
        std::cerr << fname << std::endl;
        HAVOQGT_ERROR_MSG("Failed posix_memalign");
        return;
      }

      aligned_alloc(&m_buffer, m_align_size, m_buf_size);
    }

    int m_fd;
    size_t m_offset;
    void *m_buffer;
    size_t m_buf_size;
    size_t m_align_size;
  };


  int aligned_alloc(void** actual_buffer, size_t align_size, size_t length)
  {
    int result = posix_memalign(actual_buffer, align_size, length);
    if (result != 0) {
        HAVOQGT_ERROR_MSG("Failed posix_memalign");
    }
    return result;
  }

  void aligned_free(void* ptr)
  {
    free(ptr);
  }

}
#endif // GRAPHSTORE_UTILITIES_HPP

