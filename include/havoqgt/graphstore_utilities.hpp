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

#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#if __linux__
#include <sys/sysinfo.h>
#endif

#define ENABLE_HAVOQGT_ERR_PROCEDURE 1
#if ENABLE_HAVOQGT_ERR_PROCEDURE
#include <havoqgt/utilities.hpp>
#endif

namespace graphstore {
namespace utility {


class aligned_alloc
{
 public:
    static int alloc(void** actual_buffer, size_t align_size, size_t length)
    {
        int result = ::posix_memalign(actual_buffer, align_size, length);
        if (result != 0) {
#if ENABLE_HAVOQGT_ERR_PROCEDURE
            HAVOQGT_ERROR_MSG("Failed posix_memalign");
#else
            std::cerr << "Failed posix_memalign";
            ::exit(1);
#endif
        }
        return result;
    }

    static void free(void* ptr)
    {
        ::free(ptr);
    }
};


template <size_t AlignSize>
class direct_file_reader
{
 public:
    explicit direct_file_reader()
        : m_fd(-1)
        , m_buffer(nullptr)
        , m_buf_size(0)
        , m_current_buf_pos(0)
        , m_end_buf_pos(0)
        , m_is_reached_eof(false)
    { }

    ~direct_file_reader()
    {
        close_file();
        free_buffer();
    }

    /// --- Disable copy and move assigment and constructor --- ///
    direct_file_reader& operator=(direct_file_reader const&) = delete;
    direct_file_reader& operator=(direct_file_reader const&&) = delete;
    direct_file_reader(direct_file_reader const&) = delete;
    direct_file_reader(direct_file_reader const&&) = delete;

    bool set_file(const char* const fname)
    {
        free_buffer();
        alloc_buffer(kInitBufSize);
        return open_file_with_direct_io_mode(fname);
    }


    bool getline(std::string& str)
    {

        while (1) {
            ssize_t line_start_pos;
            ssize_t line_length;
            const bool is_line_found = find_line_from_bufffer(line_start_pos, line_length);

            if (is_line_found) {
                copy_string_to_userscape_buffer(str, line_start_pos, line_length);
                return true;
            }

            if (!m_is_reached_eof) {
                const off_t length_remains = m_end_buf_pos - line_start_pos;
                seek(-length_remains, SEEK_CUR);
                aligned_read(kChunkSize);
            } else {
                if (m_end_buf_pos-line_start_pos > 0ULL) {
                    /// Now, the file is reachaed to EOF, Therefore if data is remaining in the buffer, copy it to the user space
                    copy_string_to_userscape_buffer(str, line_start_pos, m_end_buf_pos-line_start_pos);
                    return true;
                } else {
                    /// File has already reached to EOF and there is no remaining data in the buffer
                    return false;
                }
            }
        }

        /// File has already reached to EOF and there is no remaining data in the buffer
        return false;
    }


 private:
    enum {
        kChunkSize = (1ULL << 25),
        kInitBufSize = kChunkSize + 2ULL * AlignSize
    };


    void alloc_buffer(const size_t size)
    {
        aligned_alloc::alloc(&m_buffer, AlignSize, size);
        m_buf_size = size;
    }

    void free_buffer()
    {
        if (m_buffer != nullptr)
            aligned_alloc::free(m_buffer);
        m_buffer = nullptr;
        m_buf_size = 0;
        m_current_buf_pos = 0;
        m_end_buf_pos = 0;
    }

    bool open_file_with_direct_io_mode(const char* const fname)
    {
        if (m_fd != -1)
            close_file();

#ifdef O_DIRECT
        const int flags = O_RDONLY | O_DIRECT;
#else
        const int flags = O_RDONLY;
        std::cerr << "O_DIRECT is not suported\n";
#endif
        m_fd = ::open(fname, flags);
        if (m_fd == -1) {
            std::cerr << fname << std::endl;
            perror("Failed file open");
            return false;
        }

        return true;
    }

    void close_file()
    {
        m_fd = -1;
        ::close(m_fd);
    }

    off_t seek(off_t offset, int whence)
    {
        return ::lseek(m_fd, offset, whence);
    }

    /// Return: On success, the number of bytes read is returned
    ///         On error, -1 is returned
    ssize_t aligned_read(const size_t count)
    {
        /// Align the file offset of I/Os
        const off_t cur_pos = lseek(m_fd, 0, SEEK_CUR);
        const size_t off = cur_pos % AlignSize;
        seek(-off, SEEK_CUR);
        m_current_buf_pos = off;

        /// Align the lenght of read
        const size_t raw_read_size = count + off;
        const size_t aligned_read_size = calc_aligned_size(raw_read_size, AlignSize);

        /// If buffer is not allocated or current buffer is small, then allocate new
        if (m_buffer == nullptr || (m_buf_size <  aligned_read_size)) {
            free_buffer();
            alloc_buffer(aligned_read_size);
        }

        /// read data from the file using read(2)
        ssize_t read_size = ::read(m_fd, m_buffer, aligned_read_size);
        m_end_buf_pos = read_size;
        if (read_size < aligned_read_size) {
            m_is_reached_eof = true;
        }

        /// discarded unnecessary data
        if (read_size > 0) {
            read_size -= off;
        }

//        char* const buf = reinterpret_cast<char*>(m_buffer);
//        for (int i=0; i < m_end_buf_pos; ++i) {
//            std::cout << buf[i] << " ";
//        }
//        std::cout << std::endl;

        return read_size;
    }

    bool find_line_from_bufffer(ssize_t& line_start_pos, ssize_t& line_length)
    {
        skip_contiguous_charctor('\n');
        char* const buf = reinterpret_cast<char*>(m_buffer);
        line_start_pos = m_current_buf_pos;
        for (; m_current_buf_pos < m_end_buf_pos; ++m_current_buf_pos) {
            if (buf[m_current_buf_pos] == '\n') {
                line_length = m_current_buf_pos - line_start_pos;
                return true;
            }
        }
        return false;
    }

    /// Example:
    ///  input: x = 0101, align_size = 10
    ///  output: 0110
    ///  Note: this function dose not check if align_size is power of 2
    inline size_t calc_aligned_size(const size_t x, const size_t align_size)
    {
        const size_t y = (align_size - 1);
        return ((x + y) & (~y));
    }

    void copy_string_to_userscape_buffer(std::string& str, const size_t start_pos, const size_t length)
    {
        char* buf = reinterpret_cast<char*>(m_buffer);
        str.clear();
        str.append(&buf[start_pos], length);
    }

    inline void skip_contiguous_charctor(const char delim)
    {
        char* buf = reinterpret_cast<char*>(m_buffer);
        for (; m_current_buf_pos < m_end_buf_pos; ++m_current_buf_pos) {
            if (buf[m_current_buf_pos] != delim) break;
        }
    }

    int m_fd;
    void *m_buffer;
    size_t m_buf_size;
    size_t m_current_buf_pos;
    ssize_t m_end_buf_pos;
    bool m_is_reached_eof;
};


bool get_system_memory_usages( size_t* const mem_unit, size_t* const totalram, size_t* const freeram, size_t* const usedram,
                               size_t* const bufferram, size_t* const totalswap, size_t* const freeswap )
{
  bool is_succeed = false;
#if __linux__
  struct sysinfo info;
  if (sysinfo(&info) == 0) {
    *mem_unit = static_cast<size_t>(info.mem_unit);
    *totalram = static_cast<size_t>(info.totalram);
    *freeram = static_cast<size_t>(info.freeram);
    *usedram = static_cast<size_t>(info.totalram - info.freeram);
    *bufferram = static_cast<size_t>(info.bufferram);
    *totalswap = static_cast<size_t>(info.totalswap);
    *freeswap = static_cast<size_t>(info.freeswap);
    is_succeed = true;
  }
#endif
  return is_succeed;
}

bool get_my_memory_usages(size_t* const size, size_t* const resident, size_t* const share,
                          size_t* const text, size_t* const lib, size_t* const data, size_t* const dt)
{
  bool is_succeed = false;
#ifdef __linux__
  const char* statm_path = "/proc/self/statm";
  FILE *f = fopen(statm_path, "r");
  if (f) {
    if (7 == fscanf(f,"%ld %ld %ld %ld %ld %ld %ld", size, resident, share, text, lib, data, dt)) {
      is_succeed = true;
    }
  }
  fclose(f);
#endif
  return is_succeed;
}

} /// namespace utility
} /// namespace graphstore
#endif /// GRAPHSTORE_UTILITIES_HPP
