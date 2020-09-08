// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_DETAIL_NULL_OSTREAM_HPP_INCLUDED
#define HAVOQGT_DETAIL_NULL_OSTREAM_HPP_INCLUDED

#include <streambuf>
#include <ostream>

/// Null ostream, adapted from:
/// http://stackoverflow.com/questions/760301/implementing-a-no-op-stdostream

namespace havoqgt { namespace detail {
  
template <class cT, class traits = std::char_traits<cT> >
class basic_nullbuf: public std::basic_streambuf<cT, traits> {
    typename traits::int_type overflow(typename traits::int_type c)
    {
        return traits::not_eof(c); // indicate success
    }
};

template <class cT, class traits = std::char_traits<cT> >
class basic_onullstream: public std::basic_ostream<cT, traits> {
    public:
        basic_onullstream():
        std::basic_ios<cT, traits>(&m_sbuf),
        std::basic_ostream<cT, traits>(&m_sbuf)
        {
            this->init(&m_sbuf);
        }

    private:
        basic_nullbuf<cT, traits> m_sbuf;
};

inline std::ostream& get_null_ostream() {
  static basic_onullstream<char> null_ostream;
  return null_ostream;
}

//typedef basic_onullstream<char> onullstream;
//typedef basic_onullstream<wchar_t> wonullstream;  
  
}} //end namespace havoqgt::detail


#endif //HAVOQGT_DETAIL_NULL_OSTREAM_HPP_INCLUDED
