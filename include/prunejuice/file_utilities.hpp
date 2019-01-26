#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace prunejuice { namespace utilities {

template <typename FileStream>
bool is_file_empty(FileStream& file_stream) {
  return file_stream.peek() == FileStream::traits_type::eof();
}   

}} // end namespace prunejuice::utilities
