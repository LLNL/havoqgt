# Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy
# Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: MIT

#
# This function adds an mpi example.
#
function ( add_example example_name )
  set(example_source "${example_name}.cpp")
  set(example_exe    "${example_name}")
  add_executable(${example_exe} ${example_source})
  target_include_directories(${example_exe} PRIVATE ${METALL_INCLUDE_DIR})
  target_include_directories(${example_exe} PRIVATE ${NLOHMANN_JSON_INCLUDE_DIR})
  add_dependencies(${example_exe} metall)
  add_dependencies(${example_exe} nlohmann_json)
  if(UNIX AND NOT APPLE)
    target_link_libraries(${example_exe} rt)
  endif()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

add_example(grumpy)
add_example(sort_edges)
add_example(sort_strings)
add_example(howdy)
add_example(sum)
add_example(karen)
add_example(return_tuple)
add_example(st_path)
add_example(st_shortest_path)
