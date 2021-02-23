project (HAVOQGT)
cmake_minimum_required (VERSION 2.6)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

#
# Options
#
option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Use Rpath for dependencies." TRUE)
option(HAVOQGT_BUILD_TEST "Build HAVOQGT tests." TRUE)
set(METALL_ROOT "" CACHE PATH "Metall installed root directory")
option(USE_UMAP "Use Umap instead of the system mmap" OFF)
set(UMAP_ROOT "" CACHE PATH "UMap installed root directory")
option(USE_HDF5 "Build executable(s) that use HDF5" OFF)

#
# Main Inlucde directory
#
include_directories("${PROJECT_SOURCE_DIR}/include")

include(CheckSymbolExists)

#
# MPI is required
#
find_package( MPI REQUIRED)
macro(include_link_mpi target)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(${target} ${MPI_LIBRARIES})
endmacro()

#
# PThreads is required
#
find_package (Threads)
macro(include_link_threads target)
  target_link_libraries(${target} ${CMAKE_THREAD_LIBS_INIT})
endmacro()


#
# BOOST is required
#
find_package( Boost )
macro(include_link_boost target)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(${target} ${Boost_LIBRARIES})
endmacro() 

#
#  stdc++fs is required
#
macro(include_link_stdcppfs target)
  target_link_libraries(${target} stdc++fs)
endmacro()

#
#  Posix_fallocate
#
check_symbol_exists(posix_fallocate "fcntl.h" HAVE_POSIX_FALLOCATE)
if (HAVE_POSIX_FALLOCATE)
  add_definitions(-DHAVE_POSIX_FALLOCATE)
endif()

#
# Metall
#
if (METALL_ROOT)
  if (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.1)
      message(FATAL_ERROR "Metall requires GCC >= 8.1")
    endif ()
  endif ()

  if (Boost_VERSION VERSION_LESS 1.64)
    message(FATAL_ERROR "Metall requires Boost C++ Libraries >= 1.64")
  endif ()

  include_directories(${METALL_ROOT}/include)
endif()

#
# HDF5
#
if (USE_HDF5)
  find_package(HDF5 COMPONENTS CXX C HL) # TODO, add required version info
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    add_definitions(${HDF5_DEFINITIONS})
  endif()
endif()
macro(include_link_hdf5 target)
  target_link_libraries(${target} ${HDF5_LIBRARIES})
endmacro()

#
# Umap
#
if (USE_UMAP)
  if (UMAP_ROOT)
    find_library(LIBUMAP umap PATHS ${UMAP_ROOT}/lib REQUIRED)
    find_path(UMAP_INCLUDE_DIR umap/umap.h PATHS ${UMAP_ROOT}/include REQUIRED)
  else ()
    find_library(LIBUMAP umap REQUIRED)
    find_path(UMAP_INCLUDE_DIR umap/umap.h REQUIRED)
  endif ()
endif ()

macro(include_link_umap target)
  target_include_directories(${target} PUBLIC ${UMAP_INCLUDE_DIR})
  target_link_libraries(${target} ${LIBUMAP})
  target_compile_definitions(${target} PUBLIC METALL_USE_UMAP)
endmacro()

#
# Doxygen
#
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in 
               ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)


add_subdirectory (include)
add_subdirectory (src)
enable_testing()
add_subdirectory (test)

#
# Configuration Summary
#
message("== Configuration Summary ==")
message("   CMAKE_BUILD_TYPE:   ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE MATCHES None)
  message("        - CXX FLAGS:   ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  message("        - CXX FLAGS:   ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  message("        - CXX FLAGS:   ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  message("        - CXX FLAGS:   ${CMAKE_CXX_FLAGS_DEBUG}")
endif()
message("   HAVOQGT_BUILD_TEST:   ${HAVOQGT_BUILD_TEST}")
message("   Boost_FOUND:        ${Boost_FOUND}")
message("   MPI_FOUND:          ${MPI_FOUND}")
#message("   DOXYGEN_FOUND:      ${DOXYGEN_FOUND}")
#MPI params found
if(MPI_FOUND)
  message("   MPIEXEC: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS 
          ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
endif(MPI_FOUND)

if (METALL_ROOT)
  message("METALL_ROOT: ${METALL_ROOT}")
endif()
