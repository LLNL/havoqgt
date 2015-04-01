project (HAVOQGT)
cmake_minimum_required (VERSION 2.6)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

#
# Options
#
option(CMAKE_INSTALL_RPATH_USE_LINK_PATH "Use Rpath for dependencies." TRUE)
option(HAVOQGT_BUILD_TEST "Build HAVOQGT tests." TRUE)

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
#  Posix_fallocate
#
check_symbol_exists(posix_fallocate "fcntl.h" HAVE_POSIX_FALLOCATE)
if (HAVE_POSIX_FALLOCATE)
  add_definitions(-DHAVE_POSIX_FALLOCATE)
endif()

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
