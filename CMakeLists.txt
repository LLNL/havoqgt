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

#
# MPI is required
#
find_package( MPI REQUIRED)
macro(include_link_mpi target)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(${target} ${MPI_LIBRARIES})
endmacro()

#
# BOOST is required
#
find_package( Boost )
find_package (Threads)
macro(include_link_boost target)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(${target} ${Boost_LIBRARIES})
  target_link_libraries (${target} pthread)
  target_link_libraries (${target} rt)
endmacro() 

add_subdirectory (include)
add_subdirectory (src)

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