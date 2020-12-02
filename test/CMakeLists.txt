#
# Directly include GTest here, and add include.
#
add_subdirectory (include)
include_directories(include)

#
# This function adds a mpi test without Google Test.
#
function ( add_mpi_ctest_no_gtest test_name )
  set(test_source "test_${test_name}.cpp")
  set(test_exe    "test_${test_name}")
  add_executable(${test_exe} ${test_source})
  include_link_mpi(${test_exe})
  include_link_boost(${test_exe})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_link_stdcppfs(${test_exe})
  SET( num_procs 2)
  WHILE( 8 GREATER ${num_procs} )
    add_test( "${test_name}_${num_procs}" ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
            ${num_procs} ${MPIEXEC_PREFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/${test_exe}")
    MATH( EXPR num_procs "${num_procs} * 2" ) #double num_procs
  ENDWHILE()
endfunction()

add_mpi_ctest_no_gtest( bfs )
add_mpi_ctest_no_gtest( bfs_delegate )
add_mpi_ctest_no_gtest( cc )
add_mpi_ctest_no_gtest( cc_delegate )

add_subdirectory(datasets)