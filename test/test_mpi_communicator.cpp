// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include <havoqgt/mpi.hpp>

namespace havoqgt { namespace test {

TEST(my_test, test_b) {
  EXPECT_EQ(5, 5);
}

}} //end namespace havoqgt::test

//mpi main for gteset
GTEST_API_ int main(int argc, char **argv) {
  CHK_MPI( MPI_Init( &argc, &argv) );
  std::cout << "Running main() from gtest_main.cc\n";

  testing::InitGoogleTest(&argc, argv);
  int to_return RUN_ALL_TESTS();
  CHK_MPI( MPI_Finalize() );
  return to_return;
}
