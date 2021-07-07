// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <clippy/clippy.hpp>

int main(int argc, char **argv) {
  clippy::clippy clip("sum", "Sums to numbers");
  clip.add_required<double>("i", "first Number");
  clip.add_required<double>("j", "second Number");

  clip.returns<double>("i + j");
  if (clip.parse(argc, argv)) { return 0; }

  auto i = clip.get<double>("i");
  auto j = clip.get<double>("j");

  clip.to_return(i + j);
  return 0;
}
