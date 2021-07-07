// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <tuple>
#include <clippy/clippy.hpp>
#include <typeinfo>

int main(int argc, char **argv) {
  clippy::clippy clip("return_tuple", "Always returns a tuple");

  if (clip.parse(argc, argv)) { return 0; }

  clip.to_return(std::make_tuple("foo", 42, 3.24));

  return 0;
}
