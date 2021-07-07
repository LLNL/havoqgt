// Copyright 2019 Lawrence Livermore National Security, LLC and other Clippy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <clippy/clippy.hpp>

int main(int argc, char **argv)
{
  clippy::clippy clip("grumpy", "Always throws errors because he's Grumpy!");
  if (clip.parse(argc, argv))
  {
    return 0;
  }

  throw std::runtime_error("I'm Grumpy!");
  return 0;
}
