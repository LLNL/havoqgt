// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <clippy/clippy.hpp>

int main(int argc, char **argv)
{
  clippy::clippy clip("karen", "Produces stderr but exits normally.");
  if (clip.parse(argc, argv))
  {
    return 0;
  }

  std::cerr << "Can I speak to the manager!?!?!" << std::endl;
  return 0;
}
