// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <clippy/clippy.hpp>

int main(int argc, char **argv) {
  clippy::clippy clip("howdy", "Formal Texan greeting.");
  clip.add_required<clippy::string>("name", "Name to greet");
  clip.returns<clippy::string>("The greeting");
  if (clip.parse(argc, argv)) { return 0; }

  auto name = clip.get<clippy::string>("name");

  clip.to_return(std::string("Howdy, ") + name);
  return 0;
}
