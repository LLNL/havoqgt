// Copyright 2020 Lawrence Livermore National Security, LLC and other CLIPPy Project Developers.
// See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: MIT

#include <clippy/clippy.hpp>
#include <vector>
#include <algorithm>

int main(int argc, char **argv) {
  clippy::clippy clip("sort_edges", "Sorts an edgelist");
  clip.add_required<clippy::arrayintint>("edges", "Unordered edgelist");
  clip.add_optional<clippy::boolean>("reverse", "Sort in reverse order", false);

  clip.returns<clippy::arrayintint>("Sorted edgelist");
  if (clip.parse(argc, argv)) { return 0; }

  auto edges = clip.get<clippy::arrayintint>("edges");
  bool reverse = clip.get<clippy::boolean>("reverse");

  if (reverse) {
    std::sort(edges.begin(), edges.end(),
              std::greater<decltype(edges)::value_type>{});
  } else {
    std::sort(edges.begin(), edges.end());
  }

  clip.to_return(edges);
  return 0;
}
