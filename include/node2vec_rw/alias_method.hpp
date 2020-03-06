/*
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Roger Pearce <rpearce@llnl.gov>.
 * LLNL-CODE-644630.
 * All rights reserved.
 *
 * This file is part of HavoqGT, Version 0.1.
 * For details, see
 * https://computation.llnl.gov/casc/dcca-pub/dcca/Downloads.html
 *
 * Please also read this link – Our Notice and GNU Lesser General Public
 * License.
 *   http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR
 * A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC LICENSE
 *
 * Our Preamble Notice
 *
 * A. This notice is required to be provided under our contract with the
 * U.S. Department of Energy (DOE). This work was produced at the Lawrence
 * Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the
 * DOE.
 *
 * B. Neither the United States Government nor Lawrence Livermore National
 * Security, LLC nor any of their employees, makes any warranty, express or
 * implied, or assumes any liability or responsibility for the accuracy,
 * completeness, or usefulness of any information, apparatus, product, or
 * process
 * disclosed, or represents that its use would not infringe privately-owned
 * rights.
 *
 * C. Also, reference herein to any specific commercial products, process, or
 * services by trade name, trademark, manufacturer or otherwise does not
 * necessarily constitute or imply its endorsement, recommendation, or favoring
 * by
 * the United States Government or Lawrence Livermore National Security, LLC.
 * The
 * views and opinions of authors expressed herein do not necessarily state or
 * reflect those of the United States Government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product endorsement
 * purposes.
 *
 */

#ifndef HAVOQGT_NODE2VEC_RW_ALIAS_METHOD_HPP
#define HAVOQGT_NODE2VEC_RW_ALIAS_METHOD_HPP

#include <random>
#include <tuple>
#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>

namespace node2vec_rw {

/// \brief An alias method class which works like std::discrete_distribution
template <typename index_type = uint64_t, typename prob_type = double>
class alias_method {

 public:
  alias_method() = default;

  template <typename output_iterator_type>
  alias_method(const output_iterator_type first, const output_iterator_type last)
      :m_alias_table() {
    reset_table(first, last);
  }

  /// Construct a new biased distribution, reusing memory
  template <typename output_iterator_type>
  void reset_table(const output_iterator_type first, const output_iterator_type last) {

    m_alias_table.clear();

    if (first == last) return;

    const double average = static_cast<double>(std::accumulate(first, last, 0.0f)) / std::distance(first, last);

    // Move to the first under element
    auto under_itr = first;
    for (; under_itr != last; ++under_itr) {
      if (*under_itr < k_epsilon) continue; // Ignore too small values
      if (nearly_less_equal(*under_itr, average)) break;
    }

    // Move to the first over element
    auto over_itr = first;
    for (; over_itr != last; ++over_itr) {
      if (more_than(*over_itr, average)) break;
    }

    while (under_itr != last && over_itr != last) {
      const index_type under_index = std::distance(first, under_itr);
      assert((index_type)under_index < (index_type)std::distance(first, last));
      const index_type over_index = std::distance(first, over_itr);
      assert((index_type)over_index < (index_type)std::distance(first, last));

      if (nearly_equal(*under_itr, average)) { // The under element is "equal" to average
        m_alias_table.emplace_back(under_index, under_index, 1.0);
      } else {
        m_alias_table.emplace_back(under_index, over_index, *under_itr / average);
      }

      *over_itr -= std::max((double)(average - *under_itr),
                            (double)0.0f); // Get not enough weights from the over element

      *under_itr = 0; // Set to 0 so that we don't use the value again

      // the over element has just become an under element
      if (nearly_less_equal(*over_itr, average)) {

        if (over_index < under_index) {
          under_itr = over_itr;
        }

        // Move to the next over element
        for (++over_itr; over_itr != last; ++over_itr) {
          if (more_than(*over_itr, average)) break;
        }
      }

      // Move to the next under element, if needed
      for (; under_itr != last; ++under_itr) {
        if (*under_itr < k_epsilon) continue; // Ignore too small or already used values
        if (nearly_less_equal(*under_itr, average)) break;
      }
    }

    // Put unused elements into the m_alias_table
    auto unsed_itr = (under_itr != last) ? under_itr : over_itr;
    for (; unsed_itr != last; ++unsed_itr) {
      if (*unsed_itr < k_epsilon) continue; // Ignore too small or already used values
      const auto index = std::distance(first, unsed_itr);
      assert(index < std::distance(first, last));
      // Note that *unsed_itr could be somewhat different from 1.0 due to the floating operations
      m_alias_table.emplace_back(index, index, 1.0);
    }

    assert((std::size_t)m_alias_table.size() <= (std::size_t)std::distance(first, last));
  }

  template <typename generator>
  index_type operator()(generator &g) const {
    assert(m_alias_table.size() > 0);

    std::uniform_int_distribution<uint64_t> dis1(0, m_alias_table.size() - 1);
    const auto index = dis1(g);

    std::uniform_real_distribution<> dis2(k_epsilon, 1.0);
    const auto under_probability = dis2(g);
    if (under_probability <= m_alias_table[index].threshold) {
      return m_alias_table[index].under_index;
    } else {
      return m_alias_table[index].over_index;
    }

    assert(false);
  }

 private:
  static constexpr double k_epsilon = std::numeric_limits<prob_type>::epsilon();

  bool nearly_equal(const double a, const double b) {
    return (std::fabs(a - b) < k_epsilon);
  }

  bool nearly_less_equal(const double a, const double b) {
    return (a < b || nearly_equal(a, b));
  }

  bool more_than(const double a, const double b) {
    return !(nearly_less_equal(a, b));
  }

  struct alias_table_element {
    alias_table_element(index_type _under_index, index_type _over_index, prob_type _threshold)
        : under_index(_under_index),
          over_index(_over_index),
          threshold(_threshold) {}
    index_type under_index;
    index_type over_index;
    prob_type threshold; // Weights for selecting the under element
  };

  std::vector<alias_table_element> m_alias_table;
};

}

#endif //HAVOQGT_NODE2VEC_RW_ALIAS_METHOD_HPP
