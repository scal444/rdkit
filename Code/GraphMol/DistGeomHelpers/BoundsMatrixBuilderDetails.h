//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_BOUNDS_MATRIX_BUILDER_DETAILS_H
#define RD_BOUNDS_MATRIX_BUILDER_DETAILS_H

#include "BoundsMatrixBuilder.h"
#include <RDGeneral/Invariant.h>
#include <algorithm>
#include <array>
#include <optional>
#include <ostream>
#include <ranges>
#include <vector>

namespace RDKit {
namespace DGeomHelpers {
inline std::size_t getUnifiedId(const unsigned int id1, const unsigned int id2,
                                const unsigned int n) {
  // returns an id for (id1, id2) independent of order within range (0, 2*n - 1)
  // assuming id1 < n and id2 < n
  return id1 < id2 ? (static_cast<std::size_t>(id1) * n + id2)
                   : (static_cast<std::size_t>(id2) * n + id1);
}

inline std::size_t getUnifiedId(const unsigned int id1, const unsigned int id2,
                                const unsigned int id3, const unsigned int n) {
  // returns an id for (id1, id2, id3) independent of order of id1, id3 within
  // range (0, 3*(n) - 1) assuming id1 < n, id2 < n and id3 < n
  return id1 < id3 ? (static_cast<std::size_t>(id1) * n * n + id2 * n + id3)
                   : (static_cast<std::size_t>(id3) * n * n + id2 * n + id1);
}

template <unsigned int numBondIds>
auto unifiedIdToBondIds(std::size_t id, const unsigned int n) {
  // converts an id (created with getUnifiedId) back into the original bond ids
  std::array<unsigned int, numBondIds> bondIds;

  for (auto i : std::views::iota(0u, numBondIds) | std::views::reverse) {
    bondIds[i] = static_cast<unsigned int>(id % n);
    id /= n;
  }

  return bondIds;
};

struct Bounds {
  double lower{1.0}, upper{-1.0};  // we start invalid
  unsigned int aid1{0}, aid4{0};

  inline bool valid() const { return lower <= upper; }

  bool operator==(const Bounds &) const = default;

  friend std::ostream &operator<<(std::ostream &os, const Bounds &b) {
    return os << "Bounds{"
              << "lower=" << b.lower << ", upper=" << b.upper
              << ", aid1=" << b.aid1 << ", aid4=" << b.aid4 << '}';
  }
};

inline Bounds merge(std::vector<Bounds> bounds) {
  PRECONDITION(bounds.size(), "Cannot merge empty list of bounds");

  std::ranges::sort(bounds, {}, &Bounds::lower);

  Bounds current = bounds.front();
  double componentUpper = current.upper;
  std::optional<double> resultLower;

  // What we are doing here:
  // U {i'=intersection(i_j,..,i_k) | {i_j, ..., i_k}\subset(I) ^ i` !=
  // \emptyset ^ !\exists(i_l): intersection(i`, i_l) != \emptyset}
  // or in other words:
  // we aim to find the union of all intersections that are maximal in a sense
  // that adding another arbitrary bounds to it, would lead into an empty set

  // we solve this by traversing the sorted bounds in a sweep manner while
  // keeping track on the current/active non-empty intersection
  // (currentIntersection), the largest upperBound that was reached so far
  // (this is needed since the currentIntersection.upper can be smaller than
  // that, losing track of potenial overlaps/intersections).
  // To avoid storing all maximal non-overlapping intersections (only the
  // first and last one is relevant), we store the lower bound of the first
  // maximal intersection in resultLower

  for (const auto &_bound : bounds | std::views::drop(1)) {
    if (_bound.lower <= current.upper) {
      // Case 1: _bounds intersects with currentIntersection => add to current
      // intersection
      //  we know that bounds are sorted by lower bounds =>
      // _bound.lower is always greater/equal currentIntersection.lower
      current.lower = _bound.lower;
      current.upper = std::min(current.upper, _bound.upper);
    } else {
      // Case 2: _bound is not overlapping with the current intersection => we
      // know that currentIntersection is maximal

      if (!resultLower) {
        resultLower = current.lower;
      }

      current.lower = _bound.lower;
      current.upper =
          _bound.lower <= componentUpper
              ? std::min(componentUpper,
                         _bound.upper)  // there is this at least former
                                        // bounds that is overlapping and
                                        // needs to be considered
              : _bound.upper;
    }

    componentUpper = std::max(componentUpper, _bound.upper);
  }
  return Bounds{.lower = resultLower.value_or(current.lower),
                .upper = current.upper,
                .aid1 = bounds.front().aid1,
                .aid4 = bounds.front().aid4};
}

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
