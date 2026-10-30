//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Test helper for code with SIMD kernels. Tests using it need
// rdkit_link_highway(<test>) in CMake.
//
#ifndef RD_SIMDTESTUTILS_H
#define RD_SIMDTESTUTILS_H

#include <catch2/catch_all.hpp>

#ifdef RDK_USE_HIGHWAY
#include <hwy/targets.h>
#endif

namespace RDKit {

//! Runs check() once for every SIMD target that Highway compiled and the CPU
//! supports, or just once when Highway isn't in use.
template <typename Check>
void forEachSimdTarget(Check &&check) {
#ifdef RDK_USE_HIGHWAY
  for (int64_t target : hwy::SupportedAndGeneratedTargets()) {
    INFO("Highway target " << hwy::TargetName(target));
    hwy::SetSupportedTargetsForTest(target);
    check();
  }
  hwy::SetSupportedTargetsForTest(0);
#else
  check();
#endif
}

}  // namespace RDKit

#endif
