//
// Copyright (C) 2015-2018 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDTHREADS_H_2015
#define RDTHREADS_H_2015

#include <RDGeneral/Invariant.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>

namespace RDKit {
inline unsigned int getNumThreadsToUse(int target) {
  if (target >= 1) {
    return static_cast<unsigned int>(target);
  }
  unsigned int res = std::thread::hardware_concurrency();
  const auto adjusted = static_cast<long long>(res) + target;
  return adjusted > 0 ? static_cast<unsigned int>(adjusted) : 1;
}
}  // namespace RDKit

#else

namespace RDKit {
inline unsigned int getNumThreadsToUse(int target) {
  RDUNUSED_PARAM(target);
  return 1;
}
}  // namespace RDKit
#endif

#endif
