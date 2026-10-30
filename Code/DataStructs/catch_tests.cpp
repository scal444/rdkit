//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
#include "BitVects.h"
#include "BitOps.h"
#include "BitVectUtils.h"
#include "ExplicitBitVect.h"
#include "SparseIntVect.h"
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

using namespace RDKit;

TEST_CASE("special cases for the limits of sparse vectors") {
  SECTION("SparseBitVect") {
    SparseBitVect sbv(std::numeric_limits<unsigned int>::max());
    CHECK(sbv.getNumBits() == std::numeric_limits<unsigned int>::max());
    CHECK(!sbv.setBit(std::numeric_limits<unsigned int>::max()));
    CHECK(sbv.getBit(std::numeric_limits<unsigned int>::max()) == 1);
  }
}

TEST_CASE("github #9033: tversky is 1 when no bits are set") {
  ExplicitBitVect bv1(8);
  ExplicitBitVect bv2(8);
  CHECK(TverskySimilarity(bv1, bv2, 0.5, 0.5) == 0.0);
  CHECK(TanimotoSimilarity(bv1, bv2) == 0.0);
  CHECK(CosineSimilarity(bv1, bv2) == 0.0);
  CHECK(KulczynskiSimilarity(bv1, bv2) == 0.0);
  CHECK(SokalSimilarity(bv1, bv2) == 0.0);
  CHECK(McConnaugheySimilarity(bv1, bv2) == 0.0);
  CHECK(BraunBlanquetSimilarity(bv1, bv2) == 0.0);
  CHECK(RusselSimilarity(bv1, bv2) == 0.0);
  CHECK(RogotGoldbergSimilarity(bv1, bv2) == 0.0);
}
TEST_CASE("bitmap functions on unaligned bitmaps") {
  // The bitmaps needn't be aligned for the word-at-a-time popcount (run under
  // UBSan to check the loads themselves).
  std::vector<unsigned char> abuf(70), bbuf(70);
  for (unsigned int i = 0; i < abuf.size(); ++i) {
    abuf[i] = static_cast<unsigned char>(37 * i + 11);
    bbuf[i] = static_cast<unsigned char>(91 * i + 5);
  }
  for (unsigned int offset = 1; offset < 4; ++offset) {
    const unsigned char *a = abuf.data() + offset;
    const unsigned char *b = bbuf.data() + offset;
    // aligned copies of the same bitmaps
    std::vector<std::uint64_t> alignedA(8), alignedB(8);
    std::memcpy(alignedA.data(), a, 64);
    std::memcpy(alignedB.data(), b, 64);
    const auto *ca = reinterpret_cast<const unsigned char *>(alignedA.data());
    const auto *cb = reinterpret_cast<const unsigned char *>(alignedB.data());
    CHECK(CalcBitmapPopcount(a, 64) == CalcBitmapPopcount(ca, 64));
    CHECK(CalcBitmapTanimoto(a, b, 64) == CalcBitmapTanimoto(ca, cb, 64));
    CHECK(CalcBitmapDice(a, b, 64) == CalcBitmapDice(ca, cb, 64));
    CHECK(CalcBitmapTversky(a, b, 64, 0.3, 0.7) ==
          CalcBitmapTversky(ca, cb, 64, 0.3, 0.7));
    CHECK(CalcBitmapAllProbeBitsMatch(a, a, 64));
  }
}
