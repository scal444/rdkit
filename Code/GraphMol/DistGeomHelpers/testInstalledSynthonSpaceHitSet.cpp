// Copyright (C) 2026 RDKit contributors
// This file is covered by the BSD license in license.txt.

#include <GraphMol/SynthonSpaceSearch/SynthonSpaceHitSet.h>

std::size_t installedSynthonHitCount(
    const RDKit::SynthonSpaceSearch::SynthonSpaceHitSet &hits) {
  return hits.numHits;
}
