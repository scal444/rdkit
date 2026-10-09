// Copyright (C) 2026 RDKit contributors
// This file is covered by the BSD license in license.txt.

#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilder.h>

DistGeom::TorsionCandidates installedBoundsMatrixTorsions() {
  RDKit::DGeomHelpers::PATH14_VECT paths(1);
  paths.front().value.type = RDKit::DGeomHelpers::TorsionType::CISTRANS;
  paths.front().type14 = RDKit::DGeomHelpers::Type14::IN_CHAIN;
  return paths.front().toTorsionRange();
}
