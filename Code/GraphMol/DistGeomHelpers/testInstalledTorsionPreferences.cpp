// Copyright (C) 2026 RDKit contributors
// This file is covered by the BSD license in license.txt.

#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>

RDKit::DGeomHelpers::PATH14_VECT installedCrystalFFPaths() {
  ForceFields::CrystalFF::CrystalFFDetails details;
  details.path14Configs.emplace_back();
  return details.path14Configs;
}
