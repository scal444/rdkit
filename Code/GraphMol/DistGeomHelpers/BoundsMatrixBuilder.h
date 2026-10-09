//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//  Copyright (C) 2026 ETH Zurich
//  Path14 configuration types created by Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_BOUNDS_MATRIX_BUILDER_H
#define RD_BOUNDS_MATRIX_BUILDER_H

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/ZMatrixUtils.h>
#include <cstddef>
#include <optional>
#include <vector>
#include "Embedder.h"

namespace RDKit {
class ROMol;
namespace DGeomHelpers {

enum class TorsionType {
  CIS = 0,
  TRANS,
  FLEXIBLE,
  CUSTOM,
  CISTRANS
};

enum class Type14 {
  IN_CHAIN,
  IN_RING,
  TWO_IN_SAME_RING,
  TWO_IN_DIFF_RING,
  SHARE_RING_BOND,
  MACROCYCLE_TWO_IN_SAME_RING,
  MACROCYCLE_ALL_IN_SAME_RING
};

struct TorsionValue {
  TorsionType type = TorsionType::FLEXIBLE;
  std::optional<double> value = {};
  std::optional<double> extraDist = {};
  bool isForced = false;
};

inline DistGeom::TorsionRange ringTorsion(const std::size_t rSize) {
  double torsion = M_PI;
  switch (rSize) {
    case 4u:
      [[fallthrough]];
    case 5u:
      torsion = M_PI * 45.0 / 180.0;
      break;
    case 6u:
      torsion = M_PI * 60.0 / 180.0;
      break;
    case 7u:
      torsion = M_PI * 90.0 / 180.0;
      break;
    case 8u:
      torsion = M_PI * 100.0 / 180.0;
      break;
  }
  return {-torsion, torsion};
}

//! A structure used to store planar 14 paths - cis/trans
struct Path14Configuration {
  unsigned int bid1, bid2, bid3;
  unsigned int aid1, aid2, aid3, aid4;
  TorsionValue value;
  Type14 type14;
  std::size_t rSize = 0;

  DistGeom::TorsionCandidates toTorsionRange() const {
    switch (value.type) {
      case TorsionType::CIS:
        return DistGeom::TorsionValues{0.0};
      case TorsionType::TRANS:
        return DistGeom::TorsionValues{M_PI};
      case TorsionType::FLEXIBLE:
        return ringTorsion(rSize);
      case TorsionType::CISTRANS:
        return DistGeom::TorsionValues{0.0, M_PI};
      case TorsionType::CUSTOM:
        return DistGeom::TorsionValues{*value.value};
      default:
        break;
    }
    return ringTorsion(0);
  }
};

using PATH14_VECT = std::vector<Path14Configuration>;

//! Set default upper and lower distance bounds in a distance matrix
/*!
  \param mmat        pointer to the bounds matrix to be altered
  \param defaultMin  default value for the lower distance bounds
  \param defaultMax  default value for the upper distance bounds
*/
RDKIT_DISTGEOMHELPERS_EXPORT void initBoundsMat(DistGeom::BoundsMatrix *mmat,
                                                double defaultMin = 0.0,
                                                double defaultMax = 1000.0);
/*! \overload
 */
RDKIT_DISTGEOMHELPERS_EXPORT void initBoundsMat(DistGeom::BoundsMatPtr mmat,
                                                double defaultMin = 0.0,
                                                double defaultMax = 1000.0);

RDKIT_DISTGEOMHELPERS_EXPORT void setTopolBounds(
    const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    const EmbedParameters &params, bool scaleVDW = false,
    bool set15bounds = true, bool set14bounds = true, bool set13bounds = true,
    PATH14_VECT *paths14 = nullptr,
    const EmbedFF embedForceField = EmbedFF::UFF,
    InternalCoordinates *internalCoords = nullptr);

//! Set upper and lower distance bounds between atoms in a molecule based on
/// topology
/*!
  This consists of setting 1-2, 1-3 and 1-4 distance based on bond lengths,
  bond angles and torsion angle ranges. Optionally 1-5 bounds can also be set,
  in particular, for path that contain rigid 1-4 paths.
  The final step involves setting lower bound to the sum of the vdW radii for
  the remaining atom pairs.
  \param mol          The molecule of interest
  \param mmat         Bounds matrix to the bounds are written
  \param set15bounds  If true try to set 1-5 bounds also based on topology
  \param scaleVDW     Ignored.
  \param useMacrocycle14config  If 1-4 distances bound heuristics for
  macrocycles is used <b>Note</b> For some strained systems the bounds matrix
  resulting from setting 1-5 bounds may fail triangle smoothing. In these cases
  it is recommended to back out and recompute the bounds matrix with no 1-5
  bounds and with vdW scaling.
  \param forceTransAmides  If true, amide bonds are enforced to be trans
  \param set14bounds  If true, set 1-4 distance bounds based on topology
  \param set13bounds  If true, set 1-3 distance bounds based on topology
*/
inline void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                           bool set15bounds = true, bool scaleVDW = false,
                           bool useMacrocycle14config = false,
                           bool forceTransAmides = true,
                           bool set14bounds = true, bool set13bounds = true,
                           const EmbedFF embedForceField = EmbedFF::UFF) {
  EmbedParameters params{.useMacrocycle14config = useMacrocycle14config,
                         .forceTransAmides = forceTransAmides};
  setTopolBounds(mol, mmat, params, scaleVDW, set15bounds, set14bounds,
                 set13bounds, nullptr, embedForceField);
}

/* ! \overload */
RDKIT_DISTGEOMHELPERS_EXPORT void setTopolBounds(
    const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    std::vector<std::pair<int, int>> &bonds,
    std::vector<std::vector<int>> &angles, const EmbedParameters &params,
    bool scaleVDW = false, bool set15bounds = true, bool set14bounds = true,
    bool set13bounds = true, PATH14_VECT *paths14 = nullptr,
    const EmbedFF embedForceField = EmbedFF::UFF,
    InternalCoordinates *internalCoords = nullptr);
/*! \overload for experimental torsion angle preferences
 */
inline void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                           std::vector<std::pair<int, int>> &bonds,
                           std::vector<std::vector<int>> &angles,
                           bool set15bounds = true, bool scaleVDW = false,
                           bool useMacrocycle14config = false,
                           bool forceTransAmides = true,
                           bool set14bounds = true, bool set13bounds = true,
                           const EmbedFF embedForceField = EmbedFF::UFF) {
  EmbedParameters params{.useMacrocycle14config = useMacrocycle14config,
                         .forceTransAmides = forceTransAmides};
  setTopolBounds(mol, mmat, bonds, angles, params, scaleVDW, set15bounds,
                 set14bounds, set13bounds, nullptr, embedForceField);
}

//! generate the vectors of bonds and angles used by (ET)KDG
RDKIT_DISTGEOMHELPERS_EXPORT void collectBondsAndAngles(
    const ROMol &mol, std::vector<std::pair<int, int>> &bonds,
    std::vector<std::vector<int>> &angles);

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
