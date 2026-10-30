//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// Portable implementation of the Gaussian overlap volumes and gradients.
//
#include <cmath>
#include <numbers>
#include <vector>

#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlap.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>

namespace RDKit {
namespace GaussianShape {

// atoms/shape features
double calcVolAndGrads(const double *ref, const double *refAlphas,
                       const int numRefPts,
                       const boost::dynamic_bitset<> *refCarbonRadii,
                       const double *fit, const double *fitAlphas,
                       const int numFitPts,
                       const boost::dynamic_bitset<> *fitCarbonRadii,
                       std::vector<double> &gradConverters,
                       const bool useCutoff, const double distCutoff2,
                       const double *quat, double *gradients) {
  if (gradients) {
    detail::cartToQuatGrads(quat, fit, numFitPts, gradConverters, 0);
  }
  static constexpr double CARBON_A = KAPPA / (1.7 * 1.7);
  static const double CARBON_BIT =
      8.0 * pow(std::numbers::pi / (2 * CARBON_A), 1.5);
  double vol = 0.0;
  double vij;
  // If either of the carbon radii flags aren't supplied, treat them
  // both as being all carbon.  There isn't enough information to do
  // otherwise.
  const bool allCarbon = !refCarbonRadii || !fitCarbonRadii;
  for (int i = 0, i_idx = 0; i < numRefPts * 3; i += 3, i_idx++) {
    const auto ai = refAlphas[i_idx];
    if (ai < 0.0) {
      continue;
    }
    for (int j = 0, j_idx = 0, k = 0; j < numFitPts * 3;
         j += 3, j_idx++, k += 12) {
      const auto dx = ref[i] - fit[j];
      const auto dy = ref[i + 1] - fit[j + 1];
      const auto dz = ref[i + 2] - fit[j + 2];
      const auto d2 = dx * dx + dy * dy + dz * dz;
      if (useCutoff && d2 > distCutoff2) {
        continue;
      }
      const auto aj = fitAlphas[j_idx];
      if (aj < 0.0) {
        continue;
      }
      const auto mult = -(ai * aj) / (ai + aj);
      const auto kij = exp(mult * d2);
      if (allCarbon || ((*refCarbonRadii)[i_idx] && (*fitCarbonRadii)[j_idx])) {
        vij = kij * CARBON_BIT;
      } else {
        const auto pi_ai_aj = std::numbers::pi / (ai + aj);
        vij = 8 * kij * pi_ai_aj * std::sqrt(pi_ai_aj);
      }
      vol += vij;
      if (gradients) {
        const auto r = 2.0 * vij * mult;
        // Use the gradient converters to calculate the gradients in quaternion
        // space.
        // The zeroth gradient is never used, so don't waste time calculating
        // it but leave the code here for completeness and possible future use.
        // gradients[0] +=
        //     r * (dx * gradConverters[k] + dy *
        //     gradConverters[k + 1] +
        //          dz * gradConverters[k + 2]);
        gradients[1] +=
            r * (dx * gradConverters[k + 3] + dy * gradConverters[k + 4] +
                 dz * gradConverters[k + 5]);
        gradients[2] +=
            r * (dx * gradConverters[k + 6] + dy * gradConverters[k + 7] +
                 dz * gradConverters[k + 8]);
        gradients[3] +=
            r * (dx * gradConverters[k + 9] + dy * gradConverters[k + 10] +
                 dz * gradConverters[k + 11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}

// color features
double calcVolAndGrads(const double *ref, const double *refAlphas,
                       const int numRefPts, const int *refTypes,
                       const double *fit, const double *fitAlphas,
                       const int numFitPts, const int *fitTypes,
                       const int numFitShape,
                       std::vector<double> &gradConverters,
                       const bool useCutoff, const double distCutoff2,
                       const double *quat, double *gradients) {
  double vol = 0.0;
  if (gradients) {
    detail::cartToQuatGrads(quat, fit, numFitPts, gradConverters, numFitShape);
  }

  for (int i = 0, i_idx = 0; i < numRefPts * 3; i += 3, i_idx++) {
    const auto ai = refAlphas[i_idx];
    const auto aType = refTypes[i_idx];
    for (int j = 0, j_idx = 0, k = 0; j < numFitPts * 3;
         j += 3, j_idx++, k += 12) {
      const auto bType = fitTypes[j_idx];
      if (aType != bType) {
        continue;
      }
      const auto dx = ref[i] - fit[j];
      const auto dy = ref[i + 1] - fit[j + 1];
      const auto dz = ref[i + 2] - fit[j + 2];
      const auto d2 = dx * dx + dy * dy + dz * dz;
      if (useCutoff && d2 > distCutoff2) {
        continue;
      }
      const auto aj = fitAlphas[j_idx];
      const auto mult = -(ai * aj) / (ai + aj);
      const auto kij = exp(mult * d2);

      const auto pi_ai_aj = std::numbers::pi / (ai + aj);
      const auto vij = 8 * kij * pi_ai_aj * std::sqrt(pi_ai_aj);
      vol += vij;
      if (gradients) {
        const auto r = 2.0 * vij * mult;
        // Use the converters to calculate the gradients in quaternion space.
        // The zeroth gradient is never used, so don't waste time calculating
        // it but leave the code here for completeness and possible future use.
        // gradients[0] +=
        //     r * (dx * gradConverters[k + 0] + dy *
        //     gradConverters[k + 1] +
        //          dz * gradConverters[k + 2]);
        gradients[1] +=
            r * (dx * gradConverters[k + 3] + dy * gradConverters[k + 4] +
                 dz * gradConverters[k + 5]);
        gradients[2] +=
            r * (dx * gradConverters[k + 6] + dy * gradConverters[k + 7] +
                 dz * gradConverters[k + 8]);
        gradients[3] +=
            r * (dx * gradConverters[k + 9] + dy * gradConverters[k + 10] +
                 dz * gradConverters[k + 11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}
}  // namespace GaussianShape
}  // namespace RDKit
