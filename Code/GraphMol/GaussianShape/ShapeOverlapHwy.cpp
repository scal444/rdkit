//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Highway implementation of the Gaussian overlap volumes and gradients, used
// in place of ShapeOverlap.cpp when RDK_USE_HIGHWAY is enabled.
// foreach_target.h re-includes this file once per SIMD target and the best
// one supported by the CPU is chosen at runtime.
//
// The fit points are copied into structure-of-arrays scratch space and
// processed a vector at a time for each reference point; the pairs that
// ShapeOverlap.cpp skips are masked out. Results differ from ShapeOverlap.cpp
// in the last bits: the sums are accumulated in a different order and exp()
// is Highway's (max error 1 ULP).
//
#ifndef RDKIT_GAUSSIANSHAPE_SHAPEOVERLAPHWY_ONCE
#define RDKIT_GAUSSIANSHAPE_SHAPEOVERLAPHWY_ONCE

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlap.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>

namespace RDKit::GaussianShape::detail {

// The 9 gradient converters (out of 12 per point) that are used.
constexpr size_t NUM_CONVERTERS = 9;
constexpr size_t FIRST_CONVERTER = 3;

// The fit points in structure-of-arrays form.
struct FitPoints {
  size_t n = 0;
  const double *x = nullptr;
  const double *y = nullptr;
  const double *z = nullptr;
  const double *alpha = nullptr;
  // carbon flag (1 or 0) for shape points, type for color features
  const double *flag = nullptr;
  // NUM_CONVERTERS arrays of n values, or null without gradients
  const double *converters = nullptr;
};

struct RefPoints {
  const double *coords = nullptr;
  const double *alphas = nullptr;
  int n = 0;
  // shape: carbon flags (null if all carbon); color: types
  const boost::dynamic_bitset<> *carbons = nullptr;
  const int *types = nullptr;
};

struct OverlapParams {
  bool useCutoff = false;
  double distCutoff2 = 0.0;
  double carbonBit = 0.0;
  bool allCarbon = true;
};

}  // namespace RDKit::GaussianShape::detail

#endif  // RDKIT_GAUSSIANSHAPE_SHAPEOVERLAPHWY_ONCE

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "GraphMol/GaussianShape/ShapeOverlapHwy.cpp"
#include <hwy/foreach_target.h>  // must come before highway.h
#include <hwy/highway.h>
#include <hwy/contrib/math/math-inl.h>

HWY_BEFORE_NAMESPACE();
namespace RDKit::GaussianShape::detail::HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

using D = hn::ScalableTag<double>;
using V = hn::Vec<D>;

HWY_INLINE V load(D d, const double *HWY_RESTRICT p, size_t count) {
  return count == hn::Lanes(d) ? hn::LoadU(d, p) : hn::LoadN(d, p, count);
}

// Accumulates the overlap volume of every reference point with the fit
// points, and the gradients if gradients isn't null. Color features only
// overlap with features of the same type; shape points with negative alphas
// are skipped.
template <bool Color>
HWY_INLINE double overlap(const RefPoints &ref, const FitPoints &fit,
                          const OverlapParams &params, double *gradients) {
  const D d;
  const size_t N = hn::Lanes(d);
  const V zero = hn::Zero(d);
  const V two = hn::Set(d, 2.0);
  const V eight = hn::Set(d, 8.0);
  const V pi = hn::Set(d, std::numbers::pi);
  const V carbonBit = hn::Set(d, params.carbonBit);
  const V cutoff2 = hn::Set(d, params.distCutoff2);
  V vol = zero;
  std::array<V, 7> grads;
  grads.fill(zero);

  for (int i = 0; i < ref.n; ++i) {
    const double ai = ref.alphas[i];
    if (!Color && ai < 0.0) {
      continue;
    }
    const V vai = hn::Set(d, ai);
    const V rx = hn::Set(d, ref.coords[3 * i]);
    const V ry = hn::Set(d, ref.coords[3 * i + 1]);
    const V rz = hn::Set(d, ref.coords[3 * i + 2]);
    const V refFlag =
        hn::Set(d, Color ? ref.types[i]
                         : (params.allCarbon || (*ref.carbons)[i] ? 1.0 : 0.0));

    for (size_t j = 0; j < fit.n; j += N) {
      const size_t count = std::min(N, fit.n - j);
      const V dx = hn::Sub(rx, load(d, fit.x + j, count));
      const V dy = hn::Sub(ry, load(d, fit.y + j, count));
      const V dz = hn::Sub(rz, load(d, fit.z + j, count));
      const V d2 = hn::MulAdd(dz, dz, hn::MulAdd(dy, dy, hn::Mul(dx, dx)));
      const V aj = load(d, fit.alpha + j, count);
      const V fitFlag = load(d, fit.flag + j, count);

      auto valid = hn::FirstN(d, count);
      if (params.useCutoff) {
        valid = hn::And(valid, hn::Not(hn::Gt(d2, cutoff2)));
      }
      if constexpr (Color) {
        valid = hn::And(valid, hn::Eq(fitFlag, refFlag));
      } else {
        valid = hn::And(valid, hn::Not(hn::Lt(aj, zero)));
      }
      if (hn::AllFalse(d, valid)) {
        continue;
      }

      const V aiPlusAj = hn::Add(vai, aj);
      const V mult = hn::Neg(hn::Div(hn::Mul(vai, aj), aiPlusAj));
      const V kij = hn::Exp(d, hn::Mul(mult, d2));
      const V piAiAj = hn::Div(pi, aiPlusAj);
      V vij = hn::Mul(hn::Mul(hn::Mul(eight, kij), piAiAj), hn::Sqrt(piAiAj));
      if constexpr (!Color) {
        // both carbon: use the precomputed carbon volume
        const auto carbonPair =
            hn::And(hn::Ne(refFlag, zero), hn::Ne(fitFlag, zero));
        vij = hn::IfThenElse(carbonPair, hn::Mul(kij, carbonBit), vij);
      }
      vij = hn::IfThenElseZero(valid, vij);
      vol = hn::Add(vol, vij);

      if (gradients) {
        // masked again: mult may be infinite in the skipped lanes
        const V r = hn::IfThenElseZero(valid, hn::Mul(hn::Mul(two, vij), mult));
        auto converter = [&](size_t c) HWY_ATTR {
          return load(d, fit.converters + c * fit.n + j, count);
        };
        // gradient 0 is never used, so isn't calculated
        for (size_t g = 0; g < 3; ++g) {
          const V sum = hn::MulAdd(dz, converter(3 * g + 2),
                                   hn::MulAdd(dy, converter(3 * g + 1),
                                              hn::Mul(dx, converter(3 * g))));
          grads[g + 1] = hn::MulAdd(r, sum, grads[g + 1]);
        }
        grads[4] = hn::MulAdd(r, dx, grads[4]);
        grads[5] = hn::MulAdd(r, dy, grads[5]);
        grads[6] = hn::MulAdd(r, dz, grads[6]);
      }
    }
  }
  if (gradients) {
    for (size_t g = 1; g < grads.size(); ++g) {
      gradients[g] += hn::ReduceSum(d, grads[g]);
    }
  }
  return hn::ReduceSum(d, vol);
}

double ShapeOverlap(const RefPoints &ref, const FitPoints &fit,
                    const OverlapParams &params, double *gradients) {
  return overlap<false>(ref, fit, params, gradients);
}

double ColorOverlap(const RefPoints &ref, const FitPoints &fit,
                    const OverlapParams &params, double *gradients) {
  return overlap<true>(ref, fit, params, gradients);
}

}  // namespace RDKit::GaussianShape::detail::HWY_NAMESPACE
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace RDKit::GaussianShape {
namespace detail {
HWY_EXPORT(ShapeOverlap);
HWY_EXPORT(ColorOverlap);

namespace {
// Copies the fit points into structure-of-arrays scratch space, which is
// reused between calls. flag(j) gives the carbon flag or type of point j.
template <typename Flag>
FitPoints toFitPoints(const double *coords, const double *alphas, int n,
                      Flag flag, const std::vector<double> *gradConverters) {
  thread_local std::vector<double> scratch;
  const size_t numArrays = 5 + (gradConverters ? NUM_CONVERTERS : 0);
  scratch.resize(numArrays * n);
  FitPoints fit;
  fit.n = n;
  double *x = scratch.data(), *y = x + n, *z = y + n, *alpha = z + n,
         *flags = alpha + n;
  for (int j = 0; j < n; ++j) {
    x[j] = coords[3 * j];
    y[j] = coords[3 * j + 1];
    z[j] = coords[3 * j + 2];
    alpha[j] = alphas[j];
    flags[j] = flag(j);
  }
  fit.x = x;
  fit.y = y;
  fit.z = z;
  fit.alpha = alpha;
  fit.flag = flags;
  if (gradConverters) {
    double *converters = flags + n;
    for (size_t c = 0; c < NUM_CONVERTERS; ++c) {
      for (int j = 0; j < n; ++j) {
        converters[c * n + j] = (*gradConverters)[12 * j + FIRST_CONVERTER + c];
      }
    }
    fit.converters = converters;
  }
  return fit;
}
}  // namespace
}  // namespace detail

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
  // If either of the carbon radii flags aren't supplied, treat them
  // both as being all carbon.  There isn't enough information to do
  // otherwise.
  const bool allCarbon = !refCarbonRadii || !fitCarbonRadii;
  const auto fitPoints = detail::toFitPoints(
      fit, fitAlphas, numFitPts,
      [&](int j) { return allCarbon || (*fitCarbonRadii)[j] ? 1.0 : 0.0; },
      gradients ? &gradConverters : nullptr);
  const detail::RefPoints refPoints{ref, refAlphas, numRefPts, refCarbonRadii,
                                    nullptr};
  const detail::OverlapParams params{useCutoff, distCutoff2, CARBON_BIT,
                                     allCarbon};
  return HWY_DYNAMIC_DISPATCH(detail::ShapeOverlap)(refPoints, fitPoints,
                                                    params, gradients);
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
  if (gradients) {
    detail::cartToQuatGrads(quat, fit, numFitPts, gradConverters, numFitShape);
  }
  // As in ShapeOverlap.cpp, the converters are read from the start of
  // gradConverters, not from where they were just written.
  const auto fitPoints = detail::toFitPoints(
      fit, fitAlphas, numFitPts,
      [&](int j) { return static_cast<double>(fitTypes[j]); },
      gradients ? &gradConverters : nullptr);
  const detail::RefPoints refPoints{ref, refAlphas, numRefPts, nullptr,
                                    refTypes};
  const detail::OverlapParams params{useCutoff, distCutoff2, 0.0, false};
  return HWY_DYNAMIC_DISPATCH(detail::ColorOverlap)(refPoints, fitPoints,
                                                    params, gradients);
}

}  // namespace RDKit::GaussianShape
#endif  // HWY_ONCE
