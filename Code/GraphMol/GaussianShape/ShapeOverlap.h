//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Internal helpers shared by the implementations of calcVolAndGrads() in
// ShapeOverlap.cpp and ShapeOverlapHwy.cpp.
//
#ifndef RDKIT_GAUSSIANSHAPE_SHAPEOVERLAP_H
#define RDKIT_GAUSSIANSHAPE_SHAPEOVERLAP_H

#include <vector>

namespace RDKit::GaussianShape::detail {

// Set of values to convert the cartesian gradients to quaternion gradients,
// 12 per point of mol, written from gradConverters[12 * gradConvOffset].
void cartToQuatGrads(const double *quat, const double *mol, int numBPts,
                     std::vector<double> &gradConverters, int gradConvOffset);

}  // namespace RDKit::GaussianShape::detail

#endif
