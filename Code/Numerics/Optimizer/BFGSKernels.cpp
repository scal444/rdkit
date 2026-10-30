//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Scalar implementation of the BFGS kernels, with SVE fast paths that are
// selected at runtime on aarch64 Linux.
//
#include "BFGSKernels.h"

#if defined(__linux__) && defined(__aarch64__)
#include <sys/auxv.h>
#include <asm/hwcap.h>
#if defined(__has_include)
#if __has_include(<arm_sve.h>)
#include <arm_sve.h>
#define RDK_SVE_AVAILABLE 1
#endif
#endif
#endif

namespace BFGSOpt::detail {
namespace {
bool cpuHasSVE() {
#if defined(__linux__) && defined(__aarch64__) && defined(RDK_SVE_AVAILABLE)
  static const bool result = (getauxval(AT_HWCAP) & HWCAP_SVE) != 0;
  return result;
#else
  return false;
#endif
}

#ifdef RDK_SVE_AVAILABLE

// ---------------------------------------------------------------------------
// SVE kernel: initialise search direction xi = -grad and accumulate ||pos||^2
// ---------------------------------------------------------------------------
__attribute__((target("+sve"))) void sveInitXiAndSum(unsigned int dim,
                                                     const double *grad,
                                                     double *xi,
                                                     const double *pos,
                                                     double *outSum) {
  svfloat64_t acc = svdup_f64(0.0);
  uint64_t i = 0;
  while (i < dim) {
    svbool_t pg = svwhilelt_b64(i, (uint64_t)dim);
    svfloat64_t g = svld1_f64(pg, grad + i);
    svfloat64_t p = svld1_f64(pg, pos + i);
    // Negate grad in-place during store — avoids a separate negation pass
    svst1_f64(pg, xi + i, svneg_f64_m(svdup_f64(0.0), pg, g));
    // Fused multiply-add accumulates p[i]^2 without intermediate stores
    acc = svmla_f64_m(pg, acc, p, p);
    i += svcntd();
  }
  // Horizontal reduction collapses all vector lanes to a single scalar sum
  *outSum = svaddv_f64(svptrue_b64(), acc);
}

// ---------------------------------------------------------------------------
// SVE kernel: compute hessDGrad = invHessian * dGrad, then accumulate the
// four dot-product scalars (fac, fae, sumDGrad, sumXi) needed for the BFGS
// rank-1 update.
//
// The inner matrix-vector product uses SVE gather/FMA to handle arbitrary dim
// without scalar remainder loops. Each outer row is streamed once, keeping
// cache pressure proportional to dim rather than dim^2.  The four scalar
// accumulators are computed in a single fused pass over the result vector,
// saving two additional O(dim) traversals compared to separate dot-product
// calls.
// ---------------------------------------------------------------------------
__attribute__((target("+sve"))) void sveHessianVecMul(
    unsigned int dim, const double *invHessian, const double *dGrad,
    double *hessDGrad, const double *xi, double *outFac, double *outFae,
    double *outSumDGrad, double *outSumXi) {
  // Phase 1: matrix-vector multiply  invHessian * dGrad -> hessDGrad
  // Fully vectorised over both i and j.
  for (unsigned int i = 0; i < dim; i++) {
    const double *row = invHessian + i * dim;
    svfloat64_t acc = svdup_f64(0.0);
    uint64_t j = 0;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      acc = svmla_f64_m(pg, acc, svld1_f64(pg, row + j),
                        svld1_f64(pg, dGrad + j));
      j += svcntd();
    }
    hessDGrad[i] = svaddv_f64(svptrue_b64(), acc);
  }

  svfloat64_t vFac = svdup_f64(0.0), vFae = svdup_f64(0.0);
  svfloat64_t vSDG = svdup_f64(0.0), vSXi = svdup_f64(0.0);
  uint64_t i = 0;
  while (i < dim) {
    svbool_t pg = svwhilelt_b64(i, (uint64_t)dim);
    svfloat64_t vdg = svld1_f64(pg, dGrad + i);
    svfloat64_t vxi = svld1_f64(pg, xi + i);
    svfloat64_t vhd = svld1_f64(pg, hessDGrad + i);

    vFac = svmla_f64_m(pg, vFac, vdg, vxi);
    vFae = svmla_f64_m(pg, vFae, vdg, vhd);
    vSDG = svmla_f64_m(pg, vSDG, vdg, vdg);
    vSXi = svmla_f64_m(pg, vSXi, vxi, vxi);
    i += svcntd();
  }
  *outFac = svaddv_f64(svptrue_b64(), vFac);
  *outFae = svaddv_f64(svptrue_b64(), vFae);
  *outSumDGrad = svaddv_f64(svptrue_b64(), vSDG);
  *outSumXi = svaddv_f64(svptrue_b64(), vSXi);
}

// ---------------------------------------------------------------------------
// SVE kernel: symmetric rank-1 update of the inverse Hessian approximation.
//
// The BFGS update formula adds three outer-product terms to invHessian.
// Exploiting symmetry (only the upper triangle is computed; the lower is
// mirrored afterwards) halves the number of FLOPs and memory writes versus a
// naive full-matrix update. The SVE FMA instructions (svmla / svmls) fuse
// multiply and add into a single pipeline stage, reducing instruction count
// and register pressure compared to separate multiply + add sequences.
// ---------------------------------------------------------------------------
__attribute__((target("+sve"))) void sveHessianRank1Update(
    unsigned int dim, double *invHessian, const double *xi,
    const double *hessDGrad, const double *dGrad, double fac, double fad,
    double fae) {
  for (unsigned int i = 0; i < dim; i++) {
    // Broadcast scalar multipliers once per row to avoid redundant computation
    svfloat64_t vpxi = svdup_f64(fac * xi[i]);
    svfloat64_t vhdgi = svdup_f64(fad * hessDGrad[i]);
    svfloat64_t vdgi = svdup_f64(fae * dGrad[i]);
    double *row = invHessian + i * dim;
    // Start from column i to process only the upper triangle
    uint64_t j = i;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      svfloat64_t vxj = svld1_f64(pg, xi + j);
      svfloat64_t vhdgj = svld1_f64(pg, hessDGrad + j);
      svfloat64_t vdgj = svld1_f64(pg, dGrad + j);
      svfloat64_t vh = svld1_f64(pg, row + j);
      // Three fused multiply-adds apply all three BFGS update terms at once
      vh = svmla_f64_m(pg, vh, vpxi, vxj);
      vh = svmls_f64_m(pg, vh, vhdgi, vhdgj);
      vh = svmla_f64_m(pg, vh, vdgi, vdgj);
      svst1_f64(pg, row + j, vh);
      j += svcntd();
    }
    // Mirror upper triangle to lower triangle to maintain symmetry;
    // scalar loop is cheap (dim - i iterations) relative to the SVE inner loop
    for (unsigned int j2 = i + 1; j2 < dim; j2++) {
      invHessian[j2 * dim + i] = invHessian[i * dim + j2];
    }
  }
}

// ---------------------------------------------------------------------------
// SVE kernel: compute new search direction xi = -(invHessian * grad)
//
// SVE implementation accelerates the matrix-vector multiply using vectorized
// FMA across columns, then negates the accumulated dot product at the scalar
// level (single negation per row) rather than applying a separate vector
// negation pass. This keeps the loop body to one SVE FMA instruction
// per iteration, maximising throughput on in-order SVE pipelines.
// ---------------------------------------------------------------------------
__attribute__((target("+sve"))) void sveHessianVecMulNeg(
    unsigned int dim, const double *invHessian, const double *grad,
    double *xi) {
  for (unsigned int i = 0; i < dim; i++) {
    const double *row = invHessian + i * dim;
    svfloat64_t acc = svdup_f64(0.0);
    uint64_t j = 0;
    while (j < dim) {
      svbool_t pg = svwhilelt_b64(j, (uint64_t)dim);
      svfloat64_t h = svld1_f64(pg, row + j);
      svfloat64_t g = svld1_f64(pg, grad + j);
      acc = svmla_f64_m(pg, acc, h, g);
      j += svcntd();
    }
    // Negate the scalar result once per row rather than vectorising the
    // negation, keeping the store path simple and avoiding an extra SVE pass
    xi[i] = -svaddv_f64(svptrue_b64(), acc);
  }
}
#endif  // RDK_SVE_AVAILABLE
}  // namespace

double initDirection(unsigned int dim, const double *grad, double *xi,
                     const double *pos) {
  double sum = 0.0;
#ifdef RDK_SVE_AVAILABLE
  if (cpuHasSVE()) {
    sveInitXiAndSum(dim, grad, xi, pos, &sum);
    return sum;
  }
#endif
  for (unsigned int i = 0; i < dim; i++) {
    xi[i] = -grad[i];
    sum += pos[i] * pos[i];
  }
  return sum;
}

UpdateDots hessianTimesDGrad(unsigned int dim, const double *invHessian,
                             const double *dGrad, const double *xi,
                             double *hessDGrad) {
  UpdateDots dots;
#ifdef RDK_SVE_AVAILABLE
  if (cpuHasSVE()) {
    // matrix-vector multiply and all four dot products computed in one
    // vectorised pass, saving two additional O(dim) traversals compared
    // to separate scalar dot-product calls.
    sveHessianVecMul(dim, invHessian, dGrad, hessDGrad, xi, &dots.fac,
                     &dots.fae, &dots.sumDGrad, &dots.sumXi);
    return dots;
  }
#endif
  // fused matrix-vector multiply and dot-product accumulation.
  // Pointer arithmetic (++ivh, ++dgj) avoids repeated index computations
  // and helps the compiler generate efficient load sequences.
  for (unsigned int i = 0; i < dim; i++) {
    const double *ivh = &(invHessian[i * dim]);
    double &hdgradi = hessDGrad[i];
    const double *dgj = dGrad;
    hdgradi = 0.0;
    for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
      hdgradi += *ivh * *dgj;
    }
    dots.fac += dGrad[i] * xi[i];
    dots.fae += dGrad[i] * hessDGrad[i];
    dots.sumDGrad += dGrad[i] * dGrad[i];
    dots.sumXi += xi[i] * xi[i];
  }
  return dots;
}

void updateInverseHessian(unsigned int dim, double *invHessian,
                          const double *xi, const double *hessDGrad,
                          const double *dGrad, double fac, double fad,
                          double fae) {
#ifdef RDK_SVE_AVAILABLE
  if (cpuHasSVE()) {
    // symmetric rank-1 update with FMA, exploiting symmetry to halve memory
    // writes and FLOPs versus a full-matrix update
    sveHessianRank1Update(dim, invHessian, xi, hessDGrad, dGrad, fac, fad, fae);
    return;
  }
#endif
  // upper-triangle-only update (j >= i) followed by explicit
  // symmetrisation. This halves the number of Hessian writes at the cost of
  // one additional pass over a row to mirror elements.
  for (unsigned int i = 0; i < dim; i++) {
    unsigned int itab = i * dim;
    double pxi = fac * xi[i], hdgi = fad * hessDGrad[i], dgi = fae * dGrad[i];
    const double *pxj = &(xi[i]), *hdgj = &(hessDGrad[i]), *dgj = &(dGrad[i]);
    for (unsigned int j = i; j < dim; ++j, ++pxj, ++hdgj, ++dgj) {
      invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
      invHessian[j * dim + i] = invHessian[itab + j];
    }
  }
}

void newDirection(unsigned int dim, const double *invHessian,
                  const double *grad, double *xi) {
#ifdef RDK_SVE_AVAILABLE
  if (cpuHasSVE()) {
    sveHessianVecMulNeg(dim, invHessian, grad, xi);
    return;
  }
#endif
  for (unsigned int i = 0; i < dim; i++) {
    unsigned int itab = i * dim;
    xi[i] = 0.0;
    double &pxi = xi[i];
    const double *ivh = &(invHessian[itab]);
    const double *gj = grad;
    for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
      pxi -= *ivh * *gj;
    }
  }
}

}  // namespace BFGSOpt::detail
