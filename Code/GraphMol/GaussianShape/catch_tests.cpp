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
// Tests for the Roshambo2-based shape alignment code.

#include <chrono>
#include <random>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <execution>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

bool checkMolsHaveRoughlySameCoords(const ROMol &m1, const ROMol &m2,
                                    double margin = 0.005) {
  for (unsigned int i = 0; i < m1.getNumAtoms(); ++i) {
    auto pos1 = m1.getConformer().getAtomPos(i);
    auto pos2 = m2.getConformer().getAtomPos(i);
    if ((pos1 - pos2).length() > margin) {
      // So the error is printed in a relevant place.
      std::cout << i << " : " << m1.getAtomWithIdx(i)->getAtomicNum()
                << " :: " << (pos1 - pos2).length() << std::endl;
      CHECK_THAT((pos1 - pos2).length(),
                 Catch::Matchers::WithinAbs(0.0, margin));
      return false;
    }
  }
  return true;
}

TEST_CASE("basic alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/Code/GraphMol/GaussianShape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto refT = suppl[0];
  auto ref = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*refT));
  REQUIRE(ref);
  auto probeT = suppl[1];
  auto probe = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*probeT));
  REQUIRE(probe);
  GaussianShape::ShapeOverlayOptions overlayOpts;
  overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
  overlayOpts.nSteps = 50;
  GaussianShape::ShapeInputOptions shapeOpts;
  SECTION("setup") {
    auto refShape = GaussianShape::ShapeInput(*ref, -1, shapeOpts);
    CHECK_THAT(refShape.getShapeVolume(),
               Catch::Matchers::WithinAbs(591.057, 0.005));
    CHECK_THAT(refShape.getColorVolume(),
               Catch::Matchers::WithinAbs(31.935, 0.005));

    auto probeShape = GaussianShape::ShapeInput(*probe, -1, shapeOpts);
    CHECK_THAT(probeShape.getShapeVolume(),
               Catch::Matchers::WithinAbs(751.013, 0.005));
    CHECK_THAT(probeShape.getColorVolume(),
               Catch::Matchers::WithinAbs(42.530, 0.005));
  }
  SECTION("shape only") {
    ROMol cp(*probe);
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_ONLY;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;

    GaussianShape::ShapeInputOptions tShapeOpts;
    tShapeOpts.useColors = false;
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, tShapeOpts, tShapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.0, 0.005));
    // Check that a re-score gives the same answer.
    auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts, shapeOpts,
                                                 overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("shape plus color score") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.494, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.228, 0.005));
    // Check that a re-score gives the same answer.
    const auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts,
                                                       shapeOpts, overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("shape and color") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.512, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.747, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.273, 0.005));
    const auto rescores = GaussianShape::ScoreMolecule(*ref, cp, shapeOpts,
                                                       shapeOpts, overlayOpts);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
  }
  SECTION("collect transform") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::ROTATE_180;
    ROMol cp(*probe);
    RDGeom::Transform3D xform;
    const auto scores = GaussianShape::AlignMolecule(
        *ref, cp, shapeOpts, shapeOpts, &xform, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.494, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.228, 0.005));
    // Check a few values from the transform, just to be sure
    CHECK_THAT(xform.getValUnchecked(0, 0),
               Catch::Matchers::WithinAbs(-0.886, 0.005));
    CHECK_THAT(xform.getValUnchecked(1, 1),
               Catch::Matchers::WithinAbs(-0.834, 0.005));
    CHECK_THAT(xform.getValUnchecked(2, 2),
               Catch::Matchers::WithinAbs(0.827, 0.005));
    CHECK_THAT(xform.getValUnchecked(3, 3),
               Catch::Matchers::WithinAbs(1.0, 0.005));
  }
  SECTION("shape plus color score a la pubchem") {
    overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
    overlayOpts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
    GaussianShape::ShapeInputOptions shapeOpts2;
    for (const auto acr : std::vector{true, false}) {
      shapeOpts2.allCarbonRadii = acr;
      ROMol cp(*probe);
      const auto scores = GaussianShape::AlignMolecule(
          *ref, cp, shapeOpts2, shapeOpts2, nullptr, overlayOpts);
      if (acr) {
        CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.498, 0.005));
        CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.758, 0.005));
        CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.228, 0.005));
      } else {
        CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.503, 0.005));
        CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.761, 0.005));
        CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.233, 0.005));
      }
      // Check that a re-score gives the same answer.
      const auto rescores = GaussianShape::ScoreMolecule(
          *ref, cp, shapeOpts2, shapeOpts2, overlayOpts);
      CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
      CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
      CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
    }
  }
}

TEST_CASE("bulk") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";
  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/bulk.pubchem.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  std::string testout = dirName + "/bulk.pubchem_out.sdf";
  auto writer = SDWriter(testout);
  writer.write(*ref);
  GaussianShape::ShapeOverlayOptions overlayOpts;
  overlayOpts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  overlayOpts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
  GaussianShape::ShapeInputOptions shapeOpts;
  for (auto i = 1u; i < suppl.length(); ++i) {
    auto probe = suppl[1];
    REQUIRE(probe);
    auto scores = GaussianShape::AlignMolecule(*ref, *probe, shapeOpts,
                                               shapeOpts, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.575, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.824, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.332, 0.005));
    const auto rescores = GaussianShape::ScoreMolecule(*ref, *probe);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
    writer.write(*probe);
    break;
  }
  writer.close();
}

TEST_CASE("shape alignment") {
  std::string dirName = getenv("RDBASE");
  dirName += "/Code/GraphMol/GaussianShape/test_data";

  auto suppl = v2::FileParsers::SDMolSupplier(dirName + "/test1.sdf");
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  auto refShape = GaussianShape::ShapeInput(*ref, -1);
  auto probeShape = GaussianShape::ShapeInput(*probe, -1);
  const auto ovProbe =
      "FC1(F)C[C@H](C(O)=O)N(Cc2ocnc2)C1 |(-13.7583,-5.79896,4.44174;-13.5177,-6.64376,3.41248;-12.4622,-7.40849,3.74704;-13.2648,-5.86618,2.1389;-14.6452,-5.76234,1.48767;-15.1343,-4.34988,1.46304;-14.7802,-3.74082,0.305157;-15.7694,-3.81132,2.35739;-15.559,-6.51653,2.35483;-16.6354,-7.16994,1.6233;-17.7292,-6.20359,1.25389;-18.3329,-5.52754,2.27059;-19.2537,-4.75211,1.63512;-19.2878,-4.87389,0.327298;-18.3113,-5.80392,0.0765952;-14.7454,-7.46926,3.11334),wU:4.4|"_smiles;
  RDGeom::Transform3D xform;
  auto scores = GaussianShape::AlignShape(refShape, probeShape, &xform);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.498, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.760, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.228, 0.005));
  // This effectively checks that xform is correct.
  std::pair<double, double> overlapVols;
  auto rescores = GaussianShape::ScoreShape(
      refShape, probeShape, GaussianShape::ShapeOverlayOptions(), &overlapVols);
  CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.001));
  CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.001));
  CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.001));
  CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.001));
  CHECK_THAT(overlapVols.first, Catch::Matchers::WithinAbs(579.827, 0.005));
  CHECK_THAT(overlapVols.second, Catch::Matchers::WithinAbs(13.8469, 0.005));
  SmilesWriteParams params;
  params.canonical = false;
  // The input structure being from an SDF doesn't have the atoms in an order
  // that will make a SMILES string so bounce it through one for comparison.
  auto probeCp1 = v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*probe, params));
  MolTransforms::transformConformer(probeCp1->getConformer(), xform);
  CHECK(checkMolsHaveRoughlySameCoords(*ovProbe, *probeCp1));
  // And pre-normalizing the shapes
  refShape.normalizeCoords();
  probeShape.normalizeCoords();
  RDGeom::Transform3D xform1;
  auto scores1 = GaussianShape::AlignShape(refShape, probeShape, &xform1);
  CHECK_THAT(scores1[0], Catch::Matchers::WithinAbs(0.498, 0.005));
  CHECK_THAT(scores1[1], Catch::Matchers::WithinAbs(0.760, 0.005));
  CHECK_THAT(scores1[2], Catch::Matchers::WithinAbs(0.228, 0.005));
}

TEST_CASE("Overlay onto shape bug (Github8462)") {
  auto m1 =
      R"(c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|)"_smiles;
  REQUIRE(m1);
  ROMol m2(*m1);
  for (auto a : m2.atoms()) {
    auto &pos = m2.getConformer().getAtomPos(a->getIdx());
    pos.x += 3.0;
    pos.y += 2.0;
  }
  ROMol m3(m2);

  auto scores = GaussianShape::AlignMolecule(*m1, m2);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK(checkMolsHaveRoughlySameCoords(*m1, m2, 0.005));

  // Create the shape without normalization to mimic an arbitrary shape.
  auto s1 = GaussianShape::ShapeInput(*m1, -1);
  auto scores1 = AlignMolecule(s1, m3);
  CHECK_THAT(scores1[0], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores1[1], Catch::Matchers::WithinAbs(1.0, 0.005));
  CHECK_THAT(scores1[2], Catch::Matchers::WithinAbs(1.0, 0.005));
  for (unsigned int i = 0; i < m3.getNumAtoms(); ++i) {
    RDGeom::Point3D pos1(s1.getCoords()[3 * i], s1.getCoords()[3 * i + 1],
                         s1.getCoords()[3 * i + 2]);
    auto pos2 = m3.getConformer().getAtomPos(i);
    CHECK_THAT((pos1 - pos2).length(), Catch::Matchers::WithinAbs(0.0, 0.01));
  }
}

TEST_CASE("handling molecules with Hs") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  v2::FileParsers::MolFileParserParams params;
  params.removeHs = false;
  auto suppl =
      v2::FileParsers::SDMolSupplier(dirName + "/align_with_hs.sdf", params);
  auto ref = suppl[0];
  REQUIRE(ref);
  auto probe = suppl[1];
  REQUIRE(probe);
  SECTION("basics") {
    RWMol cp(*probe);
    RDGeom::Transform3D xform;
    auto scores = GaussianShape::AlignMolecule(*ref, cp);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.802, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.840, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.763, 0.005));
    for (auto i = 0u; i < cp.getNumAtoms(); ++i) {
      // the failure mode here was that Hs had HUGE coordinates
      auto pos = cp.getConformer().getAtomPos(i);
      CHECK((pos.x > -10 && pos.x < 10));
    }
    // Check the rescore
    std::pair<double, double> overlapVols;
    auto rescores = GaussianShape::ScoreMolecule(
        *ref, cp, GaussianShape::ShapeInputOptions(),
        GaussianShape::ShapeInputOptions(),
        GaussianShape::ShapeOverlayOptions(), -1, -1, &overlapVols);
    CHECK_THAT(rescores[0], Catch::Matchers::WithinAbs(scores[0], 0.005));
    CHECK_THAT(rescores[1], Catch::Matchers::WithinAbs(scores[1], 0.005));
    CHECK_THAT(rescores[2], Catch::Matchers::WithinAbs(scores[2], 0.005));
    CHECK_THAT(overlapVols.first, Catch::Matchers::WithinAbs(356.373, 0.005));
    CHECK_THAT(overlapVols.second, Catch::Matchers::WithinAbs(9.202, 0.005));
  }
}

TEST_CASE("Github #8096") {
  SECTION("as reported") {
    auto m1 =
        R"([H]c1c([H])c([H])c([2H])c([H])c1[H] |(1.55967,1.91617,0.0546381;0.885536,1.07172,0.030849;1.38172,-0.23747,0.0274262;2.44539,-0.439501,0.0483424;0.470206,-1.27516,-0.00361916;0.856925,-2.30002,-0.00633525;-0.896665,-1.07227,-0.0310991;-1.60071,-1.87642,-0.0551085;-1.36315,0.22877,-0.0271173;-2.43593,0.379132,-0.0487835;-0.479018,1.29083,0.00359778;-0.823965,2.31421,0.00720933)|)"_smiles;
    REQUIRE(m1);
    auto m2 =
        R"([H]c1c([H])c([H])c([H])c([H])c1[H] |(-2.06264,-0.844763,-0.0261403;-1.04035,-0.481453,-0.0114878;-0.00743655,-1.41861,-0.0137121;-0.215455,-2.47997,-0.0295909;1.29853,-0.949412,0.00507497;2.12524,-1.65277,0.00390664;1.58501,0.395878,0.0254188;2.61997,0.704365,0.0394811;0.550242,1.31385,0.0273741;0.783172,2.37039,0.0434262;-0.763786,0.88847,0.00908113;-1.60557,1.58532,0.0100194)|)"_smiles;
    REQUIRE(m2);
    auto scores = GaussianShape::AlignMolecule(*m1, *m2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.005));
  }
}

TEST_CASE("Hs not properly transformed when hcount = feature count") {
  std::string dirName = getenv("RDBASE");
  dirName += "/External/pubchem_shape/test_data";

  SECTION("as reported") {
    v2::FileParsers::MolFileParserParams ps;
    ps.removeHs = false;
    auto mol1 =
        v2::FileParsers::MolFromMolFile(dirName + "/hcount_ex1_1.mol", ps);
    REQUIRE(mol1);
    auto mol2 =
        v2::FileParsers::MolFromMolFile(dirName + "/hcount_ex1_2.mol", ps);
    REQUIRE(mol2);
    {
      RWMol cp(*mol2);
      auto scores = GaussianShape::AlignMolecule(*mol1, cp);
      CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.734, 0.005));
      CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.918, 0.005));
      CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.550, 0.005));
      // the bug led to H atoms in stupid positions, so we can detect it by just
      // looking at bond lengths to Hs:
      for (auto i = cp.getNumHeavyAtoms(); i < cp.getNumAtoms(); ++i) {
        INFO("checking atom " << i);
        auto at = cp.getAtomWithIdx(i);
        for (auto nbr : cp.atomNeighbors(at)) {
          auto dist = (cp.getConformer().getAtomPos(i) -
                       cp.getConformer().getAtomPos(nbr->getIdx()))
                          .length();
          CHECK(dist < 1.2);  // should be a bond to H
        }
      }
    }
  }
}

TEST_CASE("Score No Overlay") {
  // These are 2 ligands used by Andy Grant and Co in their original paper
  // https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(19961115)17:14%3C1653::AID-JCC7%3E3.0.CO;2-K
  // Ligands as extracted from PDB, with a bit of munging to get them as
  // SMILES strings (downloaded the Ideal ligand structures from RCSB
  // as SDFs and transferred the corresponding atom coords from 3tmn and 1tmn).
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  GaussianShape::ShapeInputOptions shapeOpts;
  {
    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto pdb_trp_3tmn_cp =
        R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
    RDGeom::Point3D trans{100.0, 100.0, 100.0};
    RDGeom::Transform3D transform_3d;
    transform_3d.SetTranslation(trans);
    MolTransforms::transformConformer(pdb_trp_3tmn_cp->getConformer(),
                                      transform_3d);

    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_trp_3tmn_cp, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.0, 0.001));
  }
  {
    auto scores = ScoreMolecule(*pdb_0zn_1tmn, *pdb_0zn_1tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(1.0, 0.001));
  }
  {
    auto scores = ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn, shapeOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
  }
  {
    auto shape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    std::pair<double, double> overlapVols;
    auto scores =
        ScoreMolecule(shape, *pdb_0zn_1tmn, shapeOpts,
                      GaussianShape::ShapeOverlayOptions(), -1, &overlapVols);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
    CHECK_THAT(overlapVols.first, Catch::Matchers::WithinAbs(590.831, 0.005));
    CHECK_THAT(overlapVols.second, Catch::Matchers::WithinAbs(26.736, 0.005));
  }
  {
    auto shape1 = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    auto shape2 = GaussianShape::ShapeInput(*pdb_0zn_1tmn, -1, shapeOpts);
    auto scores = GaussianShape::ScoreShape(shape1, shape2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.307, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.349, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.265, 0.005));
  }
}

TEST_CASE("Iressa onto Tagrisso") {
  // Conformations from PubChem produced by Omega. Iressa rotated and translated
  // by a random amount.  PubChem puts them both in their inertial frame which
  // makes things too easy.
  auto tagrisso =
      R"(C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C |(-0.9161,3.8415,-2.9811;0.1848,3.1933,-2.588;0.1064,1.7789,-2.12;-0.9619,1.1797,-2.0847;1.3654,1.2872,-1.7553;1.6841,0.0144,-1.273;0.6638,-0.9235,-1.1146;0.9578,-2.1997,-0.6343;-0.0813,-3.1358,-0.4783;-1.4556,-2.9979,-0.1847;-2.1716,-4.1359,-0.1085;-3.4803,-3.9673,0.173;-4.0689,-2.7353,0.3728;-3.2269,-1.647,0.2676;-3.7311,-0.317,0.4568;-5.0275,0.0291,0.153;-5.1887,1.3569,0.4454;-6.4231,2.0889,0.2595;-4.0141,1.8811,0.9361;-3.7121,3.1796,1.3615;-2.4139,3.4249,1.8179;-1.4588,2.4106,1.8467;-1.7752,1.1164,1.4181;-3.0776,0.8453,0.9537;-1.9103,-1.7423,-0.011;2.2723,-2.5382,-0.3127;2.58,-3.7798,0.1575;2.539,-3.9651,1.571;3.2927,-1.6003,-0.4713;2.9986,-0.324,-0.9514;4.0475,0.61,-1.1047;4.7738,0.6956,-2.3546;4.4021,1.497,-0.0162;5.4401,0.8254,0.8736;5.8294,1.7155,1.9601;4.8213,1.7057,3.0218;7.1361,1.3324,2.4981)|)"_smiles;
  REQUIRE(tagrisso);
  auto iressa =
      R"(COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 |(11.4672,-0.467948,5.63989;12.0133,0.532631,6.49693;11.2039,1.5801,6.81985;11.2014,2.71958,6.00975;10.3926,3.81652,6.29699;10.4038,4.90395,5.50623;9.58889,5.91871,5.85946;8.76443,5.96486,6.91838;8.77814,4.86059,7.68868;7.92337,4.86224,8.81914;7.44878,5.8925,9.64622;8.22182,7.03851,9.85619;7.75051,8.06265,10.6777;6.50441,7.94546,11.2936;6.06567,8.93802,12.0809;5.72932,6.80403,11.0875;4.19056,6.65372,11.8447;6.20047,5.78015,10.2656;9.57161,3.74547,7.43436;9.56851,2.60328,8.25407;10.3868,1.52933,7.93769;10.3797,0.419365,8.74203;11.3064,0.402096,9.81907;10.7104,-0.399165,10.9685;9.40938,0.22121,11.4678;9.64205,1.59049,11.9223;8.38006,2.22985,12.3199;8.64991,3.65266,12.8011;9.56883,3.64192,13.8942;10.8103,3.05101,13.5078;10.5931,1.61394,13.0425)|)"_smiles;
  REQUIRE(iressa);
  GaussianShape::ShapeOverlayOptions opts;
  opts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  opts.startMode = GaussianShape::StartMode::A_LA_PUBCHEM;
  opts.nSteps = 100;
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.allCarbonRadii = false;
  auto scores = GaussianShape::AlignMolecule(*tagrisso, *iressa, shapeOpts,
                                             shapeOpts, nullptr, opts);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.331, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.572, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.090, 0.005));
  auto rescores = GaussianShape::ScoreMolecule(*tagrisso, *iressa, shapeOpts,
                                               shapeOpts, opts);
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(rescores[0], 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(rescores[1], 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(rescores[2], 0.005));
  auto aligned_iressa =
      "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 |(3.07478,-4.69769,0.157642;2.27093,-4.19922,-0.909719;1.10841,-3.57503,-0.569809;-0.0646997,-4.32671,-0.453791;-1.27804,-3.73592,-0.109237;-2.39613,-4.47513,-0.00257941;-3.52557,-3.82055,0.335758;-3.65891,-2.50653,0.578133;-2.52086,-1.79632,0.462366;-2.61158,-0.403683,0.709139;-3.65256,0.522606,0.538616;-4.60037,0.337355,-0.472503;-5.63512,1.25735,-0.642125;-5.72667,2.36716,0.197977;-6.72653,3.24279,0.0217711;-4.78347,2.55667,1.2082;-4.8913,3.92427,2.24828;-3.74903,1.6364,1.37808;-1.29043,-2.35038,0.120551;-0.115093,-1.58753,0.00575868;1.07556,-2.20919,-0.339134;2.21748,-1.4591,-0.450542;2.50073,-0.908555,-1.72951;3.24779,0.405383,-1.54539;2.41621,1.41722,-0.763407;1.15935,1.68631,-1.45877;0.309613,2.58928,-0.670051;-0.990058,2.86435,-1.42102;-0.713572,3.43258,-2.70167;0.0855337,2.54759,-3.48777;1.41093,2.2636,-2.78667)|"_smiles;
  REQUIRE(aligned_iressa);
  CHECK(checkMolsHaveRoughlySameCoords(*iressa, *aligned_iressa));
}

TEST_CASE("Optimise in place") {
  // These are 2 ligands used by Andy Grant and Co in their original paper
  // https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1096-987X(19961115)17:14%3C1653::AID-JCC7%3E3.0.CO;2-K
  // Ligands as extracted from PDB, with a bit of munging to get them as
  // SMILES strings (downloaded the Ideal ligand structures from RCSB
  // as SDFs and transferred the corresponding atom coords from 3tmn and 1tmn).
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  // This is the overlay produced by the first test below, to make sure we
  // haven't broken anything.
  auto ov_pdb_0zn_1tmn =
      R"(CC(C)C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O |(38.0714,43.9721,-0.878703;37.9576,42.4574,-0.719904;37.7577,42.046,0.732339;36.7739,41.9272,-1.50675;36.9772,42.2935,-2.97812;35.7268,42.027,-3.73126;35.8612,42.1927,-5.23029;34.6327,41.8342,-6.08608;34.7832,40.375,-6.53487;33.6583,39.7599,-7.396;33.6467,39.94,-8.77358;32.6733,39.3165,-9.57847;31.7041,38.5036,-8.98298;31.7526,38.3068,-7.58829;32.7402,38.9282,-6.8042;36.2665,43.5574,-5.73451;37.1622,43.6837,-6.59647;35.6253,44.5489,-5.37542;38.2382,41.5679,-3.46678;39.3654,42.1399,-3.47633;38.0044,40.307,-3.81274;39.0541,39.3969,-4.30077;39.8703,38.7572,-3.13916;39.0028,37.7219,-2.42823;37.728,37.8163,-2.05559;37.2322,36.5297,-1.76132;38.1925,35.7015,-1.88597;38.1776,34.3299,-1.68955;39.4475,33.7105,-1.64072;40.652,34.4406,-1.79973;40.6391,35.8191,-2.08104;39.3808,36.4168,-2.13072;38.5971,38.3663,-5.32643;37.4133,38.4775,-5.77705;39.4717,37.6908,-5.97898),wD:6.6,wU:4.3,21.22|)"_smiles;
  auto initScores = GaussianShape::ScoreMolecule(*pdb_trp_3tmn, *pdb_0zn_1tmn);
  CHECK_THAT(initScores[0], Catch::Matchers::WithinAbs(0.307, 0.001));
  CHECK_THAT(initScores[1], Catch::Matchers::WithinAbs(0.349, 0.001));
  CHECK_THAT(initScores[2], Catch::Matchers::WithinAbs(0.265, 0.001));
  // The PDB atom order isn't canonical, so bounce in and out of SMILES
  // to make it easier to check.
  auto canon_probe =
      v2::SmilesParse::MolFromSmiles(MolToCXSmiles(*pdb_0zn_1tmn));
  {
    // This should just tweak the input overlay.
    GaussianShape::ShapeOverlayOptions opts;
    opts.startMode = GaussianShape::StartMode::ROTATE_0;
    opts.normalize = false;
    GaussianShape::ShapeInputOptions shapeOpts;
    ROMol cp(*canon_probe);
    RDGeom::Transform3D xform;
    auto scores = GaussianShape::AlignMolecule(*pdb_trp_3tmn, cp, shapeOpts,
                                               shapeOpts, &xform, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.337, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.395, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.278, 0.001));
    CHECK(checkMolsHaveRoughlySameCoords(cp, *ov_pdb_0zn_1tmn));
  }
  {
    // With default settings, it does a poor job.
    GaussianShape::ShapeOverlayOptions opts;
    GaussianShape::ShapeInputOptions shapeOpts;
    ROMol cp(*canon_probe);
    auto scores = GaussianShape::AlignMolecule(*pdb_trp_3tmn, cp, shapeOpts,
                                               shapeOpts, nullptr, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.206, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.371, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.040, 0.001));
  }
  {
    // And with reference as a shape the same
    GaussianShape::ShapeOverlayOptions opts;
    opts.startMode = GaussianShape::StartMode::ROTATE_0;
    opts.normalize = false;
    GaussianShape::ShapeInputOptions shapeOpts;
    auto refShape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    ROMol cp(*canon_probe);
    RDGeom::Transform3D xform;
    auto scores =
        GaussianShape::AlignMolecule(refShape, cp, shapeOpts, &xform, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.337, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.395, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.278, 0.001));
    CHECK(checkMolsHaveRoughlySameCoords(cp, *ov_pdb_0zn_1tmn));
    MolTransforms::transformConformer(cp.getConformer(), xform);
    ROMol cp1(*canon_probe);
    MolTransforms::transformConformer(cp1.getConformer(), xform);
    CHECK(checkMolsHaveRoughlySameCoords(cp1, *ov_pdb_0zn_1tmn));
  }
  {
    // And with both as shapes
    GaussianShape::ShapeOverlayOptions opts;
    opts.startMode = GaussianShape::StartMode::ROTATE_0;
    opts.normalize = false;
    GaussianShape::ShapeInputOptions shapeOpts;
    auto refShape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1, shapeOpts);
    auto fitShape = GaussianShape::ShapeInput(*canon_probe, -1, shapeOpts);
    RDGeom::Transform3D xform;
    auto scores = GaussianShape::AlignShape(refShape, fitShape, &xform, opts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.337, 0.001));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.395, 0.001));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.278, 0.001));
    auto cp = fitShape.shapeToMol(false);
    CHECK(checkMolsHaveRoughlySameCoords(*cp, *ov_pdb_0zn_1tmn));
  }
}

TEST_CASE("Fragment Mode") {
  // On the PDB overlay.
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  GaussianShape::ShapeOverlayOptions opts;
  opts.nSteps = 100;
  opts.startMode = GaussianShape::StartMode::ROTATE_180_FRAGMENT;
  opts.optimMode = GaussianShape::OptimMode::SHAPE_PLUS_COLOR_SCORE;
  auto probeShape = GaussianShape::ShapeInput(*pdb_trp_3tmn, -1);
  auto refShape = GaussianShape::ShapeInput(*pdb_0zn_1tmn, -1);
  RDGeom::Transform3D xform;
  // Use the smaller molecule as the probe
  auto scores = GaussianShape::AlignShape(refShape, probeShape, &xform, opts);
  // These are close to the values above for starting from the xtal structures.
  CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.332, 0.005));
  CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(0.413, 0.005));
  CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.249, 0.005));
}

TEST_CASE("custom feature points") {
  auto m1 =
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|"_smiles;
  SECTION("using shapes") {
    auto shape1 = GaussianShape::ShapeInput(*m1, -1);
    // each carbonyl O gets one feature:
    CHECK(shape1.getCoords().size() == 18);
    GaussianShape::ShapeInputOptions opts2;
    opts2.customFeatures = {{{1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0,
                              std::vector<unsigned int>{}},
                             {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0,
                              std::vector<unsigned int>{}}}};
    auto shape2 = GaussianShape::ShapeInput(*m1, -1, opts2);
    CHECK(shape2.getCoords().size() == 18);

    {
      // confirm that we don't add the features if not requested.
      GaussianShape::ShapeInputOptions topts;
      topts.customFeatures = {{{1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0,
                                std::vector<unsigned int>{}},
                               {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0,
                                std::vector<unsigned int>{}}}};
      topts.useColors = false;
      auto tshape = GaussianShape::ShapeInput(*m1, -1, topts);
      CHECK(tshape.getCoords().size() == 12);
    }

    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    GaussianShape::ShapeInputOptions opts3;
    opts3.customFeatures = {{{2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0,
                              std::vector<unsigned int>{}},
                             {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0,
                              std::vector<unsigned int>{}}}};

    auto m2 = ROMol(*m1);
    auto shape3 = GaussianShape::ShapeInput(m2, -1, opts3);
    CHECK(shape3.getCoords().size() == 18);
    GaussianShape::ShapeOverlayOptions overlayOpts;
    overlayOpts.optParam = 0.5;
    RDGeom::Transform3D xform;
    auto scores = AlignShape(shape2, shape3, &xform, overlayOpts);

    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.999, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.000, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.998, 0.005));
    CHECK(shape3.getCoords()[0] > 0);      // x coord of first atom
    CHECK(shape3.getCoords()[3 * 3] < 0);  // x coord of fourth atom

    auto conf = m2.getConformer(-1);
    MolTransforms::transformConformer(conf, xform);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
  SECTION("using molecules") {
    GaussianShape::ShapeInputOptions opts2;
    opts2.customFeatures = {{{1, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0,
                              std::vector<unsigned int>{}},
                             {2, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0,
                              std::vector<unsigned int>{}}}};

    auto m2 = ROMol(*m1);
    // we'll swap the features on the second shape so that the alignment has to
    // be inverted
    GaussianShape::ShapeInputOptions opts3;
    opts3.customFeatures = {{{2, RDGeom::Point3D(-1.75978, 0.148897, 0), 1.0,
                              std::vector<unsigned int>{}},
                             {1, RDGeom::Point3D(1.7571, -0.120174, 0.1), 1.0,
                              std::vector<unsigned int>{}}}};

    GaussianShape::ShapeOverlayOptions overlayOpts;
    overlayOpts.optParam = 0.5;
    std::vector<float> matrix(12, 0.0);
    auto scores = AlignMolecule(*m1, m2, opts2, opts3, nullptr, overlayOpts);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.999, 0.005));
    CHECK_THAT(scores[1], Catch::Matchers::WithinAbs(1.000, 0.005));
    CHECK_THAT(scores[2], Catch::Matchers::WithinAbs(0.998, 0.005));
    auto conf = m2.getConformer(-1);
    CHECK(conf.getAtomPos(0).x > 0);
    CHECK(conf.getAtomPos(3).x < 0);
  }
}

TEST_CASE("Non-standard radii") {
  auto m1 =
      "[Xe]c1ccccc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.useColors = false;
  shapeOpts.allCarbonRadii = false;
  auto shape1 = GaussianShape::ShapeInput(*m1, -1, shapeOpts);

  CHECK(shape1.getCoords().size() == 21);
  CHECK_THAT(shape1.getShapeVolume(),
             Catch::Matchers::WithinAbs(387.396, 0.005));
  // mol1 with atom 4 with an N radius and a bigger Xe.
  shapeOpts.atomRadii =
      std::vector<std::pair<unsigned int, double>>{{0, 2.5}, {4, 1.55}};
  shapeOpts.allCarbonRadii = false;
  auto shape2 = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  CHECK_THAT(shape2.getShapeVolume(),
             Catch::Matchers::WithinAbs(425.051, 0.005));

  // Corresponding pyridine derivative.
  auto m2 =
      "[Xe]c1ccncc1 |(0.392086,-2.22477,0.190651;0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097),atomProp:0.dummyLabel.*|"_smiles;
  auto shape3 = GaussianShape::ShapeInput(*m2, -1, shapeOpts);
  CHECK(shape3.getShapeVolume() == shape2.getShapeVolume());
}

TEST_CASE("Shape subset") {
  auto m1 =
      "c1ccc(-c2ccccc2)cc1 |(-3.26053,-0.0841607,-0.741909;-2.93383,0.123873,0.593407;-1.60713,0.377277,0.917966;-0.644758,0.654885,-0.0378428;0.743308,0.219134,0.168663;1.82376,1.0395,-0.0112769;3.01462,0.695405,0.613858;3.18783,-0.589771,1.09649;2.15761,-1.50458,1.01949;0.988307,-1.1313,0.385783;-1.1048,0.797771,-1.34022;-2.39754,0.435801,-1.69921)|"_smiles;
  REQUIRE(m1);
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.atomSubset = std::vector<unsigned int>{0, 1, 2, 3, 10, 11};
  GaussianShape::ShapeInput partShape(*m1, -1, shapeOpts);
  CHECK(partShape.getSmiles() == "c1ccccc1");
  CHECK(partShape.getNumAtoms() == 6);
  CHECK(partShape.getNumFeatures() == 1);
  CHECK(partShape.getCoords().size() == 21);
  CHECK_THAT(partShape.getShapeVolume(),
             Catch::Matchers::WithinAbs(261.166, 0.005));
  CHECK_THAT(partShape.getColorVolume(),
             Catch::Matchers::WithinAbs(5.316, 0.005));

  shapeOpts.atomSubset.clear();
  GaussianShape::ShapeInput wholeShape(*m1, -1, shapeOpts);
  CHECK(wholeShape.getSmiles() == "c1ccc(-c2ccccc2)cc1");
  CHECK(wholeShape.getNumAtoms() == 12);
  CHECK(wholeShape.getNumFeatures() == 2);
  CHECK(wholeShape.getCoords().size() == 42);
  CHECK_THAT(wholeShape.getShapeVolume(),
             Catch::Matchers::WithinAbs(556.266, 0.005));
  CHECK_THAT(wholeShape.getColorVolume(),
             Catch::Matchers::WithinAbs(10.631, 0.005));
}

// These are LOBSTER structures 437_A_355, YIW_A_1353, LSA_A_503, SU0_A_263,
// VHC_A_1, 40Z_A_301, 0J8_A_1401, 5QQ_A_1401, 054_A_578, 053_A_578
// respectively. LOBSTER is published
// https://doi.org/10.1007/s10822-024-00581-1 from the Rarey and BioSolveIT
// group.
std::
    vector<std::string>
        lobstersText =
            {
                "CC(C)(C)c1cc(NC(=O)Nc2cccc3ccccc23)n(-c2ccc(CO)cc2)n1 |(4.1858,1.2187,12.6749;3.4917,2.128,11.6409;4.5532,3.0576,11.0244;2.9098,1.2685,10.5255;2.3016,2.8612,12.2277;1.3306,2.2548,13.0356;0.4334,3.2937,13.302;-0.7316,3.3518,14.0275;-0.978,2.5581,15.0932;-0.147,1.7378,15.4901;-2.1645,2.6924,15.6938;-2.626,2.0319,16.8298;-1.8457,2.1016,18.0083;-2.2654,1.5179,19.2117;-3.4858,0.8316,19.2546;-4.2646,0.7339,18.1056;-5.4724,0.0393,18.1895;-6.2813,-0.103,17.0695;-5.8756,0.4693,15.8666;-4.6776,1.1996,15.7839;-3.8502,1.3304,16.8984;0.8824,4.3943,12.6425;0.3828,5.6491,12.5902;-0.1458,6.3636,13.6677;-0.6231,7.6713,13.5508;-0.5802,8.3258,12.3036;-1.1083,9.738,12.1288;-0.2929,10.4745,11.2057;-0.0699,7.6584,11.21;0.4028,6.3542,11.3691;2.0606,4.1572,11.9546)|",
                "CC(C)c1nnc2ccc(Sc3ccccc3CNC(=O)Nc3cc(C(C)(C)C)nn3-c3ccccc3)cn12 |(-2.677,-1.147,25.057;-1.383,-1.713,24.501;-0.937,-2.645,25.64;-1.654,-2.496,23.218;-1.533,-3.804,23.035;-1.814,-4.126,21.82;-2.125,-3.054,21.15;-2.484,-2.842,19.8;-2.737,-1.521,19.397;-2.633,-0.45,20.292;-2.934,1.2,19.817;-4.086,0.989,18.473;-5.306,0.233,18.595;-6.145,0.057,17.523;-5.802,0.63,16.3;-4.651,1.427,16.155;-3.805,1.585,17.237;-2.532,2.33,17.092;-2.317,3.124,15.849;-1.208,2.896,15.139;-0.293,2.051,15.461;-0.934,3.704,14.164;0.194,3.633,13.319;1.093,2.626,13.06;2.064,3.138,12.153;3.254,2.438,11.547;2.765,1.497,10.446;4.297,3.407,10.979;3.923,1.592,12.556;1.755,4.357,11.974;0.557,4.674,12.592;0.026,5.998,12.604;-0.011,6.769,11.405;-0.52,8.07,11.439;-0.946,8.61,12.635;-0.914,7.879,13.819;-0.396,6.572,13.782;-2.292,-0.677,21.581;-2.025,-1.986,22.016)|",
                "O=S1(=O)N=C(O)c2ccccc21 |(-5.7089,1.0252,18.4004;-6.3943,1.0779,17.1227;-7.8251,0.8757,17.0645;-5.6185,0.1418,15.9961;-4.8784,0.7874,15.1629;-4.3563,0.192,14.2273;-5.1972,2.2654,15.1438;-4.6774,3.2806,14.3416;-5.0543,4.5982,14.6171;-5.8378,4.8973,15.7426;-6.3084,3.8774,16.5767;-5.9587,2.5637,16.2652)|",
                "COc1ccc2c(CC(=O)Nc3ccc(S(N)(=O)=O)cc3)cc(=O)oc2c1 |(-2.8164,14.7062,11.3592;-3.6624,13.537,11.2643;-3.2822,12.2442,11.6717;-3.8543,11.0312,10.9854;-3.454,9.7632,11.3397;-2.4737,9.5893,12.4478;-2.0245,8.2354,12.8485;-2.4284,7.0271,12.0534;-3.4577,6.2276,12.8381;-4.6459,6.5341,12.9744;-2.9341,5.1799,13.652;-3.8043,4.2751,14.3593;-3.0353,3.1465,14.8909;-3.6696,2.1888,15.6247;-5.1279,2.3197,15.867;-5.8526,1.0165,16.8148;-5.1597,-0.2684,16.1766;-7.2644,1.0484,16.6895;-5.3596,1.2234,18.1529;-5.863,3.3573,15.3579;-5.1676,4.4008,14.5479;-1.2733,8.1678,14.0258;-0.8262,9.3446,14.6451;-0.1397,9.3025,15.6623;-1.0563,10.6342,14.1829;-1.9446,10.7231,13.0793;-2.3948,12.0761,12.7036)|",
                "Cc1cc(C)c2cc1C(=O)NCCCOc1cccc(c1)Sc1cc-2nc(N)n1 |(73.8435,34.0723,26.5156;72.3815,34.1628,26.0388;71.6766,32.9652,25.9327;70.3451,32.9485,25.4823;69.6347,31.5936,25.3834;69.7284,34.1492,25.0935;70.4497,35.3542,25.1896;71.7747,35.3826,25.6823;72.4585,36.7442,25.8391;73.1172,37.0066,26.8611;72.1633,37.6292,24.8945;72.7135,39.0047,24.8716;71.7025,40.047,24.3677;71.3796,39.8568,22.8758;70.4153,38.7939,22.9499;69.8006,38.3095,21.8404;70.2406,38.4827,20.5148;69.5441,37.8479,19.489;68.407,37.0853,19.7393;68.0155,36.8722,21.0628;68.6923,37.5061,22.0951;66.5592,35.8914,21.365;66.9463,35.055,22.8617;68.2573,34.8374,23.3077;68.4281,34.247,24.558;67.3421,33.831,25.243;66.0874,33.9792,24.7404;65.0509,33.4765,25.428;65.8935,34.6016,23.5436)|",
                "Cc1c2c(n3c1CCN(Cc1ccco1)c1cc(C(N)=O)c(Cl)cc1-3)CC(C)(C)CC2=O |(74.8244,36.0896,26.0638;73.6879,35.186,25.6743;73.8428,33.8098,25.2429;72.5848,33.3216,24.9555;71.6743,34.332,25.2034;72.3429,35.4561,25.6334;71.7112,36.7663,26.0021;70.9895,37.4573,24.8514;70.1361,36.585,24.0367;69.7698,37.0101,22.6942;70.4319,38.0048,21.8175;70.8144,39.2778,21.9818;71.3289,39.7044,20.7598;71.23,38.6778,19.9306;70.6882,37.6078,20.5486;69.5068,35.3697,24.5106;68.1241,35.239,24.3599;67.4459,34.0963,24.7908;65.9689,34.0392,24.5456;65.1976,33.5111,25.4835;65.5111,34.4884,23.4842;68.1905,33.0839,25.4012;67.4671,31.5959,25.9018;69.5536,33.1987,25.5786;70.2398,34.3166,25.103;72.2711,31.9512,24.4312;73.4524,31.2236,23.7805;73.1395,29.7332,23.6933;73.6857,31.7544,22.3618;74.7055,31.4432,24.6434;74.982,32.872,25.0594;76.1434,33.2319,25.2012)|",
                "O=[N+]([O-])c1cccc(CNc2nc(C(F)(F)F)nc3ncc(-c4cnn(C5CCNCC5)c4)cc23)c1 |(-1.438,-13.226,20.761;-2.668,-13.702,20.695;-3.225,-14.278,21.715;-3.449,-13.606,19.5;-4.865,-13.703,19.611;-5.66,-13.615,18.465;-5.057,-13.44,17.21;-3.642,-13.329,17.113;-2.985,-13.176,15.754;-3.208,-14.453,15.102;-2.183,-15.47,15.231;-1.174,-15.247,16.105;-0.207,-16.166,16.24;0.913,-15.872,17.233;1.607,-16.968,17.501;1.653,-14.932,16.679;0.431,-15.434,18.382;-0.155,-17.307,15.532;-1.134,-17.587,14.626;-1.051,-18.759,13.937;-2.004,-19.081,13.034;-3.067,-18.213,12.807;-4.134,-18.597,11.829;-4.958,-17.727,11.105;-5.807,-18.503,10.373;-5.518,-19.836,10.648;-6.22,-21.008,10.059;-7.613,-21.132,10.635;-8.385,-22.302,9.993;-8.344,-22.271,8.528;-7.465,-21.299,7.857;-6.163,-20.922,8.548;-4.513,-19.902,11.543;-3.18,-16.987,13.52;-2.183,-16.667,14.458;-2.84,-13.422,18.253)|",
                "Fc1ccc(-c2cnc3nnc(C(F)(F)c4ccc5ncccc5c4)n3n2)cc1 |(-8.9341,-13.5345,15.2941;-7.8624,-13.6169,16.1019;-8.0323,-13.604,17.4747;-6.927,-13.6949,18.3088;-5.6563,-13.7867,17.7687;-4.5363,-13.8655,18.5962;-4.6798,-14.0273,19.9861;-3.5769,-14.1,20.7695;-2.3321,-14.0007,20.1927;-1.0691,-14.049,20.7489;-0.1759,-13.9183,19.7353;-0.8597,-13.8052,18.581;-0.2877,-13.6262,17.1957;1.0077,-13.3321,17.3153;-0.932,-12.6293,16.5758;-0.3953,-14.89,16.3793;0.3679,-16.0044,16.7321;0.2944,-17.1762,15.9833;-0.5658,-17.2197,14.8472;-0.6484,-18.3768,14.1035;-1.4592,-18.4618,13.0237;-2.2424,-17.3707,12.6345;-2.1971,-16.184,13.3658;-1.3263,-16.1253,14.5037;-1.2424,-14.9308,15.2731;-2.1994,-13.8502,18.8681;-3.3091,-13.7741,18.0449;-5.4848,-13.7974,16.3944;-6.5897,-13.7169,15.5608)|",
                "Nc1cc(Cn2c(C(=O)O)c(-n3c(=O)[nH]c4cscc4c3=O)c3cc(C(F)(F)F)ccc32)ccn1 |(29.5323,45.8636,43.104;28.8655,44.7866,42.5152;27.6833,44.3191,43.0964;27.0256,43.2322,42.5026;25.769,42.5937,43.0062;25.1539,43.3311,44.0877;24.4106,44.4698,43.8927;24.1642,45.0354,42.563;23.9009,46.3564,42.5697;24.0503,44.4039,41.5301;23.9917,44.8584,45.1623;23.2031,45.9947,45.3911;21.8869,45.9724,44.834;21.3954,44.9674,44.2853;21.1644,47.1439,44.9613;21.7067,48.2779,45.5571;21.0928,49.4972,45.696;22.1115,50.6335,46.4923;23.3622,49.4671,46.6456;23.0143,48.2553,46.0988;23.8418,47.054,46.0577;24.9808,47.0096,46.5158;24.5049,43.9455,46.1213;24.4196,43.8288,47.5202;25.0722,42.7404,48.1424;25.0173,42.5333,49.6217;24.5189,41.3151,49.9335;26.2412,42.5919,50.1868;24.2862,43.4348,50.3041;25.7949,41.7994,47.3917;25.8973,41.9011,46.0107;25.2386,42.9843,45.4045;27.5907,42.6386,41.3749;28.7511,43.1735,40.8552;29.3972,44.2337,41.3948)|",
                "Nc1cc(Cn2c(C(=O)O)c(-c3ccc[nH]c3=O)c3cc(C(F)(F)F)ccc32)ccn1 |(29.4427,46.05,43.4472;28.8505,44.9996,42.7518;27.6724,44.4529,43.2669;27.0734,43.3813,42.5955;25.8189,42.6999,43.052;25.2163,43.3606,44.2071;24.4581,44.4943,44.0838;24.1931,45.0901,42.7808;24.1204,46.4333,42.7998;23.8923,44.4562,41.7912;24.0461,44.8297,45.3707;23.2097,45.9536,45.7413;23.5041,46.7963,46.7667;22.6287,47.9099,47.1092;21.5037,48.1244,46.3872;21.1624,47.2732,45.3588;21.9195,46.1834,44.9822;21.563,45.4266,44.0641;24.5817,43.8725,46.2645;24.4958,43.6833,47.6538;25.1452,42.5857,48.2292;25.0531,42.3252,49.6881;24.5934,41.0857,49.9422;26.2447,42.4118,50.3073;24.2575,43.1875,50.3428;25.8943,41.6963,47.4464;25.9994,41.8646,46.0771;25.324,42.958,45.5057;27.6791,42.8877,41.4291;28.8324,43.5019,40.9842;29.4267,44.5435,41.6155)|",
};
std::vector<std::shared_ptr<RWMol>> lobsters;
void initLobsters() {
  if (lobsters.empty()) {
    for (const auto &text : lobstersText) {
      lobsters.emplace_back(v2::SmilesParse::MolFromSmiles(text));
      REQUIRE(lobsters.back());
    }
  }
}

TEST_CASE("Tversky") {
  // Score the PDB overlay.
  auto pdb_trp_3tmn =
      R"(N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|)"_smiles;
  REQUIRE(pdb_trp_3tmn);
  auto pdb_0zn_1tmn =
      R"([C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|)"_smiles;
  REQUIRE(pdb_0zn_1tmn);
  GaussianShape::ShapeOverlayOptions ovlyOpts;
  GaussianShape::ShapeInputOptions inOpts;
  auto tan_scores = GaussianShape::ScoreMolecule(*pdb_0zn_1tmn, *pdb_trp_3tmn);
  CHECK_THAT(tan_scores[0], Catch::Matchers::WithinAbs(0.307, 0.001));
  CHECK_THAT(tan_scores[1], Catch::Matchers::WithinAbs(0.349, 0.001));
  CHECK_THAT(tan_scores[2], Catch::Matchers::WithinAbs(0.265, 0.001));

  ovlyOpts.simAlpha = 0.95;
  ovlyOpts.simBeta = 0.05;
  auto ref_tversky = GaussianShape::ScoreMolecule(*pdb_0zn_1tmn, *pdb_trp_3tmn,
                                                  inOpts, inOpts, ovlyOpts);
  CHECK_THAT(ref_tversky[0], Catch::Matchers::WithinAbs(0.362, 0.001));
  CHECK_THAT(ref_tversky[1], Catch::Matchers::WithinAbs(0.383, 0.001));
  CHECK_THAT(ref_tversky[2], Catch::Matchers::WithinAbs(0.342, 0.001));

  ovlyOpts.simAlpha = 0.05;
  ovlyOpts.simBeta = 0.95;
  auto fit_tversky = GaussianShape::ScoreMolecule(*pdb_0zn_1tmn, *pdb_trp_3tmn,
                                                  inOpts, inOpts, ovlyOpts);
  CHECK_THAT(fit_tversky[0], Catch::Matchers::WithinAbs(0.668, 0.001));
  CHECK_THAT(fit_tversky[1], Catch::Matchers::WithinAbs(0.795, 0.001));
  CHECK_THAT(fit_tversky[2], Catch::Matchers::WithinAbs(0.540, 0.001));
}

#ifdef RDK_USE_BOOST_SERIALIZATION
TEST_CASE("Serialization") {
  auto m1 =
      "[H]c1c([H])c([H])c([H])c([H])c1[H] |(-2.06264,-0.844763,-0.0261403;-1.04035,-0.481453,-0.0114878;-0.00743655,-1.41861,-0.0137121;-0.215455,-2.47997,-0.0295909;1.29853,-0.949412,0.00507497;2.12524,-1.65277,0.00390664;1.58501,0.395878,0.0254188;2.61997,0.704365,0.0394811;0.550242,1.31385,0.0273741;0.783172,2.37039,0.0434262;-0.763786,0.88847,0.00908113;-1.60557,1.58532,0.0100194)|"_smiles;
  REQUIRE(m1);
  GaussianShape::ShapeInputOptions shapeOpts;
  shapeOpts.allCarbonRadii = false;
  auto shape = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  auto istr = shape.toString();

  GaussianShape::ShapeInput shape2(istr);
  CHECK(shape2.getCoords() == shape.getCoords());
  CHECK(shape2.getFeatureTypes() == shape.getFeatureTypes());
  CHECK(shape2.getNumAtoms() == shape.getNumAtoms());
  CHECK(shape2.getNumFeatures() == shape.getNumFeatures());
  CHECK(shape2.getIsNormalized() == shape.getIsNormalized());
  CHECK(shape2.calcExtremes() == shape.calcExtremes());
  CHECK(shape2.calcCanonicalRotation() == shape.calcCanonicalRotation());
  CHECK(shape2.calcCanonicalTranslation() == shape.calcCanonicalTranslation());
  CHECK(*shape2.getCarbonRadii() == *shape.getCarbonRadii());
  CHECK_THAT(shape2.getShapeVolume(),
             Catch::Matchers::WithinAbs(261.0145, 0.005));
  CHECK_THAT(shape2.getColorVolume(), Catch::Matchers::WithinAbs(5.316, 0.005));

  // Check it handles the case of no d_carbonRadii in the ShapeInput.
  shapeOpts.allCarbonRadii = true;
  auto shape3 = GaussianShape::ShapeInput(*m1, -1, shapeOpts);
  auto istr2 = shape3.toString();
  GaussianShape::ShapeInput shape4(istr2);
  CHECK(!shape4.getCarbonRadii());
}
#endif

#ifdef RDK_TEST_MULTITHREADED
#include <thread>
#include <future>

namespace {
void runblock(
    const std::vector<std::pair<std::shared_ptr<RWMol>, std::shared_ptr<RWMol>>>
        &pairs,
    unsigned int count, unsigned int idx,
    std::vector<std::array<double, 3>> &test) {
  for (unsigned int i = idx; i < pairs.size(); i += count) {
    auto p1 = *pairs[i].first;
    auto p2 = *pairs[i].second;
    test[i] = GaussianShape::AlignMolecule(p1, p2);
  }
}
}  // namespace

TEST_CASE("multithreaded") {
  initLobsters();

  constexpr size_t numRepeats = 1000;
  std::vector<std::pair<std::shared_ptr<RWMol>, std::shared_ptr<RWMol>>> pairs;
  for (auto r = 0u; r < numRepeats; ++r) {
    for (unsigned int i = 0; i < 10; i += 2) {
      unsigned int l1 = i;
      unsigned int l2 = i + 1;
      pairs.emplace_back(lobsters[l1], lobsters[l2]);
    }
  }
  // generate reference data
  std::cerr << " generating reference data" << std::endl;
  auto start = std::chrono::steady_clock::now();
  std::vector<std::array<double, 3>> ref;
  for (auto pr : pairs) {
    auto p1 = *pr.first;
    auto p2 = *pr.second;
    auto norm_scores = GaussianShape::AlignMolecule(p1, p2);

    ref.push_back(norm_scores);
  }
  auto end = std::chrono::steady_clock::now();
  auto ref_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  std::cerr << " reference time: " << ref_time << "ms" << std::endl;

  // run the same calculations in parallel and check they match the reference
  std::cerr << " parallel loop" << std::endl;
  std::vector<std::array<double, 3>> test(pairs.size());
  std::vector<unsigned int> idx(pairs.size());
  std::iota(idx.begin(), idx.end(), 0);
  auto start2 = std::chrono::steady_clock::now();
  std::vector<std::future<void>> tg;
  unsigned int count = 4;
  for (unsigned int i = 0; i < count; ++i) {
    tg.emplace_back(std::async(std::launch::async, runblock, pairs, count, i,
                               std::ref(test)));
  }
  for (auto &fut : tg) {
    fut.get();
  }
  tg.clear();
  auto end2 = std::chrono::steady_clock::now();
  auto test_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2)
          .count();
  std::cerr << " parallel time: " << test_time << "ms" << std::endl;
  CHECK(test == ref);
}
#endif

std::unique_ptr<ROMol> loadConformers(const std::string &fileName) {
  std::ifstream ifs(fileName.c_str());
  std::string nextLine;
  std::unique_ptr<ROMol> retMol;
  while (!ifs.eof()) {
    std::getline(ifs, nextLine);
    if (!retMol) {
      retMol = v2::SmilesParse::MolFromSmiles(nextLine);
    } else {
      auto m = v2::SmilesParse::MolFromSmiles(nextLine);
      if (!m || !m->getNumConformers()) {
        continue;
      }
      retMol->addConformer(new Conformer(m->getConformer()), true);
    }
  }
  return retMol;
}

TEST_CASE("Multiple Conformers") {
  std::string dirName = getenv("RDBASE");
  dirName += "/Code/GraphMol/GaussianShape/test_data";
  auto esomeprazole = loadConformers(dirName + "/esomeprazole_multi.smi");
  CHECK(esomeprazole->getNumConformers() == 10);
  {
    GaussianShape::ShapeInput shape1(*esomeprazole);
    CHECK(shape1.getNumShapes() == 10);
    shape1.setActiveShape(0);
    auto firstVol = shape1.getShapeVolume() + shape1.getColorVolume();
    shape1.setActiveShape(9);
    auto lastVol = shape1.getShapeVolume() + shape1.getColorVolume();
    CHECK(firstVol > lastVol);

    GaussianShape::ShapeInputOptions shapeOptions;
    shapeOptions.shapePruneThreshold = 0.75;
    GaussianShape::ShapeInput shape2(*esomeprazole, -1, shapeOptions);
    CHECK(shape2.getNumShapes() == 6);
    shape2.setActiveShape(0);
    firstVol = shape2.getShapeVolume() + shape2.getColorVolume();
    shape2.setActiveShape(5);
    lastVol = shape2.getShapeVolume() + shape2.getColorVolume();
    CHECK(firstVol > lastVol);

    RDLog::InitLogs();
    {
      RDLog::CaptureLog capture{rdErrorLog};
      CHECK_THROWS(shape2.setActiveShape(8));
    }

    // The shapes in 2 are a subset of 1, so the maximum similarity
    // should be 1.0.
    auto maxSim = shape1.maxPossibleSimilarity(shape2);
    CHECK_THAT(maxSim, Catch::Matchers::WithinAbs(1.0, 0.001));
    unsigned int bestFitShape, bestRefShape;
    RDGeom::Transform3D bestXform;
    auto maxSimBF =
        shape1.bestSimilarity(shape2, bestRefShape, bestFitShape, bestXform);
    CHECK_THAT(maxSimBF[0], Catch::Matchers::WithinAbs(1.0, 0.001));
    CHECK(bestRefShape == 2);
    CHECK(bestFitShape == 0);
  }

  {
    // Single conformations
    GaussianShape::ShapeInput shape1(*esomeprazole, 1);
    CHECK(shape1.getNumShapes() == 1);
    GaussianShape::ShapeInput shape2(*esomeprazole, 8);
    CHECK(shape2.getNumShapes() == 1);
    // Demonstrate that the shapes are different:
    auto scores = GaussianShape::AlignShape(shape1, shape2);
    CHECK(scores[0] < 0.5);
  }

  {
    auto ranit = loadConformers(dirName + "/ranitidine_multi.smi");
    CHECK(ranit->getNumConformers() == 10);

    int best1, best2;
    std::vector<std::vector<double>> sims;
    RDGeom::Transform3D bestXform;
    GaussianShape::ScoreMoleculeAllConformers(
        *esomeprazole, *ranit, best1, best2, sims,
        GaussianShape::ShapeInputOptions(), GaussianShape::ShapeInputOptions(),
        GaussianShape::ShapeOverlayOptions(), &bestXform);
    auto bestRanit = ROMol(*ranit, false, best2);
    MolTransforms::transformConformer(bestRanit.getConformer(), bestXform);
    CHECK(best1 == 8);
    CHECK(best2 == 3);
    CHECK_THAT(sims[best1][best2], Catch::Matchers::WithinAbs(0.456, 0.005));
    auto bestRanitOvly =
        "CN/C(=C\\[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1 |(6.68398,0.986506,0.966042;5.9518,-0.234997,0.864737;4.78568,-0.258931,0.014363;4.85802,-0.593057,-1.26255;6.14357,-0.936503,-1.79214;6.50827,-0.586214,-2.94026;7.04388,-1.67335,-1.05581;3.56197,0.104652,0.658858;2.49061,-0.837246,0.905698;1.38577,-0.062151,1.62039;0.821805,1.2952,0.576773;-0.00932388,0.695297,-0.906402;-1.27675,0.0251991,-0.460442;-1.46957,-1.30004,-0.125472;-2.81328,-1.42941,0.220382;-3.35978,-0.17268,0.0752515;-4.75996,0.282489,0.306315;-5.50725,0.159366,-0.930455;-6.8579,0.615173,-0.676725;-4.87106,0.868053,-2.01392;-2.41969,0.660958,-0.328873)|"_smiles;
    CHECK(checkMolsHaveRoughlySameCoords(bestRanit, *bestRanitOvly, 0.02));
  }
}

namespace {
bool checkBondLengths(const ROMol &mol) {
  // DetermineBonds::connectivityVdw uses a covalent factor of 1.3.
  static constexpr double radFactor = 1.3;
  const auto conf = mol.getConformer();
  for (const auto bond : mol.bonds()) {
    if (!bond->getBeginAtom()->getAtomicNum() ||
        !bond->getEndAtom()->getAtomicNum()) {
      continue;
    }
    auto bondlen = MolTransforms::getBondLength(conf, bond->getBeginAtomIdx(),
                                                bond->getEndAtomIdx());
    auto rad1 = PeriodicTable::getTable()->getRcovalent(
        bond->getBeginAtom()->getAtomicNum());
    auto rad2 = PeriodicTable::getTable()->getRcovalent(
        bond->getEndAtom()->getAtomicNum());
    if (bondlen > radFactor * (rad1 + rad2)) {
      std::cout << bond->getIdx() << " : " << bond->getBeginAtomIdx() << " -> "
                << bond->getEndAtomIdx() << " len = " << bondlen << " vs "
                << radFactor * (rad1 + rad2) << std::endl;
      return false;
    }
  }
  return true;
}
}  // namespace

TEST_CASE("Different atom orders for ShapeInput") {
  // Make sure that different atom orders always produce a ShapeInput that gives
  // a correct molecule from shapeToMol.  This wasn't always the case.
  auto fullMol =
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(-0.346914,-0.986206,-4.28744;-0.686863,-0.0357247,-3.13265;0.429505,-0.1946,-2.14134;0.21099,0.659676,-0.907145;1.06526,0.0812473,0.104663;2.29297,0.75201,0.373712;2.5837,1.80373,-0.246936;3.23325,0.27478,1.33777;4.47647,0.99197,1.57602;4.94347,1.01294,2.99117;5.59284,-0.21541,2.82613;5.71049,0.107583,1.47157;-1.19052,0.721766,-0.417623;-2.14086,-0.0964663,-1.12904;-3.25312,-0.745252,-0.540367;-3.95877,-1.43507,-1.38825;-3.7227,-0.763533,0.803759;-4.9481,-0.204581,1.05395;-5.52492,-0.18107,2.24654;-4.86554,-0.748644,3.30451;-3.63585,-1.32731,3.13734;-3.08234,-1.32608,1.89052;-1.84839,-1.9059,1.69441;-2.02702,-0.329978,-2.57998),wD:1.0,wU:3.3|"_smiles;
  REQUIRE(fullMol);
  CHECK(checkBondLengths(*fullMol));
  std::vector<unsigned int> atomOrder(fullMol->getNumAtoms());
  std::iota(atomOrder.begin(), atomOrder.end(), 0);
  auto rng = std::default_random_engine{};
  for (unsigned int i = 0; i < 100; ++i) {
    std::ranges::shuffle(atomOrder, rng);
    std::unique_ptr<ROMol> renumMol(MolOps::renumberAtoms(*fullMol, atomOrder));
    CHECK(checkBondLengths(*renumMol));
    GaussianShape::ShapeInput shape(*renumMol);
    auto outMol = shape.shapeToMol(false);
    CHECK(checkBondLengths(*outMol));
  }

  // And the same for the shape from a subset.
  GaussianShape::ShapeInputOptions shapeOptions;
  auto bitToGo = "c1nccnc1F"_smarts;
  REQUIRE(bitToGo);
  for (unsigned int i = 0; i < 100; ++i) {
    std::ranges::shuffle(atomOrder, rng);
    std::unique_ptr<ROMol> renumMol(MolOps::renumberAtoms(*fullMol, atomOrder));
    CHECK(checkBondLengths(*renumMol));
    auto match = SubstructMatch(*renumMol, *bitToGo);
    boost::dynamic_bitset<> toGo(fullMol->getNumAtoms());
    for (auto mp : match.front()) {
      toGo[mp.second] = true;
    }
    shapeOptions.atomSubset.clear();
    shapeOptions.atomSubset.reserve(renumMol->getNumAtoms());
    for (unsigned int j = 0; j < renumMol->getNumAtoms(); ++j) {
      if (!toGo[j]) {
        shapeOptions.atomSubset.push_back(j);
      }
    }
    GaussianShape::ShapeInput shape(*renumMol, -1, shapeOptions);
    auto outMol = shape.shapeToMol(false);
    CHECK(MolToSmiles(*outMol) == "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C=O)C1");
    CHECK(checkBondLengths(*outMol));
  }
}

TEST_CASE("github #9441: problem with explicit Hs") {
  SECTION("as reported") {
    auto m1 =
        "[H]/N=C(/N)N[C@H]1Cc2ccccc2[C@@H]1NC(=O)C(=O)Nc1ccc(Cl)c(F)c1 |(2.1315,-15.781,8.6689;1.793,-16.412,9.367;0.597,-16.134,10.042;0.18,-17.047,10.999;-0.038,-15.085,9.775;-1.364,-14.651,10.415;-1.866,-13.442,9.745;-2.902,-13.917,8.784;-3.483,-13.185,7.737;-4.426,-13.794,6.949;-4.784,-15.091,7.193;-4.19,-15.827,8.264;-3.257,-15.232,9.049;-2.512,-15.785,10.184;-3.359,-15.978,11.459;-2.895,-16.13,12.693;-1.719,-16.487,12.907;-3.826,-15.86,13.835;-5.064,-15.502,13.519;-3.381,-15.995,15.067;-4.055,-15.773,16.236;-3.356,-15.637,17.417;-4.06,-15.374,18.687;-5.428,-15.251,18.717;-6.28,-14.934,20.209;-6.143,-15.369,17.54;-7.452,-15.254,17.561;-5.486,-15.623,16.293),wD:5.4,wU:13.15|"_smiles;
    auto m2 =
        "[H]/N=C(/N)N[C@@H]1Cc2ccccc2[C@H]1NC(=O)C(=O)Nc1ccc(Cl)c(F)c1 |(-1.8622,-20.4007,9.0038;-2.6909,-20.1336,8.5119;-3.0037,-18.871,8.3972;-4.0884,-18.5145,7.7538;-2.1651,-17.8781,8.9694;-2.5069,-16.4989,8.8473;-1.8532,-15.9197,7.6961;-1.6088,-14.4969,8.06;-1.3546,-13.3514,7.2821;-1.1595,-12.1122,7.9047;-1.2189,-12.0095,9.3012;-1.4695,-13.1438,10.0742;-1.6627,-14.393,9.4426;-1.9498,-15.7546,9.9967;-2.9045,-15.6225,10.9981;-2.5492,-15.8844,12.3606;-1.3948,-16.2409,12.6505;-3.6003,-15.728,13.4623;-4.7286,-15.3802,13.1585;-3.2054,-16.0128,14.8184;-4.0446,-15.9028,15.9995;-3.4411,-16.0498,17.2382;-4.1812,-15.9386,18.3961;-5.5426,-15.6693,18.3166;-6.4762,-15.526,19.7442;-6.1509,-15.5114,17.0857;-7.4532,-15.2539,17.0236;-5.4088,-15.6235,15.9289),wU:5.4,13.15|"_smiles;
    auto scores = GaussianShape::AlignMolecule(*m1, *m2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.4223, 0.0001));
  }
  SECTION("Make sure [2H] survives") {
    auto m1 =
        "[2H]/N=C(/N)N[C@H]1Cc2ccccc2[C@@H]1NC(=O)C(=O)Nc1ccc(Cl)c(F)c1 |(2.1315,-15.781,8.6689;1.793,-16.412,9.367;0.597,-16.134,10.042;0.18,-17.047,10.999;-0.038,-15.085,9.775;-1.364,-14.651,10.415;-1.866,-13.442,9.745;-2.902,-13.917,8.784;-3.483,-13.185,7.737;-4.426,-13.794,6.949;-4.784,-15.091,7.193;-4.19,-15.827,8.264;-3.257,-15.232,9.049;-2.512,-15.785,10.184;-3.359,-15.978,11.459;-2.895,-16.13,12.693;-1.719,-16.487,12.907;-3.826,-15.86,13.835;-5.064,-15.502,13.519;-3.381,-15.995,15.067;-4.055,-15.773,16.236;-3.356,-15.637,17.417;-4.06,-15.374,18.687;-5.428,-15.251,18.717;-6.28,-14.934,20.209;-6.143,-15.369,17.54;-7.452,-15.254,17.561;-5.486,-15.623,16.293),wD:5.4,wU:13.15|"_smiles;
    auto m2 =
        "[2H]/N=C(/N)N[C@@H]1Cc2ccccc2[C@H]1NC(=O)C(=O)Nc1ccc(Cl)c(F)c1 |(-1.8622,-20.4007,9.0038;-2.6909,-20.1336,8.5119;-3.0037,-18.871,8.3972;-4.0884,-18.5145,7.7538;-2.1651,-17.8781,8.9694;-2.5069,-16.4989,8.8473;-1.8532,-15.9197,7.6961;-1.6088,-14.4969,8.06;-1.3546,-13.3514,7.2821;-1.1595,-12.1122,7.9047;-1.2189,-12.0095,9.3012;-1.4695,-13.1438,10.0742;-1.6627,-14.393,9.4426;-1.9498,-15.7546,9.9967;-2.9045,-15.6225,10.9981;-2.5492,-15.8844,12.3606;-1.3948,-16.2409,12.6505;-3.6003,-15.728,13.4623;-4.7286,-15.3802,13.1585;-3.2054,-16.0128,14.8184;-4.0446,-15.9028,15.9995;-3.4411,-16.0498,17.2382;-4.1812,-15.9386,18.3961;-5.5426,-15.6693,18.3166;-6.4762,-15.526,19.7442;-6.1509,-15.5114,17.0857;-7.4532,-15.2539,17.0236;-5.4088,-15.6235,15.9289),wU:5.4,13.15|"_smiles;
    auto scores = GaussianShape::AlignMolecule(*m1, *m2);
    CHECK_THAT(scores[0], Catch::Matchers::WithinAbs(0.4188, 0.0001));
  }
}
TEST_CASE("color gradients use their own converters") {
  // The color converters are stored after the numFitShape shape points'
  // ones; where they are stored mustn't change the gradients.
  auto colorGradients = [](int numFitShape) {
    const double ref[3] = {0.0, 0.0, 0.0}, fit[3] = {1.0, 0.5, -0.3};
    const double alphas[1] = {1.0};
    const int types[1] = {1};
    const double quat[4] = {0.9, 0.1, 0.3, 0.2};
    std::vector<double> gradConverters(12 * (numFitShape + 1), 0.0);
    std::array<double, 7> grads{};
    GaussianShape::calcVolAndGrads(ref, alphas, 1, types, fit, alphas, 1, types,
                                   numFitShape, gradConverters, false, 0.0,
                                   quat, grads.data());
    return grads;
  };
  const auto atStart = colorGradients(0);
  const auto afterShape = colorGradients(2);
  for (unsigned int g = 1; g < atStart.size(); ++g) {
    CHECK(atStart[g] != 0.0);
    CHECK(afterShape[g] == atStart[g]);
  }
}
TEST_CASE("color gradients are mixed into the matching components") {
  // One atom and one color feature in each of ref and fit.
  const std::vector<double> ref{0.0, 0.0, 0.0, 0.2, -0.4, 0.1};
  const std::vector<double> fit{1.0, 0.5, -0.3, -0.6, 0.8, 0.7};
  const std::vector<double> alphas{0.8, 1.0};
  const int types[2] = {0, 1};
  const double mixingParam = 0.5;
  GaussianShape::SingleConformerAlignment sca(
      ref, alphas, types, nullptr, 1, 1, 1.0, 1.0, fit, alphas, types, nullptr,
      1, 1, 1.0, 1.0, {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      GaussianShape::OptimMode::SHAPE_PLUS_COLOR, 1.0, 0.0, mixingParam, false,
      0.0, 0.001, 100);
  const std::array<double, 7> quatTrans{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double shapeVol, colorVol;
  std::array<double, 7> mixed{};
  sca.calcVolumeAndGradients(quatTrans, shapeVol, colorVol, mixed);

  // The same thing computed piecewise.
  std::vector<double> gradConverters(24, 0.0);
  std::array<double, 7> shapeGrads{}, colorGrads{};
  GaussianShape::calcVolAndGrads(
      ref.data(), alphas.data(), 1, nullptr, fit.data(), alphas.data(), 1,
      nullptr, gradConverters, false, 0.0, quatTrans.data(), shapeGrads.data());
  GaussianShape::calcVolAndGrads(ref.data() + 3, alphas.data() + 1, 1,
                                 types + 1, fit.data() + 3, alphas.data() + 1,
                                 1, types + 1, 1, gradConverters, false, 0.0,
                                 quatTrans.data(), colorGrads.data());
  double shapeSum = 0.0, colorSum = 0.0;
  for (unsigned int g = 1; g < 7; ++g) {
    shapeSum += shapeGrads[g] * shapeGrads[g];
    colorSum += colorGrads[g] * colorGrads[g];
  }
  const auto ratio = std::sqrt(shapeSum / colorSum);
  for (unsigned int g = 1; g < 7; ++g) {
    const auto expected =
        shapeGrads[g] * (1 - mixingParam) + colorGrads[g] * ratio * mixingParam;
    CHECK_THAT(mixed[g], Catch::Matchers::WithinAbs(expected, 1.0e-10));
  }
}
TEST_CASE("rotational gradients are correct away from the identity") {
  // The optimiser applies a step as quatTrans * step, so the gradients must
  // be those of V(quatTrans * step) with respect to the step.
  using QuatTrans = std::array<double, 7>;
  // A copy of combineQuatTrans, which isn't exposed.
  auto combine = [](const QuatTrans &a, const QuatTrans &b) {
    return QuatTrans{a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                     a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                     a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                     a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
                     a[4] + b[4],
                     a[5] + b[5],
                     a[6] + b[6]};
  };
  // A small asymmetric shape, used as both ref and fit.
  const std::vector<double> pts{0.0, 0.0,  0.0, 1.5, 0.0, 0.0,  2.1, 1.3,
                                0.0, -0.8, 1.1, 0.4, 0.3, -0.9, 1.2};
  const std::vector<double> alphas(5, 0.8);
  const std::vector<int> types(5, 0);
  GaussianShape::SingleConformerAlignment sca(
      pts, alphas, types.data(), nullptr, 5, 0, 1.0, 0.0, pts, alphas,
      types.data(), nullptr, 5, 0, 1.0, 0.0,
      {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, GaussianShape::OptimMode::SHAPE_ONLY,
      1.0, 1.0, 0.5, false, 0.0, 0.001, 100);
  auto volume = [&](const QuatTrans &quatTrans) {
    double shapeVol, colorVol;
    QuatTrans grads;
    sca.calcVolumeAndGradients(quatTrans, shapeVol, colorVol, grads);
    return shapeVol;
  };
  const double h = 1.0e-5;
  for (const auto angle : {0.3, 1.0}) {
    // Rotate about (1, 2, 3) and translate a little.
    const auto s = std::sin(angle / 2) / std::sqrt(14.0);
    const QuatTrans quatTrans{
        std::cos(angle / 2), s, 2 * s, 3 * s, 0.2, -0.1, 0.3};
    double shapeVol, colorVol;
    QuatTrans grads;
    sca.calcVolumeAndGradients(quatTrans, shapeVol, colorVol, grads);
    for (unsigned int g = 1; g < 7; ++g) {
      QuatTrans stepPlus{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      if (g < 4) {
        stepPlus[0] = std::sqrt(1.0 - h * h);
      }
      auto stepMinus = stepPlus;
      stepPlus[g] = h;
      stepMinus[g] = -h;
      const auto fd = (volume(combine(quatTrans, stepPlus)) -
                       volume(combine(quatTrans, stepMinus))) /
                      (2 * h);
      INFO("angle " << angle << " gradient " << g);
      CHECK_THAT(grads[g], Catch::Matchers::WithinAbs(-fd, 1.0e-6));
    }
  }
}
