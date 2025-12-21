//
//  Copyright (C) 2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <fstream>

using namespace RDKit;
using ::RDKit::QueryOps::validateAllQueries;
using ::RDKit::v2::SmilesParse::MolFromSmarts;
TEST_CASE(
    "github #6106: Dummy atoms should not be considered to be metals for M and MH queries") {
  const auto m = "C*[Fe]"_smiles;
  REQUIRE(m);

  SECTION("M") {
    std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    REQUIRE(q);
    CHECK(!q->Match(m->getAtomWithIdx(0)));
    CHECK(!q->Match(m->getAtomWithIdx(1)));
    CHECK(q->Match(m->getAtomWithIdx(2)));
  }
  SECTION("MH") {
    std::unique_ptr<ATOM_OR_QUERY> q(makeMHAtomQuery());
    REQUIRE(q);
    CHECK(!q->Match(m->getAtomWithIdx(0)));
    CHECK(!q->Match(m->getAtomWithIdx(1)));
    CHECK(q->Match(m->getAtomWithIdx(2)));
  }
}


TEST_CASE(
  "Validate query atoms pass case"
  ) {


  SECTION("Same type different value fails") {
    const std::vector<std::string> basicPatterns = {
      "[C;c]",
      "[#7;#6]"
      "[O-;O-2]",
    };
    for (const auto& pattern : basicPatterns) {
      CAPTURE(pattern);
      const auto query = MolFromSmarts(pattern);
      CHECK(!validateAllQueries(*query.get()));
    }

    const std::vector<std::string> negatePatterns = {
      "[C;!C]",
      "[#7;!#7]"
    };
    for (const auto& pattern : negatePatterns) {
      CAPTURE(pattern);
      const auto query = MolFromSmarts(pattern);
      CHECK(!validateAllQueries(*query.get()));
    }
  }

  SECTION("Handles negation") {
    const std::vector<std::string> patterns = {
      "[C;!c]",
      "[#6;!#7]",
      "[!#1;!#5]"
    };
    for (const auto& pattern : patterns) {
      CAPTURE(pattern);
      const auto query = MolFromSmarts(pattern);
      CHECK(validateAllQueries(*query.get()));
    }
  }

  SECTION("Nesting passes regardless of validity") {
    const std::string pattern = "[C;-1,-2]";
    CAPTURE(pattern);
    const auto query = MolFromSmarts(pattern);
    REQUIRE(query);
    CHECK(validateAllQueries(*query.get()));
  }

  SECTION("MultiQuery all passing") {
    const std::string pattern = "[C;-1;!#7]";
    CAPTURE(pattern);
    const auto query = MolFromSmarts(pattern);
    REQUIRE(query);
    CHECK(validateAllQueries(*query.get()));
  }

  SECTION("MultiQuery some failing fails all") {
    const std::string pattern ="[C;#6][O-;O-2]";
    CAPTURE(pattern);
    const auto query = MolFromSmarts(pattern);
    REQUIRE(query);
    CHECK(!validateAllQueries(*query.get()));
  }
}

void validate(const std::string& inputFileName) {

  std::ifstream inputFile(inputFileName);
  REQUIRE(inputFile.is_open());
  std::string line;
  std::vector<std::string> smarts;
  std::vector<std::unique_ptr<RWMol>> queries;
  while (std::getline(inputFile, line)) {
    if (line.starts_with(' ') || line.starts_with('#')) {
      continue;
    }

    auto query = MolFromSmarts(line);
    if (!query) {
      printf("Invalid SMARTS: %s\n", line.c_str());
      continue;
    }
    smarts.push_back(line);
    queries.push_back(std::move(query));
  }

  for (int i = 0; i < queries.size(); i++) {
    if (!validateAllQueries(*queries[i].get())) {
      printf("Bad SMARTS: %s\n", smarts[i].c_str());
    }
  }
}

TEST_CASE("Validator") {
  const std::vector<std::string> patterns = {
    "/home/kevin/data/smarts/crippen.txt",
    "/home/kevin/data/smarts/functional_group_hierarchy.txt",
    "/home/kevin/data/smarts/functionalGroups.txt",
    "/home/kevin/data/smarts/patty_rules.txt",

    "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/BMS_2006_filter.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/openbabel_functional_groups.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/pwalters_alert_collection.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_fragment_descriptors.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_pattern_fingerprint.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_tautomer_transforms.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_torsionPreferences_macrocycles.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_torsionPreferences_smallrings.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/rdkit_torsionPreferences_v2.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/RLewis_smarts.txt",
  "/home/kevin/repos/nvmolkit/tests/test_data/SMARTS/wehi_pains.txt"
  };

  for (const auto& pattern : patterns) {
    printf("--------------------------\n validating %s\n\n", pattern.c_str());
    validate(pattern);
  }
}
