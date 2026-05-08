#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("SmilesToMol", "[smiles]") {
  BENCHMARK("SmilesToMol") {
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      total_atoms += mol->getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("SmilesToMol no-sanitize", "[smiles]") {
  v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  BENCHMARK("SmilesToMol no-sanitize") {
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles, params);
      REQUIRE(mol);
      total_atoms += mol->getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("MolToSmiles", "[smiles]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToSmiles") {
    auto total_length = 0;
    for (auto &mol : samples) {
      auto smiles = MolToSmiles(mol);
      total_length += smiles.size();
    }
    return total_length;
  };
}

TEST_CASE("MolToCXSmiles", "[smiles]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToCXSmiles") {
    auto total_length = 0;
    for (auto &mol : samples) {
      auto smiles = MolToCXSmiles(mol);
      total_length += smiles.size();
    }
    return total_length;
  };
}

// ---------------------------------------------------------------------------
// Per-bucket size + ring-count variants for the SMILES round-trip benchmarks.

#define BENCH_SMILES_TO_MOL(DATASET, SUFFIX, TAG)                              \
  TEST_CASE("SmilesToMol " SUFFIX, "[smiles]" TAG) {                           \
    const auto &smiles_list = bench_common::dataset_smiles(DATASET);           \
    BENCHMARK("SmilesToMol " SUFFIX) {                                         \
      auto total_atoms = 0;                                                    \
      for (const auto &smiles : smiles_list) {                                 \
        auto mol = v2::SmilesParse::MolFromSmiles(smiles);                     \
        REQUIRE(mol);                                                          \
        total_atoms += mol->getNumAtoms();                                     \
      }                                                                        \
      return total_atoms;                                                      \
    };                                                                         \
  }

#define BENCH_SMILES_TO_MOL_NOSAN(DATASET, SUFFIX, TAG)                        \
  TEST_CASE("SmilesToMol no-sanitize " SUFFIX, "[smiles]" TAG) {               \
    const auto &smiles_list = bench_common::dataset_smiles(DATASET);           \
    v2::SmilesParse::SmilesParserParams params;                                \
    params.sanitize = false;                                                   \
    params.removeHs = false;                                                   \
    BENCHMARK("SmilesToMol no-sanitize " SUFFIX) {                             \
      auto total_atoms = 0;                                                    \
      for (const auto &smiles : smiles_list) {                                 \
        auto mol = v2::SmilesParse::MolFromSmiles(smiles, params);             \
        REQUIRE(mol);                                                          \
        total_atoms += mol->getNumAtoms();                                     \
      }                                                                        \
      return total_atoms;                                                      \
    };                                                                         \
  }

#define BENCH_MOL_TO_SMILES(DATASET, SUFFIX, TAG)                              \
  TEST_CASE("MolToSmiles " SUFFIX, "[smiles]" TAG) {                           \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolToSmiles " SUFFIX) {                                         \
      auto total_length = 0;                                                   \
      for (auto &mol : samples) {                                              \
        auto smiles = MolToSmiles(mol);                                        \
        total_length += smiles.size();                                         \
      }                                                                        \
      return total_length;                                                     \
    };                                                                         \
  }

#define BENCH_MOL_TO_CXSMILES(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("MolToCXSmiles " SUFFIX, "[smiles]" TAG) {                         \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolToCXSmiles " SUFFIX) {                                       \
      auto total_length = 0;                                                   \
      for (auto &mol : samples) {                                              \
        auto smiles = MolToCXSmiles(mol);                                      \
        total_length += smiles.size();                                         \
      }                                                                        \
      return total_length;                                                     \
    };                                                                         \
  }

#define BENCH_SMILES_FOR(DATASET, SUFFIX, TAG)                                 \
  BENCH_SMILES_TO_MOL(DATASET, SUFFIX, TAG)                                    \
  BENCH_SMILES_TO_MOL_NOSAN(DATASET, SUFFIX, TAG)                              \
  BENCH_MOL_TO_SMILES(DATASET, SUFFIX, TAG)                                    \
  BENCH_MOL_TO_CXSMILES(DATASET, SUFFIX, TAG)

BENCH_SMILES_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_SMILES_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_SMILES_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_SMILES_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_SMILES_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_SMILES_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_SMILES_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_SMILES_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_SMILES_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")

// Edge-case parses (no-sanitize because the inputs are intentionally
// pre-cleanup). Round-trip through MolToSmiles is omitted: it would
// canonicalize away the pre-cleanup form we're trying to bench.
BENCH_SMILES_TO_MOL_NOSAN(Dataset::Organometallics, "organometallics",
                          "[organometallics]")

BENCH_SMILES_TO_MOL_NOSAN(Dataset::PreCanonicalNO2Azide,
                          "pre_canonical_no2_azide",
                          "[pre_canonical_no2_azide]")

// kekulize_hard SMILES are aromatic-canonical, so they round-trip through the
// default parser. Both parse and write benches are interesting since the
// writer re-perceives aromaticity / kekulizes for output.
BENCH_SMILES_TO_MOL(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")
BENCH_MOL_TO_SMILES(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")
BENCH_MOL_TO_CXSMILES(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")

// atropisomers parse via CXSMILES; the bench targets the wedge-bond tracking
// path in the parser/writer.
BENCH_SMILES_TO_MOL(Dataset::Atropisomers, "atropisomers", "[atropisomers]")
BENCH_MOL_TO_CXSMILES(Dataset::Atropisomers, "atropisomers", "[atropisomers]")
