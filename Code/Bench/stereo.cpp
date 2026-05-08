#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/test_fixtures.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("Chirality::findPotentialStereo", "[stereo]") {
  auto samples = bench_common::load_samples();

  BENCHMARK("Chirality::findPotentialStereo") {
    auto total = 0;

    for (auto &mol : samples) {
      auto stereo_infos = Chirality::findPotentialStereo(mol);

      // workaround for https://github.com/rdkit/rdkit/issues/8880
      mol.clearComputedProps();

      for (auto &info : stereo_infos) {
        total += info.controllingAtoms.size();
      }
    }

    return total;
  };
}

TEST_CASE("CIPLabeler::assignCIPLabels", "[stereo]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("CIPLabeler::assignCIPLabels") {
    for (auto &mol : samples) {
      CIPLabeler::assignCIPLabels(mol);
    }
  };
}

TEST_CASE("MolOps::clearSingleBondDirFlags", "[stereo]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::clearSingleBondDirFlags")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<ROMol> work;
    work.reserve(meter.runs() * samples.size());
    for (int run = 0; run < meter.runs(); ++run) {
      for (const auto &mol : samples) {
        work.emplace_back(mol);
      }
    }
    meter.measure([&](int i) {
      uint64_t total = 0;
      for (size_t s = 0; s < samples.size(); ++s) {
        auto &mol = work[i * samples.size() + s];
        MolOps::clearSingleBondDirFlags(mol);
        total += mol.getNumBonds();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::assignStereochemistry", "[stereo]") {
  const auto cleanIt = true;
  const auto force = true;
  const auto flagPossibleStereoCenters = true;

  auto samples = bench_common::load_samples();

  const auto legacy = GENERATE(true, false);
  UseLegacyStereoPerceptionFixture fx(legacy);
  auto str_legacy = std::string(legacy ? "true" : "false");

  BENCHMARK("MolOps::assignStereochemistry legacy=" + str_legacy) {
    auto total = 0;

    for (auto &mol : samples) {
      MolOps::assignStereochemistry(mol, cleanIt, force,
                                    flagPossibleStereoCenters);
      for (auto &atom : mol.atoms()) {
        total += atom->getChiralTag();
      }

      // workaround for https://github.com/rdkit/rdkit/issues/8880
      mol.clearComputedProps();
    }

    return total;
  };
}

TEST_CASE("Canon::rankMolAtoms", "[stereo]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("Canon::rankMolAtoms") {
    uint64_t total = 0;
    std::vector<unsigned int> ranks;
    for (auto &mol : samples) {
      Canon::rankMolAtoms(mol, ranks);
      for (auto rank : ranks) {
        total += rank;
      }
    }
    return total;
  };
}

#define BENCH_RANK_MOL_ATOMS(DATASET, SUFFIX, TAG)                             \
  TEST_CASE("Canon::rankMolAtoms " SUFFIX, "[stereo]" TAG) {                    \
    auto samples = bench_common::load_samples(DATASET);                         \
    BENCHMARK("Canon::rankMolAtoms " SUFFIX) {                                  \
      uint64_t total = 0;                                                       \
      std::vector<unsigned int> ranks;                                          \
      for (auto &mol : samples) {                                               \
        Canon::rankMolAtoms(mol, ranks);                                        \
        for (auto rank : ranks) {                                               \
          total += rank;                                                        \
        }                                                                       \
      }                                                                         \
      return total;                                                             \
    };                                                                          \
  }

BENCH_RANK_MOL_ATOMS(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_RANK_MOL_ATOMS(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_RANK_MOL_ATOMS(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_RANK_MOL_ATOMS(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_RANK_MOL_ATOMS(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_RANK_MOL_ATOMS(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_RANK_MOL_ATOMS(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_RANK_MOL_ATOMS(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_RANK_MOL_ATOMS(Dataset::Rings_6, "rings 6", "[rings_6]")
