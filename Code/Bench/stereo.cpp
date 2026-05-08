#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/test_fixtures.h>

using namespace RDKit;

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
