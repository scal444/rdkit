#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

using namespace RDKit;

TEST_CASE("MolOps::addHs", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::addHs") {
    auto total_atoms = 0;
    for (auto &mol : samples) {
      RWMol mol_copy(mol);
      MolOps::addHs(mol_copy);
      total_atoms += mol_copy.getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("MolOps::FindSSR", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::FindSSR") {
    auto total = 0;
    for (auto &mol : samples) {
      total += MolOps::findSSSR(mol);
    }
    return total;
  };
}

TEST_CASE("MolOps::getMolFrags", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::getMolFrags") {
    auto total = 0;
    for (auto &mol : samples) {
      std::vector<std::unique_ptr<ROMol>> frags;
      MolOps::getMolFrags(mol, frags);
      for (auto &frag : frags) {
        total += frag->getNumAtoms();
      }
    }
    return total;
  };
}

TEST_CASE("MolOps::setConjugation", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::setConjugation")(
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
        MolOps::setConjugation(mol);
        total += mol.getNumBonds();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::setHybridization", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::setHybridization")(
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
        MolOps::setHybridization(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::adjustHs", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::adjustHs")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
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
        MolOps::adjustHs(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::assignRadicals", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::assignRadicals")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
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
        MolOps::assignRadicals(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::cleanUp", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::cleanUp")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
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
        MolOps::cleanUp(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::cleanupChirality", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::cleanupChirality")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
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
        MolOps::cleanupChirality(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::Kekulize", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::Kekulize")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
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
        MolOps::Kekulize(mol, /*markAtomsBonds=*/true,
                         /*canonical=*/false);
        total += mol.getNumBonds();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::removeHs", "[molops]") {
  // The default-parsed samples have explicit Hs removed already; build a
  // pool that contains explicit Hs to exercise the function.
  auto baseSamples = bench_common::load_samples();
  std::vector<RWMol> samplesWithHs;
  samplesWithHs.reserve(baseSamples.size());
  for (auto &mol : baseSamples) {
    RWMol withHs(mol);
    MolOps::addHs(withHs);
    samplesWithHs.push_back(std::move(withHs));
  }

  BENCHMARK_ADVANCED("MolOps::removeHs")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> work;
    work.reserve(meter.runs() * samplesWithHs.size());
    for (int run = 0; run < meter.runs(); ++run) {
      for (const auto &mol : samplesWithHs) {
        work.emplace_back(mol);
      }
    }
    meter.measure([&](int i) {
      uint64_t total = 0;
      for (size_t s = 0; s < samplesWithHs.size(); ++s) {
        auto &mol = work[i * samplesWithHs.size() + s];
        MolOps::removeHs(mol);
        total += mol.getNumAtoms();
      }
      return total;
    });
  };
}
