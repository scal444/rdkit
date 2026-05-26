#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>

using namespace RDKit;

using bench_common::Dataset;

namespace {

template <class Mol, class Op>
auto run_per_sample(const std::vector<Mol> &samples,
                    Catch::Benchmark::Chronometer meter, Op op) {
  std::vector<Mol> work;
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
      total += op(mol);
    }
    return total;
  });
}

// Sanitize step index for primed sub-step benches.  Mirrors the enum on the
// rdmol_benchmarks branch so the same per-stage names land in the parser
// output for both legs.
enum class SanitizeStage {
  cleanUp = 0,
  cleanUpOrganometallics,
  updatePropertyCache_first,
  symmetrizeSSSR,
  Kekulize,
  assignRadicals,
  setAromaticity,
  setConjugation,
  setHybridization,
  cleanupAtropisomers,
  cleanupChirality,
  adjustHs,
};

void prime_to(RWMol &mol, SanitizeStage stage) {
  using S = SanitizeStage;
  if (stage <= S::cleanUp) return;
  MolOps::cleanUp(mol);
  if (stage <= S::cleanUpOrganometallics) return;
  MolOps::cleanUpOrganometallics(mol);
  if (stage <= S::updatePropertyCache_first) return;
  mol.updatePropertyCache(true);
  if (stage <= S::symmetrizeSSSR) return;
  (void)MolOps::symmetrizeSSSR(mol);
  if (stage <= S::Kekulize) return;
  MolOps::Kekulize(mol, /*markAtomsBonds=*/true, /*canonical=*/false);
  if (stage <= S::assignRadicals) return;
  MolOps::assignRadicals(mol);
  if (stage <= S::setAromaticity) return;
  (void)MolOps::setAromaticity(mol);
  if (stage <= S::setConjugation) return;
  MolOps::setConjugation(mol);
  if (stage <= S::setHybridization) return;
  MolOps::setHybridization(mol);
  if (stage <= S::cleanupAtropisomers) return;
  MolOps::cleanupAtropisomers(mol);
  if (stage <= S::cleanupChirality) return;
  MolOps::cleanupChirality(mol);
  if (stage <= S::adjustHs) return;
  MolOps::adjustHs(mol);
}

std::vector<RWMol> load_and_prime(Dataset dataset, SanitizeStage stage) {
  auto romol_samples = bench_common::load_samples(dataset, /*sanitize=*/false);
  std::vector<RWMol> primed;
  primed.reserve(romol_samples.size());
  for (auto &mol : romol_samples) {
    primed.emplace_back(mol);
    prime_to(primed.back(), stage);
  }
  return primed;
}

}  // namespace

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

// Primed sub-step benches: each step runs on input matching its inside-
// sanitize precondition.  Replaces the earlier hand-coded canonical TEST_CASEs
// that were running on post-sanitize input (= idempotent no-op shape).
#define BENCH_STAGE_RO(NAME, STAGE, OP, COUNT, DATASET, SUFFIX, TAG)           \
  TEST_CASE("MolOps::" NAME " " SUFFIX, "[molops]" TAG) {                       \
    auto rw_samples = load_and_prime(DATASET, STAGE);                           \
    BENCHMARK_ADVANCED("MolOps::" NAME " " SUFFIX)(                             \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                        \
        OP;                                                                     \
        return COUNT;                                                           \
      });                                                                       \
    };                                                                          \
  }

BENCH_STAGE_RO("cleanUp", SanitizeStage::cleanUp, MolOps::cleanUp(mol),
               mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("cleanUpOrganometallics",
               SanitizeStage::cleanUpOrganometallics,
               MolOps::cleanUpOrganometallics(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("updatePropertyCache_first",
               SanitizeStage::updatePropertyCache_first,
               mol.updatePropertyCache(true), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("symmetrizeSSSR", SanitizeStage::symmetrizeSSSR,
               (void)MolOps::symmetrizeSSSR(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("Kekulize", SanitizeStage::Kekulize,
               MolOps::Kekulize(mol, true, false), mol.getNumBonds(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("assignRadicals", SanitizeStage::assignRadicals,
               MolOps::assignRadicals(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("setAromaticity", SanitizeStage::setAromaticity,
               (void)MolOps::setAromaticity(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("setConjugation", SanitizeStage::setConjugation,
               MolOps::setConjugation(mol), mol.getNumBonds(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("setHybridization", SanitizeStage::setHybridization,
               MolOps::setHybridization(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("cleanupAtropisomers", SanitizeStage::cleanupAtropisomers,
               MolOps::cleanupAtropisomers(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("cleanupChirality", SanitizeStage::cleanupChirality,
               MolOps::cleanupChirality(mol), mol.getNumAtoms(),
               Dataset::Canonical, "", "[canonical]")
BENCH_STAGE_RO("adjustHs", SanitizeStage::adjustHs, MolOps::adjustHs(mol),
               mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")

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

// ---------------------------------------------------------------------------
// Per-bucket size + ring-count variants.  Mirror the macros from the
// rdmol_benchmarks branch so the master baseline measures the same shape.

#define BENCH_ADDHS(DATASET, SUFFIX, TAG)                                      \
  TEST_CASE("MolOps::addHs " SUFFIX, "[molops]" TAG) {                         \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::addHs " SUFFIX) {                                       \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        RWMol mol_copy(mol);                                                   \
        MolOps::addHs(mol_copy);                                               \
        total += mol_copy.getNumAtoms();                                       \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_FINDSSR(DATASET, SUFFIX, TAG)                                    \
  TEST_CASE("MolOps::FindSSR " SUFFIX, "[molops]" TAG) {                       \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::FindSSR " SUFFIX) {                                     \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        total += MolOps::findSSSR(mol);                                        \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_GETMOLFRAGS(DATASET, SUFFIX, TAG)                                \
  TEST_CASE("MolOps::getMolFrags " SUFFIX, "[molops]" TAG) {                   \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::getMolFrags " SUFFIX) {                                 \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        std::vector<std::unique_ptr<ROMol>> frags;                             \
        MolOps::getMolFrags(mol, frags);                                       \
        for (auto &frag : frags) {                                             \
          total += frag->getNumAtoms();                                        \
        }                                                                      \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

// Macro for ROMol-template ops that copy ROMol samples into RWMol so OP can
// mutate; otherwise identical to the rdmol_benchmarks BENCH_RO_ADV.
#define BENCH_RO_ADV(NAME, OP, COUNT, DATASET, SUFFIX, TAG)                    \
  TEST_CASE("MolOps::" NAME " " SUFFIX, "[molops]" TAG) {                      \
    auto romol_samples = bench_common::load_samples(DATASET);                  \
    std::vector<RWMol> rw_samples;                                             \
    rw_samples.reserve(romol_samples.size());                                  \
    for (auto &mol : romol_samples) {                                          \
      rw_samples.emplace_back(mol);                                            \
    }                                                                          \
    BENCHMARK_ADVANCED("MolOps::" NAME " " SUFFIX)(                            \
        Catch::Benchmark::Chronometer meter) {                                 \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                       \
        OP;                                                                    \
        return COUNT;                                                          \
      });                                                                      \
    };                                                                         \
  }

// Per-bucket variants for non-sub-step utilities only.  Sub-step benches
// come from BENCH_STAGE_FOR_BUCKETS below.
#define BENCH_OPS_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_ADDHS(DATASET, SUFFIX, TAG)                                            \
  BENCH_FINDSSR(DATASET, SUFFIX, TAG)                                          \
  BENCH_GETMOLFRAGS(DATASET, SUFFIX, TAG)

BENCH_OPS_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_OPS_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_OPS_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_OPS_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_OPS_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")

// Per-bucket primed sub-step benches.
#define BENCH_STAGE_FOR_BUCKETS(DATASET, SUFFIX, TAG)                          \
  BENCH_STAGE_RO("cleanUp", SanitizeStage::cleanUp, MolOps::cleanUp(mol),      \
                 mol.getNumAtoms(), DATASET, SUFFIX, TAG)                      \
  BENCH_STAGE_RO("cleanUpOrganometallics",                                    \
                 SanitizeStage::cleanUpOrganometallics,                        \
                 MolOps::cleanUpOrganometallics(mol), mol.getNumAtoms(),       \
                 DATASET, SUFFIX, TAG)                                         \
  BENCH_STAGE_RO("updatePropertyCache_first",                                  \
                 SanitizeStage::updatePropertyCache_first,                     \
                 mol.updatePropertyCache(true), mol.getNumAtoms(), DATASET,    \
                 SUFFIX, TAG)                                                  \
  BENCH_STAGE_RO("symmetrizeSSSR", SanitizeStage::symmetrizeSSSR,             \
                 (void)MolOps::symmetrizeSSSR(mol), mol.getNumAtoms(),         \
                 DATASET, SUFFIX, TAG)                                         \
  BENCH_STAGE_RO("Kekulize", SanitizeStage::Kekulize,                         \
                 MolOps::Kekulize(mol, true, false), mol.getNumBonds(),        \
                 DATASET, SUFFIX, TAG)                                         \
  BENCH_STAGE_RO("assignRadicals", SanitizeStage::assignRadicals,             \
                 MolOps::assignRadicals(mol), mol.getNumAtoms(), DATASET,      \
                 SUFFIX, TAG)                                                  \
  BENCH_STAGE_RO("setAromaticity", SanitizeStage::setAromaticity,             \
                 (void)MolOps::setAromaticity(mol), mol.getNumAtoms(),         \
                 DATASET, SUFFIX, TAG)                                         \
  BENCH_STAGE_RO("setConjugation", SanitizeStage::setConjugation,             \
                 MolOps::setConjugation(mol), mol.getNumBonds(), DATASET,      \
                 SUFFIX, TAG)                                                  \
  BENCH_STAGE_RO("setHybridization", SanitizeStage::setHybridization,         \
                 MolOps::setHybridization(mol), mol.getNumAtoms(), DATASET,    \
                 SUFFIX, TAG)                                                  \
  BENCH_STAGE_RO("cleanupAtropisomers", SanitizeStage::cleanupAtropisomers,   \
                 MolOps::cleanupAtropisomers(mol), mol.getNumAtoms(),          \
                 DATASET, SUFFIX, TAG)                                         \
  BENCH_STAGE_RO("cleanupChirality", SanitizeStage::cleanupChirality,         \
                 MolOps::cleanupChirality(mol), mol.getNumAtoms(), DATASET,    \
                 SUFFIX, TAG)                                                  \
  BENCH_STAGE_RO("adjustHs", SanitizeStage::adjustHs, MolOps::adjustHs(mol),  \
                 mol.getNumAtoms(), DATASET, SUFFIX, TAG)

BENCH_STAGE_FOR_BUCKETS(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Rings_4, "rings 4", "[rings_4]")

// ---------------------------------------------------------------------------
// Umbrella sanitize bench.  Reflects what SmilesToMol actually calls on the
// default parse path (removeHs(mol, sanitize=true) internally invokes
// sanitizeMol).  Use no-sanitize parser output as input so the bench input
// matches the pipeline input shape.

TEST_CASE("MolOps::sanitizeMol", "[molops]") {
  auto samples = bench_common::load_samples(Dataset::Canonical,
                                            /*sanitize=*/false);
  std::vector<RWMol> rw_samples;
  rw_samples.reserve(samples.size());
  for (auto &mol : samples) {
    rw_samples.emplace_back(mol);
  }
  BENCHMARK_ADVANCED("MolOps::sanitizeMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::sanitizeMol(mol);
      return mol.getNumAtoms();
    });
  };
}

#define BENCH_SANITIZE(DATASET, SUFFIX, TAG)                                   \
  TEST_CASE("MolOps::sanitizeMol " SUFFIX, "[molops]" TAG) {                    \
    auto samples = bench_common::load_samples(DATASET, /*sanitize=*/false);     \
    std::vector<RWMol> rw_samples;                                              \
    rw_samples.reserve(samples.size());                                         \
    for (auto &mol : samples) {                                                 \
      rw_samples.emplace_back(mol);                                             \
    }                                                                           \
    BENCHMARK_ADVANCED("MolOps::sanitizeMol " SUFFIX)(                          \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                        \
        MolOps::sanitizeMol(mol);                                               \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }

BENCH_SANITIZE(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_SANITIZE(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_SANITIZE(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_SANITIZE(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_SANITIZE(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_SANITIZE(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")

// removeHs core (sanitize=false): H-removal pass without the recursive
// sanitizeMol that SmilesToMol's default path triggers.
#define BENCH_REMOVEHS_CORE(DATASET, SUFFIX, TAG)                              \
  TEST_CASE("MolOps::removeHs core " SUFFIX, "[molops]" TAG) {                  \
    auto romol_samples =                                                        \
        bench_common::load_samples(DATASET, /*sanitize=*/false);                \
    std::vector<RWMol> rw_samples;                                              \
    rw_samples.reserve(romol_samples.size());                                   \
    for (auto &mol : romol_samples) {                                           \
      rw_samples.emplace_back(mol);                                             \
    }                                                                           \
    MolOps::RemoveHsParameters ps;                                              \
    BENCHMARK_ADVANCED("MolOps::removeHs core " SUFFIX)(                        \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [&ps](RWMol &mol) {                     \
        MolOps::removeHs(mol, ps, /*sanitize=*/false);                          \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }

BENCH_REMOVEHS_CORE(Dataset::Canonical, "", "[canonical]")
BENCH_REMOVEHS_CORE(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_REMOVEHS_CORE(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_REMOVEHS_CORE(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_REMOVEHS_CORE(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_REMOVEHS_CORE(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_REMOVEHS_CORE(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")
