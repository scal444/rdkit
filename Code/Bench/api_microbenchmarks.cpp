// ROMol-only counterpart of Code/Bench/api_microbenchmarks.cpp from the
// rdmol_benchmarks branch. Each TEST_CASE name matches its rdmol_benchmarks
// counterpart so cross-branch diffs line up cleanly. The RDMol-side test
// cases from the rdmol_benchmarks file are intentionally omitted here
// (RDMol does not exist on master).
//
// Workload patterns mirror the rdmol_benchmarks file exactly: same
// per-bucket dataset coverage, same shuffled index sequence for random
// access, same property names, same minimal-op iteration body.

#include <catch2/catch_all.hpp>

#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RingInfo.h>

#include "bench_common.hpp"

using namespace RDKit;

namespace {

using bench_common::Dataset;
using bench_common::load_samples;
using bench_common::nth_random;

template <class Mol, class Op>
auto run_per_sample_mutating(const std::vector<Mol> &samples,
                             Catch::Benchmark::Chronometer meter, Op op) {
  std::vector<Mol> work;
  work.reserve(meter.runs() * samples.size());
  for (int run = 0; run < meter.runs(); ++run) {
    for (const auto &mol : samples) {
      work.emplace_back(mol);
    }
  }
  meter.measure([&](int run) {
    uint64_t total = 0;
    for (size_t sample = 0; sample < samples.size(); ++sample) {
      auto &mol = work[run * samples.size() + sample];
      total += op(mol);
    }
    return total;
  });
}

template <class Mol, class Op>
auto run_per_sample_readonly(const std::vector<Mol> &samples,
                             Catch::Benchmark::Chronometer meter, Op op) {
  meter.measure([&](int /*run*/) {
    uint64_t total = 0;
    for (const auto &mol : samples) {
      total += op(mol);
    }
    return total;
  });
}

template <class Mol, class Prime>
std::vector<Mol> prime_samples(std::vector<Mol> samples, Prime prime) {
  for (auto &mol : samples) {
    prime(mol);
  }
  return samples;
}

std::vector<uint32_t> shuffled_indices(uint32_t n) {
  std::vector<uint32_t> indices(n);
  std::iota(indices.begin(), indices.end(), 0u);
  for (uint32_t i = 0; i + 1 < n; ++i) {
    uint32_t j = i + uint32_t(nth_random(i) % (n - i));
    std::swap(indices[i], indices[j]);
  }
  return indices;
}

template <class Mol>
std::vector<std::vector<uint32_t>> per_sample_shuffled_atom_indices(
    const std::vector<Mol> &samples) {
  std::vector<std::vector<uint32_t>> result;
  result.reserve(samples.size());
  for (const auto &mol : samples) {
    result.emplace_back(shuffled_indices(mol.getNumAtoms()));
  }
  return result;
}

// Build RWMol copies so mutating ops can run.  Mirrors how the
// rdmol_benchmarks ROMol leg constructs its working set.
std::vector<RWMol> load_rwmol_samples(Dataset dataset) {
  auto base = load_samples(dataset);
  std::vector<RWMol> out;
  out.reserve(base.size());
  for (auto &mol : base) {
    out.emplace_back(mol);
  }
  return out;
}

}  // namespace

#define BENCH_MOVE_CTOR(DATASET, SUFFIX, TAG)                                  \
  TEST_CASE("ROMol move constructor " SUFFIX, "[mol_api]" TAG) {               \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol move constructor " SUFFIX)                       \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<ROMol> sources;                                              \
      sources.reserve(meter.runs() * samples.size());                          \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
        }                                                                      \
      }                                                                        \
      std::vector<Catch::Benchmark::storage_for<ROMol>> storage(               \
          meter.runs() * samples.size());                                      \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[run * samples.size() + sample].construct(                    \
              std::move(sources[run * samples.size() + sample]));              \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_MOVE_ASSIGN(DATASET, SUFFIX, TAG)                                \
  TEST_CASE("ROMol move assign " SUFFIX, "[mol_api]" TAG) {                    \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol move assign " SUFFIX)                            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<ROMol> sources;                                              \
      std::vector<ROMol> dests;                                                \
      sources.reserve(meter.runs() * samples.size());                          \
      dests.reserve(meter.runs() * samples.size());                            \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
          dests.emplace_back();                                                \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          dests[run * samples.size() + sample] =                               \
              std::move(sources[run * samples.size() + sample]);               \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ATOM_ITER_SUM_ATOMICNUM(DATASET, SUFFIX, TAG)                    \
  TEST_CASE("ROMol atoms() sum atomicNum " SUFFIX, "[mol_api]" TAG) {          \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atoms() sum atomicNum " SUFFIX)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getAtomicNum();                                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_BOND_ITER_SUM_TWICEORDER(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol bonds() sum twiceBondType " SUFFIX, "[mol_api]" TAG) {      \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol bonds() sum twiceBondType " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto bond : mol.bonds()) {                                        \
          total += static_cast<uint32_t>(2 * bond->getBondTypeAsDouble());     \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_NEIGHBOR_WALK_SUM_DEGREE(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol atomNeighbors sum degree " SUFFIX, "[mol_api]" TAG) {       \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atomNeighbors sum degree " SUFFIX)               \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          for (auto nbr : mol.atomNeighbors(atom)) {                           \
            (void)nbr;                                                         \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_BONDS_FROM_ATOM_SUM_ORDER(DATASET, SUFFIX, TAG)                  \
  TEST_CASE("ROMol atomBonds sum twiceBondType " SUFFIX, "[mol_api]" TAG) {    \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atomBonds sum twiceBondType " SUFFIX)            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          for (auto bond : mol.atomBonds(atom)) {                              \
            total += static_cast<uint32_t>(2 * bond->getBondTypeAsDouble());   \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_SET_VARIED(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol setProp int per-atom varied " SUFFIX, "[mol_api]" TAG) {    \
    auto samples = load_rwmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("ROMol setProp int per-atom varied " SUFFIX)            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [](RWMol &mol) {                 \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.getAtomWithIdx(i)->setProp("bench_int", int(i));                 \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_SET_UNIFORM(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol setProp int per-atom uniform " SUFFIX, "[mol_api]" TAG) {   \
    auto samples = load_rwmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("ROMol setProp int per-atom uniform " SUFFIX)           \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [](RWMol &mol) {                 \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.getAtomWithIdx(i)->setProp("bench_uniform", int(42));            \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_SET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                     \
  TEST_CASE("ROMol setProp int random-access " SUFFIX, "[mol_api]" TAG) {      \
    auto raw = load_rwmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [](RWMol &mol) {              \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        mol.getAtomWithIdx(i)->setProp("bench_int", int(0));                   \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("ROMol setProp int random-access " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = const_cast<RWMol &>(samples[s]);                         \
          for (uint32_t i : orderings[s]) {                                    \
            mol.getAtomWithIdx(i)->setProp("bench_int", int(i + 1));           \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_GET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                     \
  TEST_CASE("ROMol getProp int random-access " SUFFIX, "[mol_api]" TAG) {      \
    auto raw = load_rwmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [](RWMol &mol) {              \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        mol.getAtomWithIdx(i)->setProp("bench_int", int(i));                   \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("ROMol getProp int random-access " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          const auto &mol = samples[s];                                        \
          for (uint32_t i : orderings[s]) {                                    \
            total += uint32_t(mol.getAtomWithIdx(i)->getProp<int>(             \
                "bench_int"));                                                 \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_HAS_PRESENT(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol hasProp present " SUFFIX, "[mol_api]" TAG) {                \
    auto raw = load_rwmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [](RWMol &mol) {              \
      for (auto atom : mol.atoms()) {                                          \
        atom->setProp("bench_int", int(0));                                    \
      }                                                                        \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol hasProp present " SUFFIX)                        \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RWMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->hasProp("bench_int") ? 1u : 0u;                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_HAS_ABSENT(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol hasProp absent " SUFFIX, "[mol_api]" TAG) {                 \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol hasProp absent " SUFFIX)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->hasProp("bench_absent") ? 1u : 0u;                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_CLEAR_BULK(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol clearProp per-atom " SUFFIX, "[mol_api]" TAG) {             \
    auto raw = load_rwmol_samples(DATASET);                                    \
    BENCHMARK_ADVANCED("ROMol clearProp per-atom " SUFFIX)                     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RWMol> work;                                                 \
      work.reserve(meter.runs() * raw.size());                                 \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (auto &mol : raw) {                                                \
          work.emplace_back(mol);                                              \
          for (auto atom : work.back().atoms()) {                              \
            atom->setProp("bench_int", int(0));                                \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < raw.size(); ++s) {                              \
          auto &mol = work[run * raw.size() + s];                              \
          for (auto atom : mol.atoms()) {                                      \
            atom->clearProp("bench_int");                                      \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ACCESSOR_DEGREE(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol getDegree all atoms " SUFFIX, "[mol_api]" TAG) {            \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol getDegree all atoms " SUFFIX)                    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getDegree();                                          \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ACCESSOR_NUM_IMPLICIT_HS(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol getNumImplicitHs all atoms " SUFFIX, "[mol_api]" TAG) {     \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol getNumImplicitHs all atoms " SUFFIX)             \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getNumImplicitHs();                                   \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_RING_NUM_ATOM_RINGS(DATASET, SUFFIX, TAG)                        \
  TEST_CASE("ROMol numAtomRings all atoms " SUFFIX, "[mol_api]" TAG) {         \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      MolOps::findSSSR(mol);                                                   \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol numAtomRings all atoms " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        const auto *info = mol.getRingInfo();                                  \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += info->numAtomRings(atom->getIdx());                         \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_RING_NUM_BOND_RINGS(DATASET, SUFFIX, TAG)                        \
  TEST_CASE("ROMol numBondRings all bonds " SUFFIX, "[mol_api]" TAG) {         \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      MolOps::findSSSR(mol);                                                   \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol numBondRings all bonds " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        const auto *info = mol.getRingInfo();                                  \
        uint32_t total = 0;                                                    \
        for (auto bond : mol.bonds()) {                                        \
          total += info->numBondRings(bond->getIdx());                         \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_API_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_MOVE_CTOR(DATASET, SUFFIX, TAG)                                        \
  BENCH_MOVE_ASSIGN(DATASET, SUFFIX, TAG)                                      \
  BENCH_ATOM_ITER_SUM_ATOMICNUM(DATASET, SUFFIX, TAG)                          \
  BENCH_BOND_ITER_SUM_TWICEORDER(DATASET, SUFFIX, TAG)                         \
  BENCH_NEIGHBOR_WALK_SUM_DEGREE(DATASET, SUFFIX, TAG)                         \
  BENCH_BONDS_FROM_ATOM_SUM_ORDER(DATASET, SUFFIX, TAG)                        \
  BENCH_PROP_SET_VARIED(DATASET, SUFFIX, TAG)                                  \
  BENCH_PROP_SET_UNIFORM(DATASET, SUFFIX, TAG)                                 \
  BENCH_PROP_SET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                           \
  BENCH_PROP_GET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                           \
  BENCH_PROP_HAS_PRESENT(DATASET, SUFFIX, TAG)                                 \
  BENCH_PROP_HAS_ABSENT(DATASET, SUFFIX, TAG)                                  \
  BENCH_PROP_CLEAR_BULK(DATASET, SUFFIX, TAG)                                  \
  BENCH_ACCESSOR_DEGREE(DATASET, SUFFIX, TAG)                                  \
  BENCH_ACCESSOR_NUM_IMPLICIT_HS(DATASET, SUFFIX, TAG)                         \
  BENCH_RING_NUM_ATOM_RINGS(DATASET, SUFFIX, TAG)                              \
  BENCH_RING_NUM_BOND_RINGS(DATASET, SUFFIX, TAG)

BENCH_API_FOR(Dataset::Canonical, "", "[canonical]")
BENCH_API_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_API_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_API_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_API_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_API_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
