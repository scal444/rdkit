#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <GraphMol/ROMol.h>

namespace bench_common {

constexpr const char *SAMPLES[] = {
    // highly fused/bridged ring system
    "O=C1N2[C@H]3N4CN5C(=O)N6C7C5N(C2)C(=O)N7CN2C(=O)N5CN1[C@@H]3N(CN1C5C2N(C1=O)C6)C4=O",

    // multiple stereo groups
    "OC[C@]12[C@](O)(CC[C@H]3[C@]4(O)[C@@](C)([C@@H](C5COC(=O)C5)CC4)C[C@@H](O)[C@H]13)C[C@@H](O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C[C@H]2O",

    // randomly selected
    "C[C@H](NS(/C=C/c1ccccc1)(=O)=O)C(OCC(N1CCC(C)CC1)=O)=O",
    "COc1ccc(-n2cc(CNC(C)c3ncncc3)c(-c3cc(F)ccc3)n2)cc1",
    "O=C1N=C([n+]2ccccc2)/C(=C/[O-])N1c1ccccc1",
    "c1csc(CNC(=O)C2CCC(CNS(c3ccccc3)(=O)=O)CC2)c1",
    "CC1=C2C(c3ccc(C)cc3)C3=C(N=C2NN1)CC(C)(C)CC3=O",
    "NC(N)=NN/C=C1\\C=CC=C([N+]([O-])=O)C1=O",
    "CCc1c2c(oc(=O)c1)c(CN1CCN(C)CC1)c(O)c(Cl)c2",
    "COc1ccc(C2CC(C(F)(F)F)N3NC=C(C(Nc4c(C)n(C)n(-c5ccccc5)c4=O)=O)C3=N2)cc1",
    "Cc1ccc(CC(=O)O/N=C(\\N)Cc2ccc([N+]([O-])=O)cc2)cc1",
    "Br.COc1ccc(/N=C/C=C2/OC(C)(C)OC(c3ccccc3)=C2)cc1",
};

//! Dataset selector for stress / size-scan / ring-count benches.
//!
//! `Canonical` corresponds to the in-header SAMPLES array. The remaining
//! entries are loaded lazily from `Code/Bench/data/<name>.smi` via the
//! `dataset_smiles()` accessor (RDBASE-relative, cached per process).
enum class Dataset {
  Canonical,

  // Heavy-atom size scan
  Size_00_20,
  Size_20_40,
  Size_40_60,
  Size_60_80,

  // SSSR ring-count scan
  Rings_2,
  Rings_3,
  Rings_4,
  Rings_5,
  Rings_6,

  // Edge-case datasets that target inner branches not exercised by SAMPLES
  Radicals,
  Organometallics,
  HypervalentHalogens,
  HypervalentP,
  PreCanonicalNO2Azide,
  Atropisomers,
  KekulizeHard,
};

//! Short stable name for a dataset; used as the .smi filename stem and as a
//! tag suffix in TEST_CASE titles. Do not include spaces or punctuation.
const char *dataset_name(Dataset dataset);

//! SMILES contents of a dataset. For `Dataset::Canonical` this is the
//! in-header SAMPLES list. For all other datasets the result is loaded
//! once per process from `$RDBASE/Code/Bench/data/<name>.smi` and cached.
//!
//! REQUIRE-fails if the .smi file is missing or unreadable so a
//! mis-deployed bench binary surfaces the problem immediately.
const std::vector<std::string> &dataset_smiles(Dataset dataset);

std::vector<RDKit::ROMol> load_samples();

//! Per-dataset loaders. Default parameters (sanitize=true, removeHs=true)
//! work for canonical/size/ring/radicals/kekulize_hard/atropisomers; the
//! pre-canonical edge-case datasets (organometallics, hypervalent_*,
//! pre_canonical_no2_azide) are intentionally not sanitize-clean and the
//! caller should pass `sanitize=false` to keep the inputs in their
//! pre-cleanUp form.
std::vector<RDKit::ROMol> load_samples(Dataset dataset, bool sanitize = true);

constexpr uint64_t nth_random(uint64_t n) noexcept {
  // https://xoshiro.di.unimi.it/splitmix64.c
  // inlined here becuse <boost/random/splitmix64.hpp> is not always available
  std::uint64_t seed = 1;
  std::uint64_t z = seed + (n + 1) * 0x9E3779B97F4A7C15;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EB;
  return z ^ (z >> 31);
}

static_assert(0x910a2dec89025cc1 == nth_random(0));
static_assert(0xbeeb8da1658eec67 == nth_random(1));
static_assert(0xf893a2eefb32555e == nth_random(2));
static_assert(0x71c18690ee42c90b == nth_random(3));
static_assert(0x71bb54d8d101b5b9 == nth_random(4));
static_assert(0x7760003b54a685ae == nth_random(1000));
static_assert(0x28a1928f0674d152 == nth_random(1000000000000));
static_assert(0x5692161d100b05e5 ==
              nth_random(std::numeric_limits<uint64_t>::max()));

}  // namespace bench_common
