// Copyright (c) 2026, RDKit contributors
// All rights reserved.
//
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/new_canon.h>

namespace {

using Molecules = std::vector<RDKit::ROMol>;
using Queries = std::vector<std::unique_ptr<RDKit::ROMol>>;

struct Corpus {
  Molecules molecules;
  std::size_t read = 0;
  std::size_t failures = 0;
  std::size_t atoms = 0;
};

enum class Operation {
  CanonicalSmiles,
  Morgan,
  Substructure,
  Descriptors,
  Pickle,
  RankAtoms,
  Copy,
};

constexpr std::array<std::string_view, 12> kQueries = {
    "c1ccccc1",         "[CX3](=O)[OX2H1]", "[NX3;H2,H1;!$(NC=O)]",
    "[OX2H]",           "[F,Cl,Br,I]",      "[#7]~[#6]~[#8]",
    "[R]@[R]",          "[$([N;H1](C)C)]",  "[S;X4](=O)(=O)",
    "[C;H3,H2]-[O;H1]", "n1ccccc1",         "C=C-C=C",
};

std::size_t parsePositiveInteger(const char *text, const char *name) {
  try {
    std::size_t consumed = 0;
    const auto value = std::stoull(text, &consumed);
    if (consumed != std::string(text).size() || value == 0) {
      throw std::invalid_argument("not a positive integer");
    }
    return value;
  } catch (const std::exception &) {
    throw std::invalid_argument(std::string("invalid ") + name + ": " + text);
  }
}

Operation parseOperation(std::string_view text) {
  if (text == "canonical_smiles") {
    return Operation::CanonicalSmiles;
  }
  if (text == "morgan") {
    return Operation::Morgan;
  }
  if (text == "substructure") {
    return Operation::Substructure;
  }
  if (text == "descriptors") {
    return Operation::Descriptors;
  }
  if (text == "pickle") {
    return Operation::Pickle;
  }
  if (text == "rank_atoms") {
    return Operation::RankAtoms;
  }
  if (text == "copy") {
    return Operation::Copy;
  }
  throw std::invalid_argument("unknown operation: " + std::string(text));
}

Corpus loadCorpus(const std::string &path, std::size_t requested) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open input: " + path);
  }

  Corpus corpus;
  corpus.molecules.reserve(requested);
  std::string line;
  while (corpus.read < requested && std::getline(input, line)) {
    ++corpus.read;
    const auto separator = line.find_first_of("\t ");
    const auto smiles = line.substr(0, separator);
    if (smiles.empty()) {
      ++corpus.failures;
      continue;
    }
    try {
      auto molecule = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
      if (!molecule) {
        ++corpus.failures;
        continue;
      }
      corpus.atoms += molecule->getNumAtoms();
      corpus.molecules.push_back(std::move(*molecule));
    } catch (const std::exception &) {
      ++corpus.failures;
    }
  }
  if (corpus.read != requested) {
    throw std::runtime_error("input ended before the requested molecule count");
  }
  if (corpus.molecules.empty()) {
    throw std::runtime_error("input produced no molecules");
  }
  return corpus;
}

Queries makeQueries() {
  Queries queries;
  queries.reserve(kQueries.size());
  for (const auto smarts : kQueries) {
    auto query = RDKit::v2::SmilesParse::MolFromSmarts(std::string(smarts));
    if (!query) {
      throw std::runtime_error("failed to parse benchmark SMARTS");
    }
    queries.push_back(std::move(query));
  }
  return queries;
}

std::uint64_t benchmarkCanonicalSmiles(const Molecules &molecules) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    const auto smiles = RDKit::MolToSmiles(molecule);
    checksum += smiles.size();
  }
  return checksum;
}

std::uint64_t benchmarkMorgan(
    const Molecules &molecules,
    const RDKit::FingerprintGenerator<std::uint64_t> &generator) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    std::unique_ptr<ExplicitBitVect> fingerprint(
        generator.getFingerprint(molecule));
    checksum += fingerprint->getNumOnBits();
  }
  return checksum;
}

std::uint64_t benchmarkSubstructure(const Molecules &molecules,
                                    const Queries &queries) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    for (const auto &query : queries) {
      checksum += RDKit::SubstructMatch(molecule, *query).size();
    }
  }
  return checksum;
}

std::uint64_t benchmarkDescriptors(const Molecules &molecules) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    checksum += RDKit::Descriptors::calcNumSpiroAtoms(molecule);
    checksum += RDKit::Descriptors::calcNumBridgeheadAtoms(molecule);
  }
  return checksum;
}

std::uint64_t benchmarkPickle(const Molecules &molecules) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    std::string pickle;
    RDKit::MolPickler::pickleMol(molecule, pickle);
    RDKit::ROMol restored(pickle);
    checksum += pickle.size() + restored.getNumAtoms();
  }
  return checksum;
}

std::uint64_t benchmarkRankAtoms(const Molecules &molecules) {
  std::uint64_t checksum = 0;
  std::vector<unsigned int> ranks;
  for (const auto &molecule : molecules) {
    RDKit::Canon::rankMolAtoms(molecule, ranks);
    for (const auto rank : ranks) {
      checksum += rank;
    }
  }
  return checksum;
}

std::uint64_t benchmarkCopy(const Molecules &molecules) {
  std::uint64_t checksum = 0;
  for (const auto &molecule : molecules) {
    RDKit::ROMol copy(molecule);
    checksum += copy.getNumAtoms() + copy.getNumBonds();
  }
  return checksum;
}

std::uint64_t runOperation(
    Operation operation, const Molecules &molecules, const Queries &queries,
    const RDKit::FingerprintGenerator<std::uint64_t> &generator) {
  switch (operation) {
    case Operation::CanonicalSmiles:
      return benchmarkCanonicalSmiles(molecules);
    case Operation::Morgan:
      return benchmarkMorgan(molecules, generator);
    case Operation::Substructure:
      return benchmarkSubstructure(molecules, queries);
    case Operation::Descriptors:
      return benchmarkDescriptors(molecules);
    case Operation::Pickle:
      return benchmarkPickle(molecules);
    case Operation::RankAtoms:
      return benchmarkRankAtoms(molecules);
    case Operation::Copy:
      return benchmarkCopy(molecules);
  }
  throw std::logic_error("unhandled benchmark operation");
}

void printUsage(const char *program) {
  std::cerr
      << "Usage: " << program << " OPERATION INPUT COUNT REPEATS\n"
      << "Operations: canonical_smiles, morgan, substructure, descriptors, "
         "pickle, rank_atoms, copy\n";
}

}  // namespace

int main(int argc, char *argv[]) {
  if (argc != 5) {
    printUsage(argv[0]);
    return 2;
  }

  try {
    const auto operation = parseOperation(argv[1]);
    const auto requested = parsePositiveInteger(argv[3], "count");
    const auto repeats = parsePositiveInteger(argv[4], "repeat count");

    const auto setupStart = std::chrono::steady_clock::now();
    auto corpus = loadCorpus(argv[2], requested);
    auto queries = makeQueries();
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> generator(
        RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    const auto warmupChecksum =
        runOperation(operation, corpus.molecules, queries, *generator);
    const auto setupElapsed = std::chrono::steady_clock::now() - setupStart;

    std::uint64_t checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
      checksum +=
          runOperation(operation, corpus.molecules, queries, *generator);
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;

    std::cout << "operation=" << argv[1] << " input=" << argv[2]
              << " requested=" << requested << " read=" << corpus.read
              << " parsed=" << corpus.molecules.size()
              << " failures=" << corpus.failures << " atoms=" << corpus.atoms
              << " repeats=" << repeats << " warmup_checksum=" << warmupChecksum
              << " checksum=" << checksum << " setup_s="
              << std::chrono::duration<double>(setupElapsed).count()
              << " elapsed_s=" << std::chrono::duration<double>(elapsed).count()
              << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    printUsage(argv[0]);
    return 2;
  }
}
