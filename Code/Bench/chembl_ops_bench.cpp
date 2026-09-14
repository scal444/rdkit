// Copyright (c) 2026, RDKit contributors
// All rights reserved.
//
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <DataStructs/BitOps.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#ifdef RDK_BUILD_INCHI_SUPPORT
#include <INCHI-API/inchi.h>
#endif

namespace {

using Molecules = std::vector<RDKit::ROMol>;
using Fingerprints = std::vector<std::unique_ptr<ExplicitBitVect>>;

struct Corpus {
  Molecules molecules;
  std::size_t read = 0;
  std::size_t failures = 0;
  std::size_t atoms = 0;
};

struct OperationResult {
  std::uint64_t checksum = 0;
  std::size_t attempts = 0;
  std::size_t failures = 0;
};

enum class Operation {
  CanonicalSmiles,
  Morgan,
  PainsSubstructure,
  Pickle,
  Etkdg,
  BulkSimilarity,
  InchiRoundtrip,
  HsRoundtrip,
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
  if (text == "pains_substructure") {
    return Operation::PainsSubstructure;
  }
  if (text == "pickle") {
    return Operation::Pickle;
  }
  if (text == "etkdg") {
    return Operation::Etkdg;
  }
  if (text == "bulk_similarity") {
    return Operation::BulkSimilarity;
  }
  if (text == "inchi_roundtrip") {
    return Operation::InchiRoundtrip;
  }
  if (text == "hs_roundtrip") {
    return Operation::HsRoundtrip;
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
    const auto separator = line.find_first_of("\t ");
    const auto smiles = line.substr(0, separator);
    if (corpus.read == 0 && smiles == "smiles") {
      continue;
    }
    ++corpus.read;
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

OperationResult benchmarkCanonicalSmiles(const Molecules &molecules) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    const auto smiles = RDKit::MolToSmiles(molecule);
    result.checksum += smiles.size();
    ++result.attempts;
  }
  return result;
}

OperationResult benchmarkMorgan(
    const Molecules &molecules,
    const RDKit::FingerprintGenerator<std::uint64_t> &generator) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    std::unique_ptr<ExplicitBitVect> fingerprint(
        generator.getFingerprint(molecule));
    result.checksum += fingerprint->getNumOnBits();
    ++result.attempts;
  }
  return result;
}

OperationResult benchmarkPainsSubstructure(
    const Molecules &molecules, const RDKit::FilterCatalog &catalog) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    result.checksum += catalog.hasMatch(molecule);
    ++result.attempts;
  }
  return result;
}

OperationResult benchmarkPickle(const Molecules &molecules) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    std::string pickle;
    RDKit::MolPickler::pickleMol(molecule, pickle);
    RDKit::ROMol restored(pickle);
    result.checksum += pickle.size() + restored.getNumAtoms();
    ++result.attempts;
  }
  return result;
}

Fingerprints makeFingerprints(
    const Molecules &molecules,
    const RDKit::FingerprintGenerator<std::uint64_t> &generator) {
  Fingerprints fingerprints;
  fingerprints.reserve(molecules.size());
  for (const auto &molecule : molecules) {
    fingerprints.emplace_back(generator.getFingerprint(molecule));
  }
  return fingerprints;
}

OperationResult benchmarkEtkdg(const Molecules &molecules) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    ++result.attempts;
    std::unique_ptr<RDKit::ROMol> withHs(RDKit::MolOps::addHs(molecule));
    auto parameters = RDKit::DGeomHelpers::ETKDGv3;
    parameters.numThreads = 1;
    parameters.randomSeed = 0xC0FFEE;
    const auto conformerId =
        RDKit::DGeomHelpers::EmbedMolecule(*withHs, parameters);
    if (conformerId < 0) {
      ++result.checksum;
      ++result.failures;
      continue;
    }
    std::unique_ptr<RDKit::ROMol> withoutHs(RDKit::MolOps::removeHs(*withHs));
    result.checksum +=
        2 + withoutHs->getNumAtoms() + withoutHs->getNumConformers();
  }
  return result;
}

OperationResult benchmarkBulkSimilarity(const Fingerprints &fingerprints) {
  constexpr std::size_t queryLimit = 100;
  const auto queryCount = std::min(queryLimit, fingerprints.size());
  OperationResult result;
  for (std::size_t query = 0; query < queryCount; ++query) {
    const auto queryIndex = query * fingerprints.size() / queryCount;
    for (const auto &candidate : fingerprints) {
      if (TanimotoSimilarity(*fingerprints[queryIndex], *candidate) >= 0.4) {
        ++result.checksum;
      }
      ++result.attempts;
    }
  }
  return result;
}

OperationResult benchmarkInchiRoundtrip(const Molecules &molecules) {
#ifdef RDK_BUILD_INCHI_SUPPORT
  OperationResult result;
  RDLog::LogStateSetter blockLogs;
  for (const auto &molecule : molecules) {
    ++result.attempts;
    try {
      RDKit::ExtraInchiReturnValues toInchiResult;
      const auto inchi = RDKit::MolToInchi(molecule, toInchiResult);
      RDKit::ExtraInchiReturnValues fromInchiResult;
      std::unique_ptr<RDKit::ROMol> restored(
          RDKit::InchiToMol(inchi, fromInchiResult));
      result.checksum += inchi.size();
      if (restored) {
        result.checksum += restored->getNumAtoms();
      } else {
        ++result.checksum;
        ++result.failures;
      }
    } catch (const std::exception &) {
      ++result.checksum;
      ++result.failures;
    }
  }
  return result;
#else
  throw std::runtime_error("InChI support is not enabled in this build");
#endif
}

OperationResult benchmarkHsRoundtrip(const Molecules &molecules) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    std::unique_ptr<RDKit::ROMol> withHs(RDKit::MolOps::addHs(molecule));
    std::unique_ptr<RDKit::ROMol> restored(RDKit::MolOps::removeHs(*withHs));
    result.checksum += withHs->getNumAtoms() + restored->getNumAtoms();
    ++result.attempts;
  }
  return result;
}

OperationResult runOperation(
    Operation operation, const Molecules &molecules,
    const RDKit::FilterCatalog *painsCatalog,
    const RDKit::FingerprintGenerator<std::uint64_t> *generator,
    const Fingerprints &fingerprints) {
  switch (operation) {
    case Operation::CanonicalSmiles:
      return benchmarkCanonicalSmiles(molecules);
    case Operation::Morgan:
      if (!generator) {
        throw std::logic_error("Morgan generator is not initialized");
      }
      return benchmarkMorgan(molecules, *generator);
    case Operation::PainsSubstructure:
      if (!painsCatalog) {
        throw std::logic_error("PAINS catalog is not initialized");
      }
      return benchmarkPainsSubstructure(molecules, *painsCatalog);
    case Operation::Pickle:
      return benchmarkPickle(molecules);
    case Operation::Etkdg:
      return benchmarkEtkdg(molecules);
    case Operation::BulkSimilarity:
      return benchmarkBulkSimilarity(fingerprints);
    case Operation::InchiRoundtrip:
      return benchmarkInchiRoundtrip(molecules);
    case Operation::HsRoundtrip:
      return benchmarkHsRoundtrip(molecules);
  }
  throw std::logic_error("unhandled benchmark operation");
}

void printUsage(const char *program) {
  std::cerr
      << "Usage: " << program << " OPERATION INPUT COUNT REPEATS\n"
      << "Operations: canonical_smiles, morgan, pains_substructure, pickle, "
         "etkdg, bulk_similarity, inchi_roundtrip, hs_roundtrip\n";
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
    std::unique_ptr<RDKit::FilterCatalog> painsCatalog;
    if (operation == Operation::PainsSubstructure) {
      painsCatalog = std::make_unique<RDKit::FilterCatalog>(
          RDKit::FilterCatalogParams::PAINS);
    }
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> generator;
    if (operation == Operation::Morgan ||
        operation == Operation::BulkSimilarity) {
      generator.reset(
          RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    }
    Fingerprints fingerprints;
    if (operation == Operation::BulkSimilarity) {
      fingerprints = makeFingerprints(corpus.molecules, *generator);
    }
    const auto warmup =
        runOperation(operation, corpus.molecules, painsCatalog.get(),
                     generator.get(), fingerprints);
    const auto setupElapsed = std::chrono::steady_clock::now() - setupStart;

    OperationResult result;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
      const auto iteration =
          runOperation(operation, corpus.molecules, painsCatalog.get(),
                       generator.get(), fingerprints);
      result.checksum += iteration.checksum;
      result.attempts += iteration.attempts;
      result.failures += iteration.failures;
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;

    std::cout << "operation=" << argv[1] << " input=" << argv[2]
              << " requested=" << requested << " read=" << corpus.read
              << " parsed=" << corpus.molecules.size()
              << " failures=" << corpus.failures << " atoms=" << corpus.atoms
              << " repeats=" << repeats
              << " warmup_checksum=" << warmup.checksum
              << " warmup_operation_failures=" << warmup.failures
              << " checksum=" << result.checksum
              << " operation_attempts=" << result.attempts
              << " operation_failures=" << result.failures << " setup_s="
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
