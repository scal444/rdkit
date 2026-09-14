// Copyright (c) 2026, RDKit contributors
// All rights reserved.
//
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

#include <algorithm>
#include <array>
#include <charconv>
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
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
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

struct LoadedMolecules {
  std::vector<RDKit::ROMol> molecules;
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
  Mmff,
  TanimotoSimilarity,
  InchiRoundtrip,
  HsRoundtrip,
};

struct OperationOption {
  std::string_view name;
  Operation operation;
};

constexpr std::array operationOptions{
    OperationOption{"canonical_smiles", Operation::CanonicalSmiles},
    OperationOption{"morgan", Operation::Morgan},
    OperationOption{"pains_substructure", Operation::PainsSubstructure},
    OperationOption{"pickle", Operation::Pickle},
    OperationOption{"etkdg", Operation::Etkdg},
    OperationOption{"mmff", Operation::Mmff},
    OperationOption{"tanimoto_similarity", Operation::TanimotoSimilarity},
    OperationOption{"inchi_roundtrip", Operation::InchiRoundtrip},
    OperationOption{"hs_roundtrip", Operation::HsRoundtrip},
};

std::size_t parsePositiveInteger(std::string_view text, std::string_view name) {
  std::size_t value = 0;
  const auto [end, error] =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (error != std::errc{} || end != text.data() + text.size() || value == 0) {
    throw std::invalid_argument("invalid " + std::string(name) + ": " +
                                std::string(text));
  }
  return value;
}

Operation parseOperation(std::string_view text) {
  const auto option = std::find_if(
      operationOptions.begin(), operationOptions.end(),
      [text](const auto &candidate) { return candidate.name == text; });
  if (option != operationOptions.end()) {
    return option->operation;
  }
  throw std::invalid_argument("unknown operation: " + std::string(text));
}

LoadedMolecules loadMolecules(const std::string &path, std::size_t requested) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open input: " + path);
  }

  // Parse and retain the molecules before timing an operation. This keeps file
  // I/O and SMILES parsing out of benchmarks that are intended to measure a
  // later operation. The small reader also accepts both headerless SMILES files
  // and tables whose first-column heading is "smiles".
  LoadedMolecules loaded;
  loaded.molecules.reserve(requested);
  std::string line;
  while (loaded.read < requested && std::getline(input, line)) {
    const auto separator = line.find_first_of("\t ");
    const auto smiles = line.substr(0, separator);
    if (loaded.read == 0 && smiles == "smiles") {
      continue;
    }
    ++loaded.read;
    if (smiles.empty()) {
      ++loaded.failures;
      continue;
    }
    try {
      // RDKit performs the parsing and sanitization; this function only handles
      // the input rows and keeps the resulting molecules for repeated timings.
      auto molecule = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
      if (!molecule) {
        ++loaded.failures;
        continue;
      }
      loaded.atoms += molecule->getNumAtoms();
      loaded.molecules.push_back(std::move(*molecule));
    } catch (const std::exception &) {
      ++loaded.failures;
    }
  }
  if (loaded.read != requested) {
    throw std::runtime_error("input ended before the requested molecule count");
  }
  if (loaded.molecules.empty()) {
    throw std::runtime_error("input produced no molecules");
  }
  return loaded;
}

OperationResult benchmarkCanonicalSmiles(
    const std::vector<RDKit::ROMol> &molecules) {
  OperationResult result;
  const RDKit::SmilesWriteParams parameters;  // canonical is true by default
  for (const auto &molecule : molecules) {
    const auto smiles = RDKit::MolToSmiles(molecule, parameters);
    result.checksum += smiles.size();
    ++result.attempts;
  }
  return result;
}

OperationResult benchmarkMorgan(
    const std::vector<RDKit::ROMol> &molecules,
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
    const std::vector<RDKit::ROMol> &molecules,
    const RDKit::FilterCatalog &catalog) {
  OperationResult result;
  for (const auto &molecule : molecules) {
    result.checksum += catalog.hasMatch(molecule);
    ++result.attempts;
  }
  return result;
}

OperationResult benchmarkPickle(const std::vector<RDKit::ROMol> &molecules) {
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

std::vector<std::unique_ptr<ExplicitBitVect>> makeFingerprints(
    const std::vector<RDKit::ROMol> &molecules,
    const RDKit::FingerprintGenerator<std::uint64_t> &generator) {
  std::vector<std::unique_ptr<ExplicitBitVect>> fingerprints;
  fingerprints.reserve(molecules.size());
  for (const auto &molecule : molecules) {
    fingerprints.emplace_back(generator.getFingerprint(molecule));
  }
  return fingerprints;
}

OperationResult benchmarkEtkdg(const std::vector<RDKit::ROMol> &molecules) {
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

OperationResult benchmarkMmff(const std::vector<RDKit::ROMol> &molecules,
                              std::size_t inputMoleculeCount) {
  OperationResult result;
  result.attempts = inputMoleculeCount;
  result.failures = inputMoleculeCount - molecules.size();
  for (const auto &molecule : molecules) {
    RDKit::ROMol workingMolecule(molecule);
    const auto optimization =
        RDKit::MMFF::MMFFOptimizeMolecule(workingMolecule, 200);
    if (optimization.first != 0) {
      ++result.failures;
    }
    result.checksum += workingMolecule.getNumAtoms();
    result.checksum += static_cast<std::uint64_t>(optimization.first + 2);
  }
  return result;
}

std::vector<RDKit::ROMol> prepareMmffMolecules(
    const std::vector<RDKit::ROMol> &molecules) {
  std::vector<RDKit::ROMol> prepared;
  prepared.reserve(molecules.size());
  for (const auto &molecule : molecules) {
    std::unique_ptr<RDKit::ROMol> withHs(RDKit::MolOps::addHs(molecule));
    auto parameters = RDKit::DGeomHelpers::ETKDGv3;
    parameters.numThreads = 1;
    parameters.randomSeed = 0xC0FFEE;
    if (RDKit::DGeomHelpers::EmbedMolecule(*withHs, parameters) >= 0) {
      prepared.push_back(std::move(*withHs));
    }
  }
  if (prepared.empty()) {
    throw std::runtime_error("no molecules could be prepared for MMFF");
  }
  return prepared;
}

OperationResult benchmarkTanimotoSimilarity(
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fingerprints) {
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

OperationResult benchmarkInchiRoundtrip(
    const std::vector<RDKit::ROMol> &molecules) {
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

OperationResult benchmarkHsRoundtrip(
    const std::vector<RDKit::ROMol> &molecules) {
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
    Operation operation, const std::vector<RDKit::ROMol> &molecules,
    const RDKit::FilterCatalog *painsCatalog,
    const RDKit::FingerprintGenerator<std::uint64_t> *generator,
    const std::vector<std::unique_ptr<ExplicitBitVect>> &fingerprints,
    const std::vector<RDKit::ROMol> &mmffMolecules) {
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
    case Operation::Mmff:
      return benchmarkMmff(mmffMolecules, molecules.size());
    case Operation::TanimotoSimilarity:
      return benchmarkTanimotoSimilarity(fingerprints);
    case Operation::InchiRoundtrip:
      return benchmarkInchiRoundtrip(molecules);
    case Operation::HsRoundtrip:
      return benchmarkHsRoundtrip(molecules);
  }
  throw std::logic_error("unhandled benchmark operation");
}

void printUsage(const char *program) {
  std::cerr << "Usage: " << program << " OPERATION INPUT COUNT REPEATS\n"
            << "Operations:";
  for (const auto &option : operationOptions) {
    std::cerr << ' ' << option.name;
  }
  std::cerr << '\n';
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
    auto loaded = loadMolecules(argv[2], requested);
    std::unique_ptr<RDKit::FilterCatalog> painsCatalog;
    if (operation == Operation::PainsSubstructure) {
      painsCatalog = std::make_unique<RDKit::FilterCatalog>(
          RDKit::FilterCatalogParams::PAINS);
    }
    std::unique_ptr<RDKit::FingerprintGenerator<std::uint64_t>> generator;
    if (operation == Operation::Morgan ||
        operation == Operation::TanimotoSimilarity) {
      generator.reset(
          RDKit::MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    }
    std::vector<std::unique_ptr<ExplicitBitVect>> fingerprints;
    if (operation == Operation::TanimotoSimilarity) {
      fingerprints = makeFingerprints(loaded.molecules, *generator);
    }
    std::vector<RDKit::ROMol> mmffMolecules;
    if (operation == Operation::Mmff) {
      mmffMolecules = prepareMmffMolecules(loaded.molecules);
    }
    const auto warmup =
        runOperation(operation, loaded.molecules, painsCatalog.get(),
                     generator.get(), fingerprints, mmffMolecules);
    const auto setupElapsed = std::chrono::steady_clock::now() - setupStart;

    OperationResult result;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t repeat = 0; repeat < repeats; ++repeat) {
      const auto iteration =
          runOperation(operation, loaded.molecules, painsCatalog.get(),
                       generator.get(), fingerprints, mmffMolecules);
      result.checksum += iteration.checksum;
      result.attempts += iteration.attempts;
      result.failures += iteration.failures;
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;

    std::cout << "operation=" << argv[1] << " input=" << argv[2]
              << " requested=" << requested << " read=" << loaded.read
              << " parsed=" << loaded.molecules.size()
              << " failures=" << loaded.failures << " atoms=" << loaded.atoms
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
