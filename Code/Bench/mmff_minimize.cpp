#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string_view>
#include <vector>

#include <ForceField/ForceField.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>

namespace {

struct WorkItem {
  std::unique_ptr<RDKit::ROMol> molecule;
  std::unique_ptr<ForceFields::ForceField> forceField;
  std::vector<RDGeom::Point3D> initialPositions;
  double initialEnergy = 0.0;
};

unsigned int parsePositiveInteger(std::string_view text,
                                  std::string_view argumentName) {
  unsigned int value = 0;
  const auto result =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (result.ec != std::errc{} || result.ptr != text.data() + text.size() ||
      value == 0) {
    std::cerr << argumentName << " must be a positive integer\n";
    std::exit(2);
  }
  return value;
}

std::vector<WorkItem> loadWorkItems(const char *sdfPath,
                                    unsigned int &rejectedRecords) {
  RDKit::SDMolSupplier supplier(sdfPath, true, false);
  std::vector<WorkItem> workItems;
  rejectedRecords = 0;

  while (!supplier.atEnd()) {
    std::unique_ptr<RDKit::ROMol> inputMolecule;
    try {
      inputMolecule.reset(supplier.next());
    } catch (...) {
      ++rejectedRecords;
      continue;
    }
    if (!inputMolecule) {
      ++rejectedRecords;
      continue;
    }

    std::unique_ptr<RDKit::ROMol> molecule(
        RDKit::MolOps::addHs(*inputMolecule, false, true));
    RDKit::MMFF::MMFFMolProperties properties(*molecule, "MMFF94");
    if (!properties.isValid()) {
      ++rejectedRecords;
      continue;
    }

    std::unique_ptr<ForceFields::ForceField> forceField(
        RDKit::MMFF::constructForceField(*molecule, &properties));
    forceField->initialize();

    WorkItem workItem;
    const auto &conformer = molecule->getConformer();
    workItem.initialPositions.reserve(molecule->getNumAtoms());
    for (unsigned int atomIndex = 0; atomIndex < molecule->getNumAtoms();
         ++atomIndex) {
      workItem.initialPositions.push_back(conformer.getAtomPos(atomIndex));
    }
    workItem.initialEnergy = forceField->calcEnergy();
    workItem.molecule = std::move(molecule);
    workItem.forceField = std::move(forceField);
    workItems.push_back(std::move(workItem));
  }

  return workItems;
}

void restoreInitialPositions(std::vector<WorkItem> &workItems) {
  for (auto &workItem : workItems) {
    auto &positions = workItem.forceField->positions();
    for (unsigned int atomIndex = 0;
         atomIndex < workItem.initialPositions.size(); ++atomIndex) {
      for (unsigned int coordinate = 0; coordinate < 3; ++coordinate) {
        (*positions[atomIndex])[coordinate] =
            workItem.initialPositions[atomIndex][coordinate];
      }
    }
    workItem.forceField->initialize();
  }
}

}  // namespace

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 5) {
    std::cerr << "usage: " << argv[0]
              << " <input.sdf> [repetitions=1] [max-iterations=200] "
                 "[molecules=200]\n";
    return 2;
  }

  const unsigned int repetitions =
      argc >= 3 ? parsePositiveInteger(argv[2], "repetitions") : 1;
  const unsigned int maxIterations =
      argc >= 4 ? parsePositiveInteger(argv[3], "max-iterations") : 200;
  const unsigned int moleculeLimit =
      argc >= 5 ? parsePositiveInteger(argv[4], "molecules") : 200;

  unsigned int rejectedRecords = 0;
  auto workItems = loadWorkItems(argv[1], rejectedRecords);
  if (workItems.empty()) {
    std::cerr << "no molecules with complete MMFF94 parameters were loaded\n";
    return 1;
  }
  if (moleculeLimit > workItems.size()) {
    std::cerr << "requested " << moleculeLimit << " molecules, but only "
              << workItems.size() << " have complete MMFF94 parameters\n";
    return 1;
  }
  workItems.resize(moleculeLimit);

  unsigned int explicitHydrogens = 0;
  for (const auto &workItem : workItems) {
    explicitHydrogens += workItem.molecule->getNumAtoms() -
                         workItem.molecule->getNumHeavyAtoms();
  }

  const auto atomCountBounds = std::minmax_element(
      workItems.begin(), workItems.end(),
      [](const auto &left, const auto &right) {
        return left.molecule->getNumAtoms() < right.molecule->getNumAtoms();
      });
  std::cout << "molecules=" << workItems.size()
            << " rejected_records=" << rejectedRecords
            << " explicit_hydrogens=" << explicitHydrogens
            << " atoms_min=" << atomCountBounds.first->molecule->getNumAtoms()
            << " atoms_max=" << atomCountBounds.second->molecule->getNumAtoms()
            << " max_iterations=" << maxIterations << " threads=1\n"
            << std::flush;

  const double initialEnergySum = [&workItems]() {
    double sum = 0.0;
    for (const auto &workItem : workItems) {
      sum += workItem.initialEnergy;
    }
    return sum;
  }();

  std::cout << std::setprecision(12);
  for (unsigned int repetition = 1; repetition <= repetitions; ++repetition) {
    restoreInitialPositions(workItems);

    double restoredEnergySum = 0.0;
    for (const auto &workItem : workItems) {
      restoredEnergySum += workItem.forceField->calcEnergy();
    }
    const double energyScale = std::max(std::abs(initialEnergySum), 1.0);
    if (std::abs(restoredEnergySum - initialEnergySum) >
        1.0e-12 * energyScale) {
      std::cerr << "coordinate restoration changed the starting energy in "
                   "repetition "
                << repetition << '\n';
      return 1;
    }

    const auto start = std::chrono::steady_clock::now();
    unsigned int converged = 0;
    for (auto &workItem : workItems) {
      converged += workItem.forceField->minimize(maxIterations) == 0;
    }
    const auto stop = std::chrono::steady_clock::now();

    double finalEnergySum = 0.0;
    for (const auto &workItem : workItems) {
      const double finalEnergy = workItem.forceField->calcEnergy();
      if (!std::isfinite(finalEnergy)) {
        std::cerr << "non-finite final energy in repetition " << repetition
                  << '\n';
        return 1;
      }
      finalEnergySum += finalEnergy;
    }

    const double seconds =
        std::chrono::duration<double>(stop - start).count();
    std::cout << "repetition=" << repetition << " seconds=" << seconds
              << " molecules_per_second=" << workItems.size() / seconds
              << " milliseconds_per_molecule="
              << 1000.0 * seconds / workItems.size()
              << " converged=" << converged
              << " capped=" << workItems.size() - converged
              << " initial_energy_sum=" << initialEnergySum
              << " final_energy_sum=" << finalEnergySum
              << " energy_delta_sum=" << finalEnergySum - initialEnergySum
              << '\n';
  }

  return 0;
}
