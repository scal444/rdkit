#!/usr/bin/env python3

import argparse
import time

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def positive_integer(text):
  value = int(text)
  if value <= 0:
    raise argparse.ArgumentTypeError("must be a positive integer")
  return value


def load_molecules(path, requested):
  # Parse and retain molecules before timing so the bulk-similarity measurement
  # excludes file I/O, SMILES parsing, and fingerprint generation. Accept both
  # headerless SMILES files and tables whose first-column heading is "smiles".
  molecules = []
  rows_read = 0
  failures = 0
  atom_count = 0
  with open(path, encoding="utf-8") as input_file:
    for line in input_file:
      fields = line.split(maxsplit=1)
      smiles = fields[0] if fields else ""
      if rows_read == 0 and smiles == "smiles":
        continue
      rows_read += 1
      if smiles:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
          molecules.append(molecule)
          atom_count += molecule.GetNumAtoms()
        else:
          failures += 1
      else:
        failures += 1
      if rows_read == requested:
        break

  if rows_read != requested:
    raise ValueError("input ended before the requested molecule count")
  if not molecules:
    raise ValueError("input produced no molecules")
  return molecules, rows_read, failures, atom_count


def run_bulk_similarity(fingerprints):
  query_count = min(100, len(fingerprints))
  checksum = 0
  attempts = 0
  for query in range(query_count):
    query_index = query * len(fingerprints) // query_count
    similarities = DataStructs.BulkTanimotoSimilarity(
        fingerprints[query_index], fingerprints)
    checksum += sum(similarity >= 0.4 for similarity in similarities)
    attempts += len(similarities)
  return checksum, attempts


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input")
  parser.add_argument("count", type=positive_integer)
  parser.add_argument("repeats", type=positive_integer)
  args = parser.parse_args()

  setup_start = time.perf_counter()
  molecules, rows_read, failures, atom_count = load_molecules(
      args.input, args.count)
  generator = rdFingerprintGenerator.GetMorganGenerator(radius=2)
  fingerprints = [generator.GetFingerprint(molecule) for molecule in molecules]
  warmup_checksum, _ = run_bulk_similarity(fingerprints)
  setup_seconds = time.perf_counter() - setup_start

  checksum = 0
  attempts = 0
  start = time.perf_counter()
  for _ in range(args.repeats):
    iteration_checksum, iteration_attempts = run_bulk_similarity(fingerprints)
    checksum += iteration_checksum
    attempts += iteration_attempts
  elapsed_seconds = time.perf_counter() - start

  print(
      f"operation=bulk_similarity input={args.input} requested={args.count} "
      f"read={rows_read} parsed={len(molecules)} failures={failures} "
      f"atoms={atom_count} repeats={args.repeats} "
      f"warmup_checksum={warmup_checksum} checksum={checksum} "
      f"operation_attempts={attempts} setup_s={setup_seconds} "
      f"elapsed_s={elapsed_seconds}")


if __name__ == "__main__":
  main()
