#!/usr/bin/env python3
"""Benchmark MCS and small-into-large substructure search on the same pairs.

The default pair generator uses a SMILES file and creates query/target pairs by
peeling terminal atoms from target molecules. This keeps most rings intact,
produces a smaller query, and normally makes the query a true substructure of
the larger target. For ChEMBL-scale size scans, use --pair-mode
size-combinations. That scans the input once, keeps representative molecules per
heavy-atom count, and samples small/large combinations from those examples.

The highlighted baseline is:
  - substruct_smiles: raw SMILES query, HasSubstructMatch(query, target)
  - mcs_elements_order: MCS with CompareElements and CompareOrder

Those two are the closest direct comparison for "SMILES substructure search vs
MCS": both use element and bond-order/aromaticity semantics from the molecule.
Substructure matching does not have a direct completeRingsOnly equivalent; the
ring-chain adjusted query variants are the closest ring-atom analog.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CHEMBL_INPUT = Path.home() / "data" / "chembl.smi"
DEFAULT_NCI_INPUT = REPO_ROOT / "Data" / "NCI" / "first_5K.smi"
DEFAULT_INPUT = DEFAULT_CHEMBL_INPUT if DEFAULT_CHEMBL_INPUT.exists() else DEFAULT_NCI_INPUT


@dataclass
class MolRecord:
  index: int
  smiles: str
  name: str
  mol: Chem.Mol
  atoms: int
  bonds: int
  rings: int


@dataclass
class PairRecord:
  pair_id: str
  query_smiles: str
  target_smiles: str
  query_mol: Chem.Mol
  target_mol: Chem.Mol
  query_name: str
  target_name: str
  pair_source: str


@dataclass(frozen=True)
class MethodSpec:
  family: str
  method: str
  label: str
  highlighted: bool
  notes: str
  kwargs: dict


FIELDNAMES = [
  "pair_id",
  "pair_source",
  "family",
  "method",
  "label",
  "highlighted",
  "query_name",
  "target_name",
  "query_smiles",
  "target_smiles",
  "query_atoms",
  "target_atoms",
  "total_atoms",
  "smaller_atoms",
  "larger_atoms",
  "combined_atoms",
  "atom_delta",
  "query_bonds",
  "target_bonds",
  "query_rings",
  "target_rings",
  "status",
  "elapsed_s",
  "repeats",
  "seconds_per_call",
  "ns_per_call",
  "substruct_op",
  "substruct_match",
  "substruct_match_count",
  "mcs_atoms",
  "mcs_bonds",
  "mcs_canceled",
  "mcs_smarts",
  "notes",
  "error",
]


def parse_int(value: str) -> int:
  return int(value, 0)


def open_text(path: Path):
  if path.suffix == ".gz":
    return gzip.open(path, "rt")
  return path.open("rt")


def iter_smiles(path: Path, limit: int | None = None) -> Iterable[tuple[int, str, str]]:
  with open_text(path) as handle:
    for line_no, line in enumerate(handle, start=1):
      if limit is not None and line_no > limit:
        break
      stripped = line.strip()
      if not stripped or stripped.startswith("#"):
        continue
      parts = stripped.split()
      smiles = parts[0]
      name = " ".join(parts[1:]) if len(parts) > 1 else str(line_no)
      yield line_no, smiles, name


def load_molecules(path: Path, limit: int | None) -> list[MolRecord]:
  records: list[MolRecord] = []
  for line_no, smiles, name in iter_smiles(path, limit):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
      print(f"Skipping unparsable SMILES at {path}:{line_no}: {smiles}", file=sys.stderr)
      continue
    records.append(
      MolRecord(
        index=line_no,
        smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
        name=name,
        mol=mol,
        atoms=mol.GetNumHeavyAtoms(),
        bonds=mol.GetNumBonds(),
        rings=mol.GetRingInfo().NumRings(),
      ))
  return records


def load_size_sampled_molecules(
  path: Path,
  limit: int | None,
  min_atoms: int,
  max_atoms: int,
  examples_per_size: int,
  seed: int,
) -> list[MolRecord]:
  rng = random.Random(seed)
  reservoirs: dict[int, list[MolRecord]] = {}
  seen_by_size: dict[int, int] = {}

  for line_no, smiles, name in iter_smiles(path, limit):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
      print(f"Skipping unparsable SMILES at {path}:{line_no}: {smiles}", file=sys.stderr)
      continue

    atoms = mol.GetNumHeavyAtoms()
    if atoms < min_atoms or atoms > max_atoms:
      continue

    record = MolRecord(
      index=line_no,
      smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
      name=name,
      mol=mol,
      atoms=atoms,
      bonds=mol.GetNumBonds(),
      rings=mol.GetRingInfo().NumRings(),
    )
    seen_by_size[atoms] = seen_by_size.get(atoms, 0) + 1
    bucket = reservoirs.setdefault(atoms, [])
    if len(bucket) < examples_per_size:
      bucket.append(record)
      continue

    replace_idx = rng.randrange(seen_by_size[atoms])
    if replace_idx < examples_per_size:
      bucket[replace_idx] = record

  records: list[MolRecord] = []
  for atom_count in sorted(reservoirs):
    records.extend(
      sorted(reservoirs[atom_count], key=lambda rec: (rec.bonds, rec.rings, rec.index)))
  return records


def thin_evenly(records: list[MolRecord], max_items: int) -> list[MolRecord]:
  if max_items <= 0 or len(records) <= max_items:
    return list(records)
  if max_items == 1:
    return [records[len(records) // 2]]

  selected: list[MolRecord] = []
  seen: set[int] = set()
  for i in range(max_items):
    idx = round(i * (len(records) - 1) / (max_items - 1))
    if idx in seen:
      continue
    seen.add(idx)
    selected.append(records[idx])
  return selected


def largest_fragment_atoms(mol: Chem.Mol) -> set[int]:
  fragments = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
  if not fragments:
    return set()
  return set(max(fragments, key=len))


def induced_degree(mol: Chem.Mol, atom_idx: int, atoms: set[int]) -> int:
  atom = mol.GetAtomWithIdx(atom_idx)
  return sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in atoms)


def peel_terminal_atoms(
  mol: Chem.Mol,
  desired_atoms: int,
  rng: random.Random,
) -> list[int]:
  remaining = largest_fragment_atoms(mol)
  if not remaining:
    return []

  desired_atoms = max(1, min(desired_atoms, len(remaining)))
  while len(remaining) > desired_atoms:
    leaves = [idx for idx in remaining if induced_degree(mol, idx, remaining) <= 1]
    if not leaves:
      break

    non_ring_leaves = [
      idx for idx in leaves if not mol.GetAtomWithIdx(idx).IsInRing()
    ]
    candidates = non_ring_leaves or leaves
    remaining.remove(rng.choice(candidates))

  return sorted(remaining)


def fragment_to_query_mol(mol: Chem.Mol, atoms: list[int]) -> tuple[str, Chem.Mol] | None:
  if not atoms:
    return None

  smiles_kwargs = [
    {
      "canonical": True,
      "isomericSmiles": True,
    },
    {
      "canonical": True,
      "isomericSmiles": True,
      "kekuleSmiles": True,
    },
  ]
  for kwargs in smiles_kwargs:
    try:
      smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=atoms, **kwargs)
    except Exception:
      continue
    query = Chem.MolFromSmiles(smiles)
    if query is not None:
      return Chem.MolToSmiles(query, isomericSmiles=True), query
  return None


def make_fragment_parent_pairs(
  records: list[MolRecord],
  max_pairs: int,
  min_target_atoms: int,
  max_target_atoms: int | None,
  min_query_atoms: int,
  query_fraction: float,
  seed: int,
) -> list[PairRecord]:
  rng = random.Random(seed)
  candidates = [
    rec for rec in records
    if rec.atoms >= min_target_atoms
    and (max_target_atoms is None or rec.atoms <= max_target_atoms)
  ]
  candidates.sort(key=lambda rec: (rec.atoms, rec.bonds, rec.index))
  candidates = thin_evenly(candidates, max_pairs * 4)

  pairs: list[PairRecord] = []
  for target in candidates:
    if len(pairs) >= max_pairs:
      break

    desired_atoms = max(min_query_atoms, round(target.atoms * query_fraction))
    desired_atoms = min(desired_atoms, target.atoms - 1)
    if desired_atoms < min_query_atoms:
      continue

    query_data = None
    for _ in range(32):
      atoms = peel_terminal_atoms(target.mol, desired_atoms, rng)
      if len(atoms) >= target.atoms or len(atoms) < min_query_atoms:
        continue
      query_data = fragment_to_query_mol(target.mol, atoms)
      if query_data is None:
        continue
      _, query_mol = query_data
      if query_mol.GetNumAtoms() < target.atoms and target.mol.HasSubstructMatch(query_mol):
        break
      query_data = None

    if query_data is None:
      continue

    query_smiles, query_mol = query_data
    pair_id = f"p{len(pairs):04d}"
    pairs.append(
      PairRecord(
        pair_id=pair_id,
        query_smiles=query_smiles,
        target_smiles=target.smiles,
        query_mol=query_mol,
        target_mol=target.mol,
        query_name=f"{target.name}:peeled",
        target_name=target.name,
        pair_source="fragment-parent",
      ))

  return pairs


def make_dataset_pairs(
  records: list[MolRecord],
  max_pairs: int,
  min_target_atoms: int,
  max_target_atoms: int | None,
  min_query_atoms: int,
  query_fraction: float,
) -> list[PairRecord]:
  targets = [
    rec for rec in records
    if rec.atoms >= min_target_atoms
    and (max_target_atoms is None or rec.atoms <= max_target_atoms)
  ]
  targets.sort(key=lambda rec: (rec.atoms, rec.bonds, rec.index))
  targets = thin_evenly(targets, max_pairs)

  sorted_records = sorted(records, key=lambda rec: (rec.atoms, rec.bonds, rec.index))
  pairs: list[PairRecord] = []
  for target in targets:
    max_query_atoms = max(min_query_atoms, round(target.atoms * query_fraction))
    eligible = [
      rec for rec in sorted_records
      if rec.index != target.index
      and min_query_atoms <= rec.atoms <= max_query_atoms
      and rec.atoms < target.atoms
    ]
    if not eligible:
      continue
    query = eligible[len(eligible) // 2]
    pair_id = f"p{len(pairs):04d}"
    pairs.append(
      PairRecord(
        pair_id=pair_id,
        query_smiles=query.smiles,
        target_smiles=target.smiles,
        query_mol=query.mol,
        target_mol=target.mol,
        query_name=query.name,
        target_name=target.name,
        pair_source="dataset-small-large",
      ))
  return pairs


def make_size_combination_pairs(
  records: list[MolRecord],
  max_pairs: int,
  seed: int,
) -> list[PairRecord]:
  by_size = sorted(records, key=lambda rec: (rec.atoms, rec.bonds, rec.rings, rec.index))
  combinations: list[tuple[MolRecord, MolRecord]] = []
  for small_idx, small in enumerate(by_size):
    for large in by_size[small_idx + 1:]:
      if small.atoms >= large.atoms:
        continue
      combinations.append((small, large))

  if max_pairs > 0 and len(combinations) > max_pairs:
    rng = random.Random(seed)
    combinations = rng.sample(combinations, max_pairs)
  combinations.sort(key=lambda pair: (pair[0].atoms, pair[1].atoms, pair[0].index, pair[1].index))

  pairs: list[PairRecord] = []
  for small, large in combinations:
    pair_id = f"p{len(pairs):04d}"
    pairs.append(
      PairRecord(
        pair_id=pair_id,
        query_smiles=small.smiles,
        target_smiles=large.smiles,
        query_mol=small.mol,
        target_mol=large.mol,
        query_name=small.name,
        target_name=large.name,
        pair_source="size-combinations",
      ))
  return pairs


def mcs_methods() -> list[MethodSpec]:
  return [
    MethodSpec(
      family="mcs",
      method="mcs_any_atom_bond",
      label="MCS CompareAny atoms/bonds",
      highlighted=False,
      notes="Permissive MCS: any atom type and any bond type can match.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareAny,
        "bondCompare": rdFMCS.BondCompare.CompareAny,
      },
    ),
    MethodSpec(
      family="mcs",
      method="mcs_elements_order",
      label="MCS elements/order",
      highlighted=True,
      notes="Closest MCS counterpart to a raw SMILES substructure query.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareElements,
        "bondCompare": rdFMCS.BondCompare.CompareOrder,
      },
    ),
    MethodSpec(
      family="mcs",
      method="mcs_elements_order_maximize_atoms",
      label="MCS elements/order maximize atoms",
      highlighted=False,
      notes="Same element/order comparison, but optimize atom count instead of bond count.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareElements,
        "bondCompare": rdFMCS.BondCompare.CompareOrder,
        "maximizeBonds": False,
      },
    ),
    MethodSpec(
      family="mcs",
      method="mcs_elements_order_ring_atoms",
      label="MCS elements/order ring atoms",
      highlighted=False,
      notes="Ring atoms and ring bonds only match ring atoms and ring bonds.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareElements,
        "bondCompare": rdFMCS.BondCompare.CompareOrder,
        "ringMatchesRingOnly": True,
      },
    ),
    MethodSpec(
      family="mcs",
      method="mcs_elements_order_whole_rings",
      label="MCS elements/order whole rings",
      highlighted=False,
      notes="MCS completeRingsOnly=True; partial ring MCS results are rejected.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareElements,
        "bondCompare": rdFMCS.BondCompare.CompareOrder,
        "completeRingsOnly": True,
      },
    ),
    MethodSpec(
      family="mcs",
      method="mcs_elements_order_ring_atoms_whole_rings",
      label="MCS elements/order ring atoms + whole rings",
      highlighted=False,
      notes="Both ringMatchesRingOnly and completeRingsOnly are enabled.",
      kwargs={
        "atomCompare": rdFMCS.AtomCompare.CompareElements,
        "bondCompare": rdFMCS.BondCompare.CompareOrder,
        "ringMatchesRingOnly": True,
        "completeRingsOnly": True,
      },
    ),
  ]


def substruct_methods() -> list[MethodSpec]:
  return [
    MethodSpec(
      family="substruct",
      method="substruct_smiles",
      label="Substruct raw SMILES",
      highlighted=True,
      notes=(
        "Best SMILES-substructure baseline. Substructure has no direct "
        "completeRingsOnly option."
      ),
      kwargs={
        "query_kind": "smiles",
      },
    ),
    MethodSpec(
      family="substruct",
      method="substruct_elements_any_bond",
      label="Substruct elements/any bond",
      highlighted=False,
      notes="Query keeps atom elements but uses generic bonds.",
      kwargs={
        "query_kind": "elements_any_bond",
      },
    ),
    MethodSpec(
      family="substruct",
      method="substruct_any_atom_bond",
      label="Substruct any atoms/bonds",
      highlighted=False,
      notes="Query uses generic atoms and generic bonds, similar to CompareAny topology.",
      kwargs={
        "query_kind": "any_atom_bond",
      },
    ),
    MethodSpec(
      family="substruct",
      method="substruct_ring_atoms",
      label="Substruct ring-chain atoms",
      highlighted=False,
      notes="Query adds ring/chain atom constraints; closest substruct analog to ringMatchesRingOnly.",
      kwargs={
        "query_kind": "ring_atoms",
      },
    ),
  ]


def adjusted_query(query: Chem.Mol, query_kind: str) -> Chem.Mol:
  if query_kind == "smiles":
    return Chem.Mol(query)

  params = Chem.AdjustQueryParameters.NoAdjustments()
  if query_kind == "elements_any_bond":
    params.makeBondsGeneric = True
  elif query_kind == "any_atom_bond":
    params.makeAtomsGeneric = True
    params.makeBondsGeneric = True
  elif query_kind == "ring_atoms":
    params.adjustRingChain = True
  else:
    raise ValueError(f"unknown substructure query kind: {query_kind}")
  return Chem.AdjustQueryProperties(Chem.Mol(query), params)


def substruct_parameters() -> Chem.SubstructMatchParameters:
  params = Chem.SubstructMatchParameters()
  params.numThreads = 1
  return params


def time_callable(
  func: Callable[[], object],
  repeats: int,
  warmups: int,
  min_seconds: float,
  max_repeats: int,
) -> tuple[float, int, object]:
  last_result = None
  for _ in range(warmups):
    last_result = func()

  count = 0
  start_ns = time.perf_counter_ns()
  while count < repeats:
    last_result = func()
    count += 1

  while min_seconds > 0.0 and count < max_repeats:
    elapsed_s = (time.perf_counter_ns() - start_ns) / 1.0e9
    if elapsed_s >= min_seconds:
      break
    last_result = func()
    count += 1

  elapsed_s = (time.perf_counter_ns() - start_ns) / 1.0e9
  return elapsed_s, count, last_result


def base_row(pair: PairRecord, spec: MethodSpec) -> dict[str, object]:
  query = pair.query_mol
  target = pair.target_mol
  query_atoms = query.GetNumHeavyAtoms()
  target_atoms = target.GetNumHeavyAtoms()
  smaller_atoms = min(query_atoms, target_atoms)
  larger_atoms = max(query_atoms, target_atoms)
  combined_atoms = query_atoms + target_atoms
  return {
    "pair_id": pair.pair_id,
    "pair_source": pair.pair_source,
    "family": spec.family,
    "method": spec.method,
    "label": spec.label,
    "highlighted": int(spec.highlighted),
    "query_name": pair.query_name,
    "target_name": pair.target_name,
    "query_smiles": pair.query_smiles,
    "target_smiles": pair.target_smiles,
    "query_atoms": query_atoms,
    "target_atoms": target_atoms,
    "total_atoms": combined_atoms,
    "smaller_atoms": smaller_atoms,
    "larger_atoms": larger_atoms,
    "combined_atoms": combined_atoms,
    "atom_delta": larger_atoms - smaller_atoms,
    "query_bonds": query.GetNumBonds(),
    "target_bonds": target.GetNumBonds(),
    "query_rings": query.GetRingInfo().NumRings(),
    "target_rings": target.GetRingInfo().NumRings(),
    "status": "ok",
    "elapsed_s": "",
    "repeats": "",
    "seconds_per_call": "",
    "ns_per_call": "",
    "substruct_op": "",
    "substruct_match": "",
    "substruct_match_count": "",
    "mcs_atoms": "",
    "mcs_bonds": "",
    "mcs_canceled": "",
    "mcs_smarts": "",
    "notes": spec.notes,
    "error": "",
  }


def run_mcs(
  pair: PairRecord,
  spec: MethodSpec,
  mcs_timeout: int,
  repeats: int,
  warmups: int,
  min_seconds: float,
  max_repeats: int,
) -> dict[str, object]:
  row = base_row(pair, spec)
  kwargs = dict(spec.kwargs)
  kwargs["timeout"] = mcs_timeout

  def call_mcs():
    return rdFMCS.FindMCS([pair.query_mol, pair.target_mol], **kwargs)

  try:
    elapsed_s, count, result = time_callable(
      call_mcs, repeats, warmups, min_seconds, max_repeats)
    row.update({
      "elapsed_s": f"{elapsed_s:.9g}",
      "repeats": count,
      "seconds_per_call": f"{elapsed_s / count:.9g}",
      "ns_per_call": f"{elapsed_s * 1.0e9 / count:.3f}",
      "mcs_atoms": result.numAtoms,
      "mcs_bonds": result.numBonds,
      "mcs_canceled": int(result.canceled),
      "mcs_smarts": result.smartsString,
    })
  except Exception as exc:
    row["status"] = "error"
    row["error"] = repr(exc)
  return row


def run_substruct(
  pair: PairRecord,
  spec: MethodSpec,
  substruct_op: str,
  repeats: int,
  warmups: int,
  min_seconds: float,
  max_repeats: int,
) -> dict[str, object]:
  row = base_row(pair, spec)
  row["substruct_op"] = substruct_op

  try:
    query = adjusted_query(pair.query_mol, spec.kwargs["query_kind"])
    params = substruct_parameters()

    if substruct_op == "has-match":
      def call_substruct():
        return pair.target_mol.HasSubstructMatch(query, params)
    elif substruct_op == "get-matches":
      def call_substruct():
        return pair.target_mol.GetSubstructMatches(query, params)
    else:
      raise ValueError(f"unknown substructure operation: {substruct_op}")

    elapsed_s, count, result = time_callable(
      call_substruct, repeats, warmups, min_seconds, max_repeats)

    if substruct_op == "has-match":
      matched = bool(result)
      match_count = int(matched)
    else:
      matched = bool(result)
      match_count = len(result)

    row.update({
      "elapsed_s": f"{elapsed_s:.9g}",
      "repeats": count,
      "seconds_per_call": f"{elapsed_s / count:.9g}",
      "ns_per_call": f"{elapsed_s * 1.0e9 / count:.3f}",
      "substruct_match": int(matched),
      "substruct_match_count": match_count,
    })
  except Exception as exc:
    row["status"] = "error"
    row["error"] = repr(exc)
  return row


def make_pairs(args: argparse.Namespace, records: list[MolRecord]) -> list[PairRecord]:
  max_target_atoms = args.max_target_atoms if args.max_target_atoms > 0 else None
  if args.pair_mode == "fragment-parent":
    return make_fragment_parent_pairs(
      records=records,
      max_pairs=args.max_pairs,
      min_target_atoms=args.min_target_atoms,
      max_target_atoms=max_target_atoms,
      min_query_atoms=args.min_query_atoms,
      query_fraction=args.query_fraction,
      seed=args.seed,
    )
  if args.pair_mode == "dataset":
    return make_dataset_pairs(
      records=records,
      max_pairs=args.max_pairs,
      min_target_atoms=args.min_target_atoms,
      max_target_atoms=max_target_atoms,
      min_query_atoms=args.min_query_atoms,
      query_fraction=args.query_fraction,
    )
  if args.pair_mode == "size-combinations":
    return make_size_combination_pairs(
      records=records,
      max_pairs=args.max_pairs,
      seed=args.seed,
    )
  raise ValueError(f"unknown pair mode: {args.pair_mode}")


def next_pair_index(path: Path) -> int:
  if not path.exists() or path.stat().st_size == 0:
    return 0

  max_index = -1
  unique_pair_ids: set[str] = set()
  with path.open(newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
      pair_id = row.get("pair_id", "")
      if not pair_id:
        continue
      unique_pair_ids.add(pair_id)
      match = re.fullmatch(r"p(\d+)", pair_id)
      if match:
        max_index = max(max_index, int(match.group(1)))

  if max_index >= 0:
    return max_index + 1
  return len(unique_pair_ids)


def renumber_pairs(pairs: list[PairRecord], start_index: int) -> None:
  for offset, pair in enumerate(pairs):
    pair.pair_id = f"p{start_index + offset:04d}"


def selected_methods(family: str) -> list[MethodSpec]:
  specs: list[MethodSpec] = []
  if family in {"all", "mcs"}:
    specs.extend(mcs_methods())
  if family in {"all", "substruct"}:
    specs.extend(substruct_methods())
  return specs


def write_rows(
  path: Path,
  rows: Iterable[dict[str, object]],
  flush_every: int,
  progress_every: int,
  append: bool,
) -> int:
  path.parent.mkdir(parents=True, exist_ok=True)
  count = 0
  write_header = not append or not path.exists() or path.stat().st_size == 0
  mode = "a" if append else "w"
  with path.open(mode, newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
    if write_header:
      writer.writeheader()
    for row in rows:
      writer.writerow(row)
      count += 1
      if flush_every > 0 and count % flush_every == 0:
        handle.flush()
      if progress_every > 0 and count % progress_every == 0:
        print(f"Wrote {count} benchmark rows to {path}", file=sys.stderr)
    handle.flush()
  return count


def benchmark_rows(args: argparse.Namespace, pairs: list[PairRecord]) -> Iterable[dict[str, object]]:
  specs = selected_methods(args.family)
  for pair in pairs:
    for spec in specs:
      if spec.family == "mcs":
        yield run_mcs(
          pair=pair,
          spec=spec,
          mcs_timeout=args.mcs_timeout,
          repeats=args.mcs_repeats,
          warmups=args.warmups,
          min_seconds=args.min_seconds,
          max_repeats=args.max_repeats,
        )
      else:
        yield run_substruct(
          pair=pair,
          spec=spec,
          substruct_op=args.substruct_op,
          repeats=args.substruct_repeats,
          warmups=args.warmups,
          min_seconds=args.min_seconds,
          max_repeats=args.max_repeats,
        )


def parse_args(argv: list[str]) -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Benchmark RDKit MCS variants against small-into-large substructure search.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                      help="SMILES file. The first whitespace-separated field is the SMILES.")
  parser.add_argument("--output", type=Path, default=Path("mcs_substruct_benchmark.csv"),
                      help="Output CSV path.")
  parser.add_argument("--append", action="store_true",
                      help="Append rows to an existing CSV and continue pair IDs.")
  parser.add_argument("--input-limit", type=int, default=None,
                      help="Maximum number of input lines to read.")
  parser.add_argument("--show-rdkit-parse-errors", action="store_true",
                      help="Show RDKit parser diagnostics for invalid input SMILES.")
  parser.add_argument("--pair-mode", choices=["fragment-parent", "dataset", "size-combinations"],
                      default="size-combinations",
                      help="How to construct query/target pairs.")
  parser.add_argument("--max-pairs", type=int, default=40,
                      help="Maximum number of query/target pairs.")
  parser.add_argument("--min-target-atoms", type=int, default=8,
                      help="Minimum target molecule atom count.")
  parser.add_argument("--max-target-atoms", type=int, default=0,
                      help="Maximum target molecule atom count. Use 0 for no limit.")
  parser.add_argument("--min-query-atoms", type=int, default=3,
                      help="Minimum query molecule atom count.")
  parser.add_argument("--query-fraction", type=float, default=0.45,
                      help="Approximate query atom fraction relative to the target.")
  parser.add_argument("--seed", type=parse_int, default=0x5EED,
                      help="Random seed used by fragment-parent pair generation.")
  parser.add_argument("--size-min-atoms", type=int, default=0,
                      help="Minimum heavy-atom count for --pair-mode size-combinations.")
  parser.add_argument("--size-max-atoms", type=int, default=100,
                      help="Maximum heavy-atom count for --pair-mode size-combinations.")
  parser.add_argument("--examples-per-size", type=int, default=1,
                      help="Representative molecules retained for each heavy-atom count.")
  parser.add_argument("--family", choices=["all", "mcs", "substruct"], default="all",
                      help="Which benchmark family to run.")
  parser.add_argument("--substruct-op", choices=["has-match", "get-matches"], default="has-match",
                      help="Substructure operation to benchmark.")
  parser.add_argument("--mcs-timeout", type=int, default=5,
                      help="Per-call MCS timeout in seconds. Use 0 for no timeout.")
  parser.add_argument("--mcs-repeats", type=int, default=1,
                      help="Measured repeats for each MCS pair/method.")
  parser.add_argument("--substruct-repeats", type=int, default=200,
                      help="Measured repeats for each substructure pair/method.")
  parser.add_argument("--warmups", type=int, default=1,
                      help="Warmup calls before measured repeats.")
  parser.add_argument("--min-seconds", type=float, default=0.0,
                      help="Keep repeating each pair/method until this many measured seconds.")
  parser.add_argument("--max-repeats", type=int, default=100000,
                      help="Maximum repeats when --min-seconds is active.")
  parser.add_argument("--flush-every", type=int, default=1,
                      help="Flush the CSV after this many rows so partial results remain usable.")
  parser.add_argument("--progress-every", type=int, default=25,
                      help="Write progress to stderr after this many benchmark rows. Use 0 to disable.")
  return parser.parse_args(argv)


def main(argv: list[str]) -> int:
  args = parse_args(argv)
  if not args.show_rdkit_parse_errors:
    RDLogger.DisableLog("rdApp.error")

  if not args.input.exists():
    print(f"Input SMILES file does not exist: {args.input}", file=sys.stderr)
    return 2
  if args.examples_per_size < 1:
    print("--examples-per-size must be at least 1", file=sys.stderr)
    return 2
  if args.size_min_atoms > args.size_max_atoms:
    print("--size-min-atoms must be <= --size-max-atoms", file=sys.stderr)
    return 2

  if args.pair_mode == "size-combinations":
    records = load_size_sampled_molecules(
      path=args.input,
      limit=args.input_limit,
      min_atoms=args.size_min_atoms,
      max_atoms=args.size_max_atoms,
      examples_per_size=args.examples_per_size,
      seed=args.seed,
    )
  else:
    records = load_molecules(args.input, args.input_limit)
  if not records:
    print(f"No molecules loaded from {args.input}", file=sys.stderr)
    return 2

  pairs = make_pairs(args, records)
  if not pairs:
    print("No query/target pairs were generated. Try lowering atom limits or use --pair-mode dataset.",
          file=sys.stderr)
    return 2

  start_pair_index = next_pair_index(args.output) if args.append else 0
  renumber_pairs(pairs, start_pair_index)

  row_count = write_rows(
    args.output,
    benchmark_rows(args, pairs),
    flush_every=args.flush_every,
    progress_every=args.progress_every,
    append=args.append,
  )
  print(f"Wrote {row_count} benchmark rows for {len(pairs)} pairs to {args.output}")
  print("Highlighted baseline: substruct_smiles vs mcs_elements_order")
  print("Substructure note: no direct completeRingsOnly parameter; substruct_ring_atoms is the ring-chain analog.")
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
