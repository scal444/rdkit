"""
Conformer generation stability test.

Runs the RMS pruning conformer generation from test6RmsPruning
multiple times to characterize the distribution of conformer counts.

Usage:
    python confgen_stability.py [--runs N]

Run with the development build:
    LD_LIBRARY_PATH=$RDBASE/lib PYTHONPATH=$RDBASE python confgen_stability.py

Run with conda baseline:
    conda activate baseline_rdkit
    python confgen_stability.py
"""

import argparse
import json
import sys
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

SMILES = [
    "CC(C)CC(NC(C1[N+]CCC1)=O)C([O-])=O",
    "CC(NC(CO)C(O)c1ccc([N+]([O-])=O)cc1)=O",
    "CC([N+])C(NC(C)C(N1C(C=O)CCC1)=O)=O",
    "CC(NC1C(O)C=C(C([O-])=O)OC1C(O)C(O)CO)=O",
    "CCCC=C(NC(C1CC1(C)C)=O)C([O-])=O",
    "OCC(O)C(O)C(Cn1c2c(cc(C)c(C)c2)nc-2c(=O)[nH]c(=O)nc12)O",
]


def run_block1(seed):
    counts = []
    for smi in SMILES:
        mol = Chem.MolFromSmiles(smi)
        cids = rdDistGeom.EmbedMultipleConfs(
            mol, 50, maxAttempts=30, randomSeed=seed, pruneRmsThresh=1.5
        )
        counts.append(len(cids))
    return counts


def run_block2(seed):
    params = rdDistGeom.ETKDG()
    params.randomSeed = seed
    params.maxIterations = 30
    params.pruneRmsThresh = 1.5
    params.useSymmetryForPruning = False
    counts = []
    for smi in SMILES:
        mol = Chem.MolFromSmiles(smi)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 50, params)
        counts.append(len(cids))
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--seed-start", type=int, default=100)
    args = parser.parse_args()

    from rdkit import rdBase
    print(f"RDKit version: {rdBase.rdkitVersion}")
    print(f"Runs: {args.runs}, seed range: {args.seed_start}-{args.seed_start + args.runs - 1}")
    print()

    for block_name, run_fn in [("Block1 (default)", run_block1),
                                ("Block2 (ETKDG)", run_block2)]:
        print(f"=== {block_name} ===")
        all_counts = {i: [] for i in range(len(SMILES))}

        for run in range(args.runs):
            seed = args.seed_start + run
            counts = run_fn(seed)
            for i, c in enumerate(counts):
                all_counts[i].append(c)

        for i, smi in enumerate(SMILES):
            vals = all_counts[i]
            ctr = Counter(vals)
            print(f"  mol {i} ({smi[:40]}...)")
            print(f"    min={min(vals)} max={max(vals)} mean={sum(vals)/len(vals):.1f}")
            print(f"    distribution: {dict(sorted(ctr.items()))}")
        print()

    # Also run with the fixed seed=100 (the actual test seed) once
    print("=== Single run with seed=100 (test conditions) ===")
    print(f"  Block1: {run_block1(100)}")
    print(f"  Block2: {run_block2(100)}")


if __name__ == "__main__":
    main()
