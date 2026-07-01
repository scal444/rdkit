"""
Compare conformer quality between baseline and fixed ETKDG.

Generates conformers with ETKDG for the test6RmsPruning molecules,
evaluates each conformer with both MMFF and UFF (energy + max force),
and outputs a CSV for plotting.

Usage:
    # With fixed build:
    conda run -n rdkit_bfgs_experiment python compare_conformer_quality.py --label fixed -o fixed.csv

    # With baseline:
    conda run -n baseline_rdkit python compare_conformer_quality.py --label baseline -o baseline.csv

    # Then plot:
    python compare_conformer_quality.py --plot fixed.csv baseline.csv
"""

import argparse
import csv
import math
import sys

def generate(label, outfile):
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
    RDLogger.DisableLog("rdApp.*")

    smiles = [
        "CC(C)CC(NC(C1[N+]CCC1)=O)C([O-])=O",
        "CC(NC(CO)C(O)c1ccc([N+]([O-])=O)cc1)=O",
        "CC([N+])C(NC(C)C(N1C(C=O)CCC1)=O)=O",
        "CC(NC1C(O)C=C(C([O-])=O)OC1C(O)C(O)CO)=O",
        "CCCC=C(NC(C1CC1(C)C)=O)C([O-])=O",
        "OCC(O)C(O)C(Cn1c2c(cc(C)c(C)c2)nc-2c(=O)[nH]c(=O)nc12)O",
    ]

    params = rdDistGeom.ETKDG()
    params.randomSeed = 100
    params.maxIterations = 30
    params.pruneRmsThresh = 1.5
    params.useSymmetryForPruning = False

    rows = []
    for mol_idx, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 50, params)

        for conf_idx, cid in enumerate(cids):
            row = {
                "label": label,
                "rdkit_version": rdBase.rdkitVersion,
                "mol_idx": mol_idx,
                "smiles": smi,
                "conf_idx": conf_idx,
                "conf_id": cid,
            }

            mmff_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
            for ff_name, get_ff in [
                ("MMFF", lambda m, c: rdForceFieldHelpers.MMFFGetMoleculeForceField(m, mmff_props, confId=c) if mmff_props else None),
                ("UFF", lambda m, c: rdForceFieldHelpers.UFFGetMoleculeForceField(m, confId=c)),
            ]:
                ff = get_ff(mol, cid)
                if ff:
                    energy = ff.CalcEnergy()
                    grad = ff.CalcGrad()
                    max_force = max(
                        math.sqrt(grad[3*i]**2 + grad[3*i+1]**2 + grad[3*i+2]**2)
                        for i in range(mol.GetNumAtoms())
                    )
                    row[f"{ff_name}_energy"] = f"{energy:.4f}"
                    row[f"{ff_name}_max_force"] = f"{max_force:.4f}"
                else:
                    row[f"{ff_name}_energy"] = ""
                    row[f"{ff_name}_max_force"] = ""

            rows.append(row)

    fieldnames = ["label", "rdkit_version", "mol_idx", "smiles", "conf_idx",
                  "conf_id", "MMFF_energy", "MMFF_max_force", "UFF_energy", "UFF_max_force"]

    with open(outfile, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {outfile}")
    print(f"RDKit version: {rdBase.rdkitVersion}")
    nconfs = [sum(1 for r in rows if r["mol_idx"] == i) for i in range(len(smiles))]
    print(f"Conformer counts: {nconfs}")


def plot(files):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)

    labels = sorted(df["label"].unique())
    colors = {labels[0]: "tab:blue", labels[1]: "tab:orange"}
    mol_indices = sorted(df["mol_idx"].unique())
    n_mols = len(mol_indices)
    bar_width = 0.35

    metrics = [
        ("MMFF_energy", "MMFF Energy (kcal/mol)"),
        ("UFF_energy", "UFF Energy (kcal/mol)"),
        ("MMFF_max_force", "MMFF Max Force (kcal/mol/A)"),
        ("UFF_max_force", "UFF Max Force (kcal/mol/A)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for col_idx, (metric, ylabel) in enumerate(metrics):
        ax = axes[col_idx // 2][col_idx % 2]
        positions = []
        data = []
        tick_positions = []
        tick_labels = []
        box_colors = []
        for mol_idx in mol_indices:
            base = mol_idx * 3
            tick_positions.append(base + 0.5)
            tick_labels.append(str(mol_idx))
            for i, label in enumerate(labels):
                vals = pd.to_numeric(
                    df[(df["label"] == label) & (df["mol_idx"] == mol_idx)][metric],
                    errors="coerce"
                ).dropna().tolist()
                positions.append(base + i)
                data.append(vals if vals else [0])
                box_colors.append(colors[label])
        bp = ax.boxplot(data, positions=positions, widths=0.7,
                        patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.set_xlabel("Molecule index")
        # Divider between affected (mols 0-3) and unaffected (mols 4-5)
        divider_x = 3.5 * 3 + 1  # between mol 3 and mol 4
        ax.axvline(x=divider_x, color="gray", linestyle="--", linewidth=1)
        ylim = ax.get_ylim()
        ax.text(divider_x - 4, ylim[1] * 0.95, "hit degenerate torsion path",
                ha="center", va="top", fontsize=8, fontstyle="italic", color="gray")
        ax.text(divider_x + 4, ylim[1] * 0.95, "unaffected",
                ha="center", va="top", fontsize=8, fontstyle="italic", color="gray")

        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(facecolor=colors[l], alpha=0.7, label=l) for l in labels])

    plt.suptitle("Conformer quality: baseline vs fixed ETKDG gradient")
    plt.tight_layout()
    plt.savefig("conformer_quality.png", dpi=150)
    print("Saved conformer_quality.png")

    # Summary table
    mol_indices = sorted(df["mol_idx"].unique())
    print("\nPer-molecule summary:")
    print(f"{'label':<12} {'mol':>3} {'nconf':>5} {'MMFF_E_mean':>11} {'MMFF_F_mean':>11} {'UFF_E_mean':>11} {'UFF_F_mean':>11}")
    for label in labels:
        sub = df[df["label"] == label]
        for mol_idx in mol_indices:
            msub = sub[sub["mol_idx"] == mol_idx]
            n = len(msub)
            mmff_e = pd.to_numeric(msub["MMFF_energy"], errors="coerce").mean()
            mmff_f = pd.to_numeric(msub["MMFF_max_force"], errors="coerce").mean()
            uff_e = pd.to_numeric(msub["UFF_energy"], errors="coerce").mean()
            uff_f = pd.to_numeric(msub["UFF_max_force"], errors="coerce").mean()
            print(f"{label:<12} {mol_idx:>3} {n:>5} {mmff_e:>11.2f} {mmff_f:>11.2f} {uff_e:>11.2f} {uff_f:>11.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", help="Label for this run (e.g. 'fixed' or 'baseline')")
    parser.add_argument("-o", "--output", help="Output CSV file")
    parser.add_argument("--plot", nargs="+", help="CSV files to plot")
    args = parser.parse_args()

    if args.plot:
        plot(args.plot)
    elif args.label and args.output:
        generate(args.label, args.output)
    else:
        parser.print_help()
