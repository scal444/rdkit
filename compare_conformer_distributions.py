"""
Generate 100 conformers per molecule, evaluate with MMFF and UFF,
and compare distributions between baseline and fixed ETKDG.

Usage:
    conda run -n rdkit_bfgs_experiment python compare_conformer_distributions.py --label fixed -o fixed_100.csv
    conda run -n baseline_rdkit python compare_conformer_distributions.py --label baseline -o baseline_100.csv
    conda run -n rdkit_bfgs_experiment python compare_conformer_distributions.py --plot fixed_100.csv baseline_100.csv
"""

import argparse
import csv
import math

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
    params.randomSeed = 42

    rows = []
    for mol_idx, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 100, params)
        mmff_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)

        for conf_idx, cid in enumerate(cids):
            row = {
                "label": label,
                "mol_idx": mol_idx,
                "smiles": smi,
                "conf_idx": conf_idx,
            }
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

    fieldnames = ["label", "mol_idx", "smiles", "conf_idx",
                  "MMFF_energy", "MMFF_max_force", "UFF_energy", "UFF_max_force"]
    with open(outfile, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {outfile} ({rdBase.rdkitVersion})")
    nconfs = [sum(1 for r in rows if r["mol_idx"] == i) for i in range(len(smiles))]
    print(f"Conformer counts: {nconfs}")


def plot(files):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    dfs = [pd.read_csv(f) for f in files]
    df = pd.concat(dfs, ignore_index=True)
    labels = sorted(df["label"].unique())
    mol_indices = sorted(df["mol_idx"].unique())
    metrics = [
        ("MMFF_energy", "MMFF Energy (kcal/mol)"),
        ("UFF_energy", "UFF Energy (kcal/mol)"),
        ("MMFF_max_force", "MMFF Max Force (kcal/mol/A)"),
        ("UFF_max_force", "UFF Max Force (kcal/mol/A)"),
    ]

    fig, axes = plt.subplots(len(mol_indices), len(metrics),
                             figsize=(4 * len(metrics), 3 * len(mol_indices)),
                             squeeze=False)

    colors = {"baseline": "tab:blue", "fixed": "tab:orange"}

    for row, mol_idx in enumerate(mol_indices):
        for col, (metric, ylabel) in enumerate(metrics):
            ax = axes[row][col]
            all_vals = []
            for label in labels:
                vals = pd.to_numeric(
                    df[(df["label"] == label) & (df["mol_idx"] == mol_idx)][metric],
                    errors="coerce"
                ).dropna()
                all_vals.append(vals)
            if all(len(v) > 0 for v in all_vals):
                bins = np.histogram_bin_edges(pd.concat(all_vals), bins=20)
                for label, vals in zip(labels, all_vals):
                    ax.hist(vals, bins=bins, alpha=0.7, label=label,
                            color=colors.get(label), edgecolor="white", linewidth=0.5)
            else:
                for label, vals in zip(labels, all_vals):
                    if len(vals) > 0:
                        ax.hist(vals, bins=20, alpha=0.7, label=label,
                                color=colors.get(label), edgecolor="white", linewidth=0.5)
            if row == 0:
                ax.set_title(ylabel, fontsize=9)
            if col == 0:
                ax.set_ylabel(f"Molecule {mol_idx}", fontsize=9)
            ax.legend(fontsize=7)
            ax.tick_params(labelsize=7)

    plt.suptitle("100 ETKDG conformers: baseline vs fixed gradient", fontsize=12)
    plt.tight_layout()
    plt.savefig("conformer_distributions.png", dpi=150)
    print("Saved conformer_distributions.png")

    # Print summary
    print(f"\n{'label':<10} {'mol':>3} {'n':>4} {'MMFF_E':>10} {'MMFF_F':>10} {'UFF_E':>10} {'UFF_F':>10}")
    print(f"{'':10} {'':>3} {'':>4} {'mean±std':>10} {'mean±std':>10} {'mean±std':>10} {'mean±std':>10}")
    print("-" * 72)
    for label in labels:
        for mol_idx in mol_indices:
            sub = df[(df["label"] == label) & (df["mol_idx"] == mol_idx)]
            n = len(sub)
            parts = [f"{label:<10}", f"{mol_idx:>3}", f"{n:>4}"]
            for m in ["MMFF_energy", "MMFF_max_force", "UFF_energy", "UFF_max_force"]:
                vals = pd.to_numeric(sub[m], errors="coerce").dropna()
                if len(vals):
                    parts.append(f"{vals.mean():>5.1f}±{vals.std():>4.1f}")
                else:
                    parts.append(f"{'N/A':>10}")
            print(" ".join(parts))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--label")
    parser.add_argument("-o", "--output")
    parser.add_argument("--plot", nargs="+")
    args = parser.parse_args()

    if args.plot:
        plot(args.plot)
    elif args.label and args.output:
        generate(args.label, args.output)
    else:
        parser.print_help()
