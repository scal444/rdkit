from rdkit import Chem, RDLogger
import re
RDLogger.DisableLog("rdApp.*")

smiles = [
    "CC(C)CC(NC(C1[N+]CCC1)=O)C([O-])=O",
    "CC(NC(CO)C(O)c1ccc([N+]([O-])=O)cc1)=O",
    "CC([N+])C(NC(C)C(N1C(C=O)CCC1)=O)=O",
    "CC(NC1C(O)C=C(C([O-])=O)OC1C(O)C(O)CO)=O",
    "CCCC=C(NC(C1CC1(C)C)=O)C([O-])=O",
    "OCC(O)C(O)C(Cn1c2c(cc(C)c(C)c2)nc-2c(=O)[nH]c(=O)nc12)O",
]

with open("Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionPreferences_v2.in") as f:
    raw = f.read()

v6_patterns = []
for line in raw.split("\\n"):
    line = line.strip().strip('"').strip()
    parts = line.split()
    if len(parts) >= 13:
        smarts = parts[0]
        try:
            v6 = float(parts[12])
        except ValueError:
            continue
        if abs(v6) > 1e-10:
            v6_patterns.append((smarts, v6))

print(f"Total V6 patterns: {len(v6_patterns)}")
print()

for i, smi in enumerate(smiles):
    mol = Chem.MolFromSmiles(smi)
    matches = []
    for sm, v6 in v6_patterns:
        try:
            q = Chem.MolFromSmarts(sm)
            if q and mol.HasSubstructMatch(q):
                n = len(mol.GetSubstructMatches(q))
                matches.append((sm, v6, n))
        except Exception:
            pass
    print(f"mol {i}: {len(matches)} V6 pattern matches  {smi[:50]}")
    for sm, v6, n in matches:
        print(f"  V6={v6:6.2f} x{n}  {sm}")
    print()
