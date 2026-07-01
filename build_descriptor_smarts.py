#!/usr/bin/env python3
"""Generate descriptor_smarts.py from RDKit source data (no RDKit import required)."""

from __future__ import annotations

import ast
import pprint
import re
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent
QED_PY = REPO / "rdkit/Chem/QED.py"
FRAGMENTS_CSV = REPO / "Data/FragmentDescriptors.csv"
CRIPPEN_TXT = REPO / "Data/Crippen.txt"
LIPINSKI_CPP = REPO / "Code/GraphMol/Descriptors/Lipinski.cpp"
OUTPUT = REPO / "descriptor_smarts.py"

# How each group invokes its SMARTS (from RDKit source). Keys must match output dict keys.
INVOCATION_SPECS: dict[str, dict[str, Any]] = {
    "QED_HBA": {
        "callable": "QED.properties -> HBA",
        "source": "rdkit/Chem/QED.py",
        "prefilter": "mol.HasSubstructMatch(pattern)",
        "match": "mol.GetSubstructMatches(pattern)",
        "match_kwargs": {
            "uniquify": True,
            "useChirality": False,
            "useQueryQueryMatches": False,
            "maxMatches": 1000,
        },
        "per_pattern": "0 if not prefilter else len(match)",
        "aggregate": "sum(per_pattern)",
        "result_kind": "count",
        "notes": (
            "Counts acceptor matches; overlapping patterns can contribute multiple times. "
            "Not the same SMARTS as Lipinski NumHBA / CalcNumHBA."
        ),
    },
    "QED_AROM": {
        "callable": "QED.properties -> AROM",
        "source": "rdkit/Chem/QED.py",
        "prefilter": None,
        "match": "Chem.DeleteSubstructs(mol, query_mol)",
        "match_kwargs": None,
        "per_pattern": "remove substructure matching SMARTS, then Chem.GetSSSR on remainder",
        "aggregate": "len(GetSSSR(...))",
        "result_kind": "ring_count",
        "notes": "Not HasSubstructMatch; aliphatic rings removed before SSSR aromatic-ring count.",
    },
    "QED_ALERTS": {
        "callable": "QED.properties -> ALERTS",
        "source": "rdkit/Chem/QED.py",
        "prefilter": None,
        "match": "mol.HasSubstructMatch(alert)",
        "match_kwargs": {
            "recursionPossible": True,
            "useChirality": False,
            "useQueryQueryMatches": False,
        },
        "per_pattern": "1 if match else 0",
        "aggregate": "sum(per_pattern)",
        "result_kind": "patterns_matched",
        "notes": (
            "One boolean per alert SMARTS (counts how many alert types hit, "
            "not total substructure occurrences)."
        ),
    },
    "Crippen": {
        "callable": "getCrippenAtomContribs / MolLogP / MolMR / MOE VSA",
        "source": "Code/GraphMol/Descriptors/Crippen.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern, matches, recursionPossible=False, useChirality=True)",
        "match_kwargs": {
            "recursionPossible": False,
            "useChirality": True,
            "returns_all_matches": True,
        },
        "per_pattern": "first atom index of each match; assign logP/MR to atom if still unassigned",
        "aggregate": "first matching pattern in Crippen.txt order wins per atom; stop when all atoms set",
        "result_kind": "per_atom_contrib",
        "notes": "Patterns tried in file order; not a match count descriptor.",
    },
    "Lipinski_NumHBD": {
        "callable": "CalcNumHBD / NumHDonors",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
        "notes": "C++ SMARTSCOUNTFUNC; differs from QED HBD which uses rdMolDescriptors.CalcNumHBD.",
    },
    "Lipinski_NumHBA": {
        "callable": "CalcNumHBA / NumHAcceptors",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
        "notes": "Different SMARTS and counting than QED_HBA acceptor list.",
    },
    "Lipinski_NumHeteroatoms": {
        "callable": "CalcNumHeteroatoms / NumHeteroatoms",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
    },
    "Lipinski_NumAmideBonds": {
        "callable": "CalcNumAmideBonds / NumAmideBonds",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
    },
    "Lipinski_NumRotatableBonds_nonstrict": {
        "callable": "CalcNumRotatableBonds(..., NumRotatableBondsOptions.NonStrict)",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
        "notes": "Chem.Descriptors NumRotatableBonds default may be NonStrict or Strict at build time.",
    },
    "Lipinski_NumRotatableBonds_strict": {
        "callable": (
            "CalcNumRotatableBonds(..., NumRotatableBondsOptions.Strict); "
            "QED.properties ROTB uses this mode"
        ),
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "per_pattern (single pattern)",
        "result_kind": "count",
    },
    "Lipinski_NumRotatableBonds_v3_rotBonds": {
        "callable": "CalcNumRotatableBonds(..., NumRotatableBondsOptions.StrictLinkages)",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "role": "base_count",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "res = base_count; then subtract symRings, terminalTripleBonds, distinct nonRingAmides",
        "result_kind": "count",
        "notes": "StrictLinkages mode; composite with other v3_* patterns in this file.",
    },
    "Lipinski_NumRotatableBonds_v3_symRings": {
        "callable": "CalcNumRotatableBonds(..., StrictLinkages)",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "role": "subtract_count",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "res -= per_pattern",
        "result_kind": "count_adjustment",
    },
    "Lipinski_NumRotatableBonds_v3_terminalTripleBonds": {
        "callable": "CalcNumRotatableBonds(..., StrictLinkages)",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "role": "subtract_count",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "len(matches)",
        "aggregate": "res -= per_pattern",
        "result_kind": "count_adjustment",
    },
    "Lipinski_NumRotatableBonds_v3_nonRingAmides": {
        "callable": "CalcNumRotatableBonds(..., StrictLinkages)",
        "source": "Code/GraphMol/Descriptors/Lipinski.cpp",
        "role": "subtract_distinct_matches",
        "prefilter": None,
        "match": "SubstructMatch(mol, pattern) -> all MatchVectType",
        "match_kwargs": {"returns_all_matches": True},
        "per_pattern": "for each match, if no atom overlap with prior amide matches: res -= 1",
        "aggregate": "subtract at most once per distinct non-ring amide match",
        "result_kind": "count_adjustment",
    },
}

FRAGMENT_INVOCATION: dict[str, Any] = {
    "callable": "Chem.Fragments.<name>",
    "source": "rdkit/Chem/Fragments.py",
    "prefilter": None,
    "match": "mol.GetSubstructMatches(pattern, uniquify=countUnique)",
    "match_kwargs": {
        "uniquify": True,
        "useChirality": False,
        "useQueryQueryMatches": False,
        "maxMatches": 1000,
        "countUnique_default": True,
    },
    "per_pattern": "len(match)",
    "aggregate": "per_pattern (single pattern per fr_* descriptor)",
    "result_kind": "count",
}


def string_from_ast(node: ast.AST) -> str:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return string_from_ast(node.left) + string_from_ast(node.right)
    raise ValueError(f"unsupported AST node for string: {ast.dump(node)}")


def list_of_strings_from_ast(node: ast.AST) -> list[str]:
    if not isinstance(node, ast.List):
        raise ValueError("expected list literal")
    return [string_from_ast(element) for element in node.elts]


def assign_list_from_module(path: Path, name: str) -> list[str]:
    module = ast.parse(path.read_text())
    for statement in module.body:
        if not isinstance(statement, ast.Assign):
            continue
        for target in statement.targets:
            if isinstance(target, ast.Name) and target.id == name:
                return list_of_strings_from_ast(statement.value)
    raise KeyError(f"{name} not found in {path}")


def assign_string_from_module(path: Path, name: str) -> str:
    module = ast.parse(path.read_text())
    for statement in module.body:
        if not isinstance(statement, ast.Assign):
            continue
        for target in statement.targets:
            if isinstance(target, ast.Name) and target.id == name:
                call = statement.value
                if (
                    isinstance(call, ast.Call)
                    and isinstance(call.func, ast.Attribute)
                    and call.func.attr == "MolFromSmarts"
                    and call.args
                ):
                    return string_from_ast(call.args[0])
    raise KeyError(f"{name} not found in {path}")


def extract_qed() -> dict[str, list[str]]:
    return {
        "QED_HBA": assign_list_from_module(QED_PY, "AcceptorSmarts"),
        "QED_AROM": [assign_string_from_module(QED_PY, "AliphaticRings")],
        "QED_ALERTS": assign_list_from_module(QED_PY, "StructuralAlertSmarts"),
    }


def extract_fragments() -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    with FRAGMENTS_CSV.open(newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3 or not parts[0].startswith("fr_"):
                continue
            result[parts[0]] = [parts[2].strip()]
    return result


def extract_crippen() -> dict[str, list[str]]:
    patterns: list[str] = []
    with CRIPPEN_TXT.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2 or parts[1] == "SMARTS":
                continue
            smarts = parts[1].strip()
            if smarts:
                patterns.append(smarts)
    return {"Crippen": patterns}


def join_cpp_string_literals(lines: list[str]) -> str:
    parts: list[str] = []
    for line in lines:
        while '"' in line:
            start = line.index('"')
            end = line.index('"', start + 1)
            parts.append(line[start + 1 : end])
            line = line[end + 1 :]
    return "".join(parts)


def extract_lipinski_cpp() -> dict[str, list[str]]:
    lines = LIPINSKI_CPP.read_text().splitlines()
    text = "\n".join(lines)
    result: dict[str, list[str]] = {}

    for macro_name, pattern in re.findall(
        r"SMARTSCOUNTFUNC\((\w+),\s*\"((?:[^\"\\]|\\.)*)\"",
        text,
        re.DOTALL,
    ):
        if macro_name == "NumHBA":
            continue
        result[f"Lipinski_{macro_name}"] = [pattern]

    hba_match = re.search(
        r'SMARTSCOUNTFUNC\(NumHBA,\s*("(?:[^"\\]|\\.)*"(?:\s*"[^"\\]*")*)',
        text,
        re.DOTALL,
    )
    if hba_match:
        result["Lipinski_NumHBA"] = [
            join_cpp_string_literals(hba_match.group(1).splitlines())
        ]

    non_strict = re.search(
        r'if \(strict == NonStrict\) \{\s*std::string pattern = "([^"]+)";',
        text,
    )
    if non_strict:
        result["Lipinski_NumRotatableBonds_nonstrict"] = [non_strict.group(1)]

    strict_literal_lines: list[str] = []
    capture_strict = False
    for line in lines:
        if "std::string strict_pattern =" in line:
            capture_strict = True
        if capture_strict:
            if "pattern_flyweight m(strict_pattern)" in line:
                break
            strict_literal_lines.append(line)
    if strict_literal_lines:
        result["Lipinski_NumRotatableBonds_strict"] = [
            join_cpp_string_literals(strict_literal_lines)
        ]

    for var_name, pattern in re.findall(
        r'pattern_flyweight (\w+)\("([^"]+)"\);',
        text,
    ):
        if var_name.endswith("_matcher"):
            key = f"Lipinski_NumRotatableBonds_v3_{var_name.removesuffix('_matcher')}"
            result[key] = [pattern]

    sym_match = re.search(
        r'pattern_flyweight symRings_matcher\(\s*((?:"[^"]*"\s*)+)\);',
        text,
        re.DOTALL,
    )
    if sym_match:
        result["Lipinski_NumRotatableBonds_v3_symRings"] = [
            join_cpp_string_literals(sym_match.group(1).splitlines())
        ]

    return result


def invocation_for(key: str) -> dict[str, Any]:
    if key in INVOCATION_SPECS:
        return dict(INVOCATION_SPECS[key])
    if key.startswith("fr_"):
        spec = dict(FRAGMENT_INVOCATION)
        spec["callable"] = f"Chem.Fragments.{key}"
        spec["chem_descriptor"] = key
        return spec
    raise KeyError(f"No invocation spec for {key}")


def deduplicate(
    data: dict[str, list[str]], key_order: list[str]
) -> dict[str, list[str]]:
    claimed: set[str] = set()
    result: dict[str, list[str]] = {}
    skipped = 0
    for key in key_order:
        if key not in data:
            continue
        unique: list[str] = []
        for pattern in data[key]:
            if pattern in claimed:
                skipped += 1
                continue
            claimed.add(pattern)
            unique.append(pattern)
        if unique:
            result[key] = unique
    if skipped:
        print(f"  skipped {skipped} duplicate SMARTS (already assigned to an earlier key)")
    return result


def build_key_order(raw: dict[str, list[str]]) -> list[str]:
    qed_order = ["QED_HBA", "QED_AROM", "QED_ALERTS"]
    fragment_order = sorted(key for key in raw if key.startswith("fr_"))
    lipinski_order = sorted(key for key in raw if key.startswith("Lipinski_"))
    return qed_order + fragment_order + ["Crippen"] + lipinski_order


def format_python(groups: dict[str, dict[str, Any]]) -> str:
    lines = [
        '"""SMARTS patterns used by Chem.Descriptors-related APIs.',
        "",
        "Each SMARTS string appears in exactly one group. Each group records how",
        "RDKit invokes those patterns (HasSubstructMatch vs GetSubstructMatches vs",
        "SubstructMatch vs DeleteSubstructs, and how results are aggregated).",
        "",
        "Auto-generated by build_descriptor_smarts.py.",
        '"""',
        "",
        "from __future__ import annotations",
        "",
        "from typing import Any, TypedDict",
        "",
        "",
        "class SmartsInvocation(TypedDict, total=False):",
        '    """How a descriptor group uses its SMARTS patterns."""',
        "    callable: str",
        "    source: str",
        "    chem_descriptor: str",
        "    role: str",
        "    prefilter: str | None",
        "    match: str",
        "    match_kwargs: dict[str, Any] | None",
        "    per_pattern: str",
        "    aggregate: str",
        "    result_kind: str",
        "    notes: str",
        "",
        "",
        "class SmartsGroup(TypedDict):",
        "    smarts: list[str]",
        "    invocation: SmartsInvocation",
        "",
        "",
        "DESCRIPTOR_SMARTS: dict[str, SmartsGroup] = {",
    ]
    for key, group in groups.items():
        lines.append(f'    "{key}": {{')
        lines.append('        "smarts": [')
        for pattern in group["smarts"]:
            lines.append(f"            {pattern!r},")
        lines.append("        ],")
        invocation_text = pprint.pformat(
            group["invocation"], width=96, sort_dicts=False
        )
        inv_lines = invocation_text.splitlines()
        lines.append(f'        "invocation": {inv_lines[0]}')
        for inv_line in inv_lines[1:]:
            lines.append(f"            {inv_line}")
        lines.append("    },")
    lines.append("}")
    lines.append("")
    lines.extend(
        [
            "",
            "def smarts_list(key: str) -> list[str]:",
            '    """Return SMARTS strings for a descriptor group key."""',
            "    return DESCRIPTOR_SMARTS[key][\"smarts\"]",
            "",
            "",
            "def smarts_invocation(key: str) -> SmartsInvocation:",
            '    """Return invocation metadata for a descriptor group key."""',
            "    return DESCRIPTOR_SMARTS[key][\"invocation\"]",
            "",
        ]
    )
    return "\n".join(lines)


def validate_smarts_compile(groups: dict[str, dict[str, Any]]) -> list[str]:
    """Return list of 'key: smarts' for patterns that fail MolFromSmarts."""
    try:
        from rdkit import Chem
    except ImportError:
        print("  (skip SMARTS compile check: RDKit not importable)")
        return []

    failures: list[str] = []
    for key, group in groups.items():
        for smarts in group["smarts"]:
            query = Chem.MolFromSmarts(smarts)
            if query is None or query.GetNumAtoms() == 0:
                failures.append(f"{key}: {smarts[:80]}...")
    return failures


def main() -> None:
    raw: dict[str, list[str]] = {}
    raw.update(extract_qed())
    raw.update(extract_fragments())
    raw.update(extract_crippen())
    raw.update(extract_lipinski_cpp())

    key_order = build_key_order(raw)
    smarts_only = deduplicate(raw, key_order)

    groups: dict[str, dict[str, Any]] = {}
    for key, patterns in smarts_only.items():
        groups[key] = {
            "smarts": patterns,
            "invocation": invocation_for(key),
        }

    failures = validate_smarts_compile(groups)
    if failures:
        raise SystemExit(
            f"SMARTS compile check failed ({len(failures)} patterns):\n"
            + "\n".join(failures[:10])
        )

    OUTPUT.write_text(format_python(groups))
    total = sum(len(group["smarts"]) for group in groups.values())
    print(f"Wrote {len(groups)} keys, {total} unique SMARTS to {OUTPUT}")


if __name__ == "__main__":
    main()
