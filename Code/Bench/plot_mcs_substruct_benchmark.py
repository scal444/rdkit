#!/usr/bin/env python3
"""Plot mcs_substruct_benchmark.py CSV output as paired scatter plots."""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import tempfile
import textwrap
from argparse import Namespace
from collections import OrderedDict
from pathlib import Path


SIZE_SUITE_X_KEYS = ("smaller_atoms", "larger_atoms", "combined_atoms")
MCS_SIZE_SUITE_X_KEYS = ("reference_mcs_atoms", "reference_mcs_bonds")
AXIS_LABELS = {
  "query_atoms": "query heavy atoms",
  "target_atoms": "target heavy atoms",
  "total_atoms": "combined heavy atoms",
  "smaller_atoms": "smaller molecule heavy atoms",
  "larger_atoms": "larger molecule heavy atoms",
  "combined_atoms": "combined heavy atoms",
  "atom_delta": "heavy atom count delta",
  "mcs_atoms": "MCS atoms",
  "mcs_bonds": "MCS bonds",
  "reference_mcs_atoms": "reference MCS atoms",
  "reference_mcs_bonds": "reference MCS bonds",
}
BREAKDOWNS = OrderedDict([
  (
    "baseline",
    (
      "Substruct raw SMILES vs MCS elements/order",
      ["substruct_smiles", "mcs_elements_order"],
    ),
  ),
  (
    "mcs-core",
    (
      "MCS core comparison",
      [
        "substruct_smiles",
        "mcs_any_atom_bond",
        "mcs_elements_order",
        "mcs_elements_order_maximize_atoms",
      ],
    ),
  ),
  (
    "mcs-rings",
    (
      "MCS ring comparison",
      [
        "substruct_smiles",
        "substruct_ring_atoms",
        "mcs_elements_order",
        "mcs_elements_order_ring_atoms",
        "mcs_elements_order_whole_rings",
        "mcs_elements_order_ring_atoms_whole_rings",
      ],
    ),
  ),
  (
    "substruct",
    (
      "Substructure query comparison",
      [
        "substruct_smiles",
        "substruct_elements_any_bond",
        "substruct_any_atom_bond",
        "substruct_ring_atoms",
      ],
    ),
  ),
])

BASELINE_COMPARISON_METHODS = ["substruct_smiles", "mcs_elements_order"]
BASELINE_COMPARISON_LABELS = {
  "substruct_smiles": "Substructure Search",
  "mcs_elements_order": "MCS",
}


def import_matplotlib():
  if "MPLCONFIGDIR" not in os.environ:
    mpl_cache = Path(tempfile.gettempdir()) / "rdkit_mplconfig"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_cache)
  if "XDG_CACHE_HOME" not in os.environ:
    xdg_cache = Path(tempfile.gettempdir()) / "rdkit_xdg_cache"
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ["XDG_CACHE_HOME"] = str(xdg_cache)
  try:
    import matplotlib.pyplot as plt
  except ImportError as exc:
    raise SystemExit(
      "matplotlib is required for plotting. Install it in the Python environment "
      "used to run this script.") from exc
  return plt


def as_float(row: dict[str, str], key: str) -> float | None:
  value = row.get(key, "")
  if value == "":
    return None
  try:
    return float(value)
  except ValueError:
    return None


def as_bool(row: dict[str, str], key: str) -> bool:
  return row.get(key, "").strip().lower() in {"1", "true", "yes", "y"}


def read_rows(path: Path) -> list[dict[str, str]]:
  with path.open(newline="") as handle:
    return list(csv.DictReader(handle))


def ensure_atom_columns(rows: list[dict[str, str]]) -> None:
  for row in rows:
    query_atoms = as_float(row, "query_atoms")
    target_atoms = as_float(row, "target_atoms")
    if query_atoms is None or target_atoms is None:
      continue

    smaller_atoms = min(query_atoms, target_atoms)
    larger_atoms = max(query_atoms, target_atoms)
    combined_atoms = query_atoms + target_atoms

    row.setdefault("smaller_atoms", f"{smaller_atoms:g}")
    row.setdefault("larger_atoms", f"{larger_atoms:g}")
    row.setdefault("combined_atoms", f"{combined_atoms:g}")
    if not row.get("smaller_atoms"):
      row["smaller_atoms"] = f"{smaller_atoms:g}"
    if not row.get("larger_atoms"):
      row["larger_atoms"] = f"{larger_atoms:g}"
    if not row.get("combined_atoms"):
      row["combined_atoms"] = row.get("total_atoms") or f"{combined_atoms:g}"


def ensure_reference_mcs_columns(rows: list[dict[str, str]], method: str) -> None:
  atoms_by_pair: dict[str, float] = {}
  bonds_by_pair: dict[str, float] = {}
  for row in rows:
    if row.get("method") != method:
      continue
    pair_id = row.get("pair_id")
    atoms = as_float(row, "mcs_atoms")
    bonds = as_float(row, "mcs_bonds")
    if pair_id and atoms is not None:
      atoms_by_pair[pair_id] = atoms
    if pair_id and bonds is not None:
      bonds_by_pair[pair_id] = bonds

  for row in rows:
    pair_id = row.get("pair_id")
    if pair_id in atoms_by_pair:
      row["reference_mcs_atoms"] = f"{atoms_by_pair[pair_id]:g}"
    if pair_id in bonds_by_pair:
      row["reference_mcs_bonds"] = f"{bonds_by_pair[pair_id]:g}"


def filter_rows(
  rows: list[dict[str, str]],
  family: str,
  methods: set[str] | None,
) -> list[dict[str, str]]:
  filtered = []
  for row in rows:
    if row.get("status") != "ok":
      continue
    if family != "all" and row.get("family") != family:
      continue
    if methods is not None and row.get("method") not in methods:
      continue
    if as_float(row, "seconds_per_call") is None:
      continue
    filtered.append(row)
  return filtered


def method_label_map(rows: list[dict[str, str]]) -> dict[str, str]:
  labels: dict[str, str] = {}
  for row in rows:
    method = row["method"]
    labels.setdefault(method, row.get("label") or method)
  return labels


def ordered_methods(
  rows: list[dict[str, str]],
  method_order: list[str] | None = None,
) -> OrderedDict[str, str]:
  labels = method_label_map(rows)
  methods: OrderedDict[str, str] = OrderedDict()
  if method_order is not None:
    present = {row["method"] for row in rows}
    for method in method_order:
      if method in present:
        methods[method] = labels.get(method, method)
    return methods

  for row in rows:
    method = row["method"]
    methods.setdefault(method, labels.get(method, method))
  return methods


def highlighted_methods(rows: list[dict[str, str]]) -> set[str]:
  return {row["method"] for row in rows if as_bool(row, "highlighted")}


def default_ratio_baseline(rows: list[dict[str, str]]) -> str | None:
  for row in rows:
    if row.get("method") == "substruct_smiles":
      return "substruct_smiles"
  highlighted = highlighted_methods(rows)
  if highlighted:
    return sorted(highlighted)[0]
  return None


def build_baseline_map(
  rows: list[dict[str, str]],
  baseline_method: str,
) -> dict[str, float]:
  baseline: dict[str, float] = {}
  for row in rows:
    if row.get("method") != baseline_method:
      continue
    value = as_float(row, "seconds_per_call")
    if value is not None and value > 0:
      baseline[row["pair_id"]] = value
  return baseline


def method_color_map(methods: list[str], plt):
  cmap = plt.get_cmap("tab20")
  return {method: cmap(i % cmap.N) for i, method in enumerate(methods)}


def scatter_methods(
  ax,
  rows: list[dict[str, str]],
  methods: OrderedDict[str, str],
  highlighted: set[str],
  colors: dict[str, object],
  x_key: str,
  y_key: str,
  y_transform=lambda row, value: value,
  annotate_highlight: bool = True,
):
  for method, label in methods.items():
    method_rows = [row for row in rows if row.get("method") == method]
    xs = []
    ys = []
    for row in method_rows:
      x = as_float(row, x_key)
      y = as_float(row, y_key)
      if x is None or y is None:
        continue
      transformed = y_transform(row, y)
      if transformed is None:
        continue
      xs.append(x)
      ys.append(transformed)

    if not xs:
      continue

    is_highlighted = method in highlighted
    legend_label = (f"{label} (baseline)"
                    if is_highlighted and annotate_highlight else label)
    ax.scatter(
      xs,
      ys,
      s=54 if is_highlighted else 28,
      alpha=0.9 if is_highlighted else 0.7,
      color=colors[method],
      edgecolors="black" if is_highlighted else "none",
      linewidths=0.8 if is_highlighted else 0.0,
      label=legend_label,
    )


def binned_log_stats(
  rows: list[dict[str, str]],
  x_key: str,
  y_value,
  bin_width: float,
) -> list[tuple[float, float, float, float, int]]:
  bins: dict[int, list[float]] = {}
  for row in rows:
    x = as_float(row, x_key)
    y = y_value(row)
    if x is None or y is None or y <= 0:
      continue
    bin_idx = math.floor(x / bin_width)
    bins.setdefault(bin_idx, []).append(math.log(y))

  stats: list[tuple[float, float, float, float, int]] = []
  for bin_idx in sorted(bins):
    values = bins[bin_idx]
    n = len(values)
    mean = sum(values) / n
    if n > 1:
      variance = sum((value - mean) ** 2 for value in values) / (n - 1)
      half_width = 1.96 * math.sqrt(variance) / math.sqrt(n)
    else:
      half_width = 0.0
    center = (bin_idx + 0.5) * bin_width
    stats.append((center, math.exp(mean), math.exp(mean - half_width),
                  math.exp(mean + half_width), n))
  return stats


def line_methods_with_ci(
  ax,
  rows: list[dict[str, str]],
  methods: OrderedDict[str, str],
  highlighted: set[str],
  colors: dict[str, object],
  x_key: str,
  y_value,
  bin_width: float,
):
  for method, label in methods.items():
    method_rows = [row for row in rows if row.get("method") == method]
    stats = binned_log_stats(method_rows, x_key, y_value, bin_width)
    if not stats:
      continue

    xs = [item[0] for item in stats]
    means = [item[1] for item in stats]
    lows = [item[2] for item in stats]
    highs = [item[3] for item in stats]
    is_highlighted = method in highlighted
    legend_label = f"{label} (baseline)" if is_highlighted else label

    ax.fill_between(xs, lows, highs, color=colors[method], alpha=0.16,
                    linewidth=0)
    ax.plot(xs, means, marker="o", markersize=5 if is_highlighted else 4,
            linewidth=2.2 if is_highlighted else 1.6,
            color=colors[method], label=legend_label)


def x_keys_for_args(args: argparse.Namespace) -> list[str]:
  if args.x == "size-suite":
    return list(SIZE_SUITE_X_KEYS)
  if args.x == "mcs-size-suite":
    return list(MCS_SIZE_SUITE_X_KEYS)
  return [args.x]


def axis_label(x_key: str) -> str:
  return AXIS_LABELS.get(x_key, x_key.replace("_", " "))


def scale_list(selection: str) -> list[str]:
  if selection == "both":
    return ["log", "linear"]
  return [selection]


def output_for_breakdown(args: argparse.Namespace, breakdown: str) -> Path:
  output_dir = args.output_dir or args.output.parent
  suffix = args.output.suffix or ".png"
  return output_dir / f"{args.output.stem}_{breakdown}{suffix}"


def make_plot(
  args: argparse.Namespace,
  rows: list[dict[str, str]],
  all_rows: list[dict[str, str]],
  method_order: list[str] | None = None,
) -> None:
  plt = import_matplotlib()
  methods = ordered_methods(rows, method_order=method_order)
  method_names = list(methods)
  highlighted = highlighted_methods(rows)
  colors = method_color_map(method_names, plt)

  baseline_method = args.ratio_baseline or default_ratio_baseline(rows)
  baseline = build_baseline_map(all_rows, baseline_method) if baseline_method else {}
  label_map = method_label_map(all_rows)
  show_ratio = args.ratio and bool(baseline)
  x_keys = x_keys_for_args(args)
  runtime_scales = scale_list(args.runtime_scales)
  ratio_scales = scale_list(args.ratio_scales) if show_ratio else []
  nrows = len(runtime_scales) + len(ratio_scales)
  ncols = len(x_keys)

  height_ratios = [3] * len(runtime_scales) + [2] * len(ratio_scales)
  gridspec_kw = {"height_ratios": height_ratios}
  fig, axes = plt.subplots(nrows, ncols, figsize=args.figsize, squeeze=False,
                           sharey="row", gridspec_kw=gridspec_kw)

  if args.title:
    title = args.title
  else:
    title = "MCS vs small-into-large substructure runtime"
  fig.suptitle(title)

  def ratio_transform(row: dict[str, str], value: float) -> float | None:
    denominator = baseline.get(row["pair_id"])
    if denominator is None or denominator <= 0:
      return None
    return value / denominator

  row_idx = 0
  for runtime_scale in runtime_scales:
    for col, x_key in enumerate(x_keys):
      runtime_ax = axes[row_idx][col]
      scatter_methods(
        ax=runtime_ax,
        rows=rows,
        methods=methods,
        highlighted=highlighted,
        colors=colors,
        x_key=x_key,
        y_key="seconds_per_call",
      )
      if row_idx == 0:
        runtime_ax.set_title(axis_label(x_key))
      runtime_ax.set_yscale(runtime_scale)
      runtime_ax.grid(True, which="both", alpha=0.25)
      if col == 0:
        runtime_ax.set_ylabel(f"seconds per call\n({runtime_scale})")
    row_idx += 1

  ratio_rows = rows
  if not args.show_ratio_baseline and baseline_method:
    ratio_rows = [row for row in rows if row.get("method") != baseline_method]

  for ratio_scale in ratio_scales:
    for col, x_key in enumerate(x_keys):
      ratio_ax = axes[row_idx][col]
      scatter_methods(
        ax=ratio_ax,
        rows=ratio_rows,
        methods=methods,
        highlighted=highlighted,
        colors=colors,
        x_key=x_key,
        y_key="seconds_per_call",
        y_transform=ratio_transform,
      )
      ratio_label = label_map.get(baseline_method, baseline_method)
      ratio_ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--")
      ratio_ax.set_yscale(ratio_scale)
      ratio_ax.grid(True, which="both", alpha=0.25)
      if col == 0:
        ratio_ax.set_ylabel(f"runtime / {ratio_label}\n({ratio_scale})")
    row_idx += 1

  for col, x_key in enumerate(x_keys):
    axes[-1][col].set_xlabel(axis_label(x_key))

  handles, labels = axes[0][0].get_legend_handles_labels()
  labels = [textwrap.fill(label, width=34) for label in labels]
  fig.legend(handles, labels, title="Method", loc="center left",
             bbox_to_anchor=(0.80, 0.5), frameon=True, fontsize="small",
             title_fontsize="small", borderaxespad=0.0, labelspacing=0.85)
  fig.tight_layout(rect=(0.0, 0.0, 0.78, 0.95))

  args.output.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(args.output, dpi=args.dpi)
  print(f"Wrote {args.output}")


def make_baseline_comparison_plot(
  args: argparse.Namespace,
  rows: list[dict[str, str]],
  all_rows: list[dict[str, str]],
) -> None:
  """Plot the combined-size baseline runtime and MCS/substructure ratio panels."""
  plt = import_matplotlib()
  methods = OrderedDict(
    (method, BASELINE_COMPARISON_LABELS[method])
    for method in BASELINE_COMPARISON_METHODS
  )
  # Preserve the exact colors assigned in the mcs-core breakdown, where the
  # two intervening MCS variants occupy the light-blue and light-orange slots.
  mcs_core_methods = BREAKDOWNS["mcs-core"][1]
  mcs_core_colors = method_color_map(mcs_core_methods, plt)
  colors = {method: mcs_core_colors[method] for method in methods}
  baseline = build_baseline_map(all_rows, "substruct_smiles")

  fig, axes = plt.subplots(
    2,
    1,
    figsize=args.figsize,
    sharex=True,
    gridspec_kw={"height_ratios": [3, 2]},
  )
  fig.suptitle("MCS/Substruct search time comparison")

  scatter_methods(
    ax=axes[0],
    rows=rows,
    methods=methods,
    highlighted=set(methods),
    colors=colors,
    x_key="combined_atoms",
    y_key="seconds_per_call",
    annotate_highlight=False,
  )
  axes[0].set_yscale("log")
  axes[0].set_ylabel("seconds per call")
  axes[0].grid(True, which="both", alpha=0.25)

  def ratio_transform(row: dict[str, str], value: float) -> float | None:
    denominator = baseline.get(row["pair_id"])
    if denominator is None or denominator <= 0:
      return None
    return value / denominator

  mcs_rows = [row for row in rows if row.get("method") == "mcs_elements_order"]
  scatter_methods(
    ax=axes[1],
    rows=mcs_rows,
    methods=OrderedDict([("mcs_elements_order", "MCS")]),
    highlighted={"mcs_elements_order"},
    colors=colors,
    x_key="combined_atoms",
    y_key="seconds_per_call",
    y_transform=ratio_transform,
    annotate_highlight=False,
  )
  axes[1].axhline(1.0, color="0.35", linewidth=1.0, linestyle="--")
  axes[1].set_yscale("log")
  axes[1].set_ylabel("MCS / Substructure Search\nruntime")
  axes[1].set_xlabel("combined heavy atom count")
  axes[1].grid(True, which="both", alpha=0.25)

  handles, labels = axes[0].get_legend_handles_labels()
  axes[0].legend(handles, labels, loc="upper left", frameon=True,
                 fontsize="small")
  fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))

  args.output.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(args.output, dpi=args.dpi, facecolor="white", transparent=False)
  print(f"Wrote {args.output}")


def make_average_plot(
  args: argparse.Namespace,
  rows: list[dict[str, str]],
  all_rows: list[dict[str, str]],
  method_order: list[str] | None = None,
) -> None:
  plt = import_matplotlib()
  methods = ordered_methods(rows, method_order=method_order)
  method_names = list(methods)
  highlighted = highlighted_methods(rows)
  colors = method_color_map(method_names, plt)

  baseline_method = args.ratio_baseline or default_ratio_baseline(rows)
  baseline = build_baseline_map(all_rows, baseline_method) if baseline_method else {}
  label_map = method_label_map(all_rows)
  show_ratio = args.ratio and bool(baseline)
  x_keys = x_keys_for_args(args)
  nrows = 2 if show_ratio else 1
  ncols = len(x_keys)

  fig, axes = plt.subplots(nrows, ncols, figsize=args.figsize, squeeze=False,
                           sharey="row", gridspec_kw={"height_ratios": [3, 2]}
                           if show_ratio else None)

  title = args.title or "Binned geometric mean runtime with 95% confidence intervals"
  fig.suptitle(title)

  def runtime_value(row: dict[str, str]) -> float | None:
    return as_float(row, "seconds_per_call")

  def ratio_value(row: dict[str, str]) -> float | None:
    if row.get("method") == baseline_method and not args.show_ratio_baseline:
      return None
    value = as_float(row, "seconds_per_call")
    denominator = baseline.get(row["pair_id"])
    if value is None or denominator is None or denominator <= 0:
      return None
    return value / denominator

  for col, x_key in enumerate(x_keys):
    runtime_ax = axes[0][col]
    line_methods_with_ci(
      ax=runtime_ax,
      rows=rows,
      methods=methods,
      highlighted=highlighted,
      colors=colors,
      x_key=x_key,
      y_value=runtime_value,
      bin_width=args.bin_width,
    )
    runtime_ax.set_title(axis_label(x_key))
    runtime_ax.set_yscale("log")
    runtime_ax.grid(True, which="both", alpha=0.25)
    if col == 0:
      runtime_ax.set_ylabel("geometric mean\nseconds per call")

    if show_ratio:
      ratio_ax = axes[1][col]
      line_methods_with_ci(
        ax=ratio_ax,
        rows=rows,
        methods=methods,
        highlighted=highlighted,
        colors=colors,
        x_key=x_key,
        y_value=ratio_value,
        bin_width=args.bin_width,
      )
      ratio_label = label_map.get(baseline_method, baseline_method)
      ratio_ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--")
      ratio_ax.set_yscale("log")
      ratio_ax.grid(True, which="both", alpha=0.25)
      if col == 0:
        ratio_ax.set_ylabel(f"geometric mean runtime /\n{ratio_label}")

  for col, x_key in enumerate(x_keys):
    axes[-1][col].set_xlabel(axis_label(x_key))

  handles, labels = axes[0][0].get_legend_handles_labels()
  labels = [textwrap.fill(label, width=34) for label in labels]
  fig.legend(handles, labels, title="Method", loc="center left",
             bbox_to_anchor=(0.80, 0.5), frameon=True, fontsize="small",
             title_fontsize="small", borderaxespad=0.0, labelspacing=0.85)
  fig.tight_layout(rect=(0.0, 0.0, 0.78, 0.95))

  args.output.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(args.output, dpi=args.dpi)
  print(f"Wrote {args.output}")


def parse_args(argv: list[str]) -> argparse.Namespace:
  parser = argparse.ArgumentParser(
    description="Plot the CSV produced by mcs_substruct_benchmark.py.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument("--input", type=Path, required=True, help="Benchmark CSV.")
  parser.add_argument("--output", type=Path, default=Path("mcs_substruct_benchmark.png"),
                      help="Output image path.")
  parser.add_argument("--family", choices=["all", "mcs", "substruct"], default="all",
                      help="Benchmark family to include.")
  parser.add_argument("--method", action="append", dest="methods",
                      help="Method to include. Can be supplied multiple times.")
  parser.add_argument("--x",
                      choices=[
                        "size-suite",
                        "query_atoms",
                        "target_atoms",
                        "total_atoms",
                        "smaller_atoms",
                        "larger_atoms",
                        "combined_atoms",
                        "atom_delta",
                        "mcs-size-suite",
                        "mcs_atoms",
                        "mcs_bonds",
                        "reference_mcs_atoms",
                        "reference_mcs_bonds",
                      ],
                      default="size-suite", help="X-axis size metric.")
  parser.add_argument("--mcs-size-method", default="mcs_elements_order",
                      help="MCS method used to derive reference_mcs_atoms/bonds x axes.")
  parser.add_argument("--runtime-scales", choices=["linear", "log", "both"], default="both",
                      help="Runtime y-axis scale rows.")
  parser.add_argument("--ratio", action=argparse.BooleanOptionalAction, default=True,
                      help="Add a paired ratio panel.")
  parser.add_argument("--ratio-baseline", default=None,
                      help="Method to use as the ratio denominator. Defaults to substruct_smiles.")
  parser.add_argument("--ratio-scales", choices=["linear", "log", "both"], default="both",
                      help="Ratio y-axis scale rows.")
  parser.add_argument("--show-ratio-baseline", action="store_true",
                      help="Include the ratio denominator method in ratio panels.")
  parser.add_argument("--breakdown", choices=list(BREAKDOWNS), action="append",
                      help="Write one named breakdown plot. Can be supplied multiple times.")
  parser.add_argument("--all-breakdowns", action="store_true",
                      help="Write all named breakdown plots.")
  parser.add_argument("--output-dir", type=Path, default=None,
                      help="Directory for breakdown outputs. Defaults to the --output directory.")
  parser.add_argument("--average", action="store_true",
                      help="Plot binned geometric means with shaded 95% confidence intervals.")
  parser.add_argument("--baseline-comparison", action="store_true",
                      help="Plot only the combined-size baseline runtime and ratio panels.")
  parser.add_argument("--bin-width", type=float, default=10.0,
                      help="Heavy-atom bin width for --average plots.")
  parser.add_argument("--title", default=None, help="Plot title.")
  parser.add_argument("--figsize", type=float, nargs=2, default=(16.0, 12.0),
                      metavar=("WIDTH", "HEIGHT"), help="Figure size in inches.")
  parser.add_argument("--dpi", type=int, default=160, help="Output DPI.")
  return parser.parse_args(argv)


def main(argv: list[str]) -> int:
  args = parse_args(argv)
  rows = read_rows(args.input)
  ensure_atom_columns(rows)
  ensure_reference_mcs_columns(rows, args.mcs_size_method)
  all_rows = filter_rows(rows, "all", None)
  if not all_rows:
    print("No plottable benchmark rows found.", file=sys.stderr)
    return 2

  if args.baseline_comparison:
    plot_rows = filter_rows(all_rows, "all", set(BASELINE_COMPARISON_METHODS))
    make_baseline_comparison_plot(args, plot_rows, all_rows)
    return 0

  breakdowns = []
  if args.all_breakdowns:
    breakdowns.extend(BREAKDOWNS)
  if args.breakdown:
    breakdowns.extend(args.breakdown)

  if breakdowns:
    seen_breakdowns: set[str] = set()
    for breakdown in breakdowns:
      if breakdown in seen_breakdowns:
        continue
      seen_breakdowns.add(breakdown)
      breakdown_title, method_order = BREAKDOWNS[breakdown]
      plot_rows = filter_rows(all_rows, "all", set(method_order))
      if not plot_rows:
        print(f"No plottable rows found for breakdown: {breakdown}", file=sys.stderr)
        continue
      plot_args = Namespace(**vars(args))
      plot_args.output = output_for_breakdown(args, breakdown)
      if args.title is None:
        plot_args.title = breakdown_title
      if args.average:
        make_average_plot(plot_args, plot_rows, all_rows, method_order=method_order)
      else:
        make_plot(plot_args, plot_rows, all_rows, method_order=method_order)
    return 0

  methods = set(args.methods) if args.methods else None
  plot_rows = filter_rows(all_rows, args.family, methods)
  if not plot_rows:
    print("No plottable benchmark rows found.", file=sys.stderr)
    return 2

  if args.average:
    make_average_plot(args, plot_rows, all_rows)
  else:
    make_plot(args, plot_rows, all_rows)
  return 0


if __name__ == "__main__":
  raise SystemExit(main(sys.argv[1:]))
