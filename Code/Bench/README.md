# RDKit Benchmarks

To run:

```bash
mkdir build
cd build
cmake ..
cmake --build . --target bench -j "$(nproc)"
# see `./Code/Bench/bench --help` for options
export RDBASE=".."
./Code/Bench/bench
```

## MCS vs substructure scan

The Python scripts `mcs_substruct_benchmark.py` and
`plot_mcs_substruct_benchmark.py` compare FMCS variants against small-into-large
substructure search on the same query/target pairs:

```bash
python Code/Bench/mcs_substruct_benchmark.py \
  --input ~/data/chembl.smi \
  --pair-mode size-combinations \
  --size-min-atoms 0 \
  --size-max-atoms 100 \
  --examples-per-size 1 \
  --max-pairs 200 \
  --mcs-timeout 10 \
  --output mcs_substruct_benchmark.csv
python Code/Bench/plot_mcs_substruct_benchmark.py \
  --input mcs_substruct_benchmark.csv \
  --output mcs_substruct_benchmark.png
```

`size-combinations` streams the input file, keeps representative molecules for
each heavy-atom count in the requested range, and samples small/large molecule
combinations from those representatives. The benchmark CSV is flushed as rows
are produced, so interrupted runs leave partial results that the plotting script
can use.

Use `--append` with a different `--seed` to add another random sample to an
existing CSV without overwriting it; pair IDs continue from the existing file.

The plotting script defaults to a three-panel size view with x axes for the
smaller molecule heavy-atom count, larger molecule heavy-atom count, and
combined heavy-atom count. Use `--x` to select a single size metric instead.
Use `--all-breakdowns --output-dir plots` to write separate baseline, MCS core,
MCS ring, and substructure comparison figures. Breakdown figures include both
log and linear runtime rows, and the ratio rows omit the substructure baseline
method unless `--show-ratio-baseline` is supplied.
Use `--runtime-scales log --ratio-scales log` for log-only breakdowns. Add
`--average --bin-width 10` to plot binned geometric means with shaded 95%
confidence intervals.
