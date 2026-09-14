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

## Profile-guided optimization

`RDK_PGO_MODE` supports instrumented (`GENERATE`) and optimized (`USE`)
builds with Clang or GCC. Keep the compiler, build options, and source revision
the same between the two builds. Link-time optimization may be enabled in both
builds with `RDK_LTO_MODE`.

For Clang, `RDK_PGO_PROFILE_PATH` names the raw-profile directory in the
instrumented build. Clang writes a separate `default_%m.profraw` file for each
instrumented executable or shared library, so training may exercise multiple
RDKit targets and processes without profiles overwriting one another. Start
with an empty profile directory; do not mix profiles from different builds.

```bash
cmake -S . -B build-pgo-generate \
  -DRDK_PGO_MODE=GENERATE \
  -DRDK_PGO_PROFILE_PATH="$PWD/build-pgo-generate/pgo-raw" \
  -DRDK_LTO_MODE=THIN \
  -DRDK_BUILD_INCHI_SUPPORT=ON
cmake --build build-pgo-generate --parallel 12
```

Run workloads representative of the intended use. The standalone drivers can
contribute parsing and sanitization plus a mixture of common operations to one
multi-target profile. For example:

```bash
smiles_file=/path/to/representative.smi
build-pgo-generate/Code/Bench/smiles_pipeline_bench \
  "$smiles_file" 135000
for operation in canonical_smiles morgan pains_substructure pickle \
    tanimoto_similarity inchi_roundtrip hs_roundtrip; do
  build-pgo-generate/Code/Bench/molecule_workloads_bench \
    "$operation" "$smiles_file" 10000 1
done
build-pgo-generate/Code/Bench/molecule_workloads_bench \
  etkdg "$smiles_file" 100 1
build-pgo-generate/Code/Bench/molecule_workloads_bench \
  mmff "$smiles_file" 100 1
```

The input sizes and operation frequencies determine the profile weighting and
should model the deployment workload. ETKDG will generally use a smaller,
representative molecule sample because conformer generation is substantially
more expensive than the other operations.

Merge all Clang raw profiles, then configure the optimized build with the
resulting indexed profile:

```bash
llvm-profdata merge \
  -output="$PWD/build-pgo-generate/rdkit.profdata" \
  "$PWD/build-pgo-generate/pgo-raw"
cmake -S . -B build-pgo-use \
  -DRDK_PGO_MODE=USE \
  -DRDK_PGO_PROFILE_PATH="$PWD/build-pgo-generate/rdkit.profdata" \
  -DRDK_LTO_MODE=THIN \
  -DRDK_BUILD_INCHI_SUPPORT=ON
cmake --build build-pgo-use --parallel 12
```

With GCC, `RDK_PGO_PROFILE_PATH` is a directory in both phases. Configure the
instrumented build with `GENERATE`, run the training workloads, and configure
the optimized build with `USE` and the same profile directory. GCC consumes
its generated profile files directly, so `llvm-profdata` is not used.
