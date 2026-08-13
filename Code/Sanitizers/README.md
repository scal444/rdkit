# Clang sanitizer builds

The sanitizer branches share one pinned conda-forge toolchain. The initial
configuration uses Clang/LLVM 22.1.8, Python 3.13, Ninja, and RDKit's normal C++
and Python test suites.

Create and activate the environment once:

```bash
mamba env create --file .conda/sanitizers.yml
conda activate rdkit-sanitizers
```

Configure, build, and test the sanitizer associated with the current branch:

```bash
Code/Sanitizers/run.sh asan all
Code/Sanitizers/run.sh tsan all
Code/Sanitizers/run.sh ubsan all
```

The second argument may be `configure`, `build`, or `test` when iterating. Any
remaining arguments on a `test` invocation are forwarded to CTest, for example:

```bash
Code/Sanitizers/run.sh asan test -R GraphMol -j 1
```

The runner enables symbolization, stops on the first report, and uses a small
launcher to preload the matching shared sanitizer runtime before Python starts.
ASan also enables the loose cross-translation-unit initialization-order check;
its stricter mode is deferred because shared-library initialization ordering can
produce false positives there. The build instruments pointer comparisons and
subtractions, and the runtime reports invalid pairs even when both pointers are
null.
Override `ASAN_OPTIONS`,
`TSAN_OPTIONS`, or `UBSAN_OPTIONS` in the shell when layering in more checks.
LeakSanitizer is initially disabled because it cannot run under ptrace-based
sandboxes; enable it outside such a sandbox with
`ASAN_OPTIONS=detect_leaks=1:halt_on_error=1`.

The runner also loads `asan.supp`, `tsan.supp`, or `ubsan.supp`. These files are
reserved for confirmed defects in external libraries that RDKit does not
bundle. They intentionally contain no active entries: a report should not be
suppressed merely because an external header or runtime function appears in
its stack when the invalid state is owned by RDKit. Set an explicit
`suppressions=` entry in the corresponding sanitizer options to use a different
file.
