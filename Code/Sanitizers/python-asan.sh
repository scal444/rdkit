#!/usr/bin/env bash
set -euo pipefail

runtime=$(clang --print-file-name=libclang_rt.asan-x86_64.so)
export LD_PRELOAD="$runtime${LD_PRELOAD:+:$LD_PRELOAD}"
# CMake probes this launcher while configuring. The test runner explicitly
# enables leak detection later; keeping it off for probes avoids LSan attaching
# to Conda's process-capture helper.
export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=0}"
exec "$CONDA_PREFIX/bin/python" "$@"
