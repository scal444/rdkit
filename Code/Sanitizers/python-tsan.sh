#!/usr/bin/env bash
set -euo pipefail

runtime=$(clang --print-file-name=libclang_rt.tsan-x86_64.so)
export LD_PRELOAD="$runtime${LD_PRELOAD:+:$LD_PRELOAD}"
exec "$CONDA_PREFIX/bin/python" "$@"
