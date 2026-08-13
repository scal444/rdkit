#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "usage: $0 {asan|tsan|ubsan} [configure|build|test|all] [ctest arguments...]" >&2
  exit 2
}

sanitizer=${1:-}
action=${2:-all}
if [[ $# -ge 2 ]]; then
  shift 2
else
  shift "$#"
fi

case "$sanitizer" in
  asan|tsan|ubsan) ;;
  *) usage ;;
esac

case "$action" in
  configure|build|test|all) ;;
  *) usage ;;
esac

if [[ -z ${CONDA_PREFIX:-} || ${CONDA_DEFAULT_ENV:-} != rdkit-sanitizers ]]; then
  echo "activate the rdkit-sanitizers Conda environment first" >&2
  exit 1
fi

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

configure() {
  (cd "$source_dir" && cmake --preset "$sanitizer")
}

build() {
  (cd "$source_dir" && cmake --build --preset "$sanitizer" --target install)
}

test_suite() {
  export RDBASE="$source_dir"
  export PYTHONPATH="$source_dir${PYTHONPATH:+:$PYTHONPATH}"
  export LD_LIBRARY_PATH="$source_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  export LLVM_SYMBOLIZER_PATH="${LLVM_SYMBOLIZER_PATH:-$(command -v llvm-symbolizer)}"

  case "$sanitizer" in
    asan)
      export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1}"
      ;;
    tsan)
      export TSAN_OPTIONS="${TSAN_OPTIONS:-halt_on_error=1}"
      ;;
    ubsan)
      export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}"
      ;;
  esac

  (cd "$source_dir" && ctest --preset "$sanitizer" "$@")
}

case "$action" in
  configure) configure ;;
  build) build ;;
  test) test_suite "$@" ;;
  all)
    configure
    build
    test_suite "$@"
    ;;
esac
