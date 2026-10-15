#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export ELAN_HOME="$PWD/.tools/elan"
export JULIA_DEPOT_PATH="$PWD/.tools/julia_depot"
export PATH="$PWD/.tools/elan/bin:$PATH"
JULIA_BIN="${JULIA_BIN:-$PWD/.tools/julia-1.12.0/bin/julia}"
# Run one compiler at a time; concurrent large bignum processes invalidate timings.
for name in Baseline Bench6 Bench8 Bench10 Bench12 H6 H8 H10 Matvec12 Dot12 Prepared12 PreparedH8 Gram6 Tests PlainDot12 FastDot12 Fast12 FastReverse8 PreparedMatvec12; do
  "$JULIA_BIN" --startup-file=no scripts/bench_one.jl "$name" 180
done
