#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export OPENBLAS_NUM_THREADS=1
if [[ -x "$PWD/certificates/blbq_haldane/runtime/julia-1.12.4/bin/julia" ]]; then
    export JULIA_DEPOT_PATH="$PWD/certificates/blbq_haldane/runtime/depot"
    export PATH="$PWD/certificates/blbq_haldane/runtime/julia-1.12.4/bin:$PATH"
fi
exec julia --project=julia "$@"
