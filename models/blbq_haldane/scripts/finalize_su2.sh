#!/usr/bin/env bash
# After run_su2_grid.sh: collect the L14 sectors, gzip the spectra, rebuild scaling.json and the
# RP/environment records, then run the tests.
set -euo pipefail
cd "$(dirname "$0")/../../.."
runner=models/blbq_haldane/scripts/run.sh
$runner models/blbq_haldane/scripts/collect_su2.jl 14
$runner models/blbq_haldane/scripts/archive_spectra.jl
$runner models/blbq_haldane/scripts/analyze.jl
$runner models/blbq_haldane/scripts/audit.jl
$runner models/blbq_haldane/test/runtests.jl
