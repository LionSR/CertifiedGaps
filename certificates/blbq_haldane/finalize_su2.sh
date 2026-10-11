#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
runner=certificates/blbq_haldane/run.sh
$runner certificates/blbq_haldane/collect_su2.jl 14
$runner certificates/blbq_haldane/archive_spectra.jl
$runner certificates/blbq_haldane/analyze.jl
$runner certificates/blbq_haldane/validate.jl
$runner certificates/blbq_haldane/report.jl
$runner certificates/blbq_haldane/audit.jl
