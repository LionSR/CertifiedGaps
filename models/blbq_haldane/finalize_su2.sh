#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
runner=models/blbq_haldane/run.sh
$runner models/blbq_haldane/collect_su2.jl 14
$runner models/blbq_haldane/archive_spectra.jl
$runner models/blbq_haldane/analyze.jl
$runner models/blbq_haldane/validate.jl
$runner models/blbq_haldane/report.jl
$runner models/blbq_haldane/audit.jl
