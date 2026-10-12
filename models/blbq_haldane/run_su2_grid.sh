#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
for spin in $(seq 0 14); do
    for momentum in $(seq 0 7); do
        output="models/blbq_haldane/results/su2_L14_S${spin}_k${momentum}.json"
        [[ -f "$output" ]] && continue
        timeout 2700 models/blbq_haldane/run.sh -t 4 \
            models/blbq_haldane/su2_sector.jl 14 "$spin" "$momentum"
    done
done
