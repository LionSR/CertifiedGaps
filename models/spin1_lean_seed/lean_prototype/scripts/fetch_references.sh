#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p references
for file in initialization transfer thermal72 trials; do
  curl -fsSL "https://raw.githubusercontent.com/openai/math/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/build/$file.tex" -o "references/$file.tex"
done
sha256sum -c benchmarks/paper-source-sha256.txt
