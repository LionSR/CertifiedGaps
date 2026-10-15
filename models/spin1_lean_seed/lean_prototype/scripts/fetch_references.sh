#!/usr/bin/env bash
# Download the paper excerpts cited in METHOD.md into the git-ignored references/.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p references
BASE=https://raw.githubusercontent.com/openai/math/adc7f1241b42e322a6451854ab7e4b4c146bf78a
DIR=preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/build
for file in initialization transfer thermal72 trials; do
  curl -fsSL "$BASE/$DIR/$file.tex" -o "references/$file.tex"
done
shasum -a 256 -c results/paper-source-sha256.txt
