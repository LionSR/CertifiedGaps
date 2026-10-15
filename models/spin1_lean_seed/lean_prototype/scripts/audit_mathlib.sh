#!/usr/bin/env bash
# Source-level Mathlib API audit (textual search, not a Lean proof of absence).
# MATHLIB_DIR defaults to the Mathlib checkout of the repository's lean/ project, which pins
# the same toolchain (v4.35.0-rc3) and Mathlib commit c55e6e78.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(git rev-parse --show-toplevel)"
PKG="${MATHLIB_DIR:-$REPO/lean/.lake/packages/mathlib}"
M="$PKG/Mathlib"
test -d "$M" || { echo "no Mathlib source at $M; set MATHLIB_DIR" >&2; exit 1; }
{
  git -C "$PKG" rev-parse HEAD
  rg -n --no-heading -e 'def eigenvalues|eigenvalues₀_antitone|spectral_theorem' \
    -e 'mulVec_eigenvectorBasis|trace_eq_sum_eigenvalues|eigenvalues_eq_eigenvalues_iff' \
    "$M/Analysis/Matrix/Spectrum.lean"
  rg -n --no-heading -e 'hasEigenvalue_iInf|hasEigenvalue_iSup|rayleighQuotient_le_norm' \
    -e 'norm_eq_iSup_rayleighQuotient|hasEigenvector_of_isMinOn' \
    "$M/Analysis/InnerProductSpace/Rayleigh.lean"
  rg -n --no-heading 'posSemidef_iff_eigenvalues_nonneg|eigenvalues_nonneg' \
    "$M/Analysis/Matrix/PosDef.lean"
  rg -n --no-heading 'posSemidef_iff_dotProduct_mulVec|dotProduct_mulVec_nonneg|trace_nonneg' \
    "$M/LinearAlgebra/Matrix/PosDef.lean"
  rg -n --no-heading 'map_sum_le' "$M/Analysis/Convex/Jensen.lean" | head -8
  rg -n --no-heading 'convexOn_exp' "$M/Analysis/Convex/SpecificFunctions/Basic.lean"
  rg -n --no-heading 'tracePositiveLinearMap|lemma le_iff' "$M/Analysis/Matrix/Order.lean"
} | sed "s|$M/||" > results/mathlib-audit.txt
# Negative search for min-max / Weyl / interlacing; hits are listed for manual triage.
rg -n -i --no-heading 'courant|fischer|interlac|peierls|weyl.{0,50}inequal|min.max.{0,50}eigen' \
  "$M" --glob '*.lean' --glob '!**/RootSystem/**' --glob '!**/Lie/**' \
  --glob '!**/Cartan/**' | sed "s|$M/||" > results/mathlib-negative-search.txt \
  || test "$?" -eq 1
