#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
M=.lake/packages/mathlib/Mathlib
{
  git -C .lake/packages/mathlib rev-parse HEAD
  rg -n 'def eigenvalues|eigenvalues₀_antitone|spectral_theorem|mulVec_eigenvectorBasis|trace_eq_sum_eigenvalues|eigenvalues_eq_eigenvalues_iff' "$M/Analysis/Matrix/Spectrum.lean"
  rg -n 'hasEigenvalue_iInf|hasEigenvalue_iSup|rayleighQuotient_le_norm|norm_eq_iSup_rayleighQuotient|hasEigenvector_of_isMinOn' "$M/Analysis/InnerProductSpace/Rayleigh.lean"
  rg -n 'posSemidef_iff_eigenvalues_nonneg|eigenvalues_nonneg' "$M/Analysis/Matrix/PosDef.lean"
  rg -n 'posSemidef_iff_dotProduct_mulVec|dotProduct_mulVec_nonneg|trace_nonneg' "$M/LinearAlgebra/Matrix/PosDef.lean"
  rg -n 'map_sum_le' "$M/Analysis/Convex/Jensen.lean" | head -8
  rg -n 'convexOn_exp' "$M/Analysis/Convex/SpecificFunctions/Basic.lean"
  rg -n 'tracePositiveLinearMap|lemma le_iff' "$M/Analysis/Matrix/Order.lean"
} > benchmarks/mathlib-audit.txt
# Keep the negative-search evidence separate; textual search is not a proof of absence.
rg -n -i 'courant|fischer|interlacing|peierls|weyl.{0,50}inequal|min.max.{0,50}eigen' \
  "$M" --glob '*.lean' --glob '!**/RootSystem/**' --glob '!**/Lie/**' \
  --glob '!**/Cartan/**' > benchmarks/mathlib-negative-search.txt || test "$?" -eq 1
