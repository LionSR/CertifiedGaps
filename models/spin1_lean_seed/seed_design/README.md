# Spin-1 Heisenberg: smaller finite seed inputs

A design study for the Lean spin-1 gap proof: which finite thermal traces does the seed really
need? Three seed designs are replayed in exact rational arithmetic down to
`ChainCertificate.Valid` (`lean/CertifiedGaps/ChainCheck.lean`). The recommended design A needs
nine twisted traces at ring lengths 6, 8, 10, eight times less Hilbert-space volume than the
original eighteen traces at n ≤ 12. **Every bound here is conditional**: it assumes the finite
trace tables to the stated absolute error and the analytic sector, filter and lifting lemmas
(METHOD.md), none of which is checked in Lean here.

| design | thermal inputs (lengths; twists; abs. error) | seed (n₀, β₀) | conditional gap > | lengths covered |
|---|---|---|---:|---|
| A (recommended) | 6, 8, 10; 1, P, C; 3/10000 | (60, 105/4) | 0.005966152 | every even L ≥ 60 |
| B | 6, 8, 10; 1, P, C; 4/10000 | (72, 525/16) | 0.006309519 | every even L ≥ 72 |
| C | 4–12 even; 1, P; 1/1000 | (60, 49/2) | 0.039036761 | every even L ≥ 60 |

## Model

Periodic spin-1 Heisenberg ring `H = Σ_j S_j·S_{j+1}`, exchange 1, shifted by
`a = 700741/500000` per bond as in `models/spin1`. Lengths count sites. "Twists" are the
seam rotations of `models/spin1`: identity, P (π) and C (2π/3); from them
`I = (Z + 3P + 8C)/12`, `O = (Z + 3P − 4C)/12`, `N = (Z − P)/4` (or `J = (Z + 3P)/4`, `N`
without C), with multiplicities 1, 2, 3.

## Reflection positivity

The Cartesian basis of `models/spin1`: `−h = Σ_α S^α ⊗ conj(S^α)` up to the constant shift,
with the moment inequalities of `lean/CertifiedGaps/Heisenberg/Main.lean`. Nothing new is proved here.

## Seeds and certified bounds

`results/certificate.json` stores, for each design, the rational filters `x^s q(x)^2` with their
exterior-ray caps, the sector masses, the seed witnesses (D, E, Q, the dominant-weight
refinements) and a 9-step rational chain on the grid `10^-512`:

| design | b → β₀ = r b | k; trial length | εₛ (p) | εₜ (q) | chain rate R | gap = log R / β₀ |
|---|---|---|---:|---:|---|---:|
| A | 21/2 → 105/4 | 30; 120 | 0.124048525 | 0.125580476 | 116954119/10⁸ | > 0.005966152 |
| B | 21/2 → 525/16 | 30; 144 | 0.117346685 | 0.128925715 | 123002087/10⁸ | > 0.006309519 |
| C | 49/4 → 49/2 | 30; 120 | 0.047914768 | 0.112578548 | 130115651/(5·10⁷) | > 0.039036761 |

(p, q rounded up to nine decimals; exact values are in the certificate.) B's thermal defect
exceeds the diagonal sufficient threshold 0.12731, but its unequal-defect chain is valid. All
three exceed the original target `log 20 / 784`, checked as `R^u > 20^v` for `784/β₀ = u/v`.
`results/compact_witnesses.json` (9 KB) is the minimal replay recipe for the three designs.

**Numerical Pareto front** (`results/pareto.csv`, Float64, non-certifying): 39 nondominated
rows from 94 grid representatives over (state proxy `Σ 3ⁿK`, error, gap, n₀). Representative
rows: n ≤ 10 with I/O/N at error 3e-4 gives ≈0.00674 (n₀ = 60), at 1e-4 ≈0.0119 (n₀ = 72),
at 3e-5 with odd twists ≈0.0193 (n₀ = 72); n ≤ 12 with J/N at 1e-3 gives ≈0.0470 and at
3e-5 (odd-inclusive) ≈0.0721 (n₀ = 60). **Negative results** (`results/diagnostics.json`):
length 8 fails even at zero error (best I cap 1.00693, best max defect ≈0.5288 over 132,075
lift cases); length 10 without C fails (best J cap 1.00449, best max defect ≈0.2956 over
478,980 cases); length 10 with C at error 1e-3 narrowly fails (≈0.1338 over 1,359,900
cases). These are failures of the tested filter family, not lower bounds.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`
(`D=models/spin1_lean_seed/seed_design`):

```sh
julia --project=$D -e 'using Pkg; Pkg.instantiate()'
julia --project=$D $D/check.jl          # 3 s: exact replay of all three designs
julia --project=$D $D/test/runtests.jl  # 4 s: 46 exact + forged-witness tests (+64 with cache/)
julia --project=$D $D/certify.jl        # 5 s: rewrites the two results/ witnesses
```

`check.jl` has no slow mode: it re-derives the masses from `models/spin1/results/source*.json`
(hash-audited against `verification.json`), re-runs the exact filter cap bisection and the whole
seed computation, compares them with the stored witness, and checks the chain conditions.

Discovery pipeline (Float64, non-certifying; outputs under `results/`, raw rows under the
ignored `results/raw/`, spectra under the ignored `cache/`):

```sh
julia --project=$D $D/scripts/spectra.jl 4 5 6 7 8 9 10   # 70 s; add 11 12 for the baseline
julia --project=$D $D/scripts/search_existing.jl          # results/raw/existing_search.json
julia --project=$D $D/scripts/search_new.jl 8             # results/raw/new_search_8.json
julia --project=$D $D/scripts/search_new.jl 10            # results/raw/new_search_10.json
julia --project=$D $D/scripts/refine.jl                   # results/raw/refinement.json
julia --project=$D $D/scripts/diagnostics.jl              # results/diagnostics.json
julia --project=$D $D/scripts/shift_search.jl             # 8 s; results/shift_search.json
julia --project=$D $D/scripts/fetch_paper.jl              # upstream/ (network)
julia --project=$D $D/scripts/trial_scan.jl               # results/trial_scan.json (network)
julia --project=$D $D/scripts/summarize.jl                # results/summary.json, pareto.csv
```

`summarize.jl` needs the raw searches, the n ≤ 12 spectra and `upstream/`. Only the times
shown were re-measured after the reorganisation; the full search was not rerun.

## What is rigorous

Exact: every inequality in `check.jl` (BigInt rationals; fractional powers by integer-root
enclosures; `log R` by a rational series enclosure). The chain check is a line-by-line
transcription of `ChainCertificate.Valid` (it additionally requires `0 ≤ s_j, t_j` and
`length(t) = length(s)`). Conditional: the trace tables to absolute error 3/10000 (A) are
assumed; here they are the stored `models/spin1` source-2 centres (Julia-verified to 1e-5, not
Lean). The C-symmetry sector equivalences, the two-piece capped-mass bound, the dominant-weight
isolation, temperature lifting and the MPS variational bound are analytic lemmas from the
upstream paper, not formalised. The length-144 trial (design B) is a fresh exact Nemo
contraction in `results/trial_scan.json`, stored, not recomputed by the checker. All spectra,
search rows, `pareto.csv`, `summary.json`, workload and Lean-effort estimates are Float64 or
estimates only.

## Layout

- `certify.jl`, `check.jl`: exact designs and their replay.
- `src/common.jl` (paths, Float64 search, provenance), `src/exact.jl` (exact seed and chain).
- `scripts/`: spectra, searches, diagnostics, shift scan, trial contraction, summary.
- `test/runtests.jl`: exact replay, forged chains and filters, cached-spectrum controls.
- `results/`: `certificate.json`, `compact_witnesses.json`, `trial_scan.json`, `pareto.csv`,
  `summary.json`, `diagnostics.json`, `shift_search.json`. `summary.json`, `trial_scan.json`,
  `diagnostics.json` and `shift_search.json` were compacted (deep-equal) and not regenerated, so
  their `metadata.source_sha256` names the pre-reorganisation file layout.
- `METHOD.md`: the full study (filters, error propagation, search coverage, negative results,
  twists and multiplicities, workload and Lean-effort estimates). `IDEAS.md`: route memo.

## History

PR #24 (`codex/spin1-seed-design`). Reorganised into the standard layout; the raw search rows
(≈2.7 MB) are no longer tracked. Regenerating `certificate.json` after the move changed only
provenance hashes and the trial-file path; every certified number is identical.
