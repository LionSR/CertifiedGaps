# Spin-one comparison anchors for the (36, 49/4) seed

**No gap certificate.** This is a feasibility study for one ingredient of the spin-one
Heisenberg seed at (n₀, β₀) = (36, 49/4): exact *counting anchors* on short twisted rings
(n = 8, 10, 12), obtained by comparing the ring with two open halves joined through an exact
local PSD remainder. The exact part (local PSD pivots, half-chain spectral enclosures,
comparison counts, binned-tail arithmetic) is checked in Rational/BigInt arithmetic by
`check.jl`. Ring spectra, parameter searches and the continuation are NON-CERTIFYING Float64.

| n | count-opt (u, v, c) | Sz=0 count: bare / ground-opt / count-opt | binned-tail Sz=0 (total) | certified gap |
|---:|---|---|---|---|
| 8 | (−3/10, −1/10, −5076841/2500000) | 44 / 63 / 34 | 28 (98) | none claimed |
| 10 | (−1/5, 0, −5538721/2500000) | 86 / 115 / 62 | 38 (140) | none claimed |
| 12 | (−1/5, 1/20, −48533/20000) | 191 / 212 / 121 | 49 (183) | none claimed |

## Model

Spin-one Heisenberg ring `H = Σ S_i · S_{i+1}` in the standard Sᶻ basis with the integer bond
`h = SᶻSᶻ + (S⁺S⁻ + S⁻S⁺)/2` (eigenvalues −2, −1, 1); lengths count sites. Both the untwisted
(g = I) and π-twisted (g = P) rings are treated; the π twist negates the two exchange entries
on the closing bond. Thermal conventions (β = 49/4, per-site shift a = 700741/500000) follow
`models/spin1/src/thermal.jl`.

The comparison: for a half of length l = n/2 and `B = u h + v h²`,
`A_l = H_open,l + B_12 + B_(l−1,l)`, `H0 = A_l ⊗ I + I ⊗ A_l + 2c`, and the four-site remainder
`R = h_23 − B_12 − B_34 − c`. If R ≥ 0 then `H − H0` is a sum of two embedded copies of R
(one per cut, the twisted one conjugated by a diagonal π rotation), so `H ≥ H0` and min–max
bounds the number of ring levels below any ρ by the (tensor-sum) count of H0. See
[METHOD.md](METHOD.md) for the derivation.

## Reflection positivity

Not reproved here. The seed proof uses the usual sublattice-rotated / Cartesian RP frame of
`models/spin1/`; this directory only supplies spectral anchors for the seed computation.

## Result

Findings (Float64 unless marked exact):

1. **Local comparison is exact and cheap.** For every parameter set (ground- and
   count-optimised, n = 8, 10, 12), R is PSD by exact rational LDL in all Sᶻ blocks
   (dimensions 19, 16, 10, 4, 1), with the π-twist and spin-flip identities checked exactly.
   Half-chain Sᶻ blocks have dimension ≤ 141 (l = 6); every eigenvalue is enclosed exactly
   (largest radius ≈ 2.8e−7, stored outward on a 1e−6 grid); no recursion is needed.
2. **Ground optimisation is the wrong objective.** At n = 12 it raises the H0 ground bound from
   −18.740550 (bare) to −17.527216 but increases the all-Sᶻ count below the uniform tail cutoff
   from 773 to 914; the count-optimised comparison (ground bound −17.690178) gives 475.
   Sᶻ=0 counts per n are in the table above (bare / ground-opt / count-opt). The physical
   ring has only 2, 5, 8 Sᶻ=0 levels below the cutoff (n = 8, 10, 12, untwisted).
3. **Binned tail: 121 → 49.** The uniform cutoff ρ ≈ −9.644247, −12.267846, −14.891444 spends
   the whole 3e−5 trace budget on the tail. Charging comparison levels in 32 bins of width 1/8
   (plus an overflow bin) with half the budget (1.5e−5) gives targets ρ = −9.925, −12.687,
   −15.453 and initial counts (Sᶻ = 0,1,2,3,4) of 28,22,11,2,0 / 38,32,15,4,0 / 49,40,20,6,1,
   i.e. 98 / 140 / 183 in total. Tail totals 1.48995e−5, 1.49126e−5, 1.49384e−5 (exact
   rational, Taylor-degree-256 exponential upper bounds).
4. **A one-shot count does not give ordered lower bounds.** `H = diag(−2,−1,10)` with Ritz
   vectors e₂, e₃ has zero residual and count 2 below ρ = 0, yet the lowest Ritz value −1 is
   not a lower bound on −2 (tested). Continuation (or a certified complement bound) is needed.
5. **Continuation works numerically.** Along `H(t) = H0 + t(H − H0)` with inverse-Ritz
   (Lehmann–Goerisch) lower bounds and a moving counting barrier, all 26 initially nonzero
   (n, g, Sᶻ ≥ 0) sectors completed. Stages and initial → final retained K, Sᶻ = 0:
   n=8 I 28 (28→2), P 25 (28→3); n=10 I 51 (38→2), P 42 (38→3); n=12 I 73 (49→5),
   P 72 (49→5). All sectors: 658 stages, about 40 min of Julia time (34 min at n = 12; the
   n=12 Sᶻ=0 runs took ≈ 670 s each). Uniform-tail Sᶻ=0 pilots: n=8 75 stages (34→2), n=10 73
   stages (62→7); the n=12 uniform pilot (121 vectors, ≈ 63 s per eigensolve) was stopped
   after 2 stages. Per-sector numbers: `results/homotopy_summary.csv` and METHOD.md.
6. **Other routes rejected in cheap form.** Three or more blocks increase the counts (best
   sampled Sᶻ=0: 104, 207, 407 for n = 8, 10, 12); overlapping Anderson clusters give only a
   scalar bound (no excited-level count); an edge spin-1/2 compression gives no lower comparison.
7. **Revised Lean estimate.** The earlier 4.5–7.5k LoC / 1.5–3 weeks for the anchor route is not
   supported. Planning allowance for this anchor/continuation component: 3.5–6k Lean LoC
   (sector/local-PSD lifting 0.6–1k, exact spectral enclosures 0.8–1.4k, tensor counts and
   bins 0.4–0.7k, inverse-Ritz/counting 1.2–2k, continuation bookkeeping 0.5–0.9k); about
   8–13k LoC and 4–8 weeks for the full seed route. Kernel runtime: the half-chain enclosure
   products alone are ≈ 11.9M multiply-adds (≈ 7 min at 3e4 iterations/s naively); the naive
   three-Gram continuation is 115k–140k packed dot products at n = 12, i.e. 80–97 hours at a
   hypothetical 2.5 s per product, so the "15–60 min" kernel target must become an explicit
   benchmark milestone with batching (these are estimates, not measurements).

## Reproduce

Julia 1.12, `export OPENBLAS_NUM_THREADS=1`, from this directory (times on an Apple-silicon
laptop):

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl   # 45 tests, about 10 s
julia --project=. certify.jl         # writes results/certificate.json, about 6 s
julia --project=. check.jl           # exact replay of the certificate, about 12 s
```

`check.jl` has no slow mode: it re-derives every exact record. The NON-CERTIFYING data are
regenerable (and gitignored) with:

```sh
julia --project=. scripts/scan.jl          # results/optimized.json (ground optimisation)
julia --project=. scripts/count_scan.jl    # results/count_optimized.json, scan_n*.json
julia --project=. scripts/spectra.jl 8 10 12   # ring low spectra n*_g[IP].json (minutes)
julia --project=. scripts/analyze.jl       # analysis.json, results/sectors.csv
julia --project=. scripts/tail_bins.jl     # results/binned_tails.json
export ANCHOR_FRESH=1
for n in 8 10 12; do for g in I P; do
  julia --project=. scripts/homotopy.jl $n $g binned
done; done
julia --project=. scripts/homotopy_sectors.jl      # n = 8, 10, Sz = 1..3
julia --project=. scripts/homotopy_sectors.jl I    # n = 12, Sz = 1..4
julia --project=. scripts/homotopy_sectors.jl P
julia --project=. scripts/report_tables.jl # METHOD.md tables, results/homotopy_summary.csv
julia --project=. scripts/report_metrics.jl  # plain-text summary to stdout
```

The continuation takes about 40 min in total. Without `ANCHOR_FRESH=1` a completed run is
reused and an incomplete one resumes. The recorded runs include retries and resumes.

## What is rigorous

- **Exact (Julia Rational{BigInt} / bounded Int128):** local PSD of R by LDL pivots, the twist
  and spin-flip identities, half-chain spectral enclosures (rounded Float64 eigenvectors Q to
  denominator 2³², then `δ = ‖QQᵀ − I‖∞ < 1/2`, `e = ‖A − QDQᵀ‖∞`, radius `e + 3‖D‖δ`, by polar
  decomposition and Weyl), the per-sector comparison counts and anchors, and the binned-tail
  sum. The eigenvector proposals are regenerated by LAPACK, not stored, so the JSON is not a
  standalone (Lean-checkable) certificate. The Taylor inequality `exp(−x) ≤ 1/Σ_{j≤256} x^j/j!`
  and the min–max/tensor-sum arguments are proved on paper only.
- **Floating point only:** ring spectra, physical counts and losses (`results/sectors.csv`),
  the boundary-parameter searches (no global optimality), multiblock comparisons, and the whole
  continuation (Float64 generalized eigenproblems with numerical safety padding).
- **Not claimed:** any certified ring eigenvalue, thermal trace, seed defect or gap; nothing
  under `lean/` or `models/spin1/results/` is touched.

## Layout

| path | content |
|---|---|
| `certify.jl`, `check.jl` | exact certificate writer and checker |
| `src/common.jl` | `AnchorStudy`: Sᶻ bases, integer bonds, comparison and ring matrices |
| `src/exact.jl` | exact LDL, half-chain spectral enclosure, rational `exp(−x)` bound |
| `scripts/` | parameter scans, ring spectra, analysis, binned tail, continuation, report tables |
| `test/runtests.jl` | 45 fast tests (Hamiltonian identities, checker rejection, decomposition, counterexample) |
| `results/certificate.json` | exact records read by `check.jl` |
| `results/optimized.json`, `count_optimized.json` | Float64 parameter proposals read by `certify.jl` |
| `results/binned_tails.json` | Float64 binned targets read by the continuation |
| `results/sectors.csv`, `homotopy_summary.csv` | compact per-sector and per-continuation summaries |
| `METHOD.md` | derivation, counterexample, continuation design, cost model, generated tables |

## History

PR #27 (anchor study). Sibling studies of the same seed live in the other subdirectories of
`models/spin1_lean_seed/`. The regenerable Float64 dumps (ring spectra, raw scan rows,
per-stage continuation records) were removed from git in the cleanup; they are in the PR's
first commit.
