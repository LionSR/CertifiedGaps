# Longitudinal-field TFIM

**Certified at all 30 grid points.** For every even L ≥ 6 the ring has a unique ground state
and gap at least the stored bound. The smallest bound is 2749898627347949/36028797018963968
≈ 0.0763, at g = 1, h = 1/20 (n₀ = 6, β₀ = 13/2). Nothing is claimed between grid points, for
odd L, or as h → 0.

| parameters | seed (n₀ sites, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| 30 points g ∈ {1/4, …, 2} × h ∈ {1/20, …, 1} | (6, 2 to 97/4 per point) | 0.0763 (g = 1, h = 1/20) to 3.7050 | every even L ≥ 6 |

The per-point bounds are in the table under "Seeds and certified bounds".

## Model

H_L = −Σ Z_j Z_{j+1} − g Σ X_j − h Σ Z_j on a periodic ring of L spin-1/2 sites (Pauli
matrices, one site per cell), on the grid g ∈ {1/4, 1/2, 3/4, 1, 3/2, 2} ×
h ∈ {1/20, 1/10, 1/4, 1/2, 1}. For h > 0 the model is not integrable and spin-flip symmetry is
broken. At every finite L the ground state is unique by Perron–Frobenius (g > 0 connects the
spin-flip hypercube), but that alone gives no uniform gap.

## Reflection positivity

σᶻ basis: −h_bond = Z ⊗ Z (one term, C = Z, real, coefficient 1), on-site o = −gX − hZ real
symmetric. `rp_check` verifies these identities in `Rational{BigInt}` at every grid point.

## Seeds and certified bounds

n₀ = 6 at every point; β₀ ∈ (1/4)ℤ was chosen by a numerical scan over n₀ ∈ {2, 4, 6},
β₀ ∈ {1/4, …, 64}. The rings of 6 and 12 sites are diagonalised densely in every momentum
sector (no spin-flip or reflection reduction); a polar-factor/Weyl argument with explicit
rounding allowances ([METHOD.md](METHOD.md)) encloses every eigenvalue with multiplicity;
interval partition sums give exact rational εₛ ≥ 1 − Pₛ(6, β₀) and εₜ ≥ 1 − P_τ(12, β₀);
`CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)` gives the bound, stored as an exact binary
rational. All exact values are in [`results/seeds.json`](results/seeds.json) (the Lean file
`lean/CertifiedGaps/Models/IsingLongitudinal.lean` quotes them).

Certified bound ≈ (β₀), row g, column h:

| g \ h | 1/20 | 1/10 | 1/4 | 1/2 | 1 |
|---:|---:|---:|---:|---:|---:|
| 1/4 | 0.3923 (97/4) | 0.7332 (55/4) | 1.6606 (13/2) | 1.9757 (11/2) | 2.0336 (21/4) |
| 1/2 | 0.4564 (13) | 0.9592 (8) | 2.3616 (9/2) | 2.7754 (4) | 2.7276 (4) |
| 3/4 | 0.3170 (8) | 0.7712 (21/4) | 2.1018 (7/2) | 3.0117 (7/2) | 3.1374 (7/2) |
| 1 | 0.0763 (13/2) | 0.5043 (17/4) | 1.6092 (3) | 2.7523 (3) | 3.6403 (3) |
| 3/2 | 0.1692 (7/2) | 0.3601 (13/4) | 1.0556 (11/4) | 2.1293 (5/2) | 3.7050 (5/2) |
| 2 | 0.9121 (9/4) | 0.9643 (9/4) | 1.2728 (9/4) | 2.0005 (9/4) | 3.4912 (2) |

The bounds are 9–74% of the 12-site ED gap (numerical; stored as `spectra` in each
certificate). The equal-defect threshold 0.12731 is sufficient, not necessary: at g = 1,
h = 1/20 the seed has εₛ ≈ 0.184 and the asymmetric recurrence still succeeds.

Dense ED stops well short of L = 20: one generic momentum block at L = 20 has 52,377 states
(≈ 41 GiB as a complex dense matrix; exact counts in `results/audit.json`). Smaller seeds
sufficed at every grid point.

## Reproduce

From this directory (`models/tfim_longitudinal`), Julia 1.12, `export OPENBLAS_NUM_THREADS=1`;
from the repository root use `--project=models/tfim_longitudinal` and prefix the paths. Times
measured on an Apple M-series laptop.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl                       # fast tests, ~15 s
julia --project=. check.jl                               # fast checker, ~15 s
julia --project=. check.jl --full                        # full replay, ~6 min
julia --project=. scripts/final_audit.jl                 # ~30 s
julia --project=. scripts/independent_residual_audit.jl  # ~10 s
julia --project=. scripts/audit.jl                       # ~20 s
julia --project=. certify.jl                             # regenerate all results, ~5 min
```

`certify.jl` takes optional point indices 1–30 (g-major) and `NMAX` (default 6); after a full
run it rewrites `results/seeds.json`. The audit scripts rewrite their `results/*audit*.json`.

## What is rigorous

- **Rigorous**: the RP identity (exact rationals); sector dimensions; eigenvalue enclosures for
  every eigenvalue with multiplicity; partition-function and seed enclosures; the uniform gap
  for every even L ≥ 6 at each grid point. The fast checker recomputes the partition intervals,
  rational seeds and bounds from the stored eigenvalue intervals; `--full` rebuilds every
  block, re-proves its enclosure and checks (in exact rationals) that the fresh interval of
  radius η lies inside the stored one of radius 2η. Accept regenerated certificates only after
  the deep replay.
- **Trusted assumptions**: the IEEE binary64 / MPFR / IntervalArithmetic model with gradual
  underflow (checked at load time) and single-threaded OpenBLAS without fast-math; the
  rounding allowances of [METHOD.md](METHOD.md) bound its dot-product errors but BLAS itself is
  not formally verified.
- **Numerical only**: seed ranking, 12-site ED reference gaps and bound/ED ratios, the
  finite-size diagnostics stored in each certificate, and the free-fermion cross-check.

Independent checks: an independent review (2026-10-11) accepted the enclosure theorem and its
implementation (all momentum sectors kept, no spin-flip reduction, correct Fourier
coefficients −g√(p_b/p_a) e^{2πikt/L}, rounding allowances valid). `scripts/final_audit.jl`
re-derives, without the production code, the sector dimensions (Ramanujan sums), Pauli trace
moments tr H and tr H², 256-bit interval partition ratios, the seeds and bounds of all 30
points (1592 checks). `scripts/independent_residual_audit.jl` checks the residual and
orthogonality allowances on all 180 L = 6 blocks in exact rational arithmetic (360
inequalities; largest exact/allowed ratios 0.005 and 0.017). `scripts/audit.jl` compares with
computational-basis ED and the h = 0 free-fermion log Z (max error 7e−13).

## Layout

| path | role |
|---|---|
| `certify.jl` | numerical seed scan, rigorous enclosure of the best seed, writes results |
| `check.jl` | checker entry point (`--full` for the full enclosure replay) |
| `src/LongitudinalCertificate.jl` | momentum blocks, enclosures, seed enclosure, I/O |
| `src/checks.jl` | `check_result`, `check_grid`: the checks behind `check.jl` and the tests |
| `scripts/final_audit.jl` | independent stored-artifact audit → `results/final_audit.json` |
| `scripts/independent_residual_audit.jl` | exact L = 6 residual audit → `results/independent_residual_audit.json` |
| `scripts/audit.jl` | numerical cross-checks, resource counts → `results/audit.json` |
| `test/runtests.jl` | fast tests (RP, sector counts, all 30 certificates via `src/checks.jl`, seed table) |
| `METHOD.md` | the enclosure theorem and the rounding allowances |
| `IDEAS.md` | candidate routes and the choice |
| `results/certificate_g<g>_h<h>.json.gz` | per-point certificate: exact seeds, bound, all enclosed eigenvalues |
| `results/seeds.json` | g, h, n₀, β₀, εₛ, εₜ, bound for all 30 points (exact) |

The per-point certificates are gzip-compressed JSON (≈ 260 kB each uncompressed, mostly the
4096 + 64 eigenvalue midpoints that the checker needs);
`gzip -dc results/certificate_g1_1_h1_20.json.gz` shows one. Exact numbers are
`"num//den"` strings.

## History

Moved from `certificates/tfim_longitudinal/`. The per-point `metadata.source_sha256` values
are historical: they fingerprint the generating sources (under the former core name
`ThermalMPOGap` and the former file layout) and are not checked; the deep replay re-derives
every enclosure with the current code instead. In the tidy-up the module `certify.jl` became
`src/LongitudinalCertificate.jl`, `scan.jl` became `certify.jl`, the audits moved to
`scripts/`, `seeds.jl` was folded into the module, and `REPORT.md`, `INDEPENDENT_REVIEW.md`
and the report generator `report.jl` were folded into this README. Earlier removals
(uncompressed JSONs, migration and validation records, run logs) and the removed files remain
in the git history. The three `results/*audit*.json` files are kept as produced before the
tidy-up (their source hashes name the old layout); all three audits were rerun afterwards and
passed with the same counts. In the format tidy-up the per-point files were renamed from
`results/g<g>_h<h>.json.gz` to `results/certificate_g<g>_h<h>.json.gz` (contents unchanged, so
the `certificate_sha256` values in `results/final_audit.json` still apply; the `file` fields of
`results/seeds.json` were updated), and `check.jl` gained `--full` as the name of the deep replay.
Later the checks moved from `check.jl` into `src/checks.jl`, so that `check.jl` is the entry
point and `test/runtests.jl` no longer includes it; the `DEEP_CHECK=1` alias of `--full` was
dropped. Behaviour and checks are unchanged.
