# Dimerised Heisenberg chain

**Certified** at δ = 1/2, 3/10 and 1/4: for every even L ≥ 4 cells (N = 8, 12, 16, … spins,
i.e. N divisible by 4) the ground state is unique and the gap is at least the bound below.
δ = 1/4 is the smallest certifiable point of the grid {1/2, 3/10, 1/4, 1/5, 3/20, 1/10} for
this recurrence with rings of at most 16 spins: at δ ≤ 1/5 the checker proves that **every**
β > 0 fails at both eligible seeds n₀ = 2, 4. The RP hypothesis is the general internally
reflected (conjugated-factor) form. The Lean statements in
`lean/CertifiedGaps/Models/DimerizedHeisenberg*.lean` take earlier, weaker seeds as hypotheses.

| δ | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| 1/2 | (4, 6) cells | 0.534756 | every even L ≥ 4 cells |
| 3/10 | (4, 13/2) cells | 0.173197 | every even L ≥ 4 cells |
| 1/4 | (4, 29/4) cells | 0.058225 | every even L ≥ 4 cells |
| 1/5, 3/20, 1/10 | none: all β > 0 fail at n₀ = 2, 4 | — | no certificate |

## Model

A ring of 2n spin-1/2, `H = Σ (1+δ) S·S` on the bond inside each two-spin cell
`+ Σ (1−δ) S·S` on the bond between neighbouring cells (`S = σ/2`). Lengths count cells; the
cell has one internal bond, so it is palindromic.

## Reflection positivity

Write `s = 1+δ`, `w = 1−δ`, `A_a = S_a ⊗ 1`, `B_a = 1 ⊗ S_a` on a cell. Rotate the second spin
of every cell by `J = [0 1; −1 0]` (a product of one-site unitaries; spectrum unchanged). Then

```
D = s (−A_x B_x + A_y B_y − A_z B_z),      −h = w (B_x ⊗ A_x − B_y ⊗ A_y + B_z ⊗ A_z).
```

With `F` the swap of the two spins inside a cell, `Θ = F K` and the channels
`(B_x, i B_y, B_z)`: `Θ² = 1`, `Θ(D) = D` and `−h = w Σ_a C_a ⊗ Θ(C_a)`, the antiunitary bond
form of Appendix E of `notes/thermal-mpo-gap-scope.tex`. In the cell basis `M/√2`,

```
M = [1 1 0 0; 0 0 1 i; 0 0 1 −i; 1 −1 0 0],   columns |00⟩+|11⟩, |00⟩−|11⟩, |01⟩+|10⟩, i(|01⟩−|10⟩)
```

this is plain conjugation: `−h = Σ C_a ⊗ conj(C_a)` with `C_a = √w K_a/2` (explicit complex
`K_a` in [METHOD.md](METHOD.md)) and the real diagonal on-site term `s·diag(−3, 1, 1, 1)/4`.
`src/rp.jl` verifies all identities in `Complex{Rational{BigInt}}` at the six deltas; every seed
certificate stores the model and factors as a core `RingModels` record, and the checker
verifies them with `RingModels.verify_rp` (a sign-flipped bond is rejected in the tests).

*What does not hold.* Any `Σ c C ⊗ C` commutes with the swap of two whole cells, but
`‖h − F_outer h F_outer‖²_F = 6(1−δ)² > 0`, and a common cell basis change preserves this;
blocks of `c ≥ 2` spins and alternate-cell reversal fail too (`src/reflection_conventions.jl`).
The criterion does not need the identical-factor form.

## Seeds and certified bounds

Seed n₀ = 4 cells: complete verified spectra of the 8- and 16-spin rings. The defects are
rounded up to a 10⁻⁶ grid plus one spare unit (so that an independent replay on another BLAS
proves the same stored values); `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)` gives the
bound, shown rounded down. Unequal defects are allowed: δ = 1/4 closes although εₛ exceeds the
symmetric sufficient threshold 0.12731.

| δ | β₀ | εₛ ≥ 1 − Pₛ(4,β₀) | εₜ ≥ 1 − P_τ(8,β₀) | bound | enclosed 16-spin gap | bound/gap₁₆ |
|---|---|---|---|---:|---|---:|
| 1/2 | 6 | 0.011267 | 0.012008 | 0.534756 | [1.188650, 1.188653] | 0.4499 |
| 3/10 | 13/2 | 0.095116 | 0.051608 | 0.173197 | [0.811143, 0.811146] | 0.2135 |
| 1/4 | 29/4 | 0.165679 | 0.051038 | 0.058225 | [0.709346, 0.709349] | 0.0821 |

The table's decimals are the exact stored rationals; the exact binary recurrence value, the
three log-partition enclosures and every enclosed eigenvalue are in `results/seed_delta_*.json.gz`
(indexed by `results/certificate.json`). The 16-spin gap intervals are rigorous (they are not
the thermodynamic gap). The survey's seeds are re-certified independently of the β grid:
**0.533681 at (δ, n₀, β₀) = (1/2, 4, 25/4)** and **0.173197 at (3/10, 4, 13/2)** (the
β-grid seed itself), both above the survey's floating-point values (0.533 and 0.172 in AGENTS.md). The previous certificate
of this directory (εₛ, εₜ = 0.011323, 0.012237 / 0.095159, 0.051792 / 0.165720, 0.051228;
bounds 0.532658, 0.172916, 0.057967) is superseded: the sharper two-pass rounding bound
(METHOD.md) shrinks the eigenvalue radii.

**Smaller δ: rigorous exclusion at the accessible seeds.** With the largest ring at 16 spins
(20 for diagnostics) the only even seeds are n₀ = 2, 4. For each δ and each such n₀ without a
seed, `negative_cover` proves that `certified_gap_bound` fails for every β > 0: thermal
monotonicity excludes β ≤ 1, a ground-energy defect bound excludes β ≥ B, and 6–12 interval
boxes cover [1, B], each propagating rigorous *lower* bounds on both defects to z ≥ 1/2. n₀ = 2
is excluded at all six δ, n₀ = 4 at 1/5, 3/20, 1/10. This is a statement about this recurrence
at these seeds, not about larger rings, the survey's Theorem 2.24(v) seed replacement, other
methods or the physical gap.

**Numerical diagnostics (not certified).** `scripts/numerical_scan.jl` (complete Float64
spectra up to 16 spins; `results/numerical_summary.json`) and `scripts/scan20.jl` (20 spins,
every sector, 40 Krylov Ritz values in large blocks; `results/diagnostic20_summary.json`):

| δ | Δ₂₀ (numerical) | ξ_eff (cells) | ξ_eff δ^{2/3} | first even n₀ by δ^{−2/3} scaling |
|---|---:|---:|---:|---:|
| 1/2 | 1.188631 | 0.5624 | 0.3543 | 4 |
| 3/10 | 0.810214 | 0.8448 | 0.3786 | 4 |
| 1/4 | 0.707193 | 0.9682 | 0.3842 | 4 |
| 1/5 | 0.599819 | 1.1430 | 0.3909 | 6 |
| 3/20 | 0.488096 | 1.4171 | 0.4001 | 6 |
| 1/10 | 0.374281 | 1.9130 | 0.4122 | 8 |

Here `ξ_eff = 1/log(d₄/d₅)` with `d_n = E₀(2n cells) − 2E₀(n cells)`: an energy-defect decay
length, not a correlation length. The rescaled column drifts, so these short rings do not fix
an exponent. Seeds n₀ = 6, 6, 8 need 24, 24, 32-spin doubled rings. The odd n₀ = 5 seed (outside
the even-seed theorem, partial spectrum) closes numerically at δ = 1/5 (β₀ = 8.5, about 0.086)
and fails on the β grid at 3/20 and 1/10.

**Beyond ED.** `dimer_bounds` is a rigorous norm sandwich around decoupled dimers at any ring
size (cost independent of the Hilbert dimension); the tests close a 512-spin seed at
δ = 0.999999. Its extensive error already gives εₛ > 1/2 at n₀ = 4, β = 6 for every requested δ,
so it certifies nothing there. Route assessment: [IDEAS.md](IDEAS.md).

## Reproduce

From the repository root, with `OPENBLAS_NUM_THREADS=1` (times on an Apple M-series laptop):

```sh
julia --project=models/dimerized_heisenberg -e 'using Pkg; Pkg.instantiate()'
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/check.jl         # ~15 s
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/check.jl --full  # ~3 min
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/test/runtests.jl # ~30 s
```

`check.jl` re-derives the exact RP audit and, per δ, the stored RP factors, every block's
residual/orthogonality/radius witness from its stored maxima, sector completeness and the
reused conjugate blocks, all seeds, log-partition enclosures and recurrence values, the 16-spin
gap interval, the dimer bound and every exclusion box and tail; it checks the source hashes and
that `results/certificate.json` indexes the certificates. `--full` re-solves every block from the
model and repeats all of this from the fresh spectra alone, so no stored eigenvalue or residual
maximum is trusted. The tests add exact rational residuals (all 4-spin blocks at six δ, all
8-spin blocks at three δ), the analytic four-spin spectrum, conjugate-block equality, the sector
solvers against the core dense ED, the dimer bound, and rejection of corrupted certificates.

Regeneration (not needed for checking):

```sh
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl  # ~3.5 min, all
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl seeds 3 10
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/scan20.jl 1 5
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/scan20.jl summary
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/independent_review.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/numerical_scan.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/summarize_numerics.jl
```

`certify.jl` takes a step (`rp`, `seeds [p q]`, `summary`; default all); each δ takes about half
a minute. Each `scan20.jl p q` run takes about 2 min and writes the gitignored
`results/diagnostic20_p_q.json` (~390 kB of Ritz values); `summary` condenses the six into the
tracked summary and requires each run to match the current seed certificate's hash.
`scripts/independent_review.jl` (~2 min) additionally requires every fresh eigenvalue interval
to lie inside the stored one, a stricter test than the theorem needs that may fail on a
different BLAS; it writes the gitignored `results/acceptance.json`.

## What is rigorous

- **Rigorous**: the RP identities (exact rationals); sector completeness; every eigenvalue
  enclosure; the outward partition sums, the exact rational seeds and the bounds; uniqueness
  (from purity > 1/2, not from ED); the 16-spin gap intervals; the all-β exclusions; the dimer
  perturbation bounds. The bound holds for even L ≥ 4 cells only (N divisible by 4).
- **Trusted assumptions**: residual and Gram products computed by classical binary64 CPU
  OpenBLAS GEMM (no Strassen, no reduced precision), bounded by the standard rounding model
  with explicit gamma and underflow terms; BLAS, Julia, IntervalArithmetic and LAPACK's returned
  arrays are not formally verified. LAPACK accuracy itself is not trusted. Eigensolves are not
  bitwise reproducible across machines; the 10⁻⁶ seed grid absorbs this for `--full`.
- **Numerical only**: the β grid that selects candidate seeds, the δ ≤ 1/5 Float64 scans, the
  20-spin Ritz spectra, the odd n₀ = 5 candidate, the ξ_eff fits and the δ^{−2/3} seed scales.

## Layout

| path | role |
|---|---|
| `certify.jl` | writes `results/rp_audit.json`, `seed_delta_*.json.gz`, `certificate.json` |
| `check.jl` | checker (also loaded by the tests); `--full` re-solves every spectrum |
| `src/Certification.jl` | seeds, all-β exclusions, dimer bound, stored and fresh checks |
| `src/verified_spectra.jl` | interval block construction, two-pass witnesses, Weyl, `log Z` |
| `src/rp.jl`, `src/reflection_conventions.jl` | exact RP audits |
| `src/sectors.jl` | symmetry orbits; Float64 sector solver for the numerical scans |
| `src/Diagnostics.jl` | 20-spin sparse Krylov diagnostics and their summary |
| `src/provenance.jl` | gzip reader and source-hash checks |
| `scripts/` | 20-spin scan, numerical scan and summary, independent containment review |
| `test/runtests.jl` | fast tests |
| `METHOD.md` | RP matrices, enclosure theorem, recurrence, exclusions, dimer bound |
| `IDEAS.md` | route assessments and outcome |
| `results/certificate.json` | index: per δ the certificate hash, seeds, exclusions, gap₁₆ ratio |
| `results/seed_delta_*.json.gz` | full certificates (six δ, every block, witness, eigenvalue) |
| `results/rp_audit.json` | exact RP audit at six deltas (no seeds) |
| `results/numerical_summary.json`, `results/diagnostic20_summary.json` | Float64 summaries |

Every result records the SHA-256 of its sources by path and is regenerated in minutes; there
are no layout exceptions.

## History

PR #23 re-certified the grid: exact RP factors stored as `RingModels` records, the two-pass
rounding bound, conjugate-momentum reuse, a checker that can rebuild every eigenvector and
residual, all-β exclusions establishing the grid minimum δ = 1/4, the dimer perturbation bound,
20-spin diagnostics, and certificates for all six δ (`seed_delta_*`, replacing the three
`verified_delta_*` files, whose bounds 0.532658, 0.172916, 0.057967 are superseded). Earlier
results were produced under the core's former name `ThermalMPOGap` from
`certificates/dimerized_heisenberg/`. Removed over time (all in the git history): the 10 MB
review output, uncompressed certificates, run logs, `PR_BODY.md`, `REPORT.md`,
`validate.jl`/`validation.json` (timings now in this README), `review/independent_review.md`,
and the per-δ Float64 scans `results/numerical_delta_*.json` and 20-spin runs
`results/diagnostic20_*.json` (regenerate with the scripts above).
