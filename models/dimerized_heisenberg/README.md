# Dimerised Heisenberg chain

**Certified** at δ = 1/2, 3/10 and 1/4: for every even L ≥ 4 cells (8, 12, 16, … spins) the
ground state is unique and the gap is at least the bound below. The RP hypothesis used is the
general internally reflected form; the stronger "real identical factors `C ⊗ C`" form is false
for this bond and is not claimed. No admissible seed exists for δ ≤ 1/5 at the sizes reached
(numerical scan). The Lean statements in `lean/CertifiedGaps/Models/DimerizedHeisenberg*.lean`
take these seeds as hypotheses.

## Model

A ring of 2n spin-1/2, `H = Σ (1+δ) S·S` on the bond inside each two-spin cell
`+ Σ (1−δ) S·S` on the bond between neighbouring cells (`S = σ/2`). Lengths count cells; the
cell has one internal bond, so it is palindromic.

## Reflection positivity

Write `s = 1+δ`, `w = 1−δ`, `A_a = S_a ⊗ 1`, `B_a = 1 ⊗ S_a` on a cell. Rotate the second spin
of every cell by `J = [0 1; −1 0]` (a product of one-site unitaries; spectrum unchanged). Then
the on-site term and the bond are

```
D = s (−A_x B_x + A_y B_y − A_z B_z),      −h = w (B_x ⊗ A_x − B_y ⊗ A_y + B_z ⊗ A_z).
```

With `F` the swap of the two spins *inside* a cell, `Θ = F K` (K complex conjugation) and the
real channels `(C_x, C_y, C_z) = (B_x, i B_y, B_z)`:

```
Θ² = 1,   Θ(D) = D,   Θ(C_x, C_y, C_z) = (A_x, i A_y, A_z),   −h = w Σ_a C_a ⊗ Θ(C_a).
```

This is the antiunitary bond form of Appendix E of `outputs/thermal-mpo-gap-scope.tex`
(`Θ = V K`, `Θ² = ±1`, `−h = Σ c C ⊗ Θ(C)`, `Θ(D) = D`), and here the transfer kernel is even
real symmetric. Equivalently, in the cell basis `W = M/√2`,

```
M = [1 1 0 0; 0 0 1 i; 0 0 1 −i; 1 −1 0 0],   columns |00⟩+|11⟩, |00⟩−|11⟩, |01⟩+|10⟩, i(|01⟩−|10⟩)
```

`W` turns `Θ` into plain conjugation: `−h = w Σ C'_a ⊗ conj(C'_a)` with complex `C'_a` and the
real diagonal on-site term `s·diag(−3/4, 1/4, 1/4, 1/4)`. `src/rp.jl` verifies all identities in
`Complex{Rational{BigInt}}` (using `M` and exact denominators, no `√2`) at the six deltas
1/2, 3/10, 1/4, 1/5, 3/20, 1/10.

*What does not hold.* Any `Σ c C ⊗ C` commutes with the swap `F_outer` of two whole cells, but
`‖h − F_outer h F_outer‖²_F = 6(1−δ)² > 0`, and a common cell basis change preserves this; blocks
of `c ≥ 2` spins have defect `3·4^c w²/8`; reversing alternate cells gives two different bond
operators. So the literal identical-factor form fails in every common cell basis
(`src/reflection_conventions.jl` checks the two- and four-spin blocks exactly). The criterion
does not need it.

## Seed and certified bound

Seed n₀ = 4 cells at every δ: the rings of 8 and 16 spins. Exact rational upper bounds
`εₛ ≥ 1 − Pₛ(4, β₀)`, `εₜ ≥ 1 − P_τ(8, β₀)`, and `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)`:

| δ | β₀ | εₛ (exact) | εₜ (exact) | certified bound | 16-spin ED gap |
|---|---|---|---|---|---|
| 1/2 | 6 | 50996126991899/4503599627370496 ≈ 0.011323 | 13777527984355/1125899906842624 ≈ 0.012237 | 0.5326577398571627 | 1.1887 |
| 3/10 | 13/2 | 857117450063167/9007199254740992 ≈ 0.095159 | 466505354307053/9007199254740992 ≈ 0.051792 | 0.17291646480187198 | 0.8111 |
| 1/4 | 29/4 | 373168159679129/2251799813685248 ≈ 0.165720 | 115354947892529/2251799813685248 ≈ 0.051228 | 0.057966591686137485 | 0.7093 |

The bound is a binary64 number; its exact value (`recurrence_bound_exact_binary`) and the three
log-partition enclosures are in `results/certificate.json` (plain) and in the full
certificates `results/verified_delta_*.json.gz` (with every enclosed eigenvalue). The ED gaps are
numerical, for comparison only.

**How the seeds are enclosed** (proof in [ENCLOSURE_PROOF.md](ENCLOSURE_PROOF.md)). Each ring
is split into Sz, cell-momentum and (at Sz = 0) spin-flip sectors: 21 blocks (largest 14) at 8
spins and 73 blocks (largest 1430) at 16 spins, with weighted dimensions checked to sum to 2^L.
Block matrices are built in interval arithmetic (roots of unity and orbit-size square roots
enclosed). Each block is diagonalised densely by LAPACK; the residual `HQ − QΛ` and the Gram
defect `Q*Q − 1` are then bounded rigorously with explicit binary64 rounding and underflow
allowances (γ_{4d+16}), and a polar-factor/Weyl argument encloses every eigenvalue, with
multiplicity. Largest radii: 2.7e−5 (16 spins), below 4.1e−11 (8 spins). Outward interval sums
give `log Z`; the lower purity uses the lower numerator and upper denominator endpoints.

**Smaller δ (numerical, not certified).** `scripts/numerical_scan.jl` scans n₀ ∈ {2, 3, 4}
cells, β₀ ∈ 0.25:0.25:32 on complete Float64 spectra (n₀ = 3 is diagnostic only: the seed must
be even). Best `max_β min(Pₛ(n), P_τ(2n))`:

| δ | n = 2 | n = 3 | n = 4 | best even seed | 16-spin ED gap |
|---|---|---|---|---|---|
| 1/2 | 0.7935 | 0.9468 | 0.9883 | (4, 6.00) | 1.1887 |
| 3/10 | 0.5974 | 0.7944 | 0.9126 | (4, 6.50) | 0.8111 |
| 1/4 | 0.5352 | 0.7342 | 0.8648 | (4, 7.25) | 0.7093 |
| 1/5 | 0.4880 | 0.6606 | 0.7896 | none | 0.6046 |
| 3/20 | 0.4545 | 0.5797 | 0.7017 | none | 0.4983 |
| 1/10 | 0.4295 | 0.4929 | 0.5836 | none | 0.3950 |

A two-point fit of `E₀(2n) − 2E₀(n) = A e^{−n/ℓ}` plus an independent-excitation model of
P_τ (`scripts/summarize_numerics.jl`, `results/numerical_summary.json`) predicts the first
admissible even seed at n₀ = 6 for δ = 1/5, 3/20 (24-spin doubled ring) and n₀ = 8 for δ = 1/10
(32 spins), beyond dense ED. ℓ·δ^{2/3} ≈ 0.35–0.40 is roughly constant. These are unvalidated
heuristics.

## Reproduce

From the repository root, with `OPENBLAS_NUM_THREADS=1` (times on an Apple M-series laptop):

```sh
julia --project=models/dimerized_heisenberg -e 'using Pkg; Pkg.instantiate()'
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/check.jl         # ~10 s
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/test/runtests.jl # ~20 s
```

`check.jl` re-derives the exact RP audit and, from the stored spectra, every block's residual,
orthogonality and radius witness, the three log-partition intervals, the rational seeds and the
bound, and requires exact equality; it also checks the source hashes. The tests add
independent checks: the sector solver against the core dense ED (4–10 spins), exact rational
residuals and the analytic four-spin spectrum against the enclosures, and the reflection
conventions. Regeneration (not needed for checking):

```sh
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl  # ~3 min, all
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl verified 3 10
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/independent_review.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/numerical_scan.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/summarize_numerics.jl
```

`certify.jl` takes a step (`rp`, `verified [p q]`, `summary`; default all). Each verified run
takes about 50 s. A rerun records hashes of the current sources and, since dense eigensolves need
not be bitwise reproducible, may store slightly different (equally valid) enclosures.
`scripts/independent_review.jl` (about 2.5 min) re-solves all six spectra, proves with a refined
error bound that every fresh eigenvalue interval lies inside the stored one
(`|λ_new − λ_stored| + ρ_new ≤ ρ_stored`, in exact rationals), recomputes all 63 eight-spin blocks'
residuals in exact rational arithmetic, and writes `results/acceptance.json` (10 MB,
gitignored). On the original Linux/OpenBLAS machine it passed 4012/4012 checks; the refined
bounds were 0.53478, 0.17320, 0.05823, each at least the stored one. On an Apple M-series
machine (2026-10-14) the fresh solves differ from the stored ones in most blocks; the fresh
bounds agree, but 2 of 4045 checks fail: in one 8-spin δ = 3/10 block the fresh interval
(radius about 1e−11) is not inside the stored one. Both are valid enclosures, so this only
means that the containment test is stricter than the theorem needs. It does not affect
`check.jl`, and a fresh `certify.jl verified 1 2` on that machine passes the checker's
verified-file checks.

## What is rigorous

- **Rigorous**: the RP identities and the on-site invariance (exact rationals); sector
  completeness; every eigenvalue enclosure; the outward partition sums, the exact rational
  seeds and the bound; uniqueness (from purity > 1/2, not from ED). The bound holds for even
  L ≥ 4 cells only, not for odd cell counts.
- **Trusted assumptions**: the residual and Gram products are computed by classical binary64
  CPU OpenBLAS GEMM (no Strassen, no reduced precision), whose rounding is bounded by the
  standard model with explicit gamma and underflow terms; that BLAS, Julia, IntervalArithmetic
  and LAPACK's returned arrays are not formally verified. LAPACK accuracy itself is not trusted:
  the residuals are bounded afterwards.
- **Numerical only**: the ED gaps and bound/gap ratios, the δ ≤ 1/5 scan, the ℓ fits and the
  predicted seed scales.

## Layout

| path | role |
|---|---|
| `certify.jl` | writes `results/rp_audit.json`, `verified_delta_*.json.gz`, `certificate.json` |
| `check.jl` | checker for all stored results (also loaded by the tests) |
| `src/rp.jl` | exact audit of the reflected, conjugate and identical-factor forms |
| `src/reflection_conventions.jl` | exact audit of other cell orientations and 4-spin blocks |
| `src/sectors.jl` | symmetry orbits; Float64 sector solver for the numerical scan |
| `src/verified_spectra.jl` | interval block construction, residual/Weyl enclosures, `log Z` |
| `src/provenance.jl` | gzip reader and source-hash checks (current and historical) |
| `scripts/` | independent review, numerical scan and its heuristic summary |
| `test/runtests.jl` | fast tests |
| `ENCLOSURE_PROOF.md` | the enclosure theorem and the asymmetric recurrence |
| `IDEAS.md` | route choice |
| `results/certificate.json` | the three seeds and bounds as exact rationals (plain JSON) |
| `results/verified_delta_*.json.gz` | full certificates (every block, witness, eigenvalue) |
| `results/rp_audit.json` | exact RP audit at six deltas (no seeds) |
| `results/numerical_summary.json` | Float64 scan summary and heuristic extrapolation |

## History

The verified certificates and `rp_audit.json` were produced under the core's former name
`ThermalMPOGap` (core commit `ef0a3e6`) from sources under `certificates/dimerized_heisenberg/`
(commit `fd3c3d5`). Their `source_sha256` records are left untouched; `src/provenance.jl`
compares them with the git objects of those commits (and skips with a warning if the commits
are absent from the clone). `results/certificate.json` hashes the current `certify.jl`.
Removed over time (all in the git history): the 10 MB review output, uncompressed certificates,
run logs, bitwise-reproduction probes, migration helpers, `PR_BODY.md`, `REPORT.md` and
`review/independent_review.md` (folded into this README), and the per-δ Float64 scans
`results/numerical_delta_*.json` (6 × 136 kB; regenerate with `scripts/numerical_scan.jl`).
