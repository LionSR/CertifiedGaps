# Two-leg spin-1/2 Heisenberg ladder

Certified at rung couplings J⊥ = 4, 7/2, 3, 5/2 and 163069/65536 from verified complete
spectra of rings of 2, 4 and 8 rungs: a unique ground state and Δ_L ≥ the bound below for every
even L ≥ 4 rungs. At J⊥ = 2 and 3/2 every seed with n₀ ∈ {2, 4} is rigorously rejected for every
β > 0. The certificate is [results/certificate.json](results/certificate.json); the route memo is
[IDEAS.md](IDEAS.md); the pipeline and the failed perturbation route are in
[METHOD.md](METHOD.md).

| J⊥ | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| 4 | (4, 2) | 1.102403091 | every even L ≥ 4 rungs (N = 8, 12, 16, … spins) |
| 7/2 | (4, 9/4) | 0.778648190 | every even L ≥ 4 rungs |
| 3 | (4, 5/2) | 0.454189334 | every even L ≥ 4 rungs |
| 5/2 | (4, 3) | 0.042433721 | every even L ≥ 4 rungs |
| 163069/65536 | (4, 3) | 0.001531165 | every even L ≥ 4 rungs |
| 2, 3/2 | none (n₀ ∈ {2, 4} rejected for all β) | — | no certificate |

## Model

    H_L = Σ_{j=1}^{L} Σ_{a=1,2} S_{a,j}·S_{a,j+1} + J⊥ Σ_{j=1}^{L} S_{1,j}·S_{2,j}

periodic in j, leg coupling 1, one rung (d = 4) per cell, so L rungs are N = 2L spins. At
L = 2 the leg bonds are counted twice (the repository's ring convention). Seeds are even
n₀; a certificate at n₀ = 4 covers every even L ≥ 4, i.e. **N = 8, 12, 16, …** (N divisible
by 4), not every N ≥ 8.

## Reflection positivity

With the exact real matrices X = Sˣ, A = iSʸ, Z = Sᶻ, the checkerboard rotation
u = exp(iπSʸ) = 2A on one spin of every rung (alternating legs, so L must be even) maps every
bond to a cell model with

    -h_cell = Σ_C C ⊗ C,   C ∈ {X⊗I, A⊗I, Z⊗I, I⊗X, I⊗A, I⊗Z},
    onsite  = -J⊥ (X⊗X + A⊗A + Z⊗Z)   (real symmetric).

`exact_rp` (src/rp.jl) checks the rotation identities for both orientations over
`Rational{BigInt}`; `rp_model` stores the factors and the core's `verify_rp` re-checks
them (both in `check.jl` and `check.jl --full`). The rotation is a unitary equivalence, so
spectra are computed in the original SU(2) basis.

## Seeds and certified bounds

Every positive row: unique ground state and Δ_L ≥ bound for all even L ≥ 4 rungs. The
defects are the exact rational upper bounds stored in the certificate (decimals are exact).

| J⊥ | seed (n₀, β₀) | εₛ | εₜ | certified Δ ≥ | on \|J⊥ − J*\| ≤ 10⁻⁴ |
|---|---|---:|---:|---:|---:|
| 4 | (4, 2) | 0.03508321 | 0.02462420 | 1.102403091 | 1.063075969 |
| 7/2 | (4, 9/4) | 0.05806853 | 0.02987064 | 0.778648190 | 0.751905791 |
| 3 | (4, 5/2) | 0.09896358 | 0.04497730 | 0.454189334 | 0.435856503 |
| 5/2 | (4, 3) | 0.18558434 | 0.05114736 | 0.042433721 | — |
| 163069/65536 ≈ 2.488235 | (4, 3) | 0.18772796 | 0.05269516 | 0.001531165 | — |
| 326133/131072 ≈ 2.488197 | all 128 grid seeds rejected | | | — | |
| 2 | n₀ ∈ {2, 4} rejected for every β > 0 | | | — | |
| 3/2 | n₀ ∈ {2, 4} rejected for every β > 0 | | | — | |

As exact rationals, the J⊥ = 4 seed is (n₀, β₀, εₛ, εₜ) = (4, 2, 3508321/100000000,
123121/5000000). Unequal defects beyond the diagonal sufficient cutoff 0.12731 (J⊥ ≤ 5/2) are
accepted only because the core recurrence `certified_gap_bound` closes for them.

**How.** For each coupling the suite encloses every eigenvalue of the rings of 2, 4 and 8
rungs (symmetry blocks with exact multiplicities, interval block entries, residual and
orthogonality bounds, Weyl), forms outward log Z sums, and classifies all 128 grid seeds
n₀ ∈ {2, 4}, β₀ = k/4 (k = 1, …, 64) as *certified* (upper defects close the recurrence),
*rejected_by_recurrence* (lower defects carry an exact witness that the recurrence cannot
converge) or *unresolved* (does not occur). At J⊥ = 2 and 3/2 temperature windows and a
ground-energy tail extend the rejection to every β > 0. Certified seeds are continued over
|J⊥ − J*| ≤ 10⁻⁴ with a norm bound on the rung sum. Details: METHOD.md, "Certificate pipeline".

**Finite-grid threshold.** On the grid n₀ ∈ {2, 4}, β₀ = k/4 (k = 1, …, 64), every seed is
rigorously rejected at J⊥ = 326133/131072 and a seed is certified at J⊥ = 163069/65536, so
the method's grid transition lies in (2.48819732666015625, 2.4882354736328125). This is a
statement about this seed grid only: not a continuous-β optimum, not monotone in J⊥ (not
assumed), and not the physical gap-closing coupling (the ladder is believed gapped for all
J⊥ > 0). The endpoints come from the floating-point bisection in
`scripts/threshold_bisection.jl`; the suite re-decides both rigorously.

**Negative rows.** At J⊥ = 2 and 3/2 the rejection covers every β > 0, but only for
n₀ = 2, 4. They are failures of these seeds, not evidence of a gapless ladder. The heuristic
fit in `results/scans/extrapolation.json` suggests n₀ = 6 (J⊥ = 2) and n₀ = 8 (J⊥ = 3/2),
i.e. complete traces of 24 and 32 spins; that estimate is unverified and unused.

## Reproduce

From the repository root, with one BLAS thread:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/heisenberg_ladder -e 'using Pkg; Pkg.instantiate()'
julia --project=models/heisenberg_ladder models/heisenberg_ladder/certify.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl --full
julia --project=models/heisenberg_ladder models/heisenberg_ladder/test/runtests.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/scripts/audit.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/scripts/replay_scans.jl
```

`check.jl` (fast) re-derives the source hashes, the exact RP identities and stored model,
every grid decision from its stored rational defects (recurrence bound or rejection witness),
the best seeds, the continuation and the README claims (`CLAIMED_GAPS`), without spectra.
`check.jl --full` verifies everything from scratch: fresh verified spectra, every defect and
log Z enclosure, all 1024 grid decisions and the all-β exclusions.

Run times (Julia 1.12.4, Apple M5, one BLAS thread, excluding package installation;
wall-clock including Julia start-up in parentheses): tests 6 s; `certify.jl` 87.5 s on an idle
machine (93 s; 175 s when rerun alongside another job), 10–30 s per coupling; `check.jl` 2 s
(4 s); `check.jl --full` 90 s (93 s); `scripts/audit.jl` 11 s (17 s); `scripts/replay_scans.jl`
1 s (4 s). An earlier run of the same suite on a slower Linux
machine took 367 s to certify and 308 s to check from scratch.

The scans are regenerated per coupling by `scripts/scan.jl <J>` (20–35 min each, dominated
by the ten-rung ring), then `scripts/extrapolate.jl`; `scripts/threshold_bisection.jl`
proposes the threshold endpoints. Running the six scans concurrently once exceeded 17 GiB;
run them one or two at a time.

## What is rigorous

Rigorous, and re-established by `check.jl --full` from fresh computation: the RP identities;
the exact symmetry accounting; every eigenvalue enclosure and log Z enclosure; all 1024 grid
decisions with their upper and lower defects and rejection witnesses; the core recurrence
bounds; the all-β exclusions; the continuation bounds; the README claims above
(`CLAIMED_GAPS`); and the SHA-256 of the core sources, `src/`, `certify.jl`, `check.jl` and
the model environment (tests and `scripts/` are not hashed). No stored spectrum or scan is an
input. The fast `check.jl` re-establishes the parts that need no spectrum and trusts the stored
defect enclosures. `scripts/audit.jl` independently redoes the J⊥ = 4 seed with 256-bit
directed-rounding MPFR partition sums and checks the exact trace moments Tr H = 0 and Tr H²
against the balls.

Assumed: the repository's abstract criterion (survey/notes, `lean/` formalises it with seeds
as hypotheses) and IEEE binary64 round-to-nearest arithmetic in OpenBLAS/LAPACK, for which
the error allowances in METHOD.md are stated.

Not rigorous (floating point, labelled as such, never used by a proof): everything under
`results/scans/` — the dense scans through ten rungs (20 spins, 2²⁰ states), their log Z
grids and finite-ring gaps, the threshold bisection and the seed-scale extrapolation — and all
timings.

## Layout

    certify.jl, check.jl             entry points (certificate; checker, fast or --full)
    METHOD.md                        certificate pipeline; perturbation route and continuation
    src/HeisenbergLadder.jl          module: rp.jl, spectra.jl, suite.jl
    test/runtests.jl                 fast tests
    scripts/audit.jl                 independent MPFR audit of the J⊥ = 4 seed
    scripts/                         floating-point scans and their replay (not certificates)
    results/certificate.json         the certificate (single-line JSON, exact rational strings)
    results/audit.json               output of scripts/audit.jl
    results/scan_replay.json         output of scripts/replay_scans.jl
    results/scans/                   floating-point scan evidence

Hash pinning: the certificate records the SHA-256 of the core sources, `src/`, `certify.jl`,
`check.jl`, `Project.toml`, `Manifest.toml` and `notes/thermal-mpo-gap-scope.tex` by path,
and both checkers compare them. It is regenerated by `certify.jl` in under two minutes, so
there are no layout exceptions.

## History

The directory-format tidy renamed `PERTURBATION.md` to `METHOD.md` (adding the pipeline
description from this README), moved `audit.jl` to `scripts/`, split `check.jl` into a fast
default and `--full`, and regenerated `results/certificate.json` and `results/audit.json`
because `check.jl` and `src/suite.jl` are hashed; every seed, defect, decision and bound is
identical to the previous certificate.
