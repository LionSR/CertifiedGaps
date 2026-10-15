# Two-leg spin-1/2 Heisenberg ladder

    H_L = Σ_{j=1}^{L} Σ_{a=1,2} S_{a,j}·S_{a,j+1} + J⊥ Σ_{j=1}^{L} S_{1,j}·S_{2,j}

periodic in j, leg coupling 1, one rung (d = 4) per cell, so L rungs are N = 2L spins. At
L = 2 the leg bonds are counted twice (the repository's ring convention). Seeds are even
n₀; a certificate at n₀ = 4 covers every even L ≥ 4, i.e. **N = 8, 12, 16, …** (N divisible
by 4), not every N ≥ 8.

The certificate is [results/certificate.json](results/certificate.json), produced by
`certify.jl` and verified from scratch by `check.jl`. [IDEAS.md](IDEAS.md) is the route memo;
[PERTURBATION.md](PERTURBATION.md) proves that the independent-rung perturbation route fails
at every seed size and derives the coupling continuation used below.

## Results

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

## Method

**Reflection positivity.** With the exact real matrices X = Sˣ, A = iSʸ, Z = Sᶻ, the
checkerboard rotation u = exp(iπSʸ) = 2A on one spin of every rung (alternating legs, so L
must be even) maps every bond to a cell model with

    -h_cell = Σ_C C ⊗ C,   C ∈ {X⊗I, A⊗I, Z⊗I, I⊗X, I⊗A, I⊗Z},
    onsite  = -J⊥ (X⊗X + A⊗A + Z⊗Z)   (real symmetric).

`exact_rp` (src/rp.jl) checks the rotation identities for both orientations over
`Rational{BigInt}`; `rp_model` stores the factors and the core's `verify_rp` re-checks
them. The rotation is a unitary equivalence, so spectra are computed in the original SU(2)
basis.

**Complete verified spectra** (src/spectra.jl). The basis is split by down-spin number,
rung momentum, leg exchange and (at zero magnetisation) spin inversion, from exact orbits and
stabilisers; opposite magnetisations and conjugate momenta give exact multiplicities, and the
dimension sum must be 4^L. Block entries are enclosed in interval arithmetic and centred to a
Float64 matrix H with ‖H_exact − H‖ ≤ herror (row-sum bound). LAPACK gives tentative (Q, Λ);
with interval Frobenius bounds r ≥ ‖HQ − QΛ‖ and η ≥ ‖Q*Q − I‖ < 1, both including the
explicit floating-point allowance γ_{8n+32} (unit roundoff 2⁻⁵³) and a 2⁻¹⁰⁰⁰ underflow
allowance, the polar factor satisfies ‖Q − U‖ ≤ d = η/(1 + √(1−η)) and Weyl gives the
simultaneous radius ρ = r + (‖H‖ + max|Λ|) d + herror. Every eigenvalue ball, with its
multiplicity, enters outward-rounded sums for log Z₄(β₀), log Z₈(β₀) and log Z₈(2β₀).

**Seed classification** (src/suite.jl). For each coupling the suite recomputes the spectra of
2, 4 and 8 rungs and classifies all 128 grid seeds:

- *certified*: upper defects, rounded up to multiples of 10⁻⁸, close `certified_gap_bound`;
  the gap is rounded down to a multiple of 10⁻⁹;
- *rejected_by_recurrence*: lower defects (rounded down) carry an exact-rational witness.
  With z = 1 − (1−t)(1−s)² and f(z) = z²/[2(1−z)²], once z ≥ 1/2 (or the thermal
  intermediate does) both next components are ≥ 1/2 and the order-preserving recurrence
  cannot converge to zero; intermediate lower bounds are rounded down to 10⁻²⁴;
- *unresolved*: neither (does not occur in the stored certificate).

**All-β exclusion** (J⊥ = 2, 3/2). With a_L ≤ E₀(L) ≤ b_L from the balls, set
F_L(β) = log Tr e^{−β(H_L − a_L)}, which decreases in β. On a window [l, u]:
log Pₛ ≤ F₂ₙ(l) − 2Fₙ(u) − β(a₂ₙ − 2aₙ), evaluated over the whole window (no monotonicity of
Pₛ is assumed), and P_τ is increasing in β (d log P_τ/dβ = 2(⟨H⟩_β − ⟨H⟩_{2β}) ≥ 0), so
P_τ ≤ P_τ(u). Windows are bisected from [0, B] until each has a rejection witness. For
β ≥ B, δ = a₂ₙ − 2bₙ > 0 gives log Pₛ ≤ F₂ₙ(B) − Bδ. The checker verifies contiguous coverage,
every window witness and the tail.

**Coupling continuation.** For |J⊥ − J*| ≤ δ the centred rung sum has norm ≤ L/2, so
s_δ = 1 − (1−s)e^{−2βnδ}, t_δ = 1 − (1−t)e^{−4βnδ} are valid seeds on the whole interval
(PERTURBATION.md); the suite stores the resulting bound for δ = 10⁻⁴.

## Rigour

Rigorous, and re-established by `check.jl` from fresh computation: the RP identities; the
exact symmetry accounting; every eigenvalue enclosure and log Z enclosure; all 1024 grid
decisions with their upper and lower defects and rejection witnesses; the core recurrence
bounds; the all-β exclusions; the continuation bounds; the README claims above
(`CLAIMED_GAPS`); and the SHA-256 of the core sources, `src/`, `certify.jl`, `check.jl` and
the model environment (tests and `scripts/` are not hashed). No stored spectrum or scan is an
input. `audit.jl` independently redoes the J⊥ = 4 seed with 256-bit directed-rounding MPFR
partition sums and checks the exact trace moments Tr H = 0 and Tr H² against the balls.

Assumed: the repository's abstract criterion (survey/notes, `lean/` formalises it with seeds
as hypotheses) and IEEE binary64 round-to-nearest arithmetic in OpenBLAS/LAPACK, for which
the error allowances above are stated.

Not rigorous (floating point, labelled as such, never used by a proof): everything under
`results/scans/` — the dense scans through ten rungs (20 spins, 2²⁰ states), their log Z
grids and finite-ring gaps, the threshold bisection and the seed-scale extrapolation — and all
timings.

## Layout and reproduction

    certify.jl, check.jl, audit.jl   entry points (certificate, checker, independent audit)
    src/HeisenbergLadder.jl          module: rp.jl, spectra.jl, suite.jl
    test/runtests.jl                 fast tests
    scripts/                         floating-point scans and their replay (not certificates)
    results/certificate.json         the certificate (single-line JSON, exact rational strings)
    results/audit.json               output of audit.jl
    results/scan_replay.json         output of scripts/replay_scans.jl
    results/scans/                   floating-point scan evidence

From the repository root, with one BLAS thread:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=models/heisenberg_ladder -e 'using Pkg; Pkg.instantiate()'
julia --project=models/heisenberg_ladder models/heisenberg_ladder/test/runtests.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/certify.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/audit.jl
julia --project=models/heisenberg_ladder models/heisenberg_ladder/scripts/replay_scans.jl
```

Run times (Julia 1.12.4, Apple M5, one BLAS thread, excluding package installation;
wall-clock including Julia start-up in parentheses): tests 6 s; `certify.jl` 87.5 s
recorded in the certificate (93 s), 10–14 s per coupling; `check.jl` 81 s (85 s);
`audit.jl` 8 s (13 s); `scripts/replay_scans.jl` 1 s (3 s). An earlier run of the same
suite on a slower Linux machine took 367 s to certify and 308 s to check.

The scans are regenerated per coupling by `scripts/scan.jl <J>` (20–35 min each, dominated
by the ten-rung ring), then `scripts/extrapolate.jl`; `scripts/threshold_bisection.jl`
proposes the threshold endpoints. Running the six scans concurrently once exceeded 17 GiB;
run them one or two at a time.
