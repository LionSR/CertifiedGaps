# Dimerized Heisenberg: verified gaps under general reflected RP

Model: spin-1/2 ring of 2n spins, `H = Σ (1+δ) S·S` on the bonds inside each two-spin cell
`+ Σ (1−δ) S·S` on the bonds between neighbouring cells. Lengths count cells.

| δ | seed (n₀ cells, β₀) | ε_s ≥ 1−Pₛ(4, β₀) | ε_t ≥ 1−P_τ(8, β₀) | certified bound | 16-spin ED gap (numerical) |
|---|---|---|---|---|---|
| 1/2 | (4, 6) | 0.011324 | 0.012238 | 0.5326577398571627 | 1.1887 |
| 3/10 | (4, 13/2) | 0.095160 | 0.051793 | 0.17291646480187198 | 0.8111 |
| 1/4 | (4, 29/4) | 0.165721 | 0.051229 | 0.057966591686137485 | 0.7093 |

Defects are rounded up for display. Each bound holds for **every even L ≥ 4 cells** (8, 12, 16,
… spins). The exact rational seeds, log-partition enclosures and bounds are in
[`results/certificate.json`](results/certificate.json) (a plain summary) and in the full
certificates `results/verified_delta_*.json.gz`. No admissible even seed was found on the tested
grid for δ ≤ 0.2 (numerical scan, `results/numerical_*.json`).

**RP basis.** Cells are the strong pairs. Rotate the second spin of every cell by
`J = [0 1; −1 0]`; then `−h = (1−δ) Σ_{C ∈ {Sˣ, iSʸ, Sᶻ}} C ⊗ Θ(C)` with C real, acting on the
second spin of the left cell, and `Θ(X) = F conj(X) F` with F the swap of the two spins of a cell
(the scope note's bond form with an internal reflection). The on-site term, `(1+δ)` times the
rotated `S₁·S₂`, is real and Θ-invariant. In the cell basis `M/√2` (columns `|00⟩+|11⟩`,
`|00⟩−|11⟩`, `|01⟩+|10⟩`, `i(|01⟩−|10⟩)`, so that `F = M Mᵀ/2`) this becomes
`−h = Σ C' ⊗ conj(C')` with complex `C'` and real diagonal on-site term. No common cell basis gives
a swap-symmetric `Σ C ⊗ C` (the swap defect is `6(1−δ)²`), so the stricter real identical-factor
form is not claimed. Exact checks: `rp.jl`, `reflection_conventions.jl`; details in
[REPORT.md](REPORT.md).

**Seed enclosures.** Complete dense Sz/momentum/spin-flip sectors of the 8- and 16-spin rings,
interval-built block matrices, rigorous residual and orthogonality norm bounds, a polar/Weyl
eigenvalue enclosure for every eigenvalue with multiplicity, outward partition sums, exact
rational defects and `certified_gap_bound`. Proof: [ENCLOSURE_PROOF.md](ENCLOSURE_PROOF.md);
independent review: [review/independent_review.md](review/independent_review.md).

**What is rigorous, and what is assumed.** Rigorous: the RP identity (exact rationals), the
enclosures given the arithmetic model, the seeds and the recurrence. Trusted assumption: the
residual and Gram matrices are computed with classical binary64 CPU OpenBLAS GEMM, whose rounding
is bounded by the standard model with explicit gamma/underflow terms; that BLAS is not formally
verified. Numerical only: the reference gaps, the δ ≤ 0.2 scans and the extrapolations.

## Running

From the repository root, with Julia 1.12 (either environment works):

```sh
# Fast check, about a minute: exact RP audit, small-ring exact residuals, and recomputation of
# every scalar witness, partition interval, seed and bound from the stored spectra.
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/test/runtests.jl
OPENBLAS_NUM_THREADS=1 julia --project=models/dimerized_heisenberg \
    models/dimerized_heisenberg/test/runtests.jl
```

The checker reads `results/verified_delta_*.json.gz` through the system `gzip`. The stored results
record the source hashes of their producing code, under the former core name `ThermalMPOGap` and
the former path `certificates/dimerized_heisenberg/`; `provenance.jl` checks them against the git
objects of commits `ef0a3e6` (core) and `fd3c3d5` (certificate sources), and skips with a warning
if those commits are absent from the clone.

Regeneration (not needed for checking; each verified run takes about a minute):

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/certify.jl         # rp_audit.json
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/run_verified.jl 1 2  # δ = 1/2
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/run_verified.jl 3 10
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/run_verified.jl 1 4
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/summarize_certificates.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/numerical_ed.jl     # numerical
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/summarize_numerics.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia models/dimerized_heisenberg/review/independent_review.jl
```

A rerun with the current core writes new source hashes; the eigenvector solves need not be
bitwise reproducible, which the independent review handles by interval containment.

## Removed when moving here from `certificates/`

To bring the directory from about 20 MB to about 2.4 MB, these files were deleted (they remain
in the history at commit `fd3c3d5`):

- `review/acceptance.json` (10.5 MB): output of `review/independent_review.jl`; regenerate with
  the command above (about 2.5 minutes).
- The three `results/verified_delta_*.json` (2.7 MB each) are now stored gzipped, unchanged.
- Run logs (`results/*.log`, `review/*.log`, including the 0.9 MB `bitwise_probe.log`).
- Superseded review probes `review/bitwise_probe.jl`, `review/refinement_probe.jl`.
- Migration helpers for the former core name: `replay_pinned.jl`, `validate_migration.jl`,
  `results/migration_validation.json` (the fast checker now runs on `CertifiedGaps` directly).
- `PR_BODY.md` (a copy of the pull-request text).
