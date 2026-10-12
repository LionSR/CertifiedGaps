# Two-leg transverse-field Ising ladder

**Certified at g = 3: Δ_L ≥ 0.636386863 for every even L ≥ 4 rungs** (N = 2L spin sites, N
divisible by four, N ≥ 8). The returned binary bound is `0.6363868632722607`
(= `179126977518527//281474976710656`). The 8-rung ED gap is 2.6618 (numerical), so the bound
is about 24% of it.

## Model, RP basis and seed

One cell is a rung (two spins, dimension 4). In Pauli normalization,

    H_L = −Σ_j (Z_{j,a} Z_{j+1,a} + Z_{j,b} Z_{j+1,b} + Z_{j,a} Z_{j,b}) − g Σ_j (X_{j,a} + X_{j,b})

on a ring of L rungs.

- **RP basis**: real σᶻ cell basis (↑↑, ↑↓, ↓↑, ↓↓). With C₁ = Z⊗I and C₂ = I⊗Z (real,
  diagonal, symmetric), −h_bond = C₁⊗C₁ + C₂⊗C₂ exactly. The on-site term
  o = −Z⊗Z − g(X⊗I + I⊗X) is real symmetric. `rp_identity` checks this with
  `Rational{BigInt}` matrices.
- **Seed**: (n₀, β₀) = (4, 3/2) (n₀ even). Outward-rounded rational defects
  εₛ ≥ 1 − Pₛ(4, 3/2) ≈ 0.11240836635 and εₜ ≥ 1 − P_τ(8, 3/2) ≈ 0.05161523886; the exact
  fractions are `epsilon_s` and `epsilon_t` in `results/certificate_g3_1.json.gz`.
- **Bound**: `CertifiedGaps.certified_gap_bound(εₛ, εₜ, 3//2).bound`.

## What is rigorous and what is not

Rigorous: the exact RP identity; complete dense spectra of all 32 symmetry sectors at 4 and
8 rungs (dimension sums 256 and 65,536), each with a polar/Weyl whole-spectrum enclosure that
includes algebraic matrix-construction error and explicit binary64 BLAS error terms; outward
MPFR log-partition intervals; exact rational seeds; the recurrence. No Lanczos/ARPACK spectrum
enters. The enclosures assume IEEE binary64 round-to-nearest with gradual underflow and
conventional BLAS dot products (see REPORT.md, "Implemented enclosure theorem").

Numerical only: the seed search (g ∈ {2, 5/2, 3, 4}, n₀ ∈ {2, 3, 4}, β = 1/8:1/8:32), the
g = 4 candidate (2.142), the failures at g = 2 and 5/2, the extrapolated larger seeds, ED
reference gaps, and the literature critical field g_c ≈ 1.83214(5) (not a premise).

The fast checker trusts the stored per-sector radii; it recomputes log Z, the rational seeds and
the bound from them. `replay_audit.jl` regenerates every eigenpair and radius from the model
(about 8 minutes) and checks bit-for-bit agreement with the stored certificate.

## Running

From the repository root (Julia 1.12):

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
# Fast checker: exact RP, sector counts, stored enclosures -> log Z -> seeds -> bound (~20 s)
julia --project=julia models/two_leg_tfim/test/runtests.jl
# Same, and record the run in results/tests.json
julia --project=julia models/two_leg_tfim/run_tests.jl
# Full replay from model inputs (~8 min), then its checker
julia --project=julia models/two_leg_tfim/replay_audit.jl results/replay.json
julia --project=julia models/two_leg_tfim/check_replay.jl results/replay.json results/replay_check.json
```

`models/two_leg_tfim/Project.toml` is also a self-contained environment
(`julia --project=models/two_leg_tfim ...`) with the core at `../../julia`.

Last validation (this layout, Julia 1.12.4, CertifiedGaps 0.2.0): fast checker 70/70 tests
passed in 16 s.

## Files

- `src/Ladder.jl`: symmetry blocks, block matrices, verified eigen-enclosures, log Z intervals.
- `src/Provenance.jl`: certificate loader (gzip), historical source-hash lookup.
- `scan.jl`, `extrapolate.jl`, `summarize.jl`: numerical seed scan, extrapolation, report tables.
- `certify.jl`: certifies the best stored seed (writes `certificate_g<g>_certifiedgaps.json`,
  which the fast checker also picks up).
- `replay_audit.jl`, `check_replay.jl`: full fixed-seed replay and its checker.
- `results/certificate_g3_1.json.gz`: the certificate (gzip of the original 5.3 MB JSON; the
  uncompressed bytes are unchanged, SHA-256 `aba98d98…a74928`, checked by the tests).
- `results/scan_g*.json`, `extrapolation.json`, `environment.json`, `tests.json`: small tables.
- REPORT.md (method and results), AUDIT-2026-10-11.md (independent replay), IDEAS.md (routes).

## Removed when moving from `certificates/two_leg_tfim` to `models/two_leg_tfim`

To bring the directory from 26 MB to about 2 MB, these were deleted (they remain in the Git
history of this branch):

- `results/spectra_g{2_1,5_2,3_1,4_1}.json` (4 × 2.5 MB, numerical scan spectra). Regenerate
  with `julia --project=julia models/two_leg_tfim/scan.jl 3//1` (and `2//1`, `5//2`, `4//1`);
  this also rewrites `scan_g*.json`. `extrapolate.jl` and `summarize.jl` need these spectra.
- `results/audit_2026-10-11_replay.json`, `results/migration_2026-10-11_replay.json`
  (2 × 5.3 MB, duplicate full replays agreeing exactly with the certificate) and their check
  records. Regenerate with `replay_audit.jl` / `check_replay.jl` as above.
- `results/migration_2026-10-11_*.json`, `check_migration.jl`, `environment.jl`,
  `MIGRATION-2026-10-11.md`, `PUBLICATION.md`.

Migration note (summary of the removed MIGRATION/PUBLICATION notes): the certificate was made
with the core under its old name `ThermalMPOGap`. The rename to `CertifiedGaps` changed only
the namespace: all 25 `julia/src` files and `src/Ladder.jl` were compared with the audited
snapshot `8a1aad63e1a285392c97bf5cf17e8da03f0b8314` and were identical after the substitution,
and a fresh replay on CertifiedGaps 0.2.0 reproduced every eigenvalue, radius, seed and the
bound exactly. The certificate's `source_sha256` entries refer to that snapshot (directory
`certificates/two_leg_tfim`); the fast checker verifies them against the Git blobs when the
commit is in the local object store and otherwise skips them with a warning. They are
provenance only; the mathematical checks do not depend on them.
