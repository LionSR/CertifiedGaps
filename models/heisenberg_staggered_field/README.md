# Spin-1/2 Heisenberg chain in a staggered field

**Certified at h = 1, 3/4, 1/2, 7/20, 1/4: a unique ground state and Δ_L ≥ 1.148261, 0.890782,
0.594785, 0.383508, 0.214486 for every even L ≥ 8 sites.** The seeds come from verified dense
spectra of complete 8- and 16-site rings. The method, the enclosure theorem, the sector tables,
the twenty-site exploration and the independent audit are in [METHOD.md](METHOD.md); the route
memo is [IDEAS.md](IDEAS.md).

| h | seed (n₀, β₀) | certified gap ≥ | lengths covered |
|---|---|---:|---|
| 1 | (8, 21/4) | 1.148261 | every even L ≥ 8 sites |
| 3/4 | (8, 21/4) | 0.890782 | every even L ≥ 8 sites |
| 1/2 | (8, 11/2) | 0.594785 | every even L ≥ 8 sites |
| 7/20 | (8, 23/4) | 0.383508 | every even L ≥ 8 sites |
| 1/4 | (8, 25/4) | 0.214486 | every even L ≥ 8 sites |

## Model

Periodic spin-1/2 chain `H = Σ_j S_j·S_{j+1} − h Σ_j (−1)^j S^z_j`, exchange 1. The two-site
ring keeps both periodic bonds, as in the core package. After the rotation below the model has
a one-site cell, so lengths count single spins; the rotation closes only on even rings.

## Reflection positivity

Rotate odd sites by U = 2·iS^y (a π rotation about y). With the exact rational matrices
X = S^x, J = iS^y (real antisymmetric) and Z = S^z, the rotated bond is
h = −(X⊗X + J⊗J + Z⊗Z), so −h = Σ_α C_α ⊗ C_α with real C_α ∈ {X, J, Z}, and the staggered
field becomes the uniform real on-site term o = −h Z. `exact_rp` in `src/common.jl` checks
these identities in rational arithmetic for both bond orientations.

## Seeds and certified bounds

| h | seed (n₀, β₀) | εₛ ≥ 1 − Pₛ(8, β₀) | εₜ ≥ 1 − P_τ(16, β₀) | certified bound | 16-site ED gap |
|---|---|---:|---:|---:|---:|
| 1 | (8, 21/4) | ≈ 6.937e-4 | ≈ 7.326e-4 | 1.148261 | 1.87351193 |
| 3/4 | (8, 21/4) | ≈ 2.343e-3 | ≈ 3.122e-3 | 0.890782 | 1.57089206 |
| 1/2 | (8, 11/2) | ≈ 1.074e-2 | ≈ 1.116e-2 | 0.594785 | 1.23131976 |
| 7/20 | (8, 23/4) | ≈ 3.265e-2 | ≈ 2.744e-2 | 0.383508 | 0.994380622 |
| 1/4 | (8, 25/4) | ≈ 7.987e-2 | ≈ 4.476e-2 | 0.214486 | 0.812031433 |

The exact rational defects (`"num//den"`) are `epsilon_s` and `epsilon_t` in
`results/certificate_<h>.json`; for h = 1 they are `3124087506499//4503599627370496` and
`6599110342439//9007199254740992`. Each bound is `CertifiedGaps.certified_gap_bound(εₛ, εₜ,
β₀).bound`. The published values are rounded down to six decimals from the smaller of the
original and a freshly replayed bound, so both rational seed pairs prove them (checked by the
fast checker; METHOD.md, "Independent audit"). ED gaps are numerical.

Cross-check: `scripts/core_certificate.jl` certifies the h = 1 seed independently with
`CertifiedGaps.GapCertificates` (its own RP identity, verified S^z-parity/momentum block
spectra and the exact rational chain): gap > 1.147699 for every even L ≥ 8, with
εₛ = `692723528103//10^15`, εₜ = `728796111263//10^15`, in `results/core_certificate_1_1.json`.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=models/heisenberg_staggered_field -e 'using Pkg; Pkg.instantiate()'
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/check.jl
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/check.jl --full
# Recertify one field at a given seed (rewrites results/certificate_<h>.json):
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/certify.jl \
    '1//1' 8 '21//4'
# Full spectral replay of all five certificates (rewrites results/spectral_replay_summary.json):
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/spectral_replay.jl [--full]
# Interval re-accumulation of residuals on small blocks; shared-core cross-check:
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/audit_bounds.jl
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/core_certificate.jl
# Numerical seed search (all fields to 16 sites; one field to 20 sites):
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/scan.jl 16
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/scan.jl 20 '1//4'
```

The fast checker is the test suite: `check.jl` includes `test/runtests.jl`.
`certify.jl h maxL` (two arguments) reads the seed from `results/scan_<h>_<maxL>.json`. The
scan caches its spectra in `results/spectrum_<h>_<L>.json` (gitignored). `--project=julia`
(the core environment) also works for every command.

Run times (Julia 1.12, one BLAS thread):

- fast checker (`check.jl`, = `test/runtests.jl`): 341 tests, about 45 s including compilation
  (10 s with precompiled packages on an Apple M-series laptop, 2026-10-14);
- `check.jl --full`: adds about 6 minutes (verified spectra of the cross-check; 2 min 41 s in
  total on the Apple M-series laptop);
- `certify.jl`: 44–80 s per field; `scripts/spectral_replay.jl`: about 190 s for all five;
- `scripts/audit_bounds.jl`: about 25 s; `scripts/core_certificate.jl`: about 12 minutes.

## What is rigorous

Rigorous (under IEEE binary64 round-to-nearest, correctly rounded scalar operations, BLAS
without fast-math, and IntervalArithmetic's outward rounding): the RP identity, complete
sector coverage (every 8- and 16-site momentum/S^z block with its multiplicity), the
eigenvalue enclosures, the log Z intervals, the rational seeds and the recurrence bounds. The
radius argument uses a priori rounding-error bounds on Float64 matrix products rather than
interval matrix products (METHOD.md, "Eigenvalue-enclosure theorem").

The fast checker recomputes log Z, the seeds and the bounds from the stored eigenvalues and
radii, but does not rebuild the eigenvectors; `certify.jl` and `scripts/spectral_replay.jl`
do, and `check.jl --full` re-verifies the cross-check's spectra with the shared core. It also
checks that the cross-check certificate encodes the intended Hamiltonian, not merely a
self-consistent one.

Numerical only: the (n₀, β₀) grid search, reference gaps and ratios, the 20-site partial
(Krylov) spectra and candidates, the apparent search threshold between h = 0.10 and 0.15, and
the scaling extrapolations. No physical threshold is claimed.

## Layout

- `certify.jl`, `check.jl`: certificate entry point and checker (fast; `--full` adds the
  verified spectra of the cross-check).
- `METHOD.md`: method, enclosure theorem, sector tables, twenty-site exploration, audit.
- `src/common.jl`: exact RP, momentum blocks with interval entry enclosures, verified
  eigenvalue enclosures, log Z and rational seeds.
- `scripts/`: `scan.jl` (numerical seed scan), `spectral_replay.jl` (full replay),
  `audit_bounds.jl` (independent residual accumulation), `core_certificate.jl` (shared-core
  cross-check).
- `test/runtests.jl`: the fast checker.
- `results/certificate_<h>.json`: the five certificates. `spectral_replay_summary.json`
  (fresh seeds and bounds from the replay, read by the checker), `core_certificate_1_1.json`
  (cross-check), `audit_bounds.json` (residual audit record).

No file is hash-pinned by its current path, so there are no layout exceptions (see History).

## History

The certificates were made when this directory was `certificates/heisenberg_staggered_field`
and the core was named `ThermalMPOGap`. Their `metadata.source_sha256` entries refer to the
snapshot of those sources kept in `historical_sources/` at commit
`cb635a1cdf2ab06630b5b0838fac4356d3c33509`; the fast checker verifies all of them against that
commit when it is in the local object store, and skips them with a warning otherwise.
`METHOD.md` was called `REPORT.md` until the directory-format tidy.
