# Staggered-field spin-1/2 Heisenberg certificate

Model: periodic spin-1/2 chain H = Σ_j S_j·S_{j+1} − h Σ_j (−1)^j S^z_j, exchange 1.
See [REPORT.md](REPORT.md) for proof details, sector tables and numerical exploration,
[AUDIT.md](AUDIT.md) for the review and the full spectral replay, and [IDEAS.md](IDEAS.md)
for the route choice.

## Result

| h | n0 (sites) | β0 | certified gap lower bound | numerical 16-site gap |
|---|---:|---:|---:|---:|
| 1 | 8 | 21/4 | 1.148261 | 1.87351193 |
| 3/4 | 8 | 21/4 | 0.890782 | 1.57089206 |
| 1/2 | 8 | 11/2 | 0.594785 | 1.23131976 |
| 7/20 | 8 | 23/4 | 0.383508 | 0.994380622 |
| 1/4 | 8 | 25/4 | 0.214486 | 0.812031433 |

Each bound holds, together with uniqueness of the ground state, for **every even L ≥ 8**
(lengths count single spins). Published bounds are rounded down to six decimals from the
smaller of the original and the freshly replayed outward bound; both rational seed pairs
prove them (checked in `test/runtests.jl`).

## RP basis

Rotate odd sites by U = 2·iS^y (a π rotation about y; only consistent on even rings). With
the exact rational matrices X = S^x, J = iS^y (real antisymmetric), Z = S^z, the rotated bond
is h = −(X⊗X + J⊗J + Z⊗Z), so −h = Σ_α C_α ⊗ C_α with real C_α ∈ {X, J, Z}, and the
staggered field becomes the uniform real on-site term o = −h Z. `exact_rp` in `common.jl`
checks these identities in rational arithmetic for both bond orientations.

## Seeds

`results/certificate_<h>.json` stores n0, β0 and the exact rational seed defects
ε_s ≥ 1 − P_s(n0, β0), ε_t ≥ 1 − P_τ(2n0, β0) (strings `"num//den"`), plus every 8- and
16-site momentum/S^z block with its multiplicity, all eigenvalues and a rigorous Weyl radius.
The radii come from computed residual and orthogonality norms with explicit IEEE rounding
allowances, a polar-factor bound, and interval enclosures of the exact block entries
(REPORT.md). log Z is summed in IntervalArithmetic over every level with its multiplicity.
`results/spectral_replay_summary.json` stores the seeds from a fresh, complete re-execution.

## Commands

From the repository root, Julia 1.12, `OPENBLAS_NUM_THREADS=1`. Either the core project
(`--project=julia`) or this directory's own project
(`--project=models/heisenberg_staggered_field`, which develops `../../julia`) works.

```sh
# fast checker (about 30 s): RP, small-ring ED comparison, sector counts,
# enclosure propagation from stored eigenvalues/radii, rational seeds, recurrence
julia --project=julia models/heisenberg_staggered_field/test/runtests.jl
# independent scalar-interval accumulation of residuals on small blocks (about 25 s)
julia --project=julia models/heisenberg_staggered_field/audit_bounds.jl
# full spectral replay of all 270 blocks for the five fields (about 3-4 min)
julia --project=julia models/heisenberg_staggered_field/spectral_replay.jl [--full]
# recertify one field at a given seed (about 1 min each)
julia --project=julia models/heisenberg_staggered_field/certify.jl '1//1' 8 '21//4'
# numerical seed search (writes spectrum_*/scan_* caches into results/)
julia --project=julia models/heisenberg_staggered_field/scan.jl 16
julia --project=julia models/heisenberg_staggered_field/scan.jl 20 '1//4'
```

The fast checker does not rebuild the large eigenvector matrices; `spectral_replay.jl` and
`certify.jl` do. `certify.jl h maxL` (two arguments) reads the seed from the scan cache.

## Rigorous vs numerical

Rigorous (under IEEE binary64 round-to-nearest, correctly rounded scalar operations, BLAS
without fast-math, and IntervalArithmetic's outward rounding): the RP identity, complete
sector coverage, eigenvalue enclosures, log Z intervals, rational seeds and the recurrence
bounds. The radius argument uses a priori rounding-error bounds on Float64 matrix products
rather than interval matrix products; see REPORT.md for the derivation.

Numerical only: the (n0, β0) grid search, reference gaps and ratios, the 20-site partial
(Krylov) spectra and candidates, the apparent search threshold between h = 0.10 and 0.15,
and scaling extrapolations. No physical threshold is claimed.

## Removed from the tree

To keep the directory small (39 MB in 183 files before, about 4.5 MB now), the following
were removed. All are in the branch history at commit `cb635a1` and can be regenerated:

- `results/spectrum_<h>_<L>.json`, `results/scan_<h>_<L>.json` (numerical spectra and grid
  scans): regenerate with `scan.jl 16` and `scan.jl 20 <h>`.
- `results/spectral_replay.json`, `results/migration_spectral_replay.json` (5.7 MB each,
  per-block fresh witnesses) and the two failed replay attempts with their scripts:
  `spectral_replay.jl --full` regenerates the witnesses; the compact
  `results/spectral_replay_summary.json` keeps the exact fresh seeds and bounds.
- `historical_sources/` (the pre-rename `ThermalMPOGap` source snapshot named by the hashes
  in the certificate metadata), the rename-migration scripts and records (`MIGRATION.md`,
  `environment.jl`, `historical_provenance.jl`, `migration_check.jl`, `seal.jl`, manifests),
  `report.jl` (generated REPORT.md from the caches), and all run logs.

The `metadata.source_sha256` hashes inside the certificate JSONs refer to that pre-rename
snapshot (recoverable from `cb635a1`). Since then `common.jl` and `scan.jl` changed only
by the package name `ThermalMPOGap` → `CertifiedGaps`; `certify.jl` gained the explicit-seed
option and `test/runtests.jl` gained the published-bound testset.

## Cross-check with the shared core

`core_certificate.jl` certifies the h = 1 seed independently with
`CertifiedGaps.GapCertificates`: the same rotated basis and RP factors, its own verified block
spectra (S^z parity and momentum blocks), and the exact rational chain instead of the interval
bound. `results/core_certificate_1_1.json` gives gap > 1.147699 for every even L ≥ 8, with
ε_s = 692723528103//10^15 and ε_t = 728796111263//10^15. Check it with
`check_certificate(read_exact_json(path))` (≈ 5 min; `spectra = false` checks only the exact
parts in seconds).
