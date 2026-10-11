# CertifiedGaps integration — 2026-10-11

Upstream main `22b2759e0790197b83e2571b465aa8dc69ebc0be` is integrated by merge
`ea789e76955549eba56e35b4a9dc755110879e1c`. The certificate now imports
`CertifiedGaps` 0.2.0 with its unchanged package UUID. The requested certificate
location remains `certificates/two_leg_tfim/`; this explicit task scope overrides the
new general preference for `models/`. No shared core code was customized.

## Mathematics and provenance

`check_migration.jl` compares all 25 files under `julia/src/` against the published
historical snapshot `8a1aad63e1a285392c97bf5cf17e8da03f0b8314`. Every file is identical
after replacing the namespace `ThermalMPOGap` with `CertifiedGaps` and renaming the
entry-point file. The recurrence implementation is byte-identical even without that
substitution. `src/Ladder.jl` likewise differs only in its two imports. Thus the
migration changes no model, sector construction, spectral theorem, rounding bound,
partition enclosure, seed, or recurrence mathematics. Project metadata and source
hashes do change; those differences are recorded rather than hidden.

All 14 historical JSON results are byte-identical to the pinned snapshot, including
numerical scans, the original certificate, and the first audit. Their original source
hashes remain untouched. `src/Provenance.jl` validates the certificate and replay
source manifests by hashing the actual Git blobs from that snapshot. The eight
numerical scan/spectrum artifacts retain the previously disclosed stale auxiliary
`test/runtests.jl` fingerprints; their embedded hashes are not all claimed to pass.
Preservation of all 14 JSON byte sequences is a separate check. New replay hashes instead describe current files and are
checked against current bytes. A full Git checkout containing the pinned snapshot is
required; shallow clones must fetch history. The published snapshot is tree-identical
to the original audited commit, as recorded in `PUBLICATION.md`.

The migrated replay starts with g=3, n0=4, β0=3/2 and regenerates all dense eigenpairs,
residual and orthogonality bounds, eigenvalue radii, partition intervals, and rational
seeds before reading old numerical values for comparison. Historical radii are never
inputs to that computation. The fresh artifact names begin `migration_2026-10-11_`.

Independent review confirmed the mathematical and provenance separation and caught
a generated-certificate filename that the checker would have missed. It was corrected
to `certificate_g3_1_certifiedgaps.json`, matching the checked `certificate_g` prefix.

## Validation

Julia 1.12.4 loaded CertifiedGaps 0.2.0 after successful instantiation/precompilation,
with OpenBLAS restricted to one thread. The original 69 fast checks passed in
10.799 seconds using the current package. The historical replay artifact passed all
262 comparison and recurrence checks using pinned historical source hashes.
`results/migration_2026-10-11_provenance.json` records the namespace-only comparison
and all 14 preserved result hashes.

The fresh replay completed in **495.324 seconds including serialization**. All 32
sectors were densely diagonalized and verified, with dimensions summing to 256 and
65,536 at four and eight rungs. Every eigenvalue, spectral radius, residual bound,
orthogonality bound, construction bound, and both exact rational seeds matched the
original exactly. The returned bound is unchanged: **0.6363868632722607**, supporting
**Δ_L ≥ 0.636386863 for every even L ≥ 4 rungs** (N=2L sites, N divisible by four,
N≥8). These are the same conditional IEEE/BLAS/MPFR and two-purity theorem assumptions
stated in REPORT.md; package migration introduces no additional mathematical premise.

The migrated replay checker passed **263 tests in 5.472 seconds**, including all
51 current source hashes, exact log-partition interval and recurrence equality,
and reconstruction of both rational seeds. Artifacts:

- `results/migration_2026-10-11_replay.json`: fresh spectra, error bounds, seeds,
  recurrence, current fingerprints, and explicit historical snapshot.
- `results/migration_2026-10-11_check.json`: current replay check and artifact hash.
- `results/migration_2026-10-11_historical_check.json`: historical replay check.
- `results/migration_2026-10-11_tests.json`: fast checker timing.
- `results/migration_2026-10-11_environment.json`: current runtime and dependencies.
- `results/migration_2026-10-11_provenance.json`: source comparison and preservation.

The integration blocker is resolved. No broader scan, new physical claim, CI addition,
core customization, or modification of another certificate directory was made.
