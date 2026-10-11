# Longitudinal-field TFIM certificates

See [REPORT.md](REPORT.md) for the parameter table and certified region, and
[METHOD.md](METHOD.md) for the rigorous eigenvalue and partition-function enclosure proof.
One site is one cell. RP uses the standard Z basis, coefficient 1 and C=Z.
All 30 requested (g,h) points were scanned; certificate JSONs store the exact rational
seeds, full enclosed spectra, source fingerprints and the unmodified core recurrence output.

Run from the repository root with Julia 1.12.4 and the existing `julia` project:

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia certificates/tfim_longitudinal/scan.jl
julia --project=julia certificates/tfim_longitudinal/audit.jl
julia --project=julia certificates/tfim_longitudinal/test/runtests.jl
DEEP_CHECK=1 julia --project=julia certificates/tfim_longitudinal/check.jl
julia --project=julia certificates/tfim_longitudinal/report.jl
```

In this cloud session Julia and its depot are temporary; prefix the commands with:

```sh
export PATH="/tmp/julia-1.12.4/bin:$PATH"
export JULIA_DEPOT_PATH=/tmp/tfim-julia-depot
```

`scan.jl` defaults to even n₀=2,4,6, and β₀=1/4,2/4,…,64. It uses `scan_seeds`
and `best_seed`, then rigorously recomputes only the best seed for each parameter.
Set `NMAX=8` to request larger rings, or pass parameter indices 1–30 to run a subset
(g-major, h-minor order). Such extensions can be much more expensive; the 20–24-site
full dense sectors exceed this machine's memory. Existing results are replaced per point.
The implementation uses complete momentum blocks, with both k and −k computed,
without additional reflection-parity splitting. Spin-flip parity cannot be used for h>0.

The fast checker verifies RP, dimensions, the partition-function intervals from the stored
spectral enclosures, and the bound from rational seeds. It does not merely compare decimal
seeds. The deep checker regenerates all rigorous dense eigensystem enclosures and checks
containment in the stored spectral intervals with exact rational comparisons. The numerical audit compares independent computational-basis ED,
the free-fermion formula at h=0,g=2, and exact rational classical spectra.

The bound is uniform for **every even L≥6 sites** for the default result set. A finite-ring
ED reference gap, its ratio to the bound, scan diagnostics and the free-fermion comparison
are numerical. The JSON eigenvalue enclosures, seed bounds and recurrence lower bound
are certified under IEEE/MPFR/IntervalArithmetic arithmetic as detailed in METHOD.md.
The reported parameter grid is discrete; no interpolation between parameters is certified.

No source outside this directory was changed. No changes were merged or published.

## Verified execution results

- Fast suite: **72/72 passed**, 8.73 seconds including loading and compilation
  (2.6 seconds in the test body), below the requested two-minute limit.
- Deep replay: **30/30 passed**, 248.67 seconds, independently regenerating all
  spectral enclosures and proving containment in the stored intervals.
- Provenance: **1,000 source-hash comparisons passed**.
- Independent audit: **42/42 passed**, including the h=0,g=2 free-fermion cross-check.
- Final full grid scan: 268.32 seconds summed over all parameter points.

The complete verifier is:

```sh
julia --project=julia certificates/tfim_longitudinal/verify_all.jl
```

Machine-readable verification, BLAS configuration and SHA-256 fingerprints of all 30
certificate files are in [results/validation.json](results/validation.json).

## Independent final review

[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) records the separate reviewer's acceptance
with no unresolved theorem gap. Additional reproducible audits are:

```sh
julia --project=julia certificates/tfim_longitudinal/final_audit.jl
julia --project=julia certificates/tfim_longitudinal/independent_residual_audit.jl
```

These passed 1,592 independent stored-artifact checks and 360 exact rational residual/
orthogonality checks, respectively. `results/audit_manifest.json` fingerprints the final
source, documentation, data, logs and the unchanged core/theorem sources.

## Current upstream compatibility

The PR targets upstream's `CertifiedGaps` package rename. All mathematical generation
sources and core Julia sources were checked for exact byte equality to the reviewed
version after replacing only `ThermalMPOGap` with `CertifiedGaps` (29 sources).
`results/namespace_migration.json` records that check. The original per-point certificates
are unchanged; `results/original_generation_sources.json` preserves their exact generation sources
and original hashes at commit `e67c0925f6d069d2dccf34bfbb3dbfc9971c6d0e`.
The verifier checks historical fingerprints against that archive while independently
regenerating spectral enclosures with the current core. Freshly generated results may
instead match the live source hashes. Deep replay and the independent audits were rerun
on the renamed core; current runtime evidence is in their JSON results and logs.

The task's explicit `certificates/tfim_longitudinal/` destination is retained. The PR
changes only this directory relative to current upstream; no core or CI changes are added.
The original no-publication statements describe the earlier local-only audit stage.

The pre-migration validation record is preserved byte-for-byte in
`results/validation_original.json`. `results/validation.json` records the replay on the
current `CertifiedGaps` core, with fresh runtime/source metadata. The original reviewed
commit ID identifies a local audit snapshot; its exact generation sources are portable
in the archive above. The published branch uses a GitHub noreply identity and does not
include the unpublished local commit history. `archive_generation.jl` validates the
archived sources without requiring that historical Git object.
