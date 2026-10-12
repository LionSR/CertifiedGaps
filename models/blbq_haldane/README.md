# BLBQ Haldane finite-size survey (negative result)

Model: spin-1 bilinear-biquadratic ring
`H_L = sum_j [cos(theta) S_j.S_(j+1) + sin(theta) (S_j.S_(j+1))^2]`, one site per cell,
periodic, `theta = 0, -1/10, -3/10, -1/2, -3/5, -7/10`.

**No uniform-gap certificate is claimed.** This directory is a numerical negative grid
survey, committed as a scan with clearly labelled floating-point evidence. See
[REPORT.md](REPORT.md) for the full tables, exact RP decomposition, resource limits,
conditional seed-scale forecasts and the rigor boundary, and [IDEAS.md](IDEAS.md) for routes.

## What is and is not rigorous

- **Rigorous:** local Cartesian reflection positivity. `model.jl` (`exact_rp`) verifies
  `X = -A` and `X^2 = A/2 + B` in `Rational{BigInt}` arithmetic, where `A` and `B` are built
  from the three real antisymmetric generators (coefficient `cos(theta) - sin(theta)/2`) and
  the symmetrized quadrupoles (coefficient `-sin(theta)`), plus an independent realignment
  check; `o = 0`. Valid symbolically for `-pi/4 < theta <= 0`; stored in `results/rp.json`.
- **Floating point only:** every purity, defect, energy mismatch and gap. Spectra are dense
  LAPACK solves with exact integer multiplicities, not eigenvalue enclosures.

## Findings (floating point)

- Complete spectra through L = 14 for all six angles (L14 via SU(2) highest-weight blocks:
  120 (S,k) sectors whose weights sum to 3^14); low-energy M = 0 references at L = 14, 16.
- No eligible even seed passes. Closest: theta = 0, n0 = 6, beta0 = 3.25, balanced purity
  min(Ps, Ptau) = 0.3288 (the sufficient condition needs both purities above 0.8727).
  L14 adds only the odd n0 = 7 diagnostic; the next eligible seed n0 = 8 needs complete
  L16 spectra, beyond this implementation's memory budget (26 GiB per largest block).
- The mismatch delta(n) = E0(2n) - 2E0(n) decays roughly like exp(-n/xi) on n = 4, 6, 8,
  with fitted xi between 3.7 (theta = 0) and 5.0 (theta = -0.7).
- Conditional forecasts of the first admissible even seed: n0 about 20-24 at theta = 0,
  24-32 at theta = -0.1, 30-48 at theta = -0.3; for theta <= -0.5 the two gap
  extrapolations disagree in sign and the required scale is unresolved. These are
  model-dependent extrapolations, not evidence of a seed.

## Kept files

Code (`*.jl`, `*.sh`), `REPORT.md`, `IDEAS.md`, `MIGRATION.md`, `historical/` (small
pre-rename evidence and source snapshots), and the summary results:

- `results/scaling.json`: best even seed per theta (with upward-rounded rational tokens that
  the recurrence rejects), `grid_maxima` (max over beta of min(Ps, Ptau) per n0 and theta),
  `energy_scaling` (delta(n), ground energies, gaps), `extrapolations` (seed-scale
  forecasts), and the SHA-256 of every spectrum it was computed from.
- `results/rp.json`, `reference_L14.json`, `reference_L16.json`, `preflight.json`,
  `archive.json`, `environment.json`, `migration_audit.json`, `validation*.json`, and the two
  small SU(2) benchmark files.

## Removed to keep the repository small

The directory was about 100 MB; it is now about 0.6 MB. Removed (recoverable from branch
history before the removal commit; their hashes remain in `scaling.json` and `archive.json`):

- complete spectra `results/spectrum_L{2..8,10,12,14}.json[.gz]` (L14 alone 19.6 MB gzip);
- the 120 per-sector L14 checkpoints `results/su2_L14_S*_k*.json` (about 70 MB);
- progress and PR notes (`review_notes.md`, `PR.md`, `PR_BODY.md`).

`analyze.jl` and `report.jl` need the spectra. Regenerate them from the repository root
(`run.sh` runs `julia --project=julia` with one BLAS thread):

```sh
models/blbq_haldane/run.sh -t 4 models/blbq_haldane/scan.jl 2 3 4 5 6 7 8 10 12  # ~20 min
models/blbq_haldane/run.sh models/blbq_haldane/reference.jl 14                   # ~2 min
models/blbq_haldane/run.sh models/blbq_haldane/reference.jl 16                   # ~20 min
models/blbq_haldane/run_su2_grid.sh   # 120 L14 sectors: ~4 h total, peak ~15 GiB RAM
models/blbq_haldane/finalize_su2.sh   # collect, archive, analyze, validate, report, audit
```

## Fast test

```sh
OPENBLAS_NUM_THREADS=1 julia --project=julia models/blbq_haldane/test/runtests.jl
```

About 2 minutes, 357 assertions (`check.jl` runs the same file). It checks the exact RP
identities and small Casimir-projector identities; SU(2) and symmetry-reduced spectra
against unreduced ED through L = 7; that the stored rational seed inputs are rejected by
`certified_gap_bound`; the hashes of the kept reference files; agreement of the L14 summary
rows with the independent L14 reference; and it recomputes the n0 = 2, 3 rows of
`grid_maxima` and `energy_scaling` from scratch for all six angles.
