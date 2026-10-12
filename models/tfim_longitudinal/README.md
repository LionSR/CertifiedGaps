# Longitudinal-field TFIM certificates

Model: H_L = −Σ Z_j Z_{j+1} − g Σ X_j − h Σ Z_j (Pauli matrices, periodic ring, one site per
cell), on the grid g ∈ {1/4, 1/2, 3/4, 1, 3/2, 2} × h ∈ {1/20, 1/10, 1/4, 1/2, 1}.

**All 30 points are certified**: for every even L ≥ 6 the ground state is unique and the gap
is at least the stored bound. The smallest bound is 2749898627347949/36028797018963968
≈ 0.0763 (g = 1, h = 1/20). See [REPORT.md](REPORT.md) for the table and region map,
[METHOD.md](METHOD.md) for the enclosure theorem, [IDEAS.md](IDEAS.md) for the route
choice and [INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md) for the independent review.

- **RP basis**: σᶻ basis, −h_bond = Z ⊗ Z (one term, C = Z, real), on-site
  o = −gX − hZ real symmetric. `rp_check` verifies this in `Rational{BigInt}`.
- **Seed**: n₀ = 6 at every point, β₀ ∈ (1/4)ℤ chosen per point by a numerical scan
  (n₀ ∈ {2, 4, 6}, β₀ ∈ {1/4, …, 64}). Rings of 6 and 12 sites are diagonalised densely in
  every momentum sector (no spin-flip or reflection-parity reduction); a polar-factor/Weyl
  argument with explicit rounding allowances encloses every eigenvalue with multiplicity;
  interval partition sums give exact rational εₛ ≥ 1 − Pₛ(6, β₀), εₜ ≥ 1 − P_τ(12, β₀).
- **Bound**: `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)`; exact binary value stored.

Rigorous: RP identity, sector dimensions, eigenvalue enclosures (under the IEEE binary64 /
MPFR / IntervalArithmetic model with gradual underflow and single-threaded BLAS, see
METHOD.md), partition-function and seed enclosures, the resulting uniform gap on the
discrete grid. Numerical only: seed ranking, 12-site ED reference gaps, bound/ED ratios,
finite-size diagnostics. Nothing is claimed between grid points, for odd L, or as h → 0.

## Files

| file | role |
|---|---|
| `certify.jl` | module: momentum blocks, rigorous enclosures, seed enclosure, I/O |
| `scan.jl` | numerical seed scan + rigorous best seed; writes `results/g*_h*.json.gz` |
| `check.jl` | checker; `DEEP_CHECK=1` also regenerates and contains all enclosures |
| `test/runtests.jl` | fast tests (RP, sector counts, partition intervals, recurrence, seed table) |
| `seeds.jl` | writes `results/seeds.json`, the exact seed table of all 30 points |
| `audit.jl`, `final_audit.jl`, `independent_residual_audit.jl` | independent cross-checks |
| `report.jl` | regenerates the body of `REPORT.md` |
| `results/g<g>_h<h>.json.gz` | per-point certificate: exact seeds, bound, full enclosed spectra |
| `results/seeds.json` | g, h, n₀, β₀, εₛ, εₜ, bound for every point (exact rationals) |
| `results/*audit*.json` | outputs of the audit scripts |

The per-point certificates are gzip-compressed JSON (each ~260 kB uncompressed, almost all
of it the 4096 + 64 stored eigenvalue midpoints and radii, which the checker needs to
recompute the partition intervals). `C.read_point` decompresses with the system `gzip`;
`gzip -dc results/g1_1_h1_20.json.gz` shows one. Exact numbers are `"num//den"` strings.

## Commands

From the repository root (Julia 1.12):

```sh
export OPENBLAS_NUM_THREADS=1
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia models/tfim_longitudinal/test/runtests.jl          # ~20 s
DEEP_CHECK=1 julia --project=julia models/tfim_longitudinal/check.jl     # ~6 min CPU, full replay
julia --project=julia models/tfim_longitudinal/final_audit.jl            # ~30 s
julia --project=julia models/tfim_longitudinal/independent_residual_audit.jl
julia --project=julia models/tfim_longitudinal/audit.jl
julia --project=julia models/tfim_longitudinal/scan.jl                   # ~5 min, rewrites results
julia --project=julia models/tfim_longitudinal/seeds.jl
```

The directory's own `Project.toml` points at `../../julia`. Measured on the current core:
fast tests 74/74 in about 20 s, deep replay 30/30 (330 s CPU), final artifact audit 1592/1592
in about 30 s.

`scan.jl` takes `NMAX` (default 6) and optional parameter indices 1–30 (g-major). Larger rings
are expensive: a single generic momentum block at L = 20 is 52,377-dimensional (~41 GiB as a
complex dense matrix), so the dense route stops well short of L = 20.

## Removed when moving to `models/`

The draft PR lived in `certificates/tfim_longitudinal/` and carried migration and provenance
records that are no longer needed: uncompressed per-point JSONs (now gzipped, content
unchanged), `results/original_generation_sources.json` and `archive_generation.jl` (archive of
the pre-rename `ThermalMPOGap` sources), `results/namespace_migration.json`,
`results/validation.json` and `results/validation_original.json` with `verify_all.jl`
(replaced by `DEEP_CHECK=1 check.jl`), `results/audit_manifest.json` with `seal.jl`, and run
logs. All remain in the git history of this branch. The per-point `metadata.source_sha256`
values are historical (they fingerprint the generating sources at run time); rerun `scan.jl`
to obtain fingerprints of the current sources.
