# Spin-1 bilinear–biquadratic chain: finite-size survey (no certificate)

**No uniform-gap certificate.** This directory is a negative finite-size survey. It proves the
local reflection-positivity identity exactly, then scans the two seed purities on complete
numerical spectra up to 14 sites. No even seed reachable at these sizes is admissible. All
purities, defects and gaps below are floating point; none is an eigenvalue enclosure.

## Model

Spin-1 ring, one site per cell, periodic, no on-site term and no energy shift:

    H_L = Σ_j [cos θ X_{j,j+1} + sin θ X_{j,j+1}²],   X = S·S,

at θ = 0, −1/10, −3/10, −1/2, −3/5, −7/10 (radians; θ = 0 is the Heisenberg chain, θ → −π/4 the
AKLT point). For L = 2 both directed ring bonds are included, as in the core convention.

## Reflection positivity

Cartesian basis, with generators `(G_a)_{bc} = ε_{abc}` (so `S_a = −i G_a`) and symmetrized real
quadrupoles `Q_ab = −(G_a G_b + G_b G_a)/2`. With `A = Σ_a G_a ⊗ G_a` and
`B = Σ_ab Q_ab ⊗ Q_ab`, exactly

    X = −A,   X² = A/2 + B,   −h = (cos θ − sin θ/2) A − sin θ B,   o = 0.

`src/model.jl` (`exact_rp`) proves the two coefficient identities in `Rational{BigInt}`, so the
identity holds for symbolic c = cos θ, s = sin θ (no rounding of transcendental coefficients),
and independently checks the realignment (eigenvalues 2c − s ×3, −s ×5, −4s ×1). The G's span the
skew-symmetric and six Q's the symmetric operators (exact 9×9 determinant ≠ 0 in the tests),
so the bond is RP **iff** both weights are nonnegative, that is for −π/4 < θ ≤ 0. The matrices
are stored as exact rationals in `results/rp.json`.

## Result: no certificate

**Computed (floating point, dense LAPACK solves with exact integer multiplicities).**

- Complete spectra for all six angles at L = 2–8, 10, 12 (magnetization, momentum, spin flip in
  M = 0, reflection at k = 0, π) and L = 14 (SU(2) highest weight × momentum: all 120 (S, k)
  sectors, weights summing to 3¹⁴ = 4 782 969, each dimension checked against an exact
  character formula; 3.98 h summed sector time).
- Low-energy M = 0 references at L = 14, 16 (KrylovKit, six levels per block), used only for
  ground energies and gaps, never in a partition function.
- Purity grid: β ∈ {1/8, 2/8, …, 64}, n₀ = 2…7 via `CertifiedGaps.Seeds.scan_seeds`; every even
  grid point was run through the actual asymmetric recurrence, not only the 0.12731 threshold.

**Best eligible (even) seed per angle**, ranked by min(Pₛ(n₀, β₀), P_τ(2n₀, β₀)); admissibility
needs both purities above about 0.873:

| θ | n₀ | β₀ | max min(Pₛ, P_τ) at n₀ = 2 | n₀ = 4 | n₀ = 6 | bound |
|---:|---:|---:|---:|---:|---:|:---:|
| 0 | 6 | 3.25 | 0.180 | 0.262 | **0.329** | none |
| −0.1 | 6 | 2.875 | 0.175 | 0.243 | 0.292 | none |
| −0.3 | 6 | 2.375 | 0.145 | 0.228 | 0.264 | none |
| −0.5 | 6 | 2.125 | 0.158 | 0.232 | 0.249 | none |
| −0.6 | 6 | 2 | 0.134 | 0.217 | 0.256 | none |
| −0.7 | 6 | 1.875 | 0.117 | 0.222 | 0.242 | none |

The rational defects stored for these points in `results/scaling.json` (`best_even`) are upward
decimal roundings of numerical defects, **not** enclosures; the tests confirm that
`certified_gap_bound` rejects every one of them. Odd n₀ (computed through n₀ = 7 with the full
L = 14 spectrum) are diagnostics only: the criterion needs an even n₀, and odd-n spatial ratios
need not be purities.

**Why it fails, and the estimated seed scale.** The spatial defect is controlled by the
ground-energy mismatch δ(n) = E₀(2n) − 2E₀(n) > 0, which forces Pₛ ≲ exp(−β δ(n)), while P_τ needs
β of several inverse gaps. On even n = 4, 6, 8 (L = 16 from the reference):

| θ | δ(4) | δ(6) | δ(8) | gap at L = 16 | fitted ξ (exp) |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.663 | 0.365 | 0.227 | 0.443 | 3.7 |
| −0.1 | 0.833 | 0.482 | 0.321 | 0.410 | 4.2 |
| −0.3 | 1.125 | 0.681 | 0.482 | 0.401 | 4.7 |
| −0.5 | 1.345 | 0.828 | 0.599 | 0.452 | 4.9 |
| −0.6 | 1.427 | 0.881 | 0.641 | 0.493 | 5.0 |
| −0.7 | 1.491 | 0.920 | 0.671 | 0.542 | 5.0 |

A conditional forecast (fit δ(n) to A e^{−n/ξ} or A e^{−n/ξ}/√n; extrapolate the gap in 1/L or
1/L²; independent-triplet partition ansatz; first even n₀ where the recurrence returns a positive
bound) gives n₀ ≈ 20–24 at θ = 0, 24–32 at θ = −0.1 and 30–48 at θ = −0.3. For θ ≤ −0.5 the two gap
extrapolations disagree in sign and the scale is unresolved (32–42 from the positive one). These
are model-dependent extrapolations, not evidence of a seed; the full set is `extrapolations` in
`results/scaling.json`.

**Resource limit.** The next eligible seed, n₀ = 8, needs the complete L = 16 spectrum. Its
largest SU(2) block has 59 143 states (26.1 GiB for one dense real matrix) against a 16 GiB
machine; without SU(2) the largest block (M = 1, k = 1) has 310 572 states. Exact counts are in
`results/preflight.json`. Even complete L = 16 would only reach n₀ = 8, far below the forecast
scale: the route to a certificate here is a non-ED method (see `IDEAS.md` and `AGENTS.md`).

## Reproduce

From the repository root, Julia 1.12, `OPENBLAS_NUM_THREADS=1`:

```sh
# Fast check (about 70 s, 357 assertions). check.jl runs the same test suite.
julia --project=models/blbq_haldane models/blbq_haldane/check.jl
julia --project=julia models/blbq_haldane/test/runtests.jl
```

The tests prove the RP identities and small Casimir-projector identities exactly; compare the
reduced and SU(2) spectra with unreduced ED through L = 7; confirm that the stored rational seed
inputs are rejected by `certified_gap_bound`; check the SHA-256 of the stored references
recorded in `scaling.json`; compare the L = 14 rows with the independent L = 14 reference; and
recompute the n₀ = 2, 3 rows of `grid_maxima` and `energy_scaling` from scratch.

Regenerating the spectra (not needed for checking; `run.sh` wraps `julia --project=julia` with one
BLAS thread; times from the original 4-CPU run):

```sh
models/blbq_haldane/run.sh -t 4 models/blbq_haldane/scan.jl 2 3 4 5 6 7 8 10 12  # ~20 min
models/blbq_haldane/run.sh models/blbq_haldane/scripts/reference.jl 14           # ~2 min
models/blbq_haldane/run.sh models/blbq_haldane/scripts/reference.jl 16           # ~20 min
models/blbq_haldane/run.sh models/blbq_haldane/scripts/preflight.jl              # seconds
models/blbq_haldane/scripts/run_su2_grid.sh   # 120 L14 sectors: ~4 h, peak ~15 GiB RAM
models/blbq_haldane/scripts/finalize_su2.sh   # collect, gzip, analyze, audit, test
```

## What is rigorous

- **Rigorous (exact arithmetic):** the RP coefficient identities and the iff condition on θ;
  the real (zero) on-site term; integer sector and character dimension counts.
- **Trusted assumptions:** none needed, since nothing beyond the RP identity is claimed.
- **Numerical only:** every spectrum, purity, defect, ground energy, gap, the ranking of seeds,
  the δ(n) fits and the seed-scale forecasts. No eigenvalue enclosure was implemented, because no
  candidate came close to admissible. Not established: an admissible seed, a uniform gap,
  uniqueness, or any thermodynamic extrapolation error bar.

## Layout

| path | role |
|---|---|
| `scan.jl` | complete dense spectra through L = 12 (writes `results/spectrum_L*.json`) |
| `check.jl` | runs the test suite (there is no certificate to check) |
| `run.sh` | wrapper: `julia --project=julia`, one BLAS thread |
| `src/model.jl` | exact RP identities (`exact_rp`), integer S·S bond |
| `src/sectors.jl` | translation orbits, momentum/flip/reflection blocks |
| `src/su2.jl` | numerical SU(2) highest-weight blocks (Casimir projector, CholeskyQR2) |
| `src/dimensions.jl` | exact character dimensions of magnetization/momentum sectors |
| `src/result_io.jl` | reading (gzipped) spectra, weighted spectra, moments |
| `scripts/reference.jl` | L = 14, 16 low-energy M = 0 references |
| `scripts/su2_sector.jl`, `run_su2_grid.sh` | one / all L = 14 SU(2) sectors |
| `scripts/collect_su2.jl`, `archive_spectra.jl` | assemble and gzip the L = 14 spectrum |
| `scripts/analyze.jl` | rebuild `results/scaling.json` from the spectra |
| `scripts/preflight.jl`, `su2_benchmark.jl` | exact dimension counts, block timing |
| `scripts/audit.jl` | writes `results/rp.json`, `results/environment.json` |
| `scripts/finalize_su2.sh` | post-sweep pipeline |
| `test/runtests.jl`, `test/su2_small.jl` | fast tests |
| `results/scaling.json` | grid maxima, best even seeds, δ(n), gaps, forecasts, input hashes |
| `results/reference_L14.json`, `reference_L16.json` | low-energy references (hashed in scaling) |
| `results/rp.json` | exact RP matrices |
| `results/preflight.json`, `su2_benchmark_*.json` | dimension counts, resources, benchmarks |
| `results/archive.json`, `environment.json` | gzip hashes of the spectra; run environment |

## History

Developed as `certificates/blbq_haldane/` against the core's former name `ThermalMPOGap`
(commit `ef0a3e6`), then moved here. The complete spectra (`results/spectrum_L*.json[.gz]`, about
20 MB for L = 14) and the 120 L = 14 sector files (about 70 MB) are not committed; their hashes
remain in `scaling.json` and `archive.json`, and `.gitignore` excludes them when regenerated.

The tidy-up removed (all recoverable from commit `463c4a5`): `REPORT.md` (23 KB, folded into
this README) and its generator `report.jl`; `MIGRATION.md`, `migration_audit.jl` and
`results/migration_audit.json` (old-vs-renamed-core comparison); `fix_pilot_provenance.jl`
(one-off correction of the L = 14 pilot's source hash); `validate.jl` with
`results/validation.json` and `results/validation_current_main.json` (records of a test run);
and `historical/` (the first survey's JSON and source snapshot, and the Householder SU(2) pilot
source, `su2_householder.jl`, whose hash the L = 14 benchmark records). Sources moved into
`src/` and `scripts/`; the `source_sha256` records already in `results/` name the old flat paths
and are kept unchanged (no checker verifies them; the test suite verifies only the input hashes
of the stored references). The project now has its own `Manifest.toml`.
