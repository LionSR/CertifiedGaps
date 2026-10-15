# Spin-1 bilinear–biquadratic chain

**No uniform-gap certificate ($n_0 \le 6$).** This module documents a finite-size scaling survey of the spin-1 bilinear–biquadratic (BLBQ) chain in the Haldane phase. It proves the exact reflection-positivity identity for all $-\pi/4 < \theta \le 0$ and evaluates the spatial and thermal purities on complete numerical spectra up to $L = 14$ sites. Because the correlation length is $\xi \approx 3.7\text{–}5.0$ sites, no even seed reachable by dense exact diagonalization ($n_0 \le 6$, requiring up to $12$ sites) comes close to the admissible region ($\min(P_{\mathrm{s}}, P_\tau) \gtrsim 0.873$). All purities, defects, and gaps reported below are floating-point diagnostics rather than certified enclosures.

| $\theta$ | best even seed $(n_0, \beta_0)$ | $\max \min(P_{\mathrm{s}}, P_\tau)$ | certified gap $\ge$ | lengths covered |
|---:|---|---:|---|---|
| 0 | $(6, 3.25)$ | 0.329 | no certificate | — |
| $-0.1$ | $(6, 2.875)$ | 0.292 | no certificate | — |
| $-0.3$ | $(6, 2.375)$ | 0.264 | no certificate | — |
| $-0.5$ | $(6, 2.125)$ | 0.249 | no certificate | — |
| $-0.6$ | $(6, 2)$ | 0.256 | no certificate | — |
| $-0.7$ | $(6, 1.875)$ | 0.242 | no certificate | — |

## Model

On a periodic ring of $L$ spin-1 sites (one site per unit cell, $o = 0$, $\kappa = 0$),

$$H_L = \sum_{j=1}^L \bigl(\cos\theta\, X_{j,j+1} + \sin\theta\, X_{j,j+1}^2\bigr), \qquad X_{j,j+1} = \mathbf{S}_j \cdot \mathbf{S}_{j+1},$$

evaluated at $\theta \in \{0, -1/10, -3/10, -1/2, -3/5, -7/10\}$ radians ($\theta = 0$ is the spin-1 Heisenberg chain, certified at $n_0 = 20$ in [`models/spin1/`](../spin1/README.md), and $\theta = -\pi/4$ is the gapless Takhtajan–Babujian critical point bounding the negative-$\theta$ Haldane window). At $L = 2$, both periodic bonds are included.

## Reflection positivity

In the Cartesian spin-1 basis with real skew-symmetric generators $(G_\alpha)_{bc} = \varepsilon_{\alpha bc}$ ($S^\alpha = -i G_\alpha$) and real symmetric quadrupole operators $Q_{\alpha\beta} = -\frac{1}{2}(G_\alpha G_\beta + G_\beta G_\alpha)$, define

$$A = \sum_{\alpha=1}^3 G_\alpha \otimes G_\alpha, \qquad B = \sum_{\alpha,\beta=1}^3 Q_{\alpha\beta} \otimes Q_{\alpha\beta}.$$

Then $X = -A$ and $X^2 = \frac{1}{2}A + B$, so the unshifted bond operator decomposes as

$$-h = \bigl(\cos\theta - \tfrac{1}{2}\sin\theta\bigr) A - \sin\theta\, B, \qquad o = 0.$$

The function `exact_rp` in [`src/model.jl`](src/model.jl) proves these algebraic identities over `Rational{BigInt}`, establishing that the identity holds symbolically for all $\cos\theta$ and $\sin\theta$, and verifies that the realignment $R[-h]$ has eigenvalues $2\cos\theta - \sin\theta$ ($\times 3$), $-\sin\theta$ ($\times 5$), and $-4\sin\theta$ ($\times 1$). Because the three $G_\alpha$ and six independent $Q_{\alpha\beta}$ form a linearly independent basis of $3 \times 3$ matrices, the unshifted bond is reflection positive if and only if both weights are non-negative, i.e., $-\pi + \arctan(2) \le \theta \le 0$ (which intersects the Haldane phase on $-\pi/4 < \theta \le 0$). The exact rational basis matrices are recorded in `results/rp.json`.

## Result

### Computed numerical spectra

- **Complete dense spectra ($L = 2\text{–}8, 10, 12, 14$):** computed at all six angles using magnetization, translation momentum, spin flip (in $M = 0$), and reflection (at $k = 0, \pi$) up to $L = 12$, and using $\mathrm{SU}(2)$ highest-weight $\times$ momentum blocks at $L = 14$ (all 120 $(S, k)$ sectors with multiplet weights $2S + 1$ summing to $3^{14} = 4{,}782{,}969$, each sector dimension verified against an exact character formula; 3.98 hours total dense solver time).
- **Low-energy reference levels ($L = 14, 16$):** computed in the $M = 0$ sector using `KrylovKit` (6 levels per block) to extract ground-state energies $E_0(L)$ and finite-size gaps $\Delta_L$ (never used in partition functions).
- **Seed purity grid:** scanned over $\beta_0 \in \{1/8, 2/8, \dots, 64\}$ and $n_0 = 2, \dots, 7$ via `CertifiedGaps.Seeds.scan_seeds`, testing every even grid point against the full asymmetric doubling recurrence.

### Best even seed per angle ($n_0 \le 6$)

| $\theta$ | $n_0$ | $\beta_0$ | $\max \min(P_{\mathrm{s}}, P_\tau)$ at $n_0 = 2$ | $n_0 = 4$ | $n_0 = 6$ | certified bound |
|---:|---:|---:|---:|---:|---:|:---:|
| 0 | 6 | 3.25 | 0.180 | 0.262 | **0.329** | none |
| $-0.1$ | 6 | 2.875 | 0.175 | 0.243 | 0.292 | none |
| $-0.3$ | 6 | 2.375 | 0.145 | 0.228 | 0.264 | none |
| $-0.5$ | 6 | 2.125 | 0.158 | 0.232 | 0.249 | none |
| $-0.6$ | 6 | 2 | 0.134 | 0.217 | 0.256 | none |
| $-0.7$ | 6 | 1.875 | 0.117 | 0.222 | 0.242 | none |

The rational defects stored for these points in `results/scaling.json` (`best_even`) are upward decimal roundings of floating-point defects rather than certified enclosures, and `certified_gap_bound` rejects every one of them. Odd $n_0$ values (computed up to $n_0 = 7$ from the full $L = 14$ spectrum) are recorded only as diagnostics.

### Why small seeds fail and estimated seed scale

With zero energy shift ($\kappa = 0$), the spatial purity is bounded above by the ground-energy doubling mismatch $\delta(n) = E_0(2n) - 2E_0(n) > 0$, which forces $P_{\mathrm{s}}(n, \beta) \lesssim e^{-\beta \delta(n)}$, while the thermal purity $P_\tau(2n, \beta)$ requires $\beta$ to be several inverse gaps. On even rings $n = 4, 6, 8$:

| $\theta$ | $\delta(4)$ | $\delta(6)$ | $\delta(8)$ | gap at $L = 16$ | fitted $\xi$ (exp) |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.663 | 0.365 | 0.227 | 0.443 | 3.7 |
| $-0.1$ | 0.833 | 0.482 | 0.321 | 0.410 | 4.2 |
| $-0.3$ | 1.125 | 0.681 | 0.482 | 0.401 | 4.7 |
| $-0.5$ | 1.345 | 0.828 | 0.599 | 0.452 | 4.9 |
| $-0.6$ | 1.427 | 0.881 | 0.641 | 0.493 | 5.0 |
| $-0.7$ | 1.491 | 0.920 | 0.671 | 0.542 | 5.0 |

Fitting $\delta(n)$ to an exponential decay and modeling the thermal purity with a dilute triplet gas estimates the first admissible seed at $n_0 \approx 20\text{–}24$ for $\theta = 0$ (consistent with the $n_0 = 20$ certificate in `models/spin1/`), $n_0 \approx 24\text{–}32$ for $\theta = -0.1$, and $n_0 \approx 30\text{–}48$ for $\theta = -0.3$ (`extrapolations` in `results/scaling.json`). For $\theta \le -0.5$, the $1/L$ and $1/L^2$ gap extrapolations disagree in sign ($n_0 \approx 32\text{–}42$ under the positive extrapolation).

### Resource limit for dense exact diagonalization

The next even seed size, $n_0 = 8$, requires the complete $L = 16$ spectrum ($3^{16} = 43{,}046{,}721$ states). Its largest $\mathrm{SU}(2) \times \text{momentum}$ block has dimension $59{,}143$ (requiring $26.1\text{ GiB}$ for a single dense `Float64` matrix), exceeding a 16 GiB workstation (`results/preflight.json`). Reaching the unshifted seed scale $n_0 \ge 20$ requires polynomial filtering and MPS trial bounds (as in `models/spin1/`) or tuned bond shifts $\kappa$.

## Reproduce

From the repository root, with Julia 1.12 and `OPENBLAS_NUM_THREADS=1`:

```sh
# Fast check (~19 s, 387 assertions): runs the test suite.
julia --project=models/blbq_haldane models/blbq_haldane/check.jl
# Full check (~21 s): also recomputes the n₀ = 4 rows from complete 4- and 8-site spectra.
julia --project=models/blbq_haldane models/blbq_haldane/check.jl --full
julia --project=julia models/blbq_haldane/test/runtests.jl
```

The tests prove the RP identities and small Casimir-projector identities exactly; compare the
reduced and $\mathrm{SU}(2)$ spectra with unreduced exact diagonalization through $L = 7$; confirm
that the stored rational seed inputs are rejected by `certified_gap_bound`; check the SHA-256 hashes
of the stored references recorded in `scaling.json`; compare the $L = 14$ rows with the independent
$L = 14$ reference; and recompute the $n_0 = 2, 3$ rows of `grid_maxima` and `energy_scaling` from
scratch (`--full` also recomputes $n_0 = 4$). There is no certificate, so there is no `certify.jl`;
the regeneration drivers below live in `scripts/`. Times measured on an Apple M-series laptop
(2026-10-14): `check.jl` 19 s, `check.jl --full` 21 s.

To regenerate the exploratory scan data from scratch (not required by `check.jl`; `scripts/run.jl`
runs `julia --project=julia` from the repository root with one BLAS thread; times are from the
original 4-CPU run):

```sh
julia models/blbq_haldane/scripts/run.jl -t 4 \
    models/blbq_haldane/scripts/scan.jl 2 3 4 5 6 7 8 10 12      # ~20 min
julia models/blbq_haldane/scripts/run.jl \
    models/blbq_haldane/scripts/reference.jl 14                    # ~2 min
julia models/blbq_haldane/scripts/run.jl \
    models/blbq_haldane/scripts/reference.jl 16                    # ~20 min
julia models/blbq_haldane/scripts/run.jl \
    models/blbq_haldane/scripts/preflight.jl                       # seconds
julia models/blbq_haldane/scripts/run_su2_grid.jl  # 120 L=14 SU(2) sectors: ~4 h, peak ~15 GiB RAM
julia models/blbq_haldane/scripts/finalize_su2.jl  # collect, gzip, analyze, audit, test
```

## What is rigorous

- **Exact arithmetic:** the reflection-positivity decomposition over `Rational{BigInt}`, the realignment eigenvalue formulas showing RP holds if and only if $-\pi + \arctan(2) \le \theta \le 0$ (covering the Haldane window $-\pi/4 < \theta \le 0$), and the exact character dimensions of all magnetization, momentum, and $\mathrm{SU}(2)$ sectors.
- **Numerical only:** all computed eigenvalues, purities, defects, ground-state energies, finite-size gaps, and seed-scale extrapolations.

## Layout

| path | role |
|---|---|
| `check.jl` | runs the verification and unit test suite (no certificate to check); `--full` adds the $n_0 = 4$ rows |
| `src/model.jl` | exact RP identities (`exact_rp`) and integer $\mathbf{S}_j \cdot \mathbf{S}_{j+1}$ bond matrices |
| `src/sectors.jl` | translation orbits and momentum / spin-flip / reflection blocks |
| `src/su2.jl` | numerical $\mathrm{SU}(2)$ highest-weight blocks (Casimir projector + CholeskyQR2) |
| `src/dimensions.jl` | exact character dimensions of magnetization and momentum sectors |
| `src/result_io.jl` | I/O helpers for gzipped spectra, multiplet weights, and trace moments |
| `scripts/scan.jl` | complete dense spectra through $L = 12$ (writes `results/spectrum_L*.json`) |
| `scripts/run.jl` | launcher: `julia --project=julia` from the repository root, one BLAS thread |
| `scripts/reference.jl` | $L = 14, 16$ low-energy $M = 0$ Krylov references |
| `scripts/su2_sector.jl`, `run_su2_grid.jl` | one / all $L = 14$ $\mathrm{SU}(2)$ sectors |
| `scripts/collect_su2.jl`, `archive_spectra.jl` | assemble and compress the $L = 14$ spectrum |
| `scripts/analyze.jl` | rebuild `results/scaling.json` from stored spectra |
| `scripts/preflight.jl`, `su2_benchmark.jl` | exact sector dimension counts and block timing benchmarks |
| `scripts/audit.jl` | writes `results/rp.json` and `results/environment.json` |
| `scripts/finalize_su2.jl` | post-sweep pipeline |
| `test/runtests.jl`, `test/su2_small.jl` | fast unit and consistency tests |
| `results/scaling.json` | grid maxima, best even seeds, $\delta(n)$, gaps, forecasts, and input hashes |
| `results/reference_L14.json`, `reference_L16.json` | low-energy Krylov reference levels (hashed in `scaling.json`) |
| `results/rp.json` | exact rational RP basis matrices |
| `results/preflight.json`, `su2_benchmark_*.json` | exact sector dimension counts and resource benchmarks |
| `results/archive.json`, `environment.json` | SHA-256 hashes of compressed spectra and run metadata |

Hash pinning: no checker compares a source hash by path (the tests verify only the input hashes of
the stored references in `scaling.json`), so there are no layout exceptions.

## History

This survey was originally developed under `certificates/blbq_haldane/` against the core's former
name `ThermalMPOGap` (commit `ef0a3e6`), then migrated to `models/blbq_haldane/`. The complete
spectra (`results/spectrum_L*.json[.gz]`, about 20 MB for $L = 14$) and the 120 $L = 14$ sector
files (about 70 MB) are not committed; their SHA-256 hashes remain in `results/scaling.json` and
`results/archive.json`, and `.gitignore` excludes them when regenerated.

The tidy-up removed (all recoverable from commit `463c4a5`): `REPORT.md` (23 KB, folded into this
README) and its generator `report.jl`; `MIGRATION.md`, `migration_audit.jl`, and
`results/migration_audit.json` (old-vs-renamed-core comparison); `fix_pilot_provenance.jl` (one-off
correction of the $L = 14$ pilot's source hash); `validate.jl` with `results/validation.json` and
`results/validation_current_main.json` (records of a test run); and `historical/` (the first
survey's JSON and source snapshot, and the Householder $\mathrm{SU}(2)$ pilot source,
`su2_householder.jl`, whose hash the $L = 14$ benchmark records). Sources moved into `src/` and
`scripts/`; the `source_sha256` records already in `results/` name the old flat paths and are kept
unchanged (no checker verifies them; the test suite verifies only the input hashes of the stored
references). The project now has its own `Manifest.toml`.

The directory-format tidy moved `scan.jl` and `run.sh` from the top level into `scripts/` (adjusting
the paths inside them) and added `check.jl --full`. The `"scan.jl"` and `"run.sh"` entries in
`results/scaling.json` and `results/environment.json` record earlier versions of those files at
their former top-level paths (see the git history); later regenerations hash them under `scripts/`.
The shell drivers `run.sh`, `run_su2_grid.sh`, and `finalize_su2.sh` were later ported to the Julia
scripts `run.jl`, `run_su2_grid.jl`, and `finalize_su2.jl` (same steps, same timeout); the `run.sh`
hash in `results/environment.json` names the old wrapper.
