# Transverse-field Ising model with longitudinal field

**Certified at all 30 grid points.** For every even $L \ge 6$ sites, the periodic ring has a unique ground state and spectral gap $\Delta_L \ge$ the stored bound. The smallest certified gap bound on the grid is $2749898627347949/36028797018963968 \approx 0.0763$ at $g = 1, h = 1/20$ ($n_0 = 6, \beta_0 = 13/2$).

| parameters | seed $(n_0 \text{ sites}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 30 points $g \in \{1/4, \dots, 2\} \times h \in \{1/20, \dots, 1\}$ | $(6,\; 2\text{ to }97/4\text{ per point})$ | 0.0763 ($g = 1, h = 1/20$) to 3.7050 | every even $L \ge 6$ |

## Model

In Pauli normalization on a periodic ring of $L$ spin-1/2 sites (one site per unit cell),

$$H_L = -\sum_{j=1}^L Z_j Z_{j+1} - g \sum_{j=1}^L X_j - h \sum_{j=1}^L Z_j,$$

evaluated on the 30-point coupling grid $g \in \{1/4, 1/2, 3/4, 1, 3/2, 2\} \times h \in \{1/20, 1/10, 1/4, 1/2, 1\}$. For $h > 0$, the model is non-integrable and $\mathbb{Z}_2$ spin-flip symmetry is broken. While Perron–Frobenius guarantees a unique ground state on every finite ring (since $g > 0$ connects the spin-flip hypercube with negative off-diagonal matrix elements), the two-purity doubling criterion establishes the uniform infinite-volume spectral gap and ground-state uniqueness for all even $L \ge 6$.

## Reflection positivity

In the standard $\sigma^z$ basis, the bond operator satisfies

$$-h_{\text{bond}} = Z \otimes Z \qquad (C_1 = Z,\text{ real diagonal}),$$

and the on-site term $o = -gX - hZ$ is real symmetric. The function `rp_check` verifies these identities in `Rational{BigInt}` arithmetic at every grid point.

## Seeds and certified bounds

Every grid point uses seed size $n_0 = 6$, with $\beta_0 \in \frac{1}{4}\mathbb{Z}$ selected by a numerical scan over $n_0 \in \{2, 4, 6\}$ and $\beta_0 \in \{1/4, \dots, 64\}$. The $6$- and $12$-site rings are block-diagonalized densely in every translation momentum sector $k$ (without spin-flip or reflection reduction). A polar-decomposition and Weyl perturbation bound with explicit binary64 rounding allowances ([`METHOD.md`](METHOD.md)) encloses every eigenvalue with multiplicity. Outward-rounded interval partition sums yield exact rational defect upper bounds $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(6, \beta_0)$ and $\varepsilon_\tau \ge 1 - P_\tau(12, \beta_0)$, which are passed to `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)`. All exact rational seeds and gap bounds are stored in [`results/seeds.json`](results/seeds.json) and formalized in [`lean/CertifiedGaps/Models/IsingLongitudinal.lean`](../../lean/CertifiedGaps/Models/IsingLongitudinal.lean).

Certified gap lower bound $\Delta_L \ge$ (with $\beta_0$ in parentheses), indexed by row $g$ and column $h$:

| $g \setminus h$ | 1/20 | 1/10 | 1/4 | 1/2 | 1 |
|---:|---:|---:|---:|---:|---:|
| 1/4 | 0.3923 (97/4) | 0.7332 (55/4) | 1.6606 (13/2) | 1.9757 (11/2) | 2.0336 (21/4) |
| 1/2 | 0.4564 (13) | 0.9592 (8) | 2.3616 (9/2) | 2.7754 (4) | 2.7276 (4) |
| 3/4 | 0.3170 (8) | 0.7712 (21/4) | 2.1018 (7/2) | 3.0117 (7/2) | 3.1374 (7/2) |
| 1 | 0.0763 (13/2) | 0.5043 (17/4) | 1.6092 (3) | 2.7523 (3) | 3.6403 (3) |
| 3/2 | 0.1692 (7/2) | 0.3601 (13/4) | 1.0556 (11/4) | 2.1293 (5/2) | 3.7050 (5/2) |
| 2 | 0.9121 (9/4) | 0.9643 (9/4) | 1.2728 (9/4) | 2.0005 (9/4) | 3.4912 (2) |

The certified bounds capture $9\text{–}74\%$ of the 12-site numerical gap. The symmetric diagonal threshold $\varepsilon_* \approx 0.12731$ is sufficient but not necessary: at $g = 1, h = 1/20$, the seed has $\varepsilon_{\mathrm{s}} \approx 0.184$ and $\varepsilon_\tau \approx 0.046$, which closes under the asymmetric recurrence.

## Reproduce

From this directory (`models/tfim_longitudinal`), with Julia 1.12 and
`export OPENBLAS_NUM_THREADS=1`; from the repository root, use `--project=models/tfim_longitudinal`
and prefix the paths. Times were measured on an Apple M-series laptop.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl                       # fast tests, ~15 s
julia --project=. check.jl                               # fast checker, ~15 s
julia --project=. check.jl --full                        # full replay, ~6 min
julia --project=. scripts/final_audit.jl                 # ~30 s
julia --project=. scripts/independent_residual_audit.jl  # ~10 s
julia --project=. scripts/audit.jl                       # ~20 s
julia --project=. certify.jl                             # regenerate all results, ~5 min
```

`certify.jl` accepts optional point indices `1–30` ($g$-major order) and `NMAX` (default `6`); after a full run it updates `results/seeds.json`.

## What is rigorous

- **Verified by the checker:** the exact rational RP identity; translation sector dimensions; complete eigenvalue enclosures with multiplicity; interval partition sums and exact rational seed defects; and the uniform gap bound for all even $L \ge 6$ at each of the 30 grid points. The fast checker (`check.jl`) recomputes the partition intervals, rational seeds, and recurrence bounds from the stored eigenvalue enclosures; `check.jl --full` rebuilds every momentum block from scratch, re-proves its eigenvalue enclosure of radius $\eta$, and verifies in exact rational arithmetic that each fresh interval lies inside the stored interval of radius $2\eta$.
- **Independent audits:** `scripts/final_audit.jl` independently verifies sector dimensions via Ramanujan sums, exact Pauli trace moments $\operatorname{Tr} H$ and $\operatorname{Tr}(H^2)$, 256-bit interval partition ratios, and all 30 seeds and gap bounds (1,592 checks). `scripts/independent_residual_audit.jl` checks the residual and orthogonality error bounds on all 180 $L = 6$ momentum blocks in exact rational arithmetic (360 inequalities; maximum exact-to-allowed error ratios $0.005$ and $0.017$). `scripts/audit.jl` cross-checks against full computational-basis exact diagonalization and the $h = 0$ free-fermion partition function (maximum discrepancy $7 \times 10^{-13}$).
- **Trusted base:** the IEEE binary64 / MPFR / IntervalArithmetic rounding model with gradual underflow (checked at load time) and single-threaded OpenBLAS without fast-math, under the explicit rounding-error allowances derived in [`METHOD.md`](METHOD.md).
- **Numerical only:** seed selection ranking, 12-site ED reference gaps, finite-size diagnostics, and the $h = 0$ free-fermion comparison.

## Layout

| path | role |
|---|---|
| `certify.jl` | numerical seed scan, rigorous enclosure of the best seed, writes results |
| `check.jl` | checker entry point (`--full` for the full enclosure replay) |
| `src/LongitudinalCertificate.jl` | momentum blocks, eigenvalue enclosures, seed enclosure, I/O |
| `src/checks.jl` | `check_result`, `check_grid`: the checks behind `check.jl` and the tests |
| `scripts/final_audit.jl` | independent stored-artifact audit $\to$ `results/final_audit.json` |
| `scripts/independent_residual_audit.jl` | exact $L = 6$ residual audit $\to$ `results/independent_residual_audit.json` |
| `scripts/audit.jl` | numerical cross-checks and resource counts $\to$ `results/audit.json` |
| `test/runtests.jl` | fast tests (RP, sector counts, all 30 certificates via `src/checks.jl`, seed table) |
| `METHOD.md` | eigenvalue enclosure theorem and binary64 rounding allowances |
| `IDEAS.md` | candidate routes and selection rationale |
| `results/certificate_g<g>_h<h>.json.gz` | per-point certificate: exact seeds, bound, and enclosed eigenvalues |
| `results/seeds.json` | exact $g, h, n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau$, and gap bound for all 30 points |

Each per-point certificate is gzip-compressed JSON (~260 kB uncompressed, storing the $4{,}096 + 64$ eigenvalue enclosures needed by the fast checker) with exact rational strings `"num//den"`.

## History

This model was migrated from `certificates/tfim_longitudinal/`. The per-point
`metadata.source_sha256` values are historical: they fingerprint the generating sources (under the
former core name `ThermalMPOGap` and the former file layout) and are not checked; the deep replay
re-derives every enclosure with the current code instead. In the tidy-up the module `certify.jl`
became `src/LongitudinalCertificate.jl`, `scan.jl` became `certify.jl`, the audits moved to
`scripts/`, `seeds.jl` was folded into the module, and `REPORT.md`, `INDEPENDENT_REVIEW.md`, and the
report generator `report.jl` were folded into this README. Earlier removals (uncompressed JSONs,
migration and validation records, run logs) and the removed files remain in the git history. The
three `results/*audit*.json` files are kept as produced before the tidy-up (their source hashes name
the old layout); all three audits were rerun afterward and passed with the same counts.

In the directory-format standardization the per-point files were renamed from
`results/g<g>_h<h>.json.gz` to `results/certificate_g<g>_h<h>.json.gz` (contents unchanged, so the
`certificate_sha256` values in `results/final_audit.json` still apply; the `file` fields of
`results/seeds.json` were updated), and `check.jl` gained `--full` as the name of the deep replay.
Later the checks moved from `check.jl` into `src/checks.jl`, so that `check.jl` is the entry point
and `test/runtests.jl` no longer includes it; the earlier environment-variable alias of `--full` was
dropped. Behavior and checks are unchanged. Later, duplicated helpers moved into the shared core `CertifiedGaps`: `read_point` became
`Results.read_result`, `rat` became `ExactArithmetic.parse_rational`, the Pauli matrices come from
`Spins.exact_pauli`, `certify.jl` compresses through `write_result` (same `gzip -9 -n`), and new
results record the core under `core_sha256` instead of `source_sha256`. The stored certificates
are unchanged; their recorded hashes name the earlier sources.
