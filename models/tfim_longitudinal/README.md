# Transverse-field Ising model with longitudinal field

**Certified at all 30 grid points.** For every even $L \ge 6$ sites, the periodic ring has a
unique ground state and spectral gap $\Delta_L \ge$ the interval bound `interval_bound_decimal`
(`certified_gap_bound`) of its certificate, and $\Delta_L >$ the exact rational chain bound
`gap_lower_bound`, which is at most $`0.13\%`$ smaller. Each certificate is produced end to end by
the shared core pipeline `CertifiedGaps.GapCertificates.certify_gap` (exact RP identity, verified
spectra of the 6- and 12-site rings, doubling recurrence). The smallest bound on the grid is
$\Delta_L \ge 0.07634$ (interval; chain $76277619/10^9 \approx 0.07628$) at $g = 1, h = 1/20$.

| parameters | seed $(n_0 \text{ sites}, \beta_0)$ | certified gap $\ge$ (interval bound) | lengths covered |
|---|---|---:|---|
| 30 points $`g \in \{1/4, \dots, 2\} \times h \in \{1/20, \dots, 1\}`$ | $`(6,\; 2\text{ to }97/4\text{ per point})`$ | 0.0763 ($g = 1, h = 1/20$) to 4.5202 ($g = 3/4, h = 1$) | every even $L \ge 6$ |

## Model

In Pauli normalization on a periodic ring of $L$ spin-1/2 sites (one site per unit cell, so
lengths count sites),

```math
H_L = -\sum_{j=1}^L Z_j Z_{j+1} - g \sum_{j=1}^L X_j - h \sum_{j=1}^L Z_j,
```

on the 30-point grid $`g \in \{1/4, 1/2, 3/4, 1, 3/2, 2\} \times h \in \{1/20, 1/10, 1/4, 1/2, 1\}`$.
For $h > 0$ the model is non-integrable and the $\mathbb{Z}_2$ spin-flip symmetry is broken, so no
additive charge is conserved; the core splits each ring by lattice momentum only. The model is
the core `RingModel(2, -Z⊗Z, -gX - hZ)` built by `longitudinal(g, h)` in
`src/LongitudinalCertificate.jl`.

## Reflection positivity

The certificate frame is the physical $\sigma^z$ frame: the unitary relating it to $H_L$ is the
identity. There the bond satisfies

```math
-h_{\text{bond}} = Z \otimes Z \qquad (C_1 = Z,\text{ real symmetric, weight } 1),
```

and the on-site term $o = -gX - hZ$ is real symmetric. The core `verify_rp` checks this identity
in exact rational arithmetic when a certificate is made and again in every check;
`test/runtests.jl` checks, at all 30 points and against independently written Pauli matrices, that
the stored terms equal the physical ones conjugated by the identity frame unitary, and that the
core's exact LDL† factorization also finds an RP certificate.

## Seeds and certified bounds

Every point uses $n_0 = 6$ with the $\beta_0 \in \frac14\mathbb{Z}$ of the table `BETA0` (chosen by
the floating-point scan `scripts/scan.jl`). The core encloses the complete spectra of the 6- and
12-site rings (momentum blocks $0 \le k \le L/2$; largest block 352; eigenvalue radii below
$1.4 \times 10^{-9}$), encloses $\log Z$, rounds the seed defects up to the grid $10^{-15}$, and
applies `certified_gap_bound` (interval bound) and `rational_gap_chain` (exact chain). The exact
$(n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ and bounds of all points are in
[`results/seeds.json`](results/seeds.json); every
$\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ is at most the seed defect assumed for that point in
[`lean/CertifiedGaps/Models/IsingLongitudinal.lean`](../../lean/CertifiedGaps/Models/IsingLongitudinal.lean)
(checked exactly by the tests).

Certified gap lower bounds, interval bound $\Delta_L \ge$ / exact chain bound $\Delta_L >$ (with
$\beta_0$ in parentheses):

| $g \setminus h$ | 1/20 | 1/10 | 1/4 | 1/2 | 1 |
|---:|---:|---:|---:|---:|---:|
| 1/4 | 0.5324 / 0.5319 (97/4) | 1.0743 / 1.0734 (55/4) | 2.6130 / 2.6111 (13/2) | 3.1474 / 3.1451 (11/2) | 3.2536 / 3.2512 (21/4) |
| 1/2 | 0.4572 / 0.4567 (13) | 0.9639 / 0.9631 (8) | 2.5217 / 2.5189 (9/2) | 3.6891 / 3.6860 (4) | 4.3080 / 4.3048 (4) |
| 3/4 | 0.3171 / 0.3167 (8) | 0.7714 / 0.7708 (21/4) | 2.1066 / 2.1048 (7/2) | 3.2211 / 3.2174 (7/2) | 4.5202 / 4.5166 (7/2) |
| 1 | 0.07634 / 0.07628 (13/2) | 0.5043 / 0.5039 (17/4) | 1.6097 / 1.6076 (3) | 2.7693 / 2.7672 (3) | 4.2864 / 4.2822 (3) |
| 3/2 | 0.1692 / 0.1690 (7/2) | 0.3601 / 0.3596 (13/4) | 1.0557 / 1.0546 (11/4) | 2.1303 / 2.1278 (5/2) | 3.7706 / 3.7655 (5/2) |
| 2 | 0.9121 / 0.9114 (9/4) | 0.9643 / 0.9636 (9/4) | 1.2729 / 1.2715 (9/4) | 2.0010 / 1.9996 (9/4) | 3.4989 / 3.4957 (2) |

The interval bounds capture $`9\text{–}82\%`$ of the 12-site numerical gap and are at least the
previous certificates' bounds at every point. The symmetric threshold
$\varepsilon_* \approx 0.12731$ is sufficient but not necessary: at $g = 1, h = 1/20$ the seed has
$\varepsilon_{\mathrm{s}} \approx 0.184$ and $\varepsilon_\tau \approx 0.0074$, which closes under
the asymmetric recurrence.

## Reproduce

From the repository root, in the core environment `julia/` (this model has no environment of
its own), with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1` (times on an Apple M-series
laptop):

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia models/tfim_longitudinal/certify.jl         # ~75 s
julia --project=julia models/tfim_longitudinal/check.jl           # ~16 s
julia --project=julia models/tfim_longitudinal/check.jl --full    # ~42 s
julia --project=julia models/tfim_longitudinal/test/runtests.jl   # ~30 s
julia --project=julia models/tfim_longitudinal/scripts/scan.jl    # ~27 s
```

`certify.jl` accepts optional point indices `1–30` ($g$-major order); after a run it rewrites
`results/seeds.json` from all stored certificates.

## What is rigorous

- **Checked by `check.jl`:** per certificate, the recorded model-local source hashes (enforced),
  the stored model, RP factor and seed against the grid definition, the exact RP identity and the
  exact rational chain reproducing the stored bound (core `check_certificate(cert; spectra =
  false)`), and `results/seeds.json` against the certificates. `--full` also re-encloses both
  seed spectra with the current core and requires the fresh defect enclosures not to exceed the
  stored $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$, the only spectral inputs of the chain.
- **Trusted base:** exact rational arithmetic (RP identity, chain); `IntervalArithmetic` for block
  entries, $\log Z$ and the defects; and the eigenvalue enclosure
  `CertifiedGaps.VerifiedSpectra.verified_eigen` under the IEEE binary64 + BLAS rounding model
  (round to nearest, $|\mathrm{fl}(\sum a_i b_i) - \sum a_i b_i| \le \gamma_K \sum |a_i b_i|$ for
  any summation order and FMA, with $K = N + 8$ for the real blocks $k = 0, L/2$ and $K = 2N + 8$
  for complex blocks, which use four real GEMMs; Frobenius norms; an underflow allowance), then
  Ostrowski or polar-factor and Weyl bounds. Gradual underflow is checked at load time; BLAS
  runs single-threaded.
- **Floating-point only:** the seed scan `scripts/scan.jl` and `results/scan.json` (it matches
  `BETA0` at 29 of 30 points; at $g = 1/4, h = 1$ it prefers $\beta_0 = 23/4$ because its `Float64`
  defects there are at rounding level, about $10^{-16}$, so its float bound is meaningless).

## Layout

| path | role |
|---|---|
| `certify.jl` | loops `certify_gap` over the grid, writes the certificates and `results/seeds.json` |
| `check.jl` | checker entry point (`--full` re-encloses the seed spectra) |
| `src/LongitudinalCertificate.jl` | model, RP factor, grid and seeds, provenance, per-point checks |
| `scripts/scan.jl` | floating-point seed scan $\to$ `results/scan.json` |
| `test/runtests.jl` | frame and RP, spectra vs ED, fast checks, tampering, Lean seed bounds |
| `METHOD.md` | pipeline, checker semantics and the old-versus-new comparison |
| `IDEAS.md` | candidate routes and selection rationale |
| `results/certificate_g<g>_h<h>.json.gz` | core-format certificate of one grid point (gzip JSON) |
| `results/seeds.json` | exact $g, h, n_0, \beta_0, \varepsilon_{\mathrm{s}}, \varepsilon_\tau$ and bounds of all points |
| `results/scan.json` | floating-point scan diagnostics |

No file is hash-pinned by a certificate except the model-local sources listed in each
certificate's `metadata.source_sha256` (`src/`, `certify.jl`, `check.jl`); the core sources and
environment (`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`) are recorded under
`core_sha256` as provenance only.

## History

Migrated from `certificates/tfim_longitudinal/` and later tidied to the standard layout (the
earlier audits, reports and migration records are in the git history). The certificates were
then regenerated with the shared core pipeline: the local momentum-block builder, polar/Weyl
enclosure with its own rounding constants (stored radius $2\eta$), local $\log Z$, checks and
audit scripts were replaced by `RingModel` + `certify_gap` + `check_certificate`, at the same
seeds. Every $\varepsilon_{\mathrm{s}}, \varepsilon_\tau$ decreased (e.g. $g = 1/4, h = 1/20$:
$\varepsilon_{\mathrm{s}}$ $9.06 \times 10^{-6} \to 1.21 \times 10^{-6}$), and the interval
bound grew at every point (smallest bound $0.07633 \to 0.07634$; $0.3923 \to 0.5324$ at
$g = 1/4, h = 1/20$; $2.7276 \to 4.3080$ at $g = 1/2, h = 1$). The exact chain bound is up to
$`0.13\%`$ below the old interval bound at 13 points (e.g. $0.07628$ at $g = 1, h = 1/20$).
See `METHOD.md`.
The model-local `Project.toml`/`Manifest.toml` were then removed: everything runs in the core
environment (`--project=julia`), and the certificates were regenerated with identical seeds.
