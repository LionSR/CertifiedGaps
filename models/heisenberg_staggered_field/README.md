# Spin-1/2 Heisenberg chain in a staggered field

**Certified at $h = 1, 3/4, 1/2, 7/20, 1/4$: a unique ground state and uniform spectral gap
$\Delta_L \ge 1.148261, 0.890782, 0.594782, 0.383508, 0.214486$ for every even chain length
$L \ge 8$ sites.** The seeds are certified from verified dense spectra of complete 8- and 16-site
rings decomposed into magnetization and two-site translation momentum blocks. Full mathematical
derivations, the eigenvalue-enclosure theorem, sector dimension tables, 20-site exploratory scans,
and the independent audit are documented in [`METHOD.md`](METHOD.md); the route memo is in
[`IDEAS.md`](IDEAS.md).

| Staggered field $h$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $1$ | $(8, 21/4)$ | $1.148261$ | every even $L \ge 8$ sites |
| $3/4$ | $(8, 21/4)$ | $0.890782$ | every even $L \ge 8$ sites |
| $1/2$ | $(8, 11/2)$ | $0.594785$ | every even $L \ge 8$ sites |
| $7/20$ | $(8, 23/4)$ | $0.383508$ | every even $L \ge 8$ sites |
| $1/4$ | $(8, 25/4)$ | $0.214486$ | every even $L \ge 8$ sites |

## Model

We consider the periodic spin-1/2 antiferromagnetic Heisenberg chain in a staggered longitudinal
magnetic field $h$:
```math
H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} - h \sum_{j=1}^L (-1)^j S_j^z, \qquad \mathbf{S}_{L+1} \equiv \mathbf{S}_1,
```
with unit exchange coupling ($J = 1$) and $\mathbf{S}_j = \frac{1}{2}\boldsymbol{\sigma}_j$. On a
two-site ring ($L = 2$), both periodic bonds $(1, 2)$ and $(2, 1)$ are included, matching the core
`CertifiedGaps` convention. Although the staggered field $(-1)^j S_j^z$ has a two-site spatial
period in the physical basis, applying a $\pi$-rotation about the $y$-axis on odd sites makes the
Hamiltonian translation-invariant with a **single-spin unit cell** on any even ring $L$. Lengths
$L$ and $n_0$ therefore count individual spins (sites), and the sublattice rotation closes
periodically on all even chain lengths $L$.

## Reflection positivity

We apply the real orthogonal single-site unitary $`U = 2i S^y = \begin{psmallmatrix} 0 & 1 \\ -1 & 0 \end{psmallmatrix}`$
(a $\pi$-rotation about the $y$-axis) to every odd site $j$. In terms of the exact rational
matrices $X = S^x$, $G_y = i S^y$ (denoted $J$ in `src/common.jl`, which is real and
antisymmetric), and $Z = S^z$, conjugating either bond orientation $(j, j+1)$ by $U$ on its odd
endpoint negates $S^x$ and $S^z$ while fixing $iS^y$, transforming every nearest-neighbor bond into
```math
h_{\mathrm{rot}} = -\left(X \otimes X + G_y \otimes G_y + Z \otimes Z\right), \qquad -h_{\mathrm{rot}} = \sum_{\alpha \in \{x,y,z\}} C_\alpha \otimes \overline{C_\alpha}
```
with real bond factors $`C_\alpha \in \{X, G_y, Z\}`$. Simultaneously, $U Z U^\top = -Z$ converts the
staggered field $(-1)^j S_j^z$ into a uniform, real diagonal on-site term $o = -h Z$. The function
`exact_rp` in `src/common.jl` verifies these identities in exact rational arithmetic for both bond
orientations.

## Seeds and certified bounds

| $h$ | Seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(8, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(16, \beta_0)$ | Certified gap $\Delta_L \ge$ | 16-site ED gap |
|---|---|---:|---:|---:|---:|
| $1$ | $(8, 21/4)$ | $\approx 6.937 \times 10^{-4}$ | $\approx 7.326 \times 10^{-4}$ | $1.148261$ | $1.87351193$ |
| $3/4$ | $(8, 21/4)$ | $\approx 2.343 \times 10^{-3}$ | $\approx 3.122 \times 10^{-3}$ | $0.890782$ | $1.57089206$ |
| $1/2$ | $(8, 11/2)$ | $\approx 1.074 \times 10^{-2}$ | $\approx 1.116 \times 10^{-2}$ | $0.594785$ | $1.23131976$ |
| $7/20$ | $(8, 23/4)$ | $\approx 3.265 \times 10^{-2}$ | $\approx 2.744 \times 10^{-2}$ | $0.383508$ | $0.994380622$ |
| $1/4$ | $(8, 25/4)$ | $\approx 7.987 \times 10^{-2}$ | $\approx 4.476 \times 10^{-2}$ | $0.214486$ | $0.812031433$ |

The exact rational seed defects (`"num//den"`) are stored as `epsilon_s` and `epsilon_t` in
`results/certificate_<h>.json`; for $h = 1$, they are
`3124087506499//4503599627370496` and `6599110342439//9007199254740992`. Each certified gap lower
bound is obtained from `CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀).bound`. The headline bounds
are rounded down to six decimal places from the minimum of the original certificate and an
independent full spectral replay (`results/spectral_replay_summary.json`), so both rational seed
pairs establish the stated six-decimal bounds (verified by `check.jl`; see [`METHOD.md`](METHOD.md),
"Independent audit"). The 16-site exact diagonalization (ED) gaps are shown for numerical comparison.

**Independent shared-core cross-check:** `scripts/core_certificate.jl` certifies the $h = 1$ seed
$(8, 21/4)$ a second time using `CertifiedGaps.GapCertificates` (which uses a different block
decomposition by $S^z$-flip parity and single-site momentum in the rotated frame, together with
the exact rational doubling chain), obtaining $\Delta_L > 1.147699$ for every even $L \ge 8$ with
$\varepsilon_{\mathrm{s}} = \text{692723528103}//10^{15}$ and
$\varepsilon_\tau = \text{728796111263}//10^{15}$ (`results/core_certificate_1_1.json`).

## Reproduce

From the repository root, using Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`:

```sh
julia --project=models/heisenberg_staggered_field -e 'using Pkg; Pkg.instantiate()'
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/check.jl
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/check.jl --full
# Recertify a single field at a specified seed (rewrites results/certificate_<h>.json):
julia --project=models/heisenberg_staggered_field models/heisenberg_staggered_field/certify.jl \
    '1//1' 8 '21//4'
# Full spectral replay of all five certificates (rewrites results/spectral_replay_summary.json):
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/spectral_replay.jl [--full]
# Direct interval re-accumulation of residuals on small blocks; shared-core cross-check:
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/audit_bounds.jl
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/core_certificate.jl
# Exploratory numerical seed search (all fields up to 16 sites; one field up to 20 sites):
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/scan.jl 16
julia --project=models/heisenberg_staggered_field \
    models/heisenberg_staggered_field/scripts/scan.jl 20 '1//4'
```

Running `check.jl` executes `test/runtests.jl`. Calling `certify.jl h maxL` with two arguments
reads the best candidate seed from `results/scan_<h>_<maxL>.json` (generated by `scripts/scan.jl`,
which caches raw spectra in gitignored `results/spectrum_<h>_<L>.json` files). All commands also
run under `--project=julia`.

Measured run times (Julia 1.12, 1 BLAS thread):
- `check.jl` (`test/runtests.jl`): 341 tests, ~45 s including compilation (~10 s warm on an Apple
  M-series laptop);
- `check.jl --full`: ~2 min 41 s on an Apple M-series laptop (~6 min on Intel Xeon), re-verifying
  the full spectra of the shared-core cross-check;
- `certify.jl`: 44–80 s per field; `scripts/spectral_replay.jl`: ~190 s for all five fields;
- `scripts/audit_bounds.jl`: ~25 s; `scripts/core_certificate.jl`: ~12 min.

## What is rigorous

- **Rigorous (under IEEE binary64 round-to-nearest arithmetic, standard non-fast-math BLAS, and
  outward-rounded `IntervalArithmetic`):** the exact rational reflection-positivity identity,
  complete symmetry-sector coverage (every 8- and 16-site $S^z$ and two-site translation momentum
  block with its exact multiplicity summing to $2^L$), the Weyl/polar-factor eigenvalue enclosures
  using a priori $\gamma_k$ rounding-error bounds on `Float64` residual and Gram matrix products
  (see [`METHOD.md`](METHOD.md), "Eigenvalue-enclosure theorem"), the outward-rounded $\log Z$
  intervals, the exact rational seed defects, and the doubling recurrence bounds.
- **Fast vs. full verification:** `check.jl` recomputes $\log Z$, the rational seed defects (each
  at most the stored one), and
  the gap bounds from the stored eigenvalue enclosures and checks the shared-core cross-check,
  while `certify.jl` and `scripts/spectral_replay.jl` rebuild every momentum block and re-verify
  all dense eigenvector residuals from scratch, and `check.jl --full` re-solves and re-verifies
  the full spectrum of `core_certificate_1_1.json`.
- **Numerical only:** the exploratory $(n_0, \beta_0)$ grid search, 16-site reference ED gaps,
  20-site partial Krylov spectra and candidate seeds ($n_0 = 10$), and the $h^{2/3}$ scaling
  extrapolations.

## Layout

- `certify.jl`, `check.jl`: certificate generator and verifier (fast mode by default; `--full`
  adds complete spectral re-verification of the shared-core cross-check).
- `METHOD.md`: mathematical method, eigenvalue-enclosure theorem, sector tables, 20-site
  exploration, and audit record.
- `IDEAS.md`: route evaluation and outcome.
- `src/common.jl`: exact RP check, two-site translation momentum blocks with interval entry
  enclosures, verified eigenvalue enclosures, $\log Z$ intervals, and rational seed extraction.
- `scripts/`: `scan.jl` (numerical seed scan), `spectral_replay.jl` (fresh spectral replay of all
  five certificates), `audit_bounds.jl` (direct interval residual accumulation on small blocks),
  and `core_certificate.jl` (independent certificate via `CertifiedGaps.GapCertificates`).
- `test/runtests.jl`: test suite executed by `check.jl`.
- `results/`: `certificate_<h>.json` (the five field certificates),
  `spectral_replay_summary.json` (fresh replay seeds and bounds checked by `check.jl`),
  `core_certificate_1_1.json` (shared-core cross-check at $h = 1$), and `audit_bounds.json`
  (interval residual audit).

No file is hash-pinned by its current working-tree path.

## History

These certificates were originally generated when the directory was located at
`certificates/heisenberg_staggered_field` and the shared core package was named `ThermalMPOGap`.
Their `metadata.source_sha256` entries reference the source snapshot preserved in
`historical_sources/` at commit `cb635a1cdf2ab06630b5b0838fac4356d3c33509`; `check.jl` verifies all
recorded source hashes against that git commit when present in the local repository object store.
Later, duplicated helpers moved into the shared core `CertifiedGaps` (`rat` → `ExactArithmetic.parse_rational`, the exact spin-½ matrices →
`Spins.exact_spin_half`), and new results record model-local sources under `source_sha256` and the
core under `core_sha256`. The stored certificates are unchanged. The fast checks of the seed defects are monotone: the defects recomputed from the stored
enclosures (through the core `Seeds.seed_defects`) must be at most the stored εₛ, εₜ in exact
rationals, rather than equal to them; the stored seeds and bounds are unchanged.
