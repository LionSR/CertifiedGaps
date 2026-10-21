# Spin-1/2 Heisenberg chain in a staggered field

**Certified at $h = 1, 3/4, 1/2, 7/20, 1/4$: a unique ground state and a uniform spectral gap
$\Delta_L \ge 1.148970, 0.890950, 0.594821, 0.383520, 0.214491$ for every even chain length
$L \ge 8$ sites.** Each certificate is produced by the shared core pipeline
(`CertifiedGaps.GapCertificates.certify_point`) from verified complete spectra of the 8- and
16-site rings. These spectra are split by an exact two-site blocking with a conserved U(1) cell
charge. Derivations, sector data and the 20-site exploration are in [`METHOD.md`](METHOD.md);
the route memo is in [`IDEAS.md`](IDEAS.md).

| Staggered field $h$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $1$ | $(8, 21/4)$ | $1.148970$ | every even $L \ge 8$ sites |
| $3/4$ | $(8, 21/4)$ | $0.890950$ | every even $L \ge 8$ sites |
| $1/2$ | $(8, 11/2)$ | $0.594821$ | every even $L \ge 8$ sites |
| $7/20$ | $(8, 23/4)$ | $0.383520$ | every even $L \ge 8$ sites |
| $1/4$ | $(8, 25/4)$ | $0.214491$ | every even $L \ge 8$ sites |

## Model

We consider the periodic spin-1/2 antiferromagnetic Heisenberg chain in a staggered longitudinal
magnetic field $h$:
```math
H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1} - h \sum_{j=1}^L (-1)^j S_j^z, \qquad \mathbf{S}_{L+1} \equiv \mathbf{S}_1,
```
with unit exchange coupling and $\mathbf{S}_j = \frac{1}{2}\boldsymbol{\sigma}_j$. The two-site
ring ($L = 2$) includes both bonds $(1, 2)$ and $(2, 1)$, as in the core convention. The field
has a two-site period. After the sublattice rotation below, the Hamiltonian is translation
invariant with a **one-spin unit cell** on every even ring. Lengths $L$ and $n_0$ therefore
count sites.

## Reflection positivity

**Frame.** Let $X = S^x$, $G_y = iS^y$ (real antisymmetric, `J` in `src/common.jl`) and
$Z = S^z$ be the exact rational matrices, and let $U = 2G_y = e^{i\pi S^y}$. The exact unitary
$V_L = \bigotimes_{j\ \mathrm{odd}} U_j$ (a $\pi$ rotation about $y$ on every odd site) maps
$H_L$ to the one-site ring
```math
V_L H_L V_L^\dagger = \sum_j h_{j,j+1} + \sum_j o_j, \qquad h = -(X \otimes X + G_y \otimes G_y + Z \otimes Z), \qquad o = -hZ,
```
because $U$ negates $S^x$ and $S^z$ and fixes $G_y$, and because
$S^y \otimes S^y = -G_y \otimes G_y$.
This is `site_model(h)`, the model stored in every certificate.
`test/runtests.jl` checks the rotation exactly with `RingModels.is_frame_change`: written in
two-site cells $(a, b) = (2x - 1, 2x)$, the physical chain is mapped by $U \otimes I$ onto
`cell_model(site_model(h), 2)` for every field, hence on every even ring.

**RP identity.** $`-h = \sum_{C \in \{X, G_y, Z\}} C \otimes \overline{C}`$ with real $C$, and
$o = -hZ$ is real and symmetric (`site_rp(h)`). The core verifies it exactly
(`RingModels.verify_rp`) when certifying and again in `check_certificate`.

**Lattice reflection.** The rotated bond is symmetric under the exchange of its two sites, so
`site_model` carries the lattice reflection $P = I$, which the core checks exactly. The seed
spectra do not use it. In two-site cells, every on-site map that mirrors the bond sends the cell
charge below to its negative, so the blocked model admits no reflection. Diagonalizing the site
model with its reflection and $S^z$ parity instead gives larger defects (see `METHOD.md`).

**Blocking used for the spectra.** In this frame only the $S^z$ parity is a one-site charge.
In two-site cells $(a, b)$, the rotated bond flips both spins of a pair together or neither, so
the cell charge $q = S^z_a - S^z_b$ (values $0, 1, -1, 0$ on
$\uparrow\uparrow, \uparrow\downarrow, \downarrow\uparrow, \downarrow\downarrow$) is conserved. It is minus the physical
magnetization of the cell. The core builds the exact cell model (`RingModels.cell_model`) and
diagonalizes rings of $n_0/2$ and $n_0$ cells in (cell charge, cell momentum) blocks. The
certificate stores this blocking, and `check_certificate` rebuilds it from the site model. No
on-site cell operator maps $q \to -q$, because spin flip reverses the field and the sublattice
exchange is not on-site, so no charge flip is used.

## Seeds and certified bounds

| $h$ | Seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(8, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(16, \beta_0)$ | Interval bound $\Delta_L \ge$ | Exact chain bound $\Delta_L >$ | 16-site ED gap |
|---|---|---:|---:|---:|---:|---:|
| $1$ | $(8, 21/4)$ | $6.926 \times 10^{-4}$ | $7.284 \times 10^{-4}$ | $1.148970$ | $1.147763329$ | $1.87351193$ |
| $3/4$ | $(8, 21/4)$ | $2.342 \times 10^{-3}$ | $3.119 \times 10^{-3}$ | $0.890950$ | $0.889743689$ | $1.57089206$ |
| $1/2$ | $(8, 11/2)$ | $1.074 \times 10^{-2}$ | $1.115 \times 10^{-2}$ | $0.594821$ | $0.594245333$ | $1.23131976$ |
| $7/20$ | $(8, 23/4)$ | $3.265 \times 10^{-2}$ | $2.743 \times 10^{-2}$ | $0.383520$ | $0.383244983$ | $0.994380622$ |
| $1/4$ | $(8, 25/4)$ | $7.987 \times 10^{-2}$ | $4.476 \times 10^{-2}$ | $0.214491$ | $0.214238387$ | $0.812031433$ |

The exact seed defects are stored as `seed.eps_s` and `seed.eps_t` in
`results/certificate_<h>.json`. They are upper ends of the defect intervals, rounded up to the grid
$10^{-15}$. For $h = 1$ they are `27705127843//40000000000000` and
`364207460943//500000000000000`. The headline bounds are `interval_bound_decimal`
(`certified_gap_bound(εₛ, εₜ, β₀)`, outward-rounded interval arithmetic), rounded down to six
decimals. The certificate's claim `gap_lower_bound` is the strict exact rational chain bound
(`rational_gap_chain`), which is slightly smaller. Every $\varepsilon$ is at most the seed constant
hardcoded in `code/lean/CertifiedGaps/Models/StaggeredHeisenberg.lean` (first entries of `s` and `t`),
and `check.jl` checks this in exact rationals
(`GapCertificates.lean_points`). The 16-site exact diagonalization gaps
are floating-point references.

## Reproduce

From the repository root, using Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`, in the core
environment (the model has no environment of its own):

```sh
julia --project=code/julia -e 'using Pkg; Pkg.instantiate()'
# All five certificates and results/summary.json (or one field: certify.jl '1//2'):
julia --project=code/julia code/models/heisenberg_staggered_field/certify.jl
julia --project=code/julia code/models/heisenberg_staggered_field/check.jl
julia --project=code/julia code/models/heisenberg_staggered_field/check.jl --full
julia --project=code/julia code/models/heisenberg_staggered_field/test/runtests.jl
# Exploratory floating-point seed scan (gitignored output):
julia --project=code/julia code/models/heisenberg_staggered_field/scripts/scan.jl 16
```

Measured run times (Julia 1.12, 1 BLAS thread, Apple M-series laptop under heavy concurrent
load, so wall times are about twice the CPU times):
- `certify.jl`: 14 min wall time (481 s CPU) for all five fields, 143–222 s per field;
- `check.jl`: 27 s including compilation;
- `check.jl --full`: 15 min wall time (488 s CPU), re-verifying the 8- and 16-site spectra of
  all five certificates (122–241 s per field);
- `test/runtests.jl`: 34 tests, 56 s wall time (28 s CPU).

## What is rigorous

- **Exact arithmetic:** the sublattice rotation (checked on two-site cells, hence on every
  even ring), the lattice reflection of the site model, the RP identity, the conservation of
  the cell charge, the exact cell model, the exact rational chain, and the comparison with the
  Lean seed constants.
- **IEEE binary64 + BLAS rounding model:** every block of the 8- and 16-site rings is enclosed
  by `CertifiedGaps.VerifiedSpectra.verified_eigen`. It bounds the residual and orthogonality
  defect of the LAPACK eigenvectors with a priori $\gamma_{N+8}$ (real blocks) or
  $\gamma_{2N+8}$ (complex blocks, four real GEMMs) dot-product bounds on Frobenius norms. The
  bounds hold for any summation order and with FMA, assume round-to-nearest, and cover
  gradual underflow. A Weyl, Ostrowski or polar-factor bound and the exact block-entry error
  then give the radius. The block entries (phases, orbit norms) are enclosed in
  `IntervalArithmetic`. $\log Z$ is summed in interval arithmetic
  (`log_partition_enclosure`) with exact block multiplicities summing to $2^L$. The trusted
  base is the core (`code/julia/src`), IntervalArithmetic and LAPACK/BLAS used as untrusted
  proposers under that model.
- **Fast vs. full:** `check.jl` checks the exact parts with `GapCertificates.check_point`
  (model-local source hashes, model, RP data, seed, blocking and params against `src/common.jl`,
  the Lean seed constants, the chain and the interval bound) and the summary. `check.jl --full`
  also rebuilds every spectrum from the stored model and requires the fresh defects to be at
  most the stored $\varepsilon$.
- **Floating point only:** `scripts/scan.jl`, the 16-site reference gaps and the 20-site
  candidate seeds in `METHOD.md`.

## Layout

Stored results live in [`data/models/heisenberg_staggered_field/`](../../../data/models/heisenberg_staggered_field/);
this directory's `results/` is a relative symlink to that canonical location.

- `certify.jl`: runs `certify_point` for every point in `SEEDS` and writes the certificates and
  `results/summary.json` (`certificate_row`).
- `check.jl`: `check_point` for every point (`--full`: with the spectra) and the summary.
- `src/common.jl`: the frame (`U`, physical and rotated bonds), `site_model`, `site_rp`, the
  blocking `BLOCKING`, the points `SEEDS` and their arguments `point(h)`.
- `scripts/scan.jl`: floating-point seed scan (`float_ring_spectrum` of the certificate's blocks
  up to 16 sites, the core Lanczos exact diagonalization beyond).
- `test/runtests.jl`: exact frame, RP, reflection and blocking tests, and verified spectra
  against an independent physical-frame exact diagonalization.
- `results/certificate_<h>.json`: core-format gap certificates (one per field), with
  `metadata` from `Results.model_metadata`: `source_sha256` (model-local, enforced by
  `verify_model_metadata`) and `core_sha256` (provenance only).
- `results/summary.json`: compact table of seeds and bounds, checked against the certificates.
- `METHOD.md`, `IDEAS.md`: method details and the route memo.

No file is hash-pinned by its working-tree path.

## History

- The certificates were first produced by a model-local pipeline: two-site translation momentum
  and $S^z$ blocks in the physical frame, a local residual/polar enclosure with $\gamma_{32d+32}$
  bounds, and a local $\log Z$. A shared-core cross-check at $h = 1$ used only the $S^z$ parity
  (about 12 min).
- The local pipeline was replaced by the shared core. Each point is now a core certificate on the
  one-site rotated model with a two-site U(1) blocking. All certificates were regenerated at
  the same seeds. Every $\varepsilon$ decreased (for example, at $h = 1$,
  $\varepsilon_{\mathrm{s}}$ went from $6.937 \times 10^{-4}$ to $6.926 \times 10^{-4}$ and
  $\varepsilon_\tau$ from $7.326 \times 10^{-4}$ to $7.284 \times 10^{-4}$). The interval bounds
  rose: $1.148261 \to 1.148970$, $0.890782 \to 0.890950$, $0.594785 \to 0.594821$,
  $0.383508 \to 0.383520$ and $0.214486 \to 0.214491$. The spectral replay and audit scripts
  and their summaries were retired with the local enclosure, and the model now runs in the core
  environment (`--project=code/julia`) instead of its own.
- The model was moved onto the core drivers `certify_point`, `check_point`, `certificate_row`,
  `lean_points` and `float_ring_spectrum`; `check.jl` no longer runs the tests, the frame test
  uses `is_frame_change` on two-site cells, the site model records its reflection $P = I$, and
  every certificate records `params`. The regenerated certificates have exactly the same seeds,
  defects and bounds.
