# Two-leg transverse-field Ising ladder

This module certifies uniform spectral gaps for the two-leg transverse-field Ising ladder at
transverse fields $`g \in \{11/4, 3, 4\}`$, each with a unique ground state, for every even
$L \ge 4$ rungs ($N = 2L$ spins). The certificates are produced by the shared core
(`CertifiedGaps.GapCertificates.certify_point`) from the ladder written as a
reflection-positive `RingModel` of rungs. The bounds use **coupled half-length seeds** at
$n_0 = 6$ rungs: the complete verified spectrum of the 6-rung ring and rigorous Rayleigh–Ritz
upper bounds on the low levels of the 12-rung ring (24 spins, core `RitzBounds`), plus the
verified gap of the 4-rung ring. Method, scans and diagnostics: [`METHOD.md`](METHOD.md); route
analysis: [`IDEAS.md`](IDEAS.md).

| $g$ | seed $(n_0 \text{ rungs}, \beta_0)$ | certified gap $\ge$ | lengths covered |
|---|---|---:|---|
| 11/4 | $(6, 69/32)$ | 0.371180737 | every even $L \ge 4$ rungs |
| 3 | $(6, 29/16)$ | 0.849386526 | every even $L \ge 4$ rungs |
| 4 | $(6, 11/8)$ | 2.546626505 | every even $L \ge 4$ rungs |

The $g = 3$ certificate with seed $(4, 3/2)$, $\Delta_L \ge 0.636870853$, is kept: it carries
the seeds of the Lean 4 theorem `CertifiedGaps.TwoLegIsing` (seed length 4).

## Model

Each unit cell is a two-spin rung ($d = 4$). In Pauli normalization on a periodic ring of $L$
rungs,

```math
H_L = -\sum_{j=1}^L (Z_{j,a} Z_{j+1,a} + Z_{j,b} Z_{j+1,b} + Z_{j,a} Z_{j,b}) - g \sum_{j=1}^L (X_{j,a} + X_{j,b}).
```

Lengths count rungs, and the $L = 2$ ring retains both periodic leg bonds, following the core
convention. The certificate works in the **rotated frame** obtained by a Hadamard gate on every
spin, i.e. the exact rational orthogonal involution $W = (H \otimes H)/2$ with
$`H = \begin{pmatrix}1&1\\1&-1\end{pmatrix}`$ on every rung (`Ladder.rung_hadamard`). It maps
$X \leftrightarrow Z$, so the rotated ladder (`Ladder.ladder(g)`) has bond
$h = -(X_a X_{a'} + X_b X_{b'})$ and on-site term $o = -X_a X_b - g(Z_a + Z_b)$. There the
$`\mathbb{Z}_2`$ charge (number of down spins mod 2, `charge = [0, 1, 1, 0]`) is conserved and
the leg exchange $`a \leftrightarrow b`$ is an exact on-site $`\mathbb{Z}_2`$ symmetry
(`symmetries = [LEG_SWAP]`). The lattice reflection $j \mapsto -j$ maps every rung to itself
and leaves the bond invariant under exchanging its two rungs, so it is declared with the
trivial on-site map (`reflection = I`, checked exactly by `RingModel`). `check.jl` and the
tests check exactly, with the core `is_frame_change`, that $W$ maps the physical ladder
(`Ladder.physical_ladder(g)`) onto the rotated one: $`W o_{\text{phys}} W = o_{\text{rot}}`$
and $(W \otimes W) h_{\text{phys}} (W \otimes W) = h_{\text{rot}}$. $W^{\otimes L}$ is
orthogonal, so both frames have the same spectrum and purities.

## Reflection positivity

In the rotated frame, the real bond factors $C_1 = X \otimes I$ and $C_2 = I \otimes X$ give

```math
-h = C_1 \otimes C_1 + C_2 \otimes C_2,
```

and the on-site term $o$ is real symmetric; a single rung needs no internal reflection ($F = I$).
`Ladder.ladder_rp` passes these factors to the core `rp_certificate`, whose `verify_rp` checks
the identity over `Rational{BigInt}`; the certificate stores them and `check.jl` re-verifies
them. In the physical frame the same identity holds with $C_1 = Z \otimes I$, $C_2 = I \otimes Z$
(the form used in `code/lean/CertifiedGaps/Models/TwoLegIsing.lean`), which the tests also check.

## Seeds and certified bounds

Coupled seeds at $n_0 = 6$ rungs (`RitzSeeds.RITZ_POINTS`, certificates
`results/certificate_g<g>_n6.json`). By the coupling identity (AGENTS.md; IDEAS.md, route 4)
and $P_{\mathrm{s}}(6, \beta) \le 1$ (reflection positivity, even $n_0$),

```math
\varepsilon_{\mathrm{s}} \ge 1 - \frac{Z_{12}(\beta_0)}{Z_6(\beta_0)^2}, \qquad
\varepsilon_\tau \ge 1 - \frac{Z_{12}(2\beta_0)}{Z_6(\beta_0)^4} \ge 1 - P_\tau(12, \beta_0),
```

with $Z_6$ enclosed from the complete verified 6-rung spectrum and $Z_{12}$ bounded below by
Rayleigh–Ritz upper bounds on 6 (charge 0) and 3 (charge 1) levels of each of the 28 blocks of
the 12-rung ring (core `coupled_seed_bounds`). The exact rationals are `seed.eps_s`,
`seed.eps_t` of each certificate. They are raised by $10^{-8}$ (`defect_padding`, a valid
weakening): the verified radii of LAPACK eigenpairs vary slightly from run to run, and
`check.jl` recomputes the defects from a fresh 6-rung spectrum.

| $g$ | seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}}$ | $\varepsilon_\tau$ | certified gap $\ge$ (interval) | exact chain bound | 4-rung gap $\ge$ (verified) | 12-rung gap (Lanczos) |
|---|---|---:|---:|---:|---:|---:|---:|
| 11/4 | $(6, 69/32)$ | 0.049856612220172 | 0.126479296586788 | 0.371180737 | $> 0.370813466$ | 2.2092349 | 2.1190 |
| 3 | $(6, 29/16)$ | 0.024179319356254 | 0.075932324712783 | 0.849386526 | $> 0.848512677$ | 2.7219276 | 2.6599 |
| 4 | $(6, 11/8)$ | 0.002828638097851 | 0.012860144463702 | 2.546626505 | $> 2.544322725$ | 4.7822235 | 4.7632 |

The headline is the outward-rounded interval bound `interval_bound_decimal` of
`certified_gap_bound(εₛ, εₜ, β₀)` (rounded down above); the exact rational doubling chain
`gap_lower_bound` is slightly smaller. Both are rigorous for every even $L \ge 6$. Each
certificate also stores `finite_bridge`, outward bounds on the gap of the 4-rung ring from its
complete verified spectrum; it exceeds the bound, so the claim covers every even $L \ge 4$.
The 12-rung gaps are numerical (Lanczos, `scripts/ritz_scan.jl`).

The complete-spectrum seed (`Ladder.POINTS`, `results/certificate_g3_1.json`):

| $g$ | seed $(n_0, \beta_0)$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(4, 3/2)$ | $\varepsilon_\tau \ge 1 - P_\tau(8, 3/2)$ | certified gap $\ge$ (interval) | exact chain bound |
|---|---|---:|---:|---:|---:|
| 3 | $(4, 3/2)$ | $0.112389665753866$ | $0.051536217029413$ | 0.636870853 | $> 0.636342903$ |

Its defects satisfy the seed hypotheses of the Lean theorem ($\varepsilon_{\mathrm{s}} \le 0.1125$,
$\varepsilon_\tau \le 0.05162$, read from `code/lean/CertifiedGaps/Models/TwoLegIsing.lean` by
the core `lean_points` and checked exactly by `check.jl`). $g = 5/2$ has no admissible seed at
$n_0 \le 4$ (complete spectra) nor at $n_0 = 6$ with the coupled defects (minimum coupled
$\varepsilon_\tau \approx 0.24$); $g = 2$ has none at $n_0 \le 4$.

## Reproduce

From the repository root, with Julia 1.12 and `export OPENBLAS_NUM_THREADS=1`, in the core
environment:

```sh
julia --project=code/julia -e 'using Pkg; Pkg.instantiate()'
julia --project=code/julia code/models/two_leg_tfim/certify.jl          # every certificate
julia --project=code/julia code/models/two_leg_tfim/certify.jl 3//1 --ritz      # one point
julia --project=code/julia code/models/two_leg_tfim/check.jl            # fast checker
julia --project=code/julia code/models/two_leg_tfim/check.jl --full     # + fresh Ritz bounds
julia --project=code/julia code/models/two_leg_tfim/test/runtests.jl    # model unit tests
# Floating-point scans: complete spectra up to 8 rungs, and the coupled n0 = 6 seeds
julia --project=code/julia code/models/two_leg_tfim/scripts/scan.jl 3//1
julia --project=code/julia code/models/two_leg_tfim/scripts/ritz_scan.jl 3//1 14 6
julia --project=code/julia code/models/two_leg_tfim/scripts/extrapolate.jl
julia --project=code/julia code/models/two_leg_tfim/scripts/summarize.jl
```

`certify.jl [g] [--ritz | --complete] [--reuse-partial]` certifies the points of
`Ladder.POINTS` (complete spectra, $n_0 = 4$) and `RitzSeeds.RITZ_POINTS` (coupled, $n_0 = 6$),
or those at field `g`; `--reuse-partial` reads the stored 12-rung Ritz bounds instead of
recomputing them. The stored certificates were last regenerated that way (all four in 110 s,
after the defect padding was added), from Ritz bounds computed by full runs. Measured run
times on a shared 8-core Linux machine (Julia 1.12.7, one BLAS thread, load average 8–12):
- `certify.jl g --ritz`: 35–41 min wall per point, peak memory 4.2 GB (Julia started with
  `--heap-size-hint=3G` or `2500M`); about 25 min build the 28 interval-enclosed sparse
  blocks of the 12-rung ring (dimension about $3.5 \times 10^5$ each), 12–13 min run Lanczos and
  the Ritz bounds, and the 6-rung spectrum and the chain take seconds;
- `certify.jl 3//1 --complete` ($n_0 = 4$): 119 s wall (92 s in `certify_point`), peak memory
  1.3 GB;
- `check.jl`: 33 s including compilation, peak memory 1.2 GB (recomputes the 6- and 4-rung
  spectra and reads the stored 12-rung bounds);
- `check.jl --full`: not run as one command in this session (about 2 h: the fresh 12-rung
Ritz bounds take 37–40 min per coupled point, plus about 4 min for the $n_0 = 4$ spectra).
Its steps for each coupled point (fresh Ritz bounds and 6-rung spectrum, then `check_point`
with these spectra and the bridge) passed for all three points: $g = 4$ in 39 min,
$g = 3$ in 37 min, $g = 11/4$ in 37 min, peak memory 4.2 GB each; the fresh defects equal
the unpadded certified ones to all printed digits;
- `test/runtests.jl`: 34 tests, 47 s including compilation;
- `scripts/ritz_scan.jl`: 15–60 min per field (Lanczos on 28 complex sectors of the
  physical-frame ring), `scripts/scan.jl`: 2–4 min per field.

## What is rigorous

- **Rigorous:** the exact RP identity and frame change (`Rational{BigInt}`); the complete
  spectra of the rings of 4, 6 and 8 rungs, split exactly by the core
  `VerifiedSpectra.ring_blocks` into $`\mathbb{Z}_2`$ charge, lattice momentum (±k merged with
  multiplicity 2) and leg-exchange parity, with block entries enclosed in interval arithmetic
  and every block real symmetric (the lattice reflection makes the momentum blocks real); each
  block enclosed by `CertifiedGaps.VerifiedSpectra.verified_eigen` (residual and Gram matrices
  in `Float64` with the a priori $`\gamma_{N+8}`$ rounding bound and Frobenius norms, then
  Ostrowski/polar-factor and Weyl bounds, plus the exact construction error). For the 12-rung
  ring, the same sparse interval-enclosed blocks and `RitzBounds.ritz_upper_bounds`: from any
  `Float64` basis $V$ (Lanczos vectors, never trusted), a priori $\gamma_n$ bounds on
  $\mathrm{fl}(MV)$, $V^{\mathsf T} W$ and $V^{\mathsf T} V$, a Gershgorin shift, Ostrowski
  for the pencil and `verified_eigen` on the small Rayleigh-quotient matrix give rigorous upper
  bounds on the lowest eigenvalues of each block (Courant–Fischer), hence a rigorous lower
  bound on $Z_{12}$. $\log Z$ in interval arithmetic; seed defects rounded up to a $10^{-15}$
  grid; the coupling identity with $P_{\mathrm{s}}(6, \beta) \le 1$ (reflection positivity at
  even $n_0$); the exact rational doubling chain and the interval `certified_gap_bound`; the
  4-rung gap bounds of `finite_bridge` (`spectral_gap_enclosure`, rounded outward). Trusted
  base: IEEE binary64 round-to-nearest with gradual underflow, standard (non-fast-math) BLAS
  and LAPACK producing finite output, `IntervalArithmetic`, and the repository's two-purity
  theorem.
- **Fast vs. full verification:** `check.jl` (core `check_point`) requires the model-local
  source hashes (`RitzSeeds.SOURCES`), the stored model, RP factors, seed and params to be
  exactly those of the point and of `Ladder.ladder(g)`, the Lean seed hypotheses (seed
  $(4, 3/2)$) and the Hadamard frame change, and re-verifies the RP identity, the exact chain
  and the interval bound. For the coupled seeds it also recomputes the defects from a freshly
  verified 6-rung spectrum and the stored Ritz bounds `results/partial_g*_L12.json` (whose
  SHA-256 the certificate pins) and recomputes the 4-rung bridge. **The fast check trusts the
  stored 12-rung Ritz bounds as written by `certify.jl`**; `--full` recomputes them (and every
  complete spectrum of the $n_0 = 4$ point) and requires the fresh defects to be at most the
  stored ones. Lanczos starts from deterministic vectors; on another platform the fresh Ritz
  bounds may differ in the last bits, so `--full` can in principle fail by $10^{-15}$ without
  any error in the certificate.
- **Numerical only (not used in the proof):** the seed searches (`scripts/scan.jl`:
  $`g \in \{2, 5/2, 3, 4\}`$, $`n_0 \in \{2, 3, 4\}`$; `scripts/ritz_scan.jl`:
  $`g \in \{5/2, 11/4, 3, 4\}`$ at $n_0 = 6$ with Lanczos values of the 12-rung ring), the
  choice of $\beta_0$ and of the number of bounded levels, the negative results at $g = 2$
  and $5/2$, the extrapolated seed sizes, finite-ring ED and Lanczos reference gaps, and the
  literature critical field $g_c \approx 1.83214(5)$.

## Layout

Stored results live in [`data/models/two_leg_tfim/`](../../../data/models/two_leg_tfim/);
this directory's `results/` is a relative symlink to that canonical location.

- `certify.jl`, `check.jl`: certificate generator (core `certify_point`, with
  `coupled_seed_bounds` for the coupled points) and checker (core `check_point`, with the Lean
  seeds from `lean_points`).
- `src/Ladder.jl`: the model: physical and rotated terms and models, the frame change, the core
  `RingModel` with its charge, leg-exchange symmetry and lattice reflection, the RP factors,
  the complete-spectrum points `POINTS` and their certificate file names.
- `src/RitzSeeds.jl`: the coupled seeds: Ritz levels per block, `doubled_ring_bounds`,
  `finite_bridge`, file names and the pinned `SOURCES`; `src/ritz_points.jl`: the list
  `RITZ_POINTS` (not pinned, so that adding a point leaves the other certificates valid).
- `scripts/`: `scan.jl` (numerical seed scan, $n_0 \le 4$), `ritz_scan.jl` and `sectors.jl`
  (floating-point screening of the coupled $n_0 = 6$ seeds with a sparse sector builder of the
  physical-frame ring), `extrapolate.jl`, `summarize.jl` (rewrites the `## Results` section of
  `METHOD.md`).
- `test/runtests.jl`: exact frame and RP tests, core blocks against an independent full
  Hamiltonian, Ritz bounds and coupled defects on small rings against complete spectra, and
  the 4-rung bridge (the stored certificates are checked by `check.jl`).
- `results/certificate_g3_1.json` ($n_0 = 4$) and `results/certificate_g<g>_n6.json`
  (coupled): core-format certificates (model, RP data, seed, log Z enclosures, spectrum
  summaries, exact chain, `params`, `metadata` with the model-local `source_sha256` and the
  core `core_sha256`; the coupled ones also `partial_spectrum` and `finite_bridge`);
  `results/partial_g<g>_L12.json`: the Ritz upper bounds of the 12-rung ring (hash-pinned by
  the coupled certificates); `scan_g*.json`, `ritz_scan_g*.json`, `ritz_spectra_g*.json`
  (numerical scans) and `extrapolation.json`.

## History

The original certificate (`results/certificate_g3_1.json.gz`, generated when the directory was
`certificates/two_leg_tfim`) used local symmetry blocks (parity, leg exchange, reflection, half
translation), a local polar/Weyl enclosure and MPFR log Z. It was retired together with its
replay scripts, `environment.json` and the model's own Julia environment. The model now runs
entirely on the shared core (`certify_gap` with the leg-exchange symmetry) and was regenerated
with the same seed: $\varepsilon_{\mathrm{s}}$ $0.1124084 \to 0.1123897$,
$\varepsilon_\tau$ $0.0516152 \to 0.0515362$, interval bound $0.636387 \to 0.636871$ (exact
chain bound $0.636343$). It was then regenerated through the core drivers `certify_point` and
`check_point`, with the lattice reflection that makes every block real (same seed):
$\varepsilon_{\mathrm{s}}$ $0.11238966699 \to 0.11238966593$, $\varepsilon_\tau$
$0.05153622006 \to 0.05153621778$, interval bound $0.6368708320 \to 0.6368708493$, exact chain
bound $0.636342881 \to 0.636342898$, largest 8-rung radius
$1.5 \times 10^{-7} \to 5.7 \times 10^{-8}$.

Coupled half-length seeds at $n_0 = 6$ (core `RitzBounds`) then raised the $g = 3$ bound from
$0.636870849$ to $0.849386526$ and added the points $g = 11/4$ and $g = 4$. `check.jl` and
`certify.jl` were extended for them, and every certificate now pins `RitzSeeds.SOURCES`
instead of the whole `src/`, so the $n_0 = 4$ certificate was regenerated with the new
Coupled half-length seeds at $n_0 = 6$ (core `RitzBounds`) then raised the $g = 3$ bound from
$0.636870849$ to $0.849386526$ and added the points $g = 11/4$ and $g = 4$. `check.jl` and
`certify.jl` were extended for them, and every certificate now pins `RitzSeeds.SOURCES`
instead of the whole `src/`, so the $n_0 = 4$ certificate was regenerated with the new drivers
(same seed). The enclosures differ between platforms in the last digits, and on the Linux
machine of the coupled runs they came out slightly tighter: $\varepsilon_{\mathrm{s}}$
$0.112389665926876 \to 0.112389665753866$, $\varepsilon_\tau$
$0.051536217782645 \to 0.051536217029413$, interval bound $0.6368708493 \to 0.6368708539$,
exact chain bound $0.636342898 \to 0.636342903$; the Lean seed hypotheses still hold.
