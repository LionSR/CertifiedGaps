# Dimerized spin-1/2 Heisenberg chain

**Certified at $\delta = 1/2, 3/10,$ and $1/4$:** for every even ring length $L \ge 4$ two-spin
unit cells ($N = 2L = 8, 12, 16, \dots$ spins, i.e., $N$ divisible by $4$), the ground state is
unique and the spectral gap satisfies $\Delta_L \ge 0.534794, 0.173204, 0.058232$ (the interval
bound `certified_gap_bound`; the exact rational chain gives
$\Delta_L > 0.534267, 0.173082, 0.058177$). Each certificate is a core
`GapCertificates.certify_point` certificate of the chain in a sublattice-rotated frame (exact RP
identity, verified spectra, both bounds). On the grid
$`\delta \in \{1/2, 3/10, 1/4, 1/5, 3/20, 1/10\}`$, $\delta = 1/4$ is the smallest dimerization
certifiable from rings of at most $16$ spins: at $`\delta \in \{1/5, 3/20, 1/10\}`$,
continuous-$`\beta`$ interval covers (core `Exclusions`) prove that **every** $\beta_0 > 0$
fails the doubling recurrence at both even seed lengths $n_0 = 2$ and $n_0 = 4$ cells.

| Dimerization $\delta$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $1/2$ | $(4, 6)$ cells | $0.534794$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $3/10$ | $(4, 13/2)$ cells | $0.173204$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $1/4$ | $(4, 29/4)$ cells | $0.058232$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $1/5, 3/20, 1/10$ | none: all $\beta_0 > 0$ excluded at $n_0 = 2, 4$ | — | no certificate at $\le 16$ spins |

## Model

A periodic ring of $N = 2L$ spin-1/2 sites grouped into $L$ two-spin unit cells (strong dimers),
with bond strength $1 + \delta$ inside each cell and $1 - \delta$ between adjacent cells:
```math
H_L = (1 + \delta) \sum_{j=1}^L \mathbf{S}_{2j-1} \cdot \mathbf{S}_{2j}
      + (1 - \delta) \sum_{j=1}^L \mathbf{S}_{2j} \cdot \mathbf{S}_{2j+1}, \qquad
\mathbf{S}_m = \frac{1}{2}\boldsymbol{\sigma}_m.
```
Lengths $L$ and $n_0$ count two-spin cells.

**Certified frame.** Rotating the second spin of every cell by
$`J = 2iS^y = \begin{psmallmatrix} 0 & 1 \\ -1 & 0 \end{psmallmatrix}`$, i.e. $U = (I \otimes J)^{\otimes L}$,
gives $U^\top H_L U = \sum_j h_{j,j+1} + \sum_j o_j$ (the same spectrum) with
```math
o = (1 + \delta)(-S^xS^x + S^yS^y - S^zS^z), \qquad
h = (1 - \delta)\, I \otimes (-S^xS^x + S^yS^y - S^zS^z) \otimes I,
```
both real rational $4 \times 4$ and $16 \times 16$ matrices on the cell basis

$|{\uparrow\uparrow}\rangle, |{\uparrow\downarrow}\rangle, |{\downarrow\uparrow}\rangle, |{\downarrow\downarrow}\rangle$. In this frame the cell charge $z_a - z_b = [0, 1, -1, 0]$ is
conserved and reversed by the cell spin flip $\sigma^x \otimes \sigma^x$. The lattice reflection
exchanges the two spins of every cell, which also reverses the charge, so its on-site part is the
cell swap composed with the spin flip (it fixes $|{\uparrow\downarrow}\rangle$ and $|{\downarrow\uparrow}\rangle$). `rotated_model(δ)` in
`src/DimerizedHeisenberg.jl` is this core `RingModel(4, h, o; charge, flip, reflection)`, and
`frame_identity(δ)` checks exactly (core `is_frame_change`) that $W = I \otimes J$ is orthogonal
and maps the cell term and bond of `physical_model(δ)` onto $o$ and $h$, which implies the
unitary equivalence for every $L$ (also cross-checked against the `Float64` spectrum of the
physical chain in the tests).

## Reflection positivity

With $F$ the swap of the two spins of a cell and $`P_\alpha \in \{S^x, iS^y, S^z\}`$ (real),
```math
-h = (1 - \delta) \sum_\alpha C_\alpha \otimes F\, \overline{C_\alpha}\, F, \qquad
C_\alpha = I \otimes P_\alpha, \qquad F o F = o,
```
the reflected form of `paper/thermal-mpo-gap-scope.tex` (Appendix E) and of the Lean structure
`RPRing` (this is the frame of `code/lean/CertifiedGaps/Models/DimerizedHeisenberg.lean`).
`rotated_rp(δ)` builds it with the core `rp_certificate(...; reflection = F)`; every certificate
stores factors, weights, and $F$, and `RingModels.verify_rp` checks the identity exactly. Without
$F$ the bond is not RP in any common cell basis (
$`\|h - F_{\mathrm{outer}} h F_{\mathrm{outer}}\|_F^2 = \frac{3}{8}(1 - \delta)^2 > 0`$ for the swap of two whole cells,
`src/reflection_conventions.jl`). `src/rp.jl` audits an equivalent conjugate form in the
singlet/triplet basis exactly (`results/rp_audit.json`); see [`METHOD.md`](METHOD.md).

## Seeds and certified bounds

For each seed, `certify.jl` computes the verified spectra of the rings of $n_0 = 4$ and
$2n_0 = 8$ cells ($8$ and $16$ spins) with the core `verified_ring_spectrum` (charge, momentum,
and charge-flip blocks, the self-conjugate charge split by flip parity, momenta $\pm k$ merged by
the lattice reflection and built as real symmetric blocks) and runs `certify_point`:
seed defects rounded up to $10^{-15}$, the interval bound `certified_gap_bound`
($\Delta_L \ge$, the headline), and the exact rational doubling chain ($\Delta_L >$, the form of
the Lean theorem). Both are rigorous and stored in each certificate.

| $\delta$ | $\beta_0$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(4, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(8, \beta_0)$ | Interval $\Delta_L \ge$ | Chain $\Delta_L >$ | Enclosed $\Delta_{16}$ | Bound / $\Delta_{16}$ |
|---|---|---:|---:|---:|---:|---|---:|
| $1/2$ | $6$ | $0.011265192$ | $0.012004580$ | $0.534794$ | $0.534267$ | $[1.188651, 1.188652]$ | $0.4499$ |
| $1/2$ | $25/4$ | $0.011707776$ | $0.008704102$ | $0.533722$ | $0.533215$ | (survey $\beta_0$) | |
| $3/10$ | $13/2$ | $0.095113917$ | $0.051604791$ | $0.173204$ | $0.173082$ | $[0.811144, 0.811145]$ | $0.2135$ |
| $1/4$ | $29/4$ | $0.165676942$ | $0.051034421$ | $0.058232$ | $0.058177$ | $[0.709347, 0.709348]$ | $0.0821$ |

($\varepsilon$ rounded up to nine decimals here; the certificates store them exactly.) At
$\delta = 1/4$ the spatial defect exceeds the diagonal threshold $\varepsilon_* \approx 0.12731$,
but the unequal-defect recurrence closes thanks to the smaller thermal defect. The Lean theorems in
`code/lean/CertifiedGaps/Models/DimerizedHeisenberg*.lean` take

$(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) \le (0.01133, 0.01224), (0.09516, 0.0518), (0.1658, 0.05123)$ at the same $(4, \beta_0)$ as hypotheses; the certified defects satisfy them
(checked exactly by `check.jl` through the core `check_point` with the Lean file's `points`,
read by `lean_points`).

### Rigorous all-$`\beta`$ exclusions at $n_0 = 2, 4$ for smaller $\delta$

With the largest ring capped at $16$ spins, the only even seed lengths are $n_0 = 2$ and
$n_0 = 4$ cells. For $n_0 = 2$ at all six $\delta$ values and for $n_0 = 4$ at
$`\delta \in \{1/5, 3/20, 1/10\}`$, the core `Exclusions.exclude_all_beta` proves from the
verified spectra that the doubling recurrence fails for **every** $\beta_0 > 0$:
1. a positive ground-energy mismatch $`d_E = \underline{E_0(2n_0)} - 2\,\overline{E_0(n_0)} > 0`$
   rejects every $\beta_0 \ge B$ ($B = 4$ at $n_0 = 2$, $B = 8$ at $n_0 = 4$);
2. $[0, B]$ is covered by $1\text{--}4$ dyadic windows (core `shifted_log_partition`), on each
   of which lower bounds on both defects leave the admissible region of the monotone doubling
   map (`Exclusions.rejection`).

Certifying $\delta \le 1/5$ therefore requires larger seed rings ($n_0 \ge 6$ cells, i.e.,
$\ge 24$ spins in the doubled ring) or a different method.

### Numerical diagnostics up to 20 spins and dimer perturbation

The scripts `scripts/numerical_scan.jl` (complete `Float64` spectra up to 16 spins,
`results/numerical_summary.json`) and `scripts/scan20.jl` (20-spin rings across all symmetry
sectors using 40 Krylov Ritz values in large blocks, `results/diagnostic20_summary.json`) provide
numerical diagnostics across all six $\delta$ values:

| $\delta$ | $\Delta_{20}$ (numerical) | $\xi_{\mathrm{eff}}$ (cells) | $`\xi_{\mathrm{eff}}\,\delta^{2/3}`$ | Estimated even $n_0$ by $\delta^{-2/3}$ scaling |
|---|---:|---:|---:|---:|
| $1/2$ | $1.188631$ | $0.5624$ | $0.3543$ | $4$ |
| $3/10$ | $0.810214$ | $0.8448$ | $0.3786$ | $4$ |
| $1/4$ | $0.707193$ | $0.9682$ | $0.3842$ | $4$ |
| $1/5$ | $0.599819$ | $1.1430$ | $0.3909$ | $6$ |
| $3/20$ | $0.488096$ | $1.4171$ | $0.4001$ | $6$ |
| $1/10$ | $0.374281$ | $1.9130$ | $0.4122$ | $8$ |

Here $\xi_{\mathrm{eff}} = 1 / \log(d_4 / d_5)$ with $d_n = E_0(2n\text{ cells}) - 2E_0(n\text{ cells})$
measures the spatial energy-defect decay length in cells. At $\delta = 1/5$, the odd seed
$n_0 = 5$ cells ($10$ and $20$ spins) enters the admissible region numerically near $\beta_0 = 8.5$
(yielding $\approx 0.086$), while the next even seeds $n_0 = 6, 6, 8$ cells would require $24$- and
$32$-spin doubled rings.

Finally, `dimer_bounds` implements an exact, arbitrary-$`n`$ operator-norm perturbation sandwich
around the decoupled dimer limit $\delta = 1$ (closing a 512-spin ring seed at $\delta = 0.999999$
in the unit tests in $O(\log n)$ operations); because unlinked norm bounds lose bulk free-energy
cancellation, their extensive error $O(\beta n (1-\delta))$ is too large to certify
$\delta \le 1/2$ (recorded per $\delta$ in `results/summary.json`; see [`METHOD.md`](METHOD.md)).

## Reproduce

From the repository root, in the core environment, with `export OPENBLAS_NUM_THREADS=1` (wall
times on an Apple M-series laptop that was heavily loaded by other jobs; CPU times in brackets):

```sh
julia --project=code/julia -e 'using Pkg; Pkg.instantiate()'
julia --project=code/julia code/models/dimerized_heisenberg/check.jl         # 1.5 min [28 s]
julia --project=code/julia code/models/dimerized_heisenberg/check.jl --full  # 4.2 min [2.4 min]
julia --project=code/julia code/models/dimerized_heisenberg/test/runtests.jl # 2.5 min [58 s]
```

- **Fast check (`check.jl`):** re-derives the exact RP audit; for every $\delta$ checks the frame
  identity, the model-local source hashes, the SHA-256 and summary row of every certificate, and
  runs the core `check_point` (the certificate's model, RP data, seed and params are exactly
  `rotated_model(δ)`, `rotated_rp(δ)`, `SEEDS` and $\delta$, its defects are at most those of the
  Lean point, and the exact RP identity, rational chain and interval bound are re-derived);
  checks every exclusion with the core `check_exclusion` (tiling of $[0, B]$, the tail, every
  window's rejection witness), the dimer perturbation bound, and the bound-to-gap ratio.
- **Full check (`check.jl --full`):** additionally recomputes the verified spectra of the rings of
  $2, 4, 8$ cells from the model and re-derives every seed defect, every exclusion bound, and the
  16-spin gap from them alone.
- **Test suite (`test/runtests.jl`):** the fast check, the exact frame, RP and lattice-reflection
  identities, the verified 4-spin spectrum and $\log Z$ against the closed form, the verified
  8-spin spectrum against the physical-frame `float_ring_spectrum`, the sparse 20-spin
  diagnostic blocks against complete spectra, the Lean points, the dimer bound, and
  corrupted-record rejection.

To regenerate (not needed for checking):

```sh
julia --project=code/julia code/models/dimerized_heisenberg/certify.jl  # 3.7 min [2.4 min], all steps
julia --project=code/julia code/models/dimerized_heisenberg/certify.jl delta 3 10
julia --project=code/julia code/models/dimerized_heisenberg/scripts/scan20.jl 1 5
julia --project=code/julia code/models/dimerized_heisenberg/scripts/scan20.jl summary
julia --project=code/julia code/models/dimerized_heisenberg/scripts/numerical_scan.jl
julia --project=code/julia code/models/dimerized_heisenberg/scripts/summarize_numerics.jl
```

## What is rigorous

- **Rigorous:** the exact frame identity and RP identity; the complete verified spectra
  (core `VerifiedSpectra`: exact charge, momentum, and charge-flip reduction; every block
  enclosed by `verified_eigen`); outward-rounded $\log Z$ enclosures, rational seed defects, and
  the exact rational gap chain (core `GapCertificates`); the 16-spin gap enclosures; the
  continuous-$`\beta`$ exclusions (core `Exclusions`); and the dimer perturbation bounds.
- **Trusted base:** exact rational arithmetic, `IntervalArithmetic` outward rounding, and the IEEE
  binary64 + BLAS rounding model of `CertifiedGaps.VerifiedSpectra.verified_eigen`
  (round-to-nearest, gradual underflow, the $\gamma_{N+8}$ / $\gamma_{2N+8}$ dot-product bounds
  for real blocks and for complex blocks evaluated as four real GEMMs, valid for any summation
  order and FMA, with Frobenius norms). LAPACK eigensolver accuracy is never trusted.
- **Numerical only:** the 20-spin partial Krylov spectra, the `Float64` scans
  (`float_ring_spectrum`), the odd $n_0 = 5$ candidate at $\delta = 1/5$, and
  the $\xi_{\mathrm{eff}} \sim \delta^{-2/3}$ scaling fits.

## Layout

Stored results live in [`data/models/dimerized_heisenberg/`](../../../data/models/dimerized_heisenberg/);
this directory's `results/` is a relative symlink to that canonical location.

| Path | Role |
|---|---|
| `certify.jl` | Generates `results/rp_audit.json`, the core certificates, and `results/summary.json` |
| `check.jl` | Certificate verifier (also loaded by the tests; `--full` recomputes every spectrum) |
| `src/DimerizedHeisenberg.jl` | `rotated_model`, `physical_model`, `rotated_rp`, `frame_identity`, `SEEDS`, dimer bound |
| `src/rp.jl`, `src/reflection_conventions.jl` | Exact RP audit (conjugate form) and basis-obstruction checks |
| `src/Diagnostics.jl` | Sparse 20-spin Krylov diagnostics in the physical $S^z$ basis (numerical only) |
| `scripts/` | 20-spin scan, numerical scan, and summary |
| `test/runtests.jl` | Fast unit and regression tests |
| `METHOD.md` | Frame, RP, lattice reflection, core pipeline, exclusions, and dimer bound |
| `IDEAS.md` | Evaluation of candidate routes and outcome |
| `results/certificate_delta_*_beta_*.json` | Core gap certificates, one per seed |
| `results/summary.json` | Per $\delta$: certificate rows and hashes, exclusions, 16-spin gap, ratio, dimer bound |
| `results/rp_audit.json` | Exact rational RP audit at all six $\delta$ values (no seeds) |
| `results/numerical_summary.json`, `results/diagnostic20_summary.json` | Exploratory `Float64` summaries |

The model has no Julia environment of its own; it runs in the core's (`--project=code/julia`). Every
result records (`Results.model_metadata`) the SHA-256 hashes of its model-local sources
(`metadata.source_sha256`: `src/`, `certify.jl`, `check.jl`), which `check.jl` enforces with
`Results.verify_model_metadata`, and of the shared core including its environment
(`metadata.core_sha256`), which is provenance only; `paper/` and `METHOD.md` are not hashed. There
are no hash-pinned layout exceptions.

## History

PR #23 re-certified the six-point $\delta$ grid (exact RP factors, a local two-pass rounding
bound, continuous-$`\beta`$ exclusion covers establishing the grid minimum $\delta = 1/4$, the dimer
perturbation bound, 20-spin diagnostics), superseding the first pass whose bounds
`0.532658, 0.172916, 0.057967` the Lean statements still reflect. Earlier results were produced
under the core's former name `ThermalMPOGap` from `certificates/dimerized_heisenberg/`; later,
shared helpers moved into the core (see the git history).

The certificates are now produced by the shared core pipeline and the model's own Julia
environment is gone: the model is a core `RingModel` in the sublattice-rotated frame with charge
$z_a - z_b$ and charge flip, and `certify_gap` replaces the local block builder, the local
enclosure (component-maximum $\gamma_{4d+16}$ witnesses), the local $\log Z$, and the local
record format; the six `seed_delta_*.json.gz` files and `certificate.json` were replaced by four
core certificates and `summary.json`. Old → new, at the same seeds: $\varepsilon_{\mathrm{s}}$ and
$\varepsilon_\tau$ decreased at every seed (e.g. $\delta = 1/2$:
$0.011267, 0.012008 \to 0.0112652, 0.0120046$) and the interval bound rose (
$0.534756, 0.173197, 0.058225 \to 0.534794, 0.173204, 0.058232$; survey $0.533681 \to 0.533722$); the exact rational chain now stored as
well gives $0.534267, 0.173082, 0.058177$ (survey $0.533215$). The exclusion covers have the
same box counts. The 20-spin diagnostic now takes its 8/16-spin ground energies from the
`Float64` solver instead of the old seed files.

The model now runs through the core drivers: `certify_point` / `check_point` (with `params`
and the Lean seeds from `lean_points`), `Exclusions` in place of the local all-$`\beta`$ covers
(windows from $\beta = 0$ and a ground-energy-mismatch tail replace the
$P_\tau(2n_0, 1) \le 1/2$ argument and the $2^{4n_0}$ tail), `float_ring_spectrum` in place of
the local `Float64` sector solver, and `swap_matrix` / `is_frame_change`. The rotated model gained
the lattice reflection (cell swap times spin flip), so every block at $0 < k < \pi$ is real; all
certificates were regenerated and every defect decreased and every bound rose (e.g.
$\delta = 1/4$: $\varepsilon_{\mathrm{s}}$ $0.165676942264 \to 0.165676941955$, interval bound
$0.0582321291 \to 0.0582321303$).
