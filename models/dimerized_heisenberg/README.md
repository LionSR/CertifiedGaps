# Dimerized spin-1/2 Heisenberg chain

**Certified at $\delta = 1/2, 3/10,$ and $1/4$:** for every even chain length $L \ge 4$ two-spin
unit cells ($N = 2L = 8, 12, 16, \dots$ spins, i.e., $N$ divisible by $4$), the ground state is
unique and the spectral gap satisfies $\Delta_L \ge 0.534756, 0.173197, 0.058225$. On the six-point
parameter grid $\delta \in \{1/2, 3/10, 1/4, 1/5, 3/20, 1/10\}$, $\delta = 1/4$ is the smallest
dimerization certifiable by the two-purity doubling recurrence from rings of at most $16$ spins: at
$\delta \in \{1/5, 3/20, 1/10\}$, the checker rigorously proves via continuous-$\beta$ interval
covers that **every** $\beta_0 > 0$ fails to enter the admissible region at both eligible even seed
lengths $n_0 = 2$ and $n_0 = 4$ cells. The RP hypothesis is the general internally reflected
(conjugated-factor) form. The Lean statements in
`lean/CertifiedGaps/Models/DimerizedHeisenberg*.lean` take earlier, weaker seeds as hypotheses.

| Dimerization $\delta$ | Seed $(n_0, \beta_0)$ | Certified gap $\Delta_L \ge$ | Lengths covered |
|---|---|---:|---|
| $1/2$ | $(4, 6)$ cells | $0.534756$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $3/10$ | $(4, 13/2)$ cells | $0.173197$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $1/4$ | $(4, 29/4)$ cells | $0.058225$ | every even $L \ge 4$ cells ($N \ge 8$ spins, $4 \mid N$) |
| $1/5, 3/20, 1/10$ | none: all $\beta_0 > 0$ excluded at $n_0 = 2, 4$ | — | no certificate at $\le 16$ spins |

## Model

We consider a periodic ring of $N = 2L$ spin-1/2 sites grouped into $L$ two-spin unit cells
(strong dimers), with alternating bond strengths $1 + \delta$ inside each cell and $1 - \delta$
between adjacent cells:
$$H_L = (1 + \delta) \sum_{j=1}^L \mathbf{S}_{2j-1} \cdot \mathbf{S}_{2j}
      + (1 - \delta) \sum_{j=1}^L \mathbf{S}_{2j} \cdot \mathbf{S}_{2j+1}, \qquad
\mathbf{S}_m = \frac{1}{2}\boldsymbol{\sigma}_m.$$
Lengths $L$ and $n_0$ count two-spin unit cells. Each unit cell contains a single internal bond,
so the palindromic unit-cell condition is automatically satisfied.

## Reflection positivity

Let $s = 1 + \delta$ and $w = 1 - \delta$, and write $A_\alpha = S^\alpha \otimes I$ and
$B_\alpha = I \otimes S^\alpha$ for the two spins inside a unit cell. First, apply the single-site
orthogonal rotation $J = 2i S^y = \begin{psmallmatrix} 0 & 1 \\ -1 & 0 \end{psmallmatrix}$ to the
second spin of every cell ($R = I \otimes J$). Under this local unitary transformation, the
intra-cell on-site operator $D$ and the inter-cell bond operator $h$ take the form
$$D_{\mathrm{rot}} = s(-A_x B_x + A_y B_y - A_z B_z), \qquad
-h_{\mathrm{rot}} = w(B_x \otimes A_x - B_y \otimes A_y + B_z \otimes A_z).$$
Let $F$ be the two-spin swap operator inside a cell and $\Theta(\cdot) = F \overline{(\cdot)} F$
the internal cell reflection composed with entrywise complex conjugation. Then $\Theta^2 = I$,
$\Theta(D_{\mathrm{rot}}) = D_{\mathrm{rot}}$, and $-h_{\mathrm{rot}} = w \sum_\alpha \tilde{C}_\alpha \otimes \Theta(\tilde{C}_\alpha)$
with real factors $\tilde{C}_\alpha \in \{B_x, i B_y, B_z\}$—matching the general reflected
unit-cell criterion of `notes/thermal-mpo-gap-scope.tex` (Appendix E) and `RPRing` in `lean/`.

Transforming further into the unitary cell basis $U_{\mathrm{cell}} = M / \sqrt{2}$ with columns
$\frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$, $\frac{1}{\sqrt{2}}(|00\rangle - |11\rangle)$,
$\frac{1}{\sqrt{2}}(|01\rangle + |10\rangle)$, and $\frac{i}{\sqrt{2}}(|01\rangle - |10\rangle)$,
$$M = \begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & i \\ 0 & 0 & 1 & -i \\ 1 & -1 & 0 & 0 \end{pmatrix}, \qquad M^\dagger M = 2I, \qquad \frac{1}{2} M M^\top = F,$$
diagonalizes the internal reflection $F$ so that $\Theta$ becomes plain entrywise conjugation:
$$-h = \sum_{\alpha \in \{x,y,z\}} C_\alpha \otimes \overline{C_\alpha}, \qquad C_\alpha = \frac{\sqrt{w}}{2} K_\alpha, \qquad o = \frac{s}{4}\operatorname{diag}(-3, 1, 1, 1)$$
(with explicit complex matrices $K_\alpha$ given in [`METHOD.md`](METHOD.md)). The module `src/rp.jl`
verifies these identities in exact `Complex{Rational{BigInt}}` arithmetic at all six $\delta$
values, and each seed certificate stores the model and factors as a core `RingModels` record
checked by `RingModels.verify_rp`.

*Why the internal cell swap $F$ is essential:* Any bond decomposition of the unreflected real form
$\sum_\alpha C_\alpha \otimes C_\alpha$ ($F = I$) must commute with the outer swap $F_{\mathrm{outer}}$
of two whole unit cells. However, the inter-cell bond $h$ links spin 2 of cell $j$ with spin 1 of
cell $j+1$, so $\|h - F_{\mathrm{outer}} h F_{\mathrm{outer}}\|_F^2 = \frac{3}{8}(1 - \delta)^2 > 0$
(or $6(1-\delta)^2$ in Pauli normalization), which is invariant under any common single-cell basis
change (`src/reflection_conventions.jl`). Including the internal cell reflection $F$ resolves this
obstruction cleanly.

## Seeds and certified bounds

At seed length $n_0 = 4$ cells, we compute complete verified spectra of the $4$-cell ($8$-spin) and
$8$-cell ($16$-spin) rings. The upper endpoints of the seed defect intervals are rounded up to a
$10^{-6}$ decimal grid plus one spare grid unit (so that independent replays across different BLAS
implementations verify the exact same stored rational seed defects), and passed to
`CertifiedGaps.certified_gap_bound(εₛ, εₜ, β₀)`:

| $\delta$ | $\beta_0$ | $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(4, \beta_0)$ | $\varepsilon_\tau \ge 1 - P_\tau(8, \beta_0)$ | Certified gap $\Delta_L \ge$ | Enclosed 16-spin gap $\Delta_{16}$ | Bound / $\Delta_{16}$ |
|---|---|---:|---:|---:|---|---:|
| $1/2$ | $6$ | $0.011267$ | $0.012008$ | $0.534756$ | $[1.188650, 1.188653]$ | $0.4499$ |
| $3/10$ | $13/2$ | $0.095116$ | $0.051608$ | $0.173197$ | $[0.811143, 0.811146]$ | $0.2135$ |
| $1/4$ | $29/4$ | $0.165679$ | $0.051038$ | $0.058225$ | $[0.709346, 0.709349]$ | $0.0821$ |

Notice that at $\delta = 1/4$, the spatial defect $\varepsilon_{\mathrm{s}} = 0.165679$ exceeds the
diagonal threshold $\varepsilon_* \approx 0.12731$, but the unequal-defect recurrence closes thanks
to the smaller thermal defect $\varepsilon_\tau = 0.051038$. All exact rational seed defects,
$\log Z$ enclosures, and block eigenvalue enclosures are stored in `results/seed_delta_*.json.gz`
(indexed by `results/certificate.json`). We also certify the seed at the prior survey's
inverse temperature $(\delta, n_0, \beta_0) = (1/2, 4, 25/4)$, obtaining $\Delta_L \ge 0.533681$.
Note that the Lean theorems in `lean/CertifiedGaps/Models/DimerizedHeisenberg*.lean` take the
earlier, slightly looser rational seeds from the first certificate pass (`0.532658`, `0.172916`,
`0.057967`) as hypotheses; the tighter two-pass Cauchy–Schwarz rounding bound in
[`METHOD.md`](METHOD.md) shrinks the eigenvalue radii and yields the stronger bounds above.

### Rigorous all-$\beta$ exclusions at $n_0 = 2, 4$ for smaller $\delta$

With the largest ring capped at $16$ spins ($8$ cells), the only even seed lengths are $n_0 = 2$
and $n_0 = 4$ cells. For $n_0 = 2$ at all six $\delta$ values and for $n_0 = 4$ at
$\delta \in \{1/5, 3/20, 1/10\}$, the function `negative_cover` proves rigorously that
`certified_gap_bound` fails for **every** $\beta_0 > 0$:
1. Monotonicity of $P_\tau(2n_0, \beta)$ excludes $0 < \beta_0 \le 1$ once $P_\tau(2n_0, 1) \le 1/2$.
2. A positive ground-state energy defect $d_E = E_0(2n_0)_{\mathrm{lower}} - 2E_0(n_0)_{\mathrm{upper}} > 0$
   implies $P_{\mathrm{s}}(n_0, \beta) \le 2^{4n_0} e^{-\beta d_E} \le 1/2$ for all $\beta_0 \ge B$.
3. The compact interval $[1, B]$ is covered by $6\text{--}12$ rigorous interval boxes, on each of
   which lower bounds on $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$ are propagated through the
   monotone doubling recurrence until $z \ge 1/2$.

This establishes that no $n_0 \in \{2, 4\}$ seed exists at $\delta \le 1/5$ for this recurrence;
certifying $\delta \le 1/5$ requires larger seed rings ($n_0 \ge 6$ cells, i.e., $\ge 24$ spins in
the doubled ring) or a different method.

### Numerical diagnostics up to 20 spins and dimer perturbation

The scripts `scripts/numerical_scan.jl` (complete `Float64` spectra up to 16 spins,
`results/numerical_summary.json`) and `scripts/scan20.jl` (20-spin rings across all symmetry
sectors using 40 Krylov Ritz values in large blocks, `results/diagnostic20_summary.json`) provide
numerical diagnostics across all six $\delta$ values:

| $\delta$ | $\Delta_{20}$ (numerical) | $\xi_{\mathrm{eff}}$ (cells) | $\xi_{\mathrm{eff}}\,\delta^{2/3}$ | Estimated even $n_0$ by $\delta^{-2/3}$ scaling |
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

Finally, `dimer_bounds` in `src/Certification.jl` implements an exact, arbitrary-$n$ operator-norm
perturbation sandwich around the decoupled dimer limit $\delta = 1$ (closing a 512-spin ring seed at
$\delta = 0.999999$ in the unit tests in $O(\log n)$ operations); however, because unlinked norm
bounds lose bulk free-energy cancellation, their extensive error $O(\beta n (1-\delta))$ is too
large to certify $\delta \le 1/2$ (see [`IDEAS.md`](IDEAS.md) and [`METHOD.md`](METHOD.md)).

## Reproduce

From the repository root, with `export OPENBLAS_NUM_THREADS=1` (run times measured on an Apple
M-series laptop):

```sh
julia --project=models/dimerized_heisenberg -e 'using Pkg; Pkg.instantiate()'
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/check.jl         # ~15 s
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/check.jl --full  # ~3 min
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/test/runtests.jl # ~30 s
```

- **Fast check (`check.jl`):** Re-verifies the exact RP audit and, for every $\delta$, checks the
  stored `RingModels` RP factors, re-derives every symmetry block's residual, orthogonality, and
  eigenvalue-radius witness from its stored component maxima, checks sector dimension completeness
  and conjugate-momentum block reuse, recomputes all $\log Z$ intervals, rational seeds, recurrence
  bounds, 16-spin gap intervals, dimer perturbation bounds, and continuous-$\beta$ exclusion boxes
  and tails, checks the model-local source SHA-256 hashes, and checks that
  `results/certificate.json` indexes the certificates.
- **Full check (`check.jl --full`):** Rebuilds every symmetry block from the model definition,
  re-solves all dense eigenproblems from scratch, and repeats all positive and negative certificate
  checks from the fresh spectra alone without trusting any stored eigenvalue or residual maximum.
- **Test suite (`test/runtests.jl`):** Adds exact rational residual recomputations (all 4-spin
  blocks at six $\delta$ values and all 8-spin blocks at three $\delta$ values), closed-form 4-spin
  spectrum checks, conjugate-block equality tests, sector-solver cross-checks against the core dense
  ED, the dimer perturbation bound, and corrupted-certificate rejection tests.

To regenerate certificates or numerical diagnostics from scratch (not needed for checking):

```sh
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl  # ~3.5 min, all steps
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/certify.jl seeds 3 10
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/scan20.jl 1 5
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/scan20.jl summary
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/independent_review.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/numerical_scan.jl
julia --project=models/dimerized_heisenberg models/dimerized_heisenberg/scripts/summarize_numerics.jl
```

## What is rigorous

- **Rigorous:** the exact rational reflection-positivity identities in the two-spin cell basis;
  complete $(S^z, k, p)$ symmetry-sector coverage; every eigenvalue enclosure via the two-pass
  Cauchy–Schwarz rounding bound and polar-factor/Weyl theorem; outward-rounded $\log Z$ sums, exact
  rational seed defects, and uniform gap lower bounds (including ground-state uniqueness from
  $P_\tau > 1/2$); the 16-spin gap enclosures; the continuous-$\beta$ exclusion covers; and the
  dimer perturbation bounds.
- **Trusted base:** IEEE binary64 round-to-nearest arithmetic, standard single-threaded CPU
  OpenBLAS real `GEMM` (`OpenBLAS_jll 0.3.29+0`, with no fast-math, reduced precision, or Strassen
  algorithms), bounded by explicit $\gamma_{4d+16}$ and gradual-underflow allowances, plus
  `IntervalArithmetic` outward rounding. LAPACK eigensolver accuracy is never trusted; the $10^{-6}$
  rational seed grid absorbs last-bit platform variations during `check.jl --full`.
- **Numerical only:** the initial $\beta_0$ grid scan, 20-spin partial Krylov spectra, the odd
  $n_0 = 5$ candidate at $\delta = 1/5$, and the $\xi_{\mathrm{eff}} \sim \delta^{-2/3}$ scaling fits.

## Layout

| Path | Role |
|---|---|
| `certify.jl` | Generates `results/rp_audit.json`, `seed_delta_*.json.gz`, and `certificate.json` |
| `check.jl` | Certificate verifier (also loaded by the tests; fast mode by default, `--full` re-solves every spectrum) |
| `src/Certification.jl` | Seed verification, continuous-$\beta$ exclusions, dimer perturbation bound, stored/fresh checks |
| `src/verified_spectra.jl` | Interval block construction, two-pass rounding witnesses, Weyl enclosures, and $\log Z$ |
| `src/rp.jl`, `src/reflection_conventions.jl` | Exact rational reflection-positivity verifiers and basis-obstruction checks |
| `src/sectors.jl` | Symmetry orbits and `Float64` sector eigensolver for exploratory scans |
| `src/Diagnostics.jl` | 20-spin sparse Krylov diagnostics and summary generator |
| `src/provenance.jl` | Model-local source-hash check of stored results (via the core `verify_source_hashes`) |
| `scripts/` | 20-spin scan, numerical scan and summary, and strict interval-containment review |
| `test/runtests.jl` | Fast unit and regression test suite |
| `METHOD.md` | Detailed mathematical derivations (RP matrices, enclosure theorem, recurrence, exclusions, dimer bound) |
| `IDEAS.md` | Evaluation of candidate routes and outcome |
| `results/certificate.json` | Master index recording, per $\delta$, the certificate hash, seeds, exclusions, and $\mathrm{gap}_{16}$ ratio |
| `results/seed_delta_*.json.gz` | Complete compressed certificates (all six $\delta$ values, every block, witness, and eigenvalue) |
| `results/rp_audit.json` | Standalone exact rational RP audit at all six $\delta$ values (no seeds) |
| `results/numerical_summary.json`, `results/diagnostic20_summary.json` | Exploratory `Float64` summaries |

Every result records, by repository path, the SHA-256 hashes of its model-local sources
(`metadata.source_sha256`: the producing `src/` files, `certify.jl`, `check.jl`, `Project.toml`,
`Manifest.toml`), which `check.jl` enforces, and of the shared core (`metadata.core_sha256`:
`julia/src`, `julia/Project.toml`, `julia/Manifest.toml`), which is provenance only and not compared
with the checkout; `notes/` and `METHOD.md` are not hashed. Certificates regenerate in minutes;
there are no hash-pinned layout exceptions.

## History

PR #23 re-certified the six-point $\delta$ grid: exact RP factors stored as `RingModels` records,
the two-pass Cauchy–Schwarz rounding bound, conjugate-momentum block reuse, a checker that can
rebuild every eigenvector and residual, continuous-$\beta$ exclusion covers for $n_0 = 2, 4$
establishing the grid minimum $\delta = 1/4$, the arbitrary-$n$ dimer perturbation bound, 20-spin
Krylov diagnostics, and certificates for all six $\delta$ (`seed_delta_*`, replacing the three
`verified_delta_*` files, whose bounds `0.532658, 0.172916, 0.057967` are superseded by
`0.534756, 0.173197, 0.058225`). Earlier results were produced under the core's former name
`ThermalMPOGap` from `certificates/dimerized_heisenberg/`. Removed over time (all in the git
history): the 10 MB review output, uncompressed certificates, run logs, `PR_BODY.md`, `REPORT.md`,
`validate.jl` / `validation.json` (timings now in this README), `review/independent_review.md`, the
per-$\delta$ `Float64` scans `results/numerical_delta_*.json`, and the 20-spin runs
`results/diagnostic20_*.json` (regenerate them with the scripts above).

The source-hash policy change (the checker enforces only model-local sources; core hashes are
provenance; notes and `METHOD.md` are no longer hashed) regenerated `rp_audit.json`, all six
`seed_delta_*.json.gz`, and `certificate.json` with `certify.jl`. A deep comparison of old and new
JSON found every eigenvalue, radius, witness, seed, $\log Z$ enclosure, exclusion box, and bound
identical; only hashes, timestamps, and timings changed. The seed hashes recorded in
`diagnostic20_summary.json` were updated by hand to the regenerated, content-identical seeds.

Later, duplicated helpers moved into the shared core `CertifiedGaps`: `load_record` and
`read_json_text` became `Results.read_result` (`.json` and `.json.gz`), `exact` and `parseq` became
`ExactArithmetic.parse_rational`, the seed files are compressed by `Results.write_result` (same
`gzip -n`), the source-hash check uses `Results.verify_source_hashes`, and the exact spin-½
matrices come from `Spins.exact_spin_half`/`exact_pauli`. Since these edits change hashed sources,
`certify.jl` regenerated `rp_audit.json`, the six seed files, and `certificate.json`; a deep
comparison again found every certified number identical (only hashes, timestamps, and timings
changed), and the seed hashes in `diagnostic20_summary.json` were updated to the regenerated files.
