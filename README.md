# A diagrammatic proof of the Haldane gap

The spin-1 antiferromagnetic Heisenberg chain $H_L = \sum_j \mathbf{S}_j \cdot \mathbf{S}_{j+1}$ has a
spectral gap $\Delta_L = E_1(L) - E_0(L)$ that stays open as the ring grows (Haldane's conjecture;
numerically $\Delta \approx 0.4105$). This repository proves a uniform lower bound on that gap using
spacetime tensor networks drawn as diagrams, certifies the finite numerical inputs by computer, and
applies the same method to a family of one-dimensional quantum spin chains and ladders.

The method, the **spacetime bootstrap** (or *two-purity doubling criterion*), has two parts:

1. **An analytic criterion.** Write the partition function $Z_L(\beta) = \operatorname{Tr} e^{-\beta H_L}$
   of a periodic ring as $\operatorname{Tr} X_\beta^L$, where $X_\beta$ is the spatial transfer operator of an
   exact continuous-time thermal matrix product operator (MPO). When the bond interaction is reflection
   positive, $X_\beta$ is self-adjoint, so $Z_L(\beta)$ is a sum of $L$-th powers of real transfer
   eigenvalues just as it is a sum of exponentials of energy eigenvalues. Consequently, both the
   **spatial purity** $P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$ and the **thermal purity**
   $P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$ improve quadratically when ring length and inverse
   temperature are doubled together. Rigorous lower bounds on $P_{\mathrm{s}}(n_0, \beta_0)$ and
   $P_\tau(2n_0, \beta_0)$ at a single finite scale $(n_0, \beta_0)$ (the *seed*) then propagate across all
   dyadic scales to yield a uniform lower bound on $\Delta_L$ for every even ring length $L \ge n_0$.
2. **Certified finite inputs.** The seed purities are bounded by computer: for the spin-1 Heisenberg
   chain using exact integer and rational arithmetic, reproducing the certificates of OpenAI's
   preprint [*The periodic spin-one Haldane gap*](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026);
   and for the other models using verified eigenvalue enclosures of complete seed-ring spectra.

## Spin-1 Heisenberg chain

| statement | status | where |
|---|---|---|
| $\Delta_L > (\log 20)/784 \approx 0.00382$ for every even $L \ge 2304$ | proved in the main note, given the seed purity bounds at $(n_0, \beta_0) = (2304, 784)$ | [`notes/thermal-mpo-gap-note.pdf`](notes/thermal-mpo-gap-note.pdf) |
| seed purity bounds $P_{\mathrm{s}}(2304, 784) \ge 0.99$ and $P_\tau(4608, 784) \ge 0.99$ | reproduced in exact arithmetic (status `passed`) | [`models/spin1/`](models/spin1/) |
| the two-purity criterion and the spin-1 theorem above, given the two seed bounds | formalized in Lean 4 | [`lean/`](lean/) |
| [`explainer/`](explainer/) | an illustrated, undergraduate-level exposition of the proof (HTML page and two animated films); floating-point illustrations only, not part of any certificate |
| $\Delta_L > 0.0729$ for every even $L \ge 36$ | exact-arithmetic computation combining temperature-lifting with a tighter seed chain; not yet independently reviewed | [`models/spin1/`](models/spin1/) (`experiments/`) |
| open chains with endpoint fields $h = 3/5$: $\Delta_n > (\log 10)/392$ for odd $n \ge 1921$ | proved in Appendices A–C, given four finite estimates from the boundary-field paper | [`notes/thermal-mpo-gap-note.pdf`](notes/thermal-mpo-gap-note.pdf) |
| general scope of the method (bond reflection positivity and multi-site unit cells) | proved in Appendices D–H | [`notes/thermal-mpo-gap-scope.pdf`](notes/thermal-mpo-gap-scope.pdf) |
| diagrammatic exposition of the spin-1 seed certificate at $(2304, 784)$ | Appendices I–M | [`notes/thermal-mpo-gap-certificate.pdf`](notes/thermal-mpo-gap-certificate.pdf) |
| smaller seed route for Lean verification: $\Delta_L > 0.039$ for even $L \ge 60$ from twisted traces at even $n \le 12$ using real $\mathrm{SU}(2)$ twists only | exact rational doubling chain, conditional on the finite trace tables and analytic seed lemmas | [`models/lean_seeds/`](models/lean_seeds/) (`seed_design`) |
| kernel-checked seed arithmetic in Lean (packed-integer operations and counting anchors) | **one unconditional Lean gap theorem**: the longitudinal-field Ising chain at $g = 1/4$, $h = 1/20$, with $\Delta_L > 0.392$ for every even $L \ge 6$ (`SeedKernel.Ising.ising_gap`, no `sorry`); spin-1 seed: feasibility studies | [`models/lean_seeds/`](models/lean_seeds/), [`lean/`](lean/) (`SeedKernel`) |

## Other models

Each certified bound establishes a unique ground state and a uniform spectral gap lower bound for
every even ring length $L \ge n_0$ (where lengths count unit cells). Bounds in the table below are
rounded down for display; each model's `README.md` gives the exact rational inputs, the reflection-positive
basis, the exact seed defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$, run times, and the verification
scope. [`lean/`](lean/) formalizes each certified model's gap theorem with its two seed inequalities as the
only open hypotheses.

| model | seed $(n_0, \beta_0)$ | certified uniform gap bound | directory |
|---|---|---|---|
| spin-½ Heisenberg + staggered field, $h = 1, 3/4, 1/2, 7/20, 1/4$ | $(8, 21/4 \dots 25/4)$ | $1.148, 0.890, 0.594, 0.383, 0.214$ | [`heisenberg_staggered_field`](models/heisenberg_staggered_field/) |
| spin-½ dimerized Heisenberg, $\delta = 1/2, 3/10, 1/4$ (2-site cells; $\delta = 1/4$ is the smallest certifiable grid point with rings $\le 16$ spins) | $(4, 6 \dots 29/4)$ | $0.534, 0.173, 0.058$ | [`dimerized_heisenberg`](models/dimerized_heisenberg/) |
| spin-1 dimerized Heisenberg, $\delta = 9/10, 7/10$ (2-site cells) | $(2, 7/2), (3, 17/4)$ | $0.528, 0.218$ | [`spin1_dimerized`](models/spin1_dimerized/) |
| spin-1 Heisenberg $+ D \sum (S^z)^2$, $D = 4$ and $D \in [3.2059998, 3.206]$ | $(12, 29/10), (6, 381/125)$ | $0.899, 0.0010$ (every even $L \ge 6$) | [`spin1_largeD`](models/spin1_largeD/) |
| two-leg Heisenberg ladder, $J_\perp = 4, 7/2, 3, 5/2$ (rung cells) | $(4, 2 \dots 3)$ | $1.102, 0.778, 0.454, 0.042$ | [`heisenberg_ladder`](models/heisenberg_ladder/) |
| two-leg transverse-field Ising ladder, $g = 3$ (rung cells) | $(4, 3/2)$ | $0.636$ | [`two_leg_tfim`](models/two_leg_tfim/) |
| Ising chain $-ZZ - gX - hZ$, 30 points $g \in [1/4, 2]$, $h \in [1/20, 1]$ | $(6, \text{per point})$ | $\ge 0.0763$ at every point (smallest at $g = 1, h = 1/20$) | [`tfim_longitudinal`](models/tfim_longitudinal/) |
| spin-1 bilinear–biquadratic, $\theta \in [-7/10, 0]$ | no admissible even seed with $n_0 \le 6$ | negative result; required seed scale estimated | [`blbq_haldane`](models/blbq_haldane/) |
| spin-2 Heisenberg chain | — | numerical and certified-MPO investigation; no certificate | [`spin2`](models/spin2/) |
| hexagonal spin-3/2 AKLT ($2+1$ dimensions) | — | numerical pilot of the $2+1$D thermal-purity recursion | [`aklt`](models/aklt/) |

For all certified models outside `spin1/`, the seed defects are computed from complete spectra of the
seed rings at lengths $n_0$, $2n_0$ (at $\beta_0$ and $2\beta_0$), where every eigenvalue is enclosed by a
residual and Weyl/Ostrowski bound whose floating-point rounding error is bounded a priori under the
IEEE binary64 model (including BLAS matrix products).

## Layout

| directory | content |
|---|---|
| [`notes/`](notes/) | LaTeX notes and compiled PDFs, the shared diagram package `tnnotes.sty`, and canonical notation and diagram conventions |
| [`julia/`](julia/) | `CertifiedGaps`, the shared Julia core package: the doubling recurrence with outward-rounded and exact-rational bounds, exact arithmetic, results I/O, reflection-positivity tests and exact RP certificates, exact diagonalization, verified ring spectra, certified MPO arithmetic (Nemo/Arb), and thermal purification (ITensors) |
| [`models/`](models/) | one directory per model in a standardized layout ([`AGENTS.md`](AGENTS.md), "Directory format"): `README.md`, `IDEAS.md`, certificate JSON, `certify.jl`, `check.jl`, and tests; documented layout exceptions are listed in each README's `## Layout`. `spin1/` reproduces the exact-arithmetic spin-1 certificate; `lean_seeds/` develops kernel-checked seed witnesses |
| [`lean/`](lean/) | Lean 4 + Mathlib formalization of the two-purity criterion and each certified model's gap theorem (with the seeds as hypotheses), plus the opt-in library `SeedKernel` for kernel-checked seed arithmetic |
| [`tools/`](tools/) | `fetch_upstream.jl`: downloads the OpenAI preprint's evidence files (for `spin1`) or, with `--paper`, its LaTeX sources (for `lean_seeds`) at a pinned commit and verifies their SHA-256 hashes |

See [`AGENTS.md`](AGENTS.md) for the exact specification of a gap certificate and instructions on adding new models.

## Reproduce

Requirements: Julia 1.12, `elan` (Lean 4), and a LaTeX installation with REVTeX 4.2 and `latexmk`.

```sh
# Julia core and its tests
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# spin-1 certificates: fetch the preprint's evidence (pinned commit, hash-checked), then test
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1 models/spin1/test/runtests.jl
julia --project=models/spin1 models/spin1/check.jl                      # stored-output audit
julia --project=models/spin1 --threads=4 models/spin1/check.jl --full   # fresh run, ≈20–45 min

# any other model runs in the core environment: fast checker (add --full for the full
# spectral replay), and its tests
OPENBLAS_NUM_THREADS=1 julia --project=julia models/heisenberg_ladder/check.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia models/heisenberg_ladder/test/runtests.jl

# Lean formalization
cd lean && lake exe cache get && lake build

# LaTeX notes
latexmk -pdf -outdir=work/latex-build notes/thermal-mpo-gap-note.tex
```

Each model directory's `README.md` provides its exact reproduction commands, measured run times, and
trusted arithmetic base.

## Status of the claims

* **Spin-1 Heisenberg chain:** The analytic proof is complete on paper given the seed purity bounds at
  $(n_0, \beta_0) = (2304, 784)$. Those bounds are reproduced in [`models/spin1/`](models/spin1/) in
  exact integer and rational arithmetic, and the deduction of the uniform spectral gap from the two seed
  bounds is formalized in [`lean/`](lean/). The finite seed computations themselves are not
  formalized in `CertifiedGaps`: there each model's two seed bounds appear as named `sorry` lemmas.
  Kernel-checked seed arithmetic proves one seed outright: `SeedKernel.Ising.ising_gap`
  ([`models/lean_seeds/`](models/lean_seeds/)) is an unconditional, fully kernel-checked gap theorem
  for the longitudinal-field Ising chain at $g = 1/4$, $h = 1/20$.
* **Other certified models:** The seed enclosures rest on verified block-diagonal eigenvalue enclosures
  whose rounding analysis assumes standard IEEE binary64 arithmetic and conventional BLAS dot-product
  error bounds; each was independently replayed or reviewed (see the model's `README.md` and
  `METHOD.md`).
* **Exploratory numerics:** Quantities labeled as numerical diagnostics (DMRG, Lanczos, unverified
  floating-point exact diagonalization, thermal MPO purification, and finite-size extrapolations) are
  used only for screening and scale estimation, never as part of a gap certificate.
