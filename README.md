# Certified Spectral Gaps of Quantum Spin Chains

[![Check certificates](https://github.com/LionSR/CertifiedGaps/actions/workflows/certificates.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/certificates.yml)
[![Kernel re-check](https://github.com/LionSR/CertifiedGaps/actions/workflows/kernel-check.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/kernel-check.yml)
[![Compile blueprint](https://github.com/LionSR/CertifiedGaps/actions/workflows/blueprint.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/blueprint.yml)
[![Website](https://img.shields.io/badge/Website-CertifiedGaps-0f4c81)](https://sirui-lu.com/CertifiedGaps/)
[![Leaderboard](https://img.shields.io/badge/Leaderboard-52_records-1f7a3f)](https://sirui-lu.com/CertifiedGaps/leaderboard.html)
[![Blueprint](https://img.shields.io/badge/Blueprint-HTML_%7C_PDF-2b5b84)](https://sirui-lu.com/CertifiedGaps/blueprint/)
[![Explainer](https://img.shields.io/badge/Explainer-Interactive-b84a1b)](https://sirui-lu.com/CertifiedGaps/explainer/)
[![Lean 4](https://img.shields.io/badge/Lean_4-v4.35.0--rc3-blue)](https://github.com/leanprover/lean4)
[![Julia](https://img.shields.io/badge/Julia-1.12-9558B2)](https://julialang.org/)

In 1983, Haldane conjectured that the spin-1 antiferromagnetic Heisenberg chain
$`H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1}`$ has a uniform spectral gap
$\Delta_L = E_1(L) - E_0(L) > 0$ in the thermodynamic limit (numerically $\Delta \approx 0.4105$).
In October 2026, OpenAI released two preprints proving uniform gap bounds for the
[periodic spin-1 ring](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026)
and for
[open spin-1 chains with endpoint fields](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026)
by combining finite-scale numerical certificates with a coupled spatial-and-thermal doubling
argument. This repository reformulates that mechanism as a general **spacetime bootstrap** (*two-purity
doubling criterion*) using continuous-time tensor networks, reproduces the spin-1 finite certificates
in exact integer and rational arithmetic, formalizes the analytic proofs in Lean 4, and applies the
method to certify uniform spectral gaps across eight families of one-dimensional quantum spin chains
and ladders.

The **spacetime bootstrap** combines two components:

1. **An analytic criterion (formalized in Lean 4).** Write the partition function
   $`Z_L(\beta) = \mathrm{Tr}\, e^{-\beta H_L}`$ of a periodic ring of $L$ unit cells as
   $`\mathrm{Tr}\, X_\beta^L`$, where $X_\beta$ is the spatial transfer operator of an exact
   continuous-time thermal matrix product operator (MPO). When the nearest-neighbor bond interaction
   is reflection positive, $X_\beta$ is self-adjoint, so $Z_L(\beta)$ is a sum of $L$-th powers of
   real transfer eigenvalues just as it is a sum of exponentials of energy eigenvalues. Consequently,
   both the **spatial purity** $P_{\mathrm{s}}(n, \beta) = Z_{2n}(\beta)/Z_n(\beta)^2$ and the
   **thermal purity** $P_\tau(n, \beta) = Z_n(2\beta)/Z_n(\beta)^2$ contract quadratically toward $1$
   when ring length $n$ and inverse temperature $\beta$ are doubled together. Rigorous upper bounds
   on the seed purity defects $\varepsilon_{\mathrm{s}} \ge 1 - P_{\mathrm{s}}(n_0, \beta_0)$ and
   $\varepsilon_\tau \ge 1 - P_\tau(2n_0, \beta_0)$ at a single even seed scale $(n_0, \beta_0)$ in
   the admissible region ($\max(\varepsilon_{\mathrm{s}}, \varepsilon_\tau) < 0.12731$) propagate
   across all dyadic scales to prove that every even ring length $L \ge n_0$ has a unique ground state
   and an explicit lower bound on $\Delta_L$.
2. **Computer-certified finite seeds.** The seed purity inequalities at $(n_0, \beta_0)$ are certified
   by computer: for the spin-1 Heisenberg chain in exact integer and rational arithmetic (reproducing
   and sharpening the finite certificates of the OpenAI periodic preprint), and for the remaining
   models via verified eigenvalue enclosures of complete symmetry-reduced seed-ring spectra (with one
   seed proved outright inside the Lean 4 kernel).

**Website & resources:** [Home page](https://sirui-lu.com/CertifiedGaps/) ·
[Leaderboard](https://sirui-lu.com/CertifiedGaps/leaderboard.html) ·
[Interactive explainer](https://sirui-lu.com/CertifiedGaps/explainer/) ([open chains](https://sirui-lu.com/CertifiedGaps/explainer/open.html)) ·
[Lean 4 blueprint](https://sirui-lu.com/CertifiedGaps/blueprint/) ([PDF](https://sirui-lu.com/CertifiedGaps/blueprint.pdf)) ·
[Lean API docs](https://sirui-lu.com/CertifiedGaps/docs/) ·
[Contribution guide](CONTRIBUTING.md)

## Spin-1 Heisenberg chain

| statement | status | where |
|---|---|---|
| $\Delta_L > 0.0729$ for every even $L \ge 36$ | seed certificates at $(60, 833/32)$ and $(36, 49/4)$ certified in exact rational arithmetic and independently audited (#49); two-seed doubling chain formalized in Lean 4 (`Leaderboard.Records.SpinOneHeisenberg.temperature_lifting`) | [`models/spin1/experiments/`](models/spin1/experiments/), [`lean/`](lean/), [`leaderboard/`](leaderboard/) |
| $\Delta_L > \log(41.196029)/784 \approx 0.00474$ (and the preprint's rounded bound $(\log 20)/784 \approx 0.00382$) for every even $L \ge 2304$; $\Delta_L > \tfrac{4}{105}\log\tfrac{80}{79} \approx 0.000479$ for every even $L \ge 60$ | seed purity bounds at $(2304, 784)$, $(36, 49/4)$, and $(60, 105/4)$ reproduced in exact integer and rational arithmetic; doubling chains formalized in Lean 4 (`SpinOneHeisenberg`) | [`models/spin1/`](models/spin1/), [`lean/`](lean/), [`leaderboard/`](leaderboard/) |
| open chains with endpoint fields $h = 3/5$: $\Delta_n > (\log 10)/392 \approx 0.00587$ for every odd $n \ge 1921$ | formalized in Lean 4 (`SpinOneHeisenberg.spin_one_boundary_field_gap`), conditional on the four finite seed bounds at $\beta_2 = 49$ from the boundary-field preprint | [`lean/CertifiedGaps/Heisenberg/OpenChain.lean`](lean/CertifiedGaps/Heisenberg/OpenChain.lean), [`explainer/open.html`](explainer/open.html) |
| general reflection-positive ring criterion (bond factorization and multi-site unit cells with internal reflection) | formalized in Lean 4 (`RPRing`, `PurityCriterion`, `Chain`) and implemented in `CertifiedGaps` | [`lean/CertifiedGaps/Ring/`](lean/CertifiedGaps/Ring/), [`julia/`](julia/), [`blueprint/`](blueprint/) |
| smaller seed route for Lean verification: $\Delta_L > 0.039$ for even $L \ge 60$ from twisted traces at even $n \le 12$ using real $\mathrm{SU}(2)$ twists only | exact rational doubling chain, conditional on the finite trace tables and analytic seed lemmas | [`models/lean_seeds/`](models/lean_seeds/) (`seed_design`) |
| kernel-checked spin-1 seed arithmetic (base-3 packed-integer operations, Ritz/residual witnesses for $H_6$, $H_8$, $H_{10}$, and local PSD counting anchors) | kernel-checked arithmetic witnesses and feasibility studies | [`models/lean_seeds/`](models/lean_seeds/), [`lean/`](lean/) (`SeedKernel`) |

## Other models

Each certified bound establishes a unique ground state and a uniform spectral gap lower bound for
every even ring length $L \ge n_0$ (where lengths count unit cells). Bounds in the table below are
rounded down for display; each model's `README.md` gives the exact rational inputs, the reflection-positive
basis, the exact seed defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$, run times, and the verification
scope. [`lean/`](lean/) and [`leaderboard/`](leaderboard/) formalize 51 of the 52 certified records in Lean 4:
- **Unconditional Lean 4 theorem (no `sorry`):** the longitudinal-field Ising chain at $g = 1/4, h = 1/20$ has
  $\Delta_L > \log(13480)/(97/4) \approx 0.3921$ for every even $L \ge 6$ (`SeedKernel.Ising.ising_gap`,
  with the 6- and 12-site seed inequalities proved inside the Lean kernel).
- **Lean 4 theorems from certified seeds:** all other even-seed certified points are proved in Lean 4 from their
  computer-certified seed inequalities (`Record.Seeded`).

| model | seed $(n_0, \beta_0)$ | certified uniform gap bound | directory |
|---|---|---|---|
| spin-½ Heisenberg + staggered field, $h = 1, 3/4, 1/2, 7/20, 1/4$ | $(8, 21/4 \dots 25/4)$ | $1.148, 0.890, 0.594, 0.383, 0.214$ | [`heisenberg_staggered_field`](models/heisenberg_staggered_field/) |
| spin-½ dimerized Heisenberg, $\delta = 1/2, 3/10, 1/4$ (2-site cells; $\delta = 1/4$ is the smallest certifiable grid point with rings $\le 16$ spins) | $(4, 6 \dots 29/4)$ | $0.534, 0.173, 0.058$ | [`dimerized_heisenberg`](models/dimerized_heisenberg/) |
| spin-1 dimerized Heisenberg, $\delta = 9/10, 7/10$ (2-site cells) | $(2, 7/2), (3, 17/4)$ | $0.528, 0.218$ | [`spin1_dimerized`](models/spin1_dimerized/) |
| spin-1 Heisenberg $+ D \sum (S^z)^2$, $D = 4$ and $D \in [3.2059998, 3.206]$ | $(12, 29/10), (6, 381/125)$ | $0.899, 0.0018$ (every even $L \ge 6$) | [`spin1_largeD`](models/spin1_largeD/) |
| two-leg Heisenberg ladder, $J_\perp = 4, 7/2, 3, 5/2, 163069/65536$ (rung cells) | $(4, 2 \dots 3)$ | $1.102, 0.778, 0.454, 0.042, 0.0015$ | [`heisenberg_ladder`](models/heisenberg_ladder/) |
| two-leg transverse-field Ising ladder, $g = 3$ (rung cells) | $(4, 3/2)$ | $0.636$ | [`two_leg_tfim`](models/two_leg_tfim/) |
| Ising chain $-ZZ - gX - hZ$, 30 points $g \in [1/4, 2]$, $h \in [1/20, 1]$ | $(6, 2 \dots 97/4)$ | $0.0763$ to $3.705$ (unconditional in Lean at $g = 1/4, h = 1/20$: $0.392$) | [`tfim_longitudinal`](models/tfim_longitudinal/) |
| spin-1 bilinear–biquadratic, $\theta \in [-7/10, 0]$ | no admissible even seed with $n_0 \le 6$ | negative result; required seed scale estimated | [`blbq_haldane`](models/blbq_haldane/) |
| spin-2 Heisenberg chain | — | certified MPS upper bounds ($N \le 1024$) and certified-MPO study; no gap certificate | [`spin2`](models/spin2/) |
| hexagonal spin-3/2 AKLT ($2+1$ dimensions) | — | numerical pilot of the $`2+1`$D thermal-purity recursion | [`aklt`](models/aklt/) |

All certified models outside `spin1/` share one pipeline in the core (`CertifiedGaps.GapCertificates`):
each model is an exact `RingModel` with exact reflection-positivity data, and `certify_gap` computes
the complete spectra of the seed rings at lengths $n_0$, $2n_0$ by exact symmetry blocks (charge,
momentum, spin flip, on-site parities, cell blocking), encloses every eigenvalue by one residual,
polar-factor and Weyl/Ostrowski bound (`verified_eigen`) whose floating-point rounding error is bounded
a priori under the IEEE binary64 model (including BLAS matrix products), and runs the exact rational
doubling chain; `check_certificate` re-derives it. The bounds in the table are the outward-rounded
interval bounds; each certificate also states the exact rational chain bound.

## Layout

| directory | content |
|---|---|
| [`julia/`](julia/) | `CertifiedGaps`, the shared Julia core package: the doubling recurrence with outward-rounded and exact-rational bounds, exact arithmetic, results I/O, reflection-positivity tests and exact RP certificates, exact diagonalization, verified ring spectra, certified MPO arithmetic (Nemo/Arb), and thermal purification (ITensors) |
| [`models/`](models/) | one directory per model in a standardized layout ([`AGENTS.md`](AGENTS.md), "Directory format"): `README.md`, `IDEAS.md`, certificate JSON, `certify.jl`, `check.jl`, and tests; documented layout exceptions are listed in each README's `## Layout`. `spin1/` reproduces the exact-arithmetic spin-1 certificate; `lean_seeds/` develops kernel-checked seed witnesses |
| [`lean/`](lean/) | Lean 4 + Mathlib formalization of the two-purity criterion and each certified model's gap theorem (with the seeds as hypotheses), plus the opt-in library `SeedKernel` for kernel-checked seed arithmetic |
| [`leaderboard/`](leaderboard/) | Lean library `Leaderboard` packaging every Lean-verified gap record in a uniform `Record.Unconditional` or `Record.Seeded` statement, plus the kernel axiom/parameter auditor `leaderboard_audit` |
| [`blueprint/`](blueprint/), [`home_page/`](home_page/) | the `leanblueprint` blueprint of the Lean formalization (web, PDF, dependency graph) and the project home page and leaderboard (`home_page/_data/records.yml`), built from the root `lakefile.toml` and deployed to GitHub Pages |
| [`explainer/`](explainer/) | interactive visual exposition of the periodic ring and open-chain proofs (two HTML pages and two narrated films; floating-point illustrations only, not part of any certificate) |
| [`tools/`](tools/) | `fetch_upstream.jl`: downloads the OpenAI preprint's evidence files (for `spin1`) or, with `--paper`, its LaTeX sources (for `lean_seeds`) at a pinned commit and verifies their SHA-256 hashes; `leaderboard/`: the validator of `home_page/_data/records.yml` against the certificates and Lean audit |
| [`notes/`](notes/) | working LaTeX notes and diagram macros (`tnnotes.sty`) |

See [`AGENTS.md`](AGENTS.md) for the exact specification of a gap certificate and instructions on adding new models,
and [`CONTRIBUTING.md`](CONTRIBUTING.md) for submitting a record to the leaderboard.

## Reproduce

Requirements: Julia 1.12, `elan` (Lean 4), and `leanblueprint` (with a LaTeX installation) for the
blueprint.

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

# Lean formalization, kernel-checked seeds, and leaderboard audit (from the repository root)
lake exe cache get && lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl \
  --lean-records home_page/_data/lean_records.json

# Blueprint of the Lean formalization
leanblueprint pdf && leanblueprint web
lake exe checkdecls blueprint/lean_decls
```

The root `lake build` compiles `CertifiedGaps`, `SeedKernel`, and `Leaderboard` into `.lake/`
(`lean/` also contains a standalone `lakefile.toml` used by `models/lean_seeds/check.jl`). Each
model directory's `README.md` provides its exact reproduction commands, measured run times, and
trusted arithmetic base.

## Status of the claims

* **Spin-1 Heisenberg chain:** For periodic rings, the seed purity bounds at $(60, 833/32)$,
  $(36, 49/4)$, $(60, 105/4)$, and $(2304, 784)$ are certified in [`models/spin1/`](models/spin1/) in
  exact integer and rational arithmetic, and the analytic deductions of the uniform spectral gap from
  those seeds ($\Delta_L > 0.0729$ for even $L \ge 36$, $\Delta_L > \log(41.196029)/784$ and
  $(\log 20)/784$ for even $L \ge 2304$, and $\Delta_L > \tfrac{4}{105}\log\tfrac{80}{79}$ for even
  $L \ge 60$) are formalized in [`lean/`](lean/). For open chains with endpoint fields $h = 3/5$, the
  analytic deduction of $\Delta_n > (\log 10)/392$ for odd $n \ge 1921$ from the four finite seed bounds
  at $\beta_2 = 49$ of the OpenAI boundary-field preprint is formalized in
  [`lean/CertifiedGaps/Heisenberg/OpenChain.lean`](lean/CertifiedGaps/Heisenberg/OpenChain.lean). In `CertifiedGaps`, each
  seeded theorem takes the finite seed inequalities as named `sorry` lemmas; in `SeedKernel`,
  `SeedKernel.Ising.ising_gap` ([`models/lean_seeds/`](models/lean_seeds/)) proves both the seed
  inequalities and the doubling recurrence inside the Lean kernel for the longitudinal-field Ising chain
  at $g = 1/4$, $h = 1/20$.
* **Other certified models:** The seed enclosures rest on verified block-diagonal eigenvalue enclosures
  whose rounding analysis assumes standard IEEE binary64 arithmetic and conventional BLAS dot-product
  error bounds; each was independently replayed or reviewed (see the model's `README.md` and
  `METHOD.md`).
* **Exploratory numerics:** Quantities labeled as numerical diagnostics (DMRG, Lanczos, unverified
  floating-point exact diagonalization, thermal MPO purification, and finite-size extrapolations) are
  used only for screening and scale estimation, never as part of a gap certificate.
