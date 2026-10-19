# Certified Spectral Gaps of Quantum Spin Chains

[![Check certificates](https://github.com/LionSR/CertifiedGaps/actions/workflows/certificates.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/certificates.yml)
[![Kernel re-check](https://github.com/LionSR/CertifiedGaps/actions/workflows/kernel-check.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/kernel-check.yml)
[![Compile blueprint](https://github.com/LionSR/CertifiedGaps/actions/workflows/blueprint.yml/badge.svg)](https://github.com/LionSR/CertifiedGaps/actions/workflows/blueprint.yml)

[Website](https://sirui-lu.com/CertifiedGaps/) ·
[Leaderboard](https://sirui-lu.com/CertifiedGaps/leaderboard.html) ·
[Interactive explainer](https://sirui-lu.com/CertifiedGaps/explainer/) ([open chains](https://sirui-lu.com/CertifiedGaps/explainer/open.html)) ·
[Lean blueprint](https://sirui-lu.com/CertifiedGaps/blueprint/) ([PDF](https://sirui-lu.com/CertifiedGaps/blueprint.pdf)) ·
[Lean API docs](https://sirui-lu.com/CertifiedGaps/docs/) ·
[Contributing](CONTRIBUTING.md)

This repository proves lower bounds on the spectral gap $\Delta_L = E_1(L) - E_0(L)$ of quantum
spin chains. Each bound holds for all even ring lengths $L \ge L_0$. A computer checks one finite
calculation for each bound, in Julia. Lean 4 checks the remaining proof.

OpenAI proved the spin-1 Haldane gap in two preprints in October 2026, for
[periodic rings](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026)
and for
[open chains](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026).
Their proof doubles the ring length and the inverse temperature together. This repository states
that argument as a general criterion, the **spacetime bootstrap**, and applies it to more models.

Let $`Z_L(\beta) = \mathrm{Tr}\, e^{-\beta H_L}`$ be the partition function of a ring of $L$ unit
cells. The method uses two ratios:

```math
P_{\mathrm S}(n, \beta) = \frac{Z_{2n}(\beta)}{Z_n(\beta)^2}, \qquad
P_{\mathrm T}(n, \beta) = \frac{Z_n(2\beta)}{Z_n(\beta)^2}.
```

1. Show that the interaction between neighboring cells is reflection positive.
2. At one seed $(n_0, \beta_0)$, calculate upper bounds on the defects
   $p = 1 - P_{\mathrm S}(n_0, \beta_0)$ and $q = 1 - P_{\mathrm T}(2n_0, \beta_0)$. The method
   works if $p$ and $q$ are both less than $`\varepsilon_* \approx 0.127`$.
3. Double $n$ and $\beta$ together. Each doubling makes $p$ and $q$ smaller. The result is a
   unique ground state and a lower bound on the gap for all even $L \ge n_0$.

## Models

The table rounds bounds down. The [leaderboard](https://sirui-lu.com/CertifiedGaps/leaderboard.html)
and the `README.md` of each model give the exact values, seeds and trusted arithmetic. Each bound
with an even seed has a Lean 4 theorem in [`leaderboard/`](leaderboard/).

| model | seed $(n_0, \beta_0)$ | certified gap $\Delta_L \ge$ | lengths | directory |
|---|---|---|---|---|
| spin-1 Heisenberg chain | $(60, 833/32)$ and $(36, 49/4)$ | $0.0729$ | even $L \ge 36$ | [`models/spin1`](models/spin1/) |
| spin-1 Heisenberg chain (preprint seed) | $(2304, 784)$ from $(36, 49/4)$ | $0.00474$ ($(\log 20)/784 \approx 0.00382$) | even $L \ge 2304$ | [`models/spin1`](models/spin1/) |
| spin-1 open chain, endpoint fields $h = 3/5$ | $\beta_2 = 49$ (Lean theorem from preprint seeds) | $(\log 10)/392 \approx 0.00587$ | odd $n \ge 1921$ | [`lean/.../OpenChain.lean`](lean/CertifiedGaps/Heisenberg/OpenChain.lean) |
| Ising chain $-ZZ - gX - hZ$, 30 points ($g \in [1/4, 2]$, $h \in [1/20, 1]$) | $(6, 2 \dots 97/4)$ | $0.0763$ to $4.520$ ($0.392$ unconditional in Lean at $g = 1/4, h = 1/20$) | even $L \ge 6$ | [`models/tfim_longitudinal`](models/tfim_longitudinal/) |
| spin-½ Heisenberg + staggered field ($h = 1, 3/4, 1/2, 7/20, 1/4$) | $(8, 21/4 \dots 25/4)$ | $1.148, 0.890, 0.594, 0.383, 0.214$ | even $L \ge 8$ | [`models/heisenberg_staggered_field`](models/heisenberg_staggered_field/) |
| spin-½ dimerized Heisenberg ($\delta = 1/2, 3/10, 1/4$, 2-site cells) | $(4, 6 \dots 29/4)$ | $0.534, 0.173, 0.058$ | even $L \ge 4$ | [`models/dimerized_heisenberg`](models/dimerized_heisenberg/) |
| spin-1 dimerized Heisenberg ($\delta = 9/10, 7/10$, 2-site cells) | $(2, 7/2), (3, 17/4)$ | $0.528, 0.218$ | even $L \ge 2, 4$ | [`models/spin1_dimerized`](models/spin1_dimerized/) |
| spin-1 Heisenberg $+ D \sum (S^z)^2$ ($D = 4$, $1603/500$, and intervals near $3.206$ and $4$) | $(12, 29/10), (6, 381/125), (6, 23/10)$ | $0.899, 0.00275, 0.00182, 0.00259$ | even $L \ge 6$ | [`models/spin1_largeD`](models/spin1_largeD/) |
| two-leg Heisenberg ladder ($J_\perp = 4, 7/2, 3, 5/2, 163069/65536$) | $(4, 2 \dots 3)$ | $1.102, 0.778, 0.454, 0.042, 0.0015$ | even $L \ge 4$ | [`models/heisenberg_ladder`](models/heisenberg_ladder/) |
| two-leg transverse-field Ising ladder ($g = 3$) | $(4, 3/2)$ | $0.636$ | even $L \ge 4$ | [`models/two_leg_tfim`](models/two_leg_tfim/) |
| spin-1 bilinear–biquadratic ($\theta \in [-7/10, 0]$) | no seed with $n_0 \le 6$ | negative result; seed scale estimated | — | [`models/blbq_haldane`](models/blbq_haldane/) |
| spin-2 Heisenberg chain | — | certified MPS upper bounds ($N \le 1024$) and MPO study | — | [`models/spin2`](models/spin2/) |
| honeycomb spin-3/2 AKLT ($2+1$D) | — | finite-torus pilot of the $`2+1`$D purity recursion | — | [`models/aklt`](models/aklt/) |

## Status

* Lean 4 proves step 3 for each bound with an even seed. The Lean theorems use the seed inequalities as
  hypotheses. For the Ising chain at $g = 1/4$, $h = 1/20$, the Lean kernel also proves the seed
  (`SeedKernel.Ising.ising_gap`).
* The spin-1 seeds use exact integer and rational arithmetic on the finite data of the periodic
  preprint. An independent audit (#49) repeats them in exact rationals and in Arb ball arithmetic.
* The seeds of the other models use verified eigenvalue enclosures. Their rounding analysis
  assumes IEEE binary64 arithmetic and the standard error bounds for BLAS dot products.
* The open-chain theorem uses the four finite bounds at $\beta_2 = 49$ of the boundary-field
  preprint as hypotheses. This repository does not certify those four bounds yet.

## Layout

| directory | content |
|---|---|
| [`julia/`](julia/) | `CertifiedGaps` core package: doubling recurrence, exact and interval arithmetic, RP certificates, verified symmetry-sector spectra (`VerifiedSpectra`), and certified MPO tools |
| [`models/`](models/) | one directory per model (`README.md`, `IDEAS.md`, `certify.jl`, `check.jl`, `results/`); `lean_seeds/` generates kernel-checked seed witnesses for `lean/SeedKernel` |
| [`lean/`](lean/) | Lean 4 + Mathlib formalization of the doubling criterion (`CertifiedGaps`) and kernel-checked seed arithmetic (`SeedKernel`) |
| [`leaderboard/`](leaderboard/) | Lean library `Leaderboard` packaging every Lean-verified gap record (`Record.Unconditional`, `Record.Seeded`) and the `leaderboard_audit` executable |
| [`blueprint/`](blueprint/), [`home_page/`](home_page/) | `leanblueprint` dependency graph/PDF and Jekyll website (`home_page/_data/records.yml`) |
| [`explainer/`](explainer/) | interactive visual exposition of the periodic ring and open-chain proofs |
| [`tools/`](tools/) | `fetch_upstream.jl` (downloads hash-pinned preprint data) and `tools/leaderboard/` (`check_records.jl` validator) |
| [`notes/`](notes/) | LaTeX notes on the criterion, scope, and tensor-network diagrams |

See [`AGENTS.md`](AGENTS.md) for the specification of a gap certificate and [`CONTRIBUTING.md`](CONTRIBUTING.md)
for submitting a record.

## Reproduce

Requirements: Julia 1.12 and `elan` (Lean 4).

```sh
# Julia core tests
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Spin-1 certificate (fetches hash-pinned upstream evidence)
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1 models/spin1/check.jl                      # fast audit (add --full to recompute from scratch)

# Any other model runs in the core environment (add --full for full spectral replay)
OPENBLAS_NUM_THREADS=1 julia --project=julia models/heisenberg_ladder/check.jl
OPENBLAS_NUM_THREADS=1 julia --project=julia models/heisenberg_ladder/test/runtests.jl

# Lean 4 proofs and leaderboard validation
lake exe cache get && lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl \
  --lean-records home_page/_data/lean_records.json
```
