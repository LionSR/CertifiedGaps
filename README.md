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
spin chains. A computer certifies a finite seed calculation in Julia, and Lean 4 proves the
doubling step to all larger even lengths.

The method, the **spacetime bootstrap**, generalizes the two-purity doubling argument from OpenAI's
October 2026 spin-1 Haldane gap preprints for
[periodic rings](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026)
and
[open chains](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026).
See the [website](https://sirui-lu.com/CertifiedGaps/) and
[interactive explainer](https://sirui-lu.com/CertifiedGaps/explainer/) for an overview.

## Models

All values are rounded down. The [leaderboard](https://sirui-lu.com/CertifiedGaps/leaderboard.html)
and each model's `README.md` list every parameter point, seed, and Lean 4 theorem.

| Model | Points | Largest certified gap | Lengths | Directory |
|---|---:|---|---|---|
| Spin-1 Heisenberg chain | 6 | $0.3350$ | even $L \ge 12$ | [`models/spin1`](models/spin1/) |
| Spin-1 Heisenberg open chain | 1 | $0.0058$, at $h = 3/5$ | odd $n \ge 1921$ | [`lean/.../OpenChain.lean`](lean/CertifiedGaps/Heisenberg/OpenChain.lean) |
| Transverse- and longitudinal-field Ising chain | 30 | $4.5206$, at $g = 2, h = 1$ | even $L \ge 6$ | [`models/tfim_longitudinal`](models/tfim_longitudinal/) |
| Spin-½ Heisenberg chain in a staggered field | 5 | $1.1489$, at $h = 1$ | even $L \ge 8$ | [`models/heisenberg_staggered_field`](models/heisenberg_staggered_field/) |
| Spin-½ dimerized Heisenberg chain | 3 | $0.5349$, at $\delta = 1/2$ | even $L \ge 4$ | [`models/dimerized_heisenberg`](models/dimerized_heisenberg/) |
| Spin-1 dimerized Heisenberg chain | 2 | $0.5284$, at $\delta = 9/10$ | even $L \ge 2$ | [`models/spin1_dimerized`](models/spin1_dimerized/) |
| Spin-1 Heisenberg chain with single-ion anisotropy | 4 | $0.8995$, at $D = 4$ | even $L \ge 12$ | [`models/spin1_largeD`](models/spin1_largeD/) |
| Two-leg spin-½ Heisenberg ladder | 5 | $1.1027$, at $J_\perp = 4$ | even $L \ge 4$ | [`models/heisenberg_ladder`](models/heisenberg_ladder/) |
| Two-leg transverse-field Ising ladder | 1 | $0.6366$, at $g = 3$ | even $L \ge 4$ | [`models/two_leg_tfim`](models/two_leg_tfim/) |

Pilot studies without a gap certificate: [`models/blbq_haldane`](models/blbq_haldane/)
(bilinear–biquadratic chain), [`models/spin2`](models/spin2/) (spin-2 chain), and
[`models/aklt`](models/aklt/) (honeycomb AKLT).

## Layout

| Directory | Content |
|---|---|
| [`julia/`](julia/) | `CertifiedGaps` core package (recurrence, exact/interval arithmetic, verified spectra) |
| [`models/`](models/) | Model certificates (`certify.jl`, `check.jl`, `results/`) |
| [`lean/`](lean/) | Lean 4 formalization (`CertifiedGaps`) and kernel-checked seeds (`SeedKernel`) |
| [`leaderboard/`](leaderboard/) | Lean record table (`Leaderboard`) and `leaderboard_audit` |
| [`blueprint/`](blueprint/), [`home_page/`](home_page/) | Lean blueprint and Jekyll website |
| [`explainer/`](explainer/) | Interactive visual explainer |
| [`tools/`](tools/) | Upstream data fetcher and leaderboard validator |
| [`notes/`](notes/) | LaTeX notes |

See [`CONTRIBUTING.md`](CONTRIBUTING.md) to submit a record and [`AGENTS.md`](AGENTS.md) for the
certificate specification.

## Reproduce

Requirements: Julia 1.12 and `elan` (Lean 4).

```sh
# Julia core tests
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Check a model certificate (add --full to recompute spectra)
OPENBLAS_NUM_THREADS=1 julia --project=julia models/heisenberg_ladder/check.jl

# Check the spin-1 certificate
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1 models/spin1/check.jl

# Build Lean 4 proofs and validate the leaderboard
lake exe cache get && lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
julia --project=tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/leaderboard tools/leaderboard/check_records.jl \
  --lean-records home_page/_data/lean_records.json
```
