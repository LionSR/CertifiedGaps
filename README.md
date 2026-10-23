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
October 2026 spin-1 Haldane gap preprints
([result 268](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/CONTENTS.md?plain=1#L6747) of the OpenAI Math Release) for
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
| Spin-1 Heisenberg chain | 2 | $0.0729$ | even $L \ge 36$ | [`code/models/spin1`](code/models/spin1/) |
| Spin-1 Heisenberg open chain | 1 | $0.0058$, at $h = 3/5$ | odd $n \ge 1921$ | [`code/lean/.../OpenChain.lean`](code/lean/CertifiedGaps/Heisenberg/OpenChain.lean) |
| Transverse- and longitudinal-field Ising chain | 30 | $4.5202$, at $g = 3/4, h = 1$ | even $L \ge 6$ | [`code/models/tfim_longitudinal`](code/models/tfim_longitudinal/) |
| Spin-½ Heisenberg chain in a staggered field | 5 | $1.1489$, at $h = 1$ | even $L \ge 8$ | [`code/models/heisenberg_staggered_field`](code/models/heisenberg_staggered_field/) |
| Spin-½ dimerized Heisenberg chain | 3 | $0.5347$, at $\delta = 1/2$ | even $L \ge 4$ | [`code/models/dimerized_heisenberg`](code/models/dimerized_heisenberg/) |
| Spin-1 dimerized Heisenberg chain | 2 | $0.5288$, at $\delta = 9/10$ | even $L \ge 2$ | [`code/models/spin1_dimerized`](code/models/spin1_dimerized/) |
| Spin-1 Heisenberg chain with single-ion anisotropy | 4 | $0.8994$, at $D = 4$ | even $L \ge 6$ | [`code/models/spin1_largeD`](code/models/spin1_largeD/) |
| Two-leg spin-½ Heisenberg ladder | 5 | $1.1024$, at $J_\perp = 4$ | even $L \ge 4$ | [`code/models/heisenberg_ladder`](code/models/heisenberg_ladder/) |
| Two-leg transverse-field Ising ladder | 3 | $2.5466$, at $g = 4$ | even $L \ge 4$ | [`code/models/two_leg_tfim`](code/models/two_leg_tfim/) |

Pilot studies without a gap certificate: [`code/models/blbq_haldane`](code/models/blbq_haldane/)
(bilinear–biquadratic chain), [`code/models/spin2`](code/models/spin2/) (spin-2 chain), and
[`code/models/aklt`](code/models/aklt/) (honeycomb AKLT).

## Layout

The folders follow the
[Agentic Publication Protocol](https://github.com/LionSR/AgenticPublicationProtocol/blob/main/PROTOCOL.md)
layout. This is a repository organization change, not a validated APP release.

| Directory | Content |
|---|---|
| [`paper/`](paper/) | Main diagrammatic Haldane-gap note, companion notes, and PDFs |
| [`code/`](code/) | Julia package, model generators/checkers, Lean proofs, leaderboard, and tools |
| [`data/`](data/) | Stored certificates, witnesses, numerical summaries, and external-data instructions |
| [`environment/`](environment/) | Setup commands and canonical dependency specifications |
| [`supplementary/`](supplementary/) | Interactive visual explainer and narrated-film instructions |
| [`blueprint/`](blueprint/), [`home_page/`](home_page/) | Lean blueprint and Jekyll website |

The paper entry point is
[*A Diagrammatic Proof of the Haldane Gap via Spacetime Tensor Networks*](paper/thermal-mpo-gap-note.pdf).
The root Lake files, `blueprint/`, and `home_page/` remain where the Lean documentation build
expects them. Each model's `results/` is a relative symlink to its canonical data directory;
certificate bytes and model-local source hashes are preserved.

See [`CONTRIBUTING.md`](CONTRIBUTING.md) to submit a record and [`AGENTS.md`](AGENTS.md) for the
certificate specification.

## Reproduce

Requirements: Julia 1.12 and `elan` (Lean 4). See
[`environment/README.md`](environment/README.md) for setup and dependency locations.

```sh
# Julia core tests
julia --project=code/julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Check a model certificate (add --full to recompute spectra)
OPENBLAS_NUM_THREADS=1 julia --project=code/julia code/models/heisenberg_ladder/check.jl

# Check the spin-1 certificate
julia code/tools/fetch_upstream.jl
julia --project=code/models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=code/models/spin1 code/models/spin1/check.jl

# Build Lean 4 proofs and validate the leaderboard
lake exe cache get && lake build
lake exe leaderboard_audit home_page/_data/lean_records.json
julia --project=code/tools/leaderboard -e 'using Pkg; Pkg.instantiate()'
julia --project=code/tools/leaderboard code/tools/leaderboard/check_records.jl \
  --lean-records home_page/_data/lean_records.json
```
