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

## Paper

**Bridging Gaps with Tensor Networks: From the Haldane Gap and Beyond**, Sirui Lu (2026).
[PDF](paper/certified-gap.pdf) · [LaTeX source](paper/certified-gap.tex) ·
[PDF on the website](https://sirui-lu.com/CertifiedGaps/notes/certified-gap.pdf). The arXiv link
goes here once the listing appears. [`paper/README.md`](paper/README.md) has the build command
and the earlier companion notes.

### Talk to this paper

The repository follows the
[Agentic Publication Protocol](https://github.com/LionSR/AgenticPublicationProtocol). Clone it and
open it in an AI coding agent that reads [`AGENTS.md`](AGENTS.md); in Claude Code,
`/load-paper https://github.com/LionSR/CertifiedGaps`.

## Models

All values are rounded down. The [leaderboard](https://sirui-lu.com/CertifiedGaps/leaderboard.html)
and each model's `README.md` list every parameter point, seed, and Lean 4 theorem.

| Model | Points | Largest certified gap | Lengths | Directory |
|---|---:|---|---|---|
| Spin-1 Heisenberg chain | 3 | $0.0729$ | even $L \ge 36$ | [`code/models/spin1`](code/models/spin1/) |
| Spin-1 Heisenberg open chain | 1 | $0.0058$, at $h = 3/5$ | odd $n \ge 1921$ | [`code/lean/.../OpenChain.lean`](code/lean/CertifiedGaps/Heisenberg/OpenChain.lean) |
| Transverse- and longitudinal-field Ising chain | 30 | $4.5202$, at $g = 3/4, h = 1$ | even $L \ge 6$ | [`code/models/tfim_longitudinal`](code/models/tfim_longitudinal/) |
| Transverse-field Ising chain (free fermions, $h = 0$) | 3 | $1.9189$, at $g = 2$ | even $L \ge 128$ (at $g = 2$) | [`code/models/tfim_free_fermions`](code/models/tfim_free_fermions/) |
| XY chain in a field ($J_x J_y < 0$, free fermions) | 3 | $1.8649$, at $J_x = -1, J_y = 1, h = 1$ | even $L \ge 128$ | [`code/models/xy_field`](code/models/xy_field/) |
| Cluster chain in a transverse field (free fermions) | 1 | $0.9515$, at $g = 1/2$ | even $L \ge 128$ cells | [`code/models/cluster_chain`](code/models/cluster_chain/) |
| Spin-½ Heisenberg chain in a staggered field | 5 | $1.1489$, at $h = 1$ | even $L \ge 8$ | [`code/models/heisenberg_staggered_field`](code/models/heisenberg_staggered_field/) |
| Spin-½ dimerized Heisenberg chain | 4 | $0.6572$, at $\delta = 1/2$ | even $L \ge 4$ | [`code/models/dimerized_heisenberg`](code/models/dimerized_heisenberg/) |
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
layout. [`AGENTS.md`](AGENTS.md) has the APP frontmatter and the paper-agent sections. There is no
validated release yet; a tagged release with an `APP_PUBLICATION.json` manifest will follow.

| Directory | Content |
|---|---|
| [`paper/`](paper/) | Main diagrammatic Haldane-gap note, companion notes, and PDFs |
| [`code/`](code/) | Julia package, model generators/checkers, Lean proofs, leaderboard, and tools |
| [`data/`](data/) | Stored certificates, witnesses, numerical summaries, and external-data instructions |
| [`environment/`](environment/) | Setup commands and canonical dependency specifications |
| [`supplementary/`](supplementary/) | Interactive visual explainer and narrated-film instructions |
| [`blueprint/`](blueprint/), [`home_page/`](home_page/) | Lean blueprint and Jekyll website |

The paper is [`paper/certified-gap.pdf`](paper/certified-gap.pdf); the earlier full-length note is
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

## Citation

```bibtex
@article{Lu2026BridgingGaps,
  title  = {Bridging Gaps with Tensor Networks: From the {Haldane} Gap and Beyond},
  author = {Lu, Sirui},
  year   = {2026},
  note   = {Companion repository: \url{https://github.com/LionSR/CertifiedGaps}}
}
```

## License

Manuscript and notes: CC BY 4.0. Code and data: MIT. See [`LICENSE`](LICENSE).
