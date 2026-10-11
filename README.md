# A diagrammatic proof of the Haldane gap

The spin-1 antiferromagnetic Heisenberg chain H_L = Σⱼ **S**ⱼ·**S**ⱼ₊₁ has a spectral gap that
stays open as the ring grows (Haldane's conjecture; numerically Δ ≈ 0.4105). This repository
proves a uniform lower bound on that gap with spacetime tensor networks, drawn as diagrams,
and checks every finite input by computer.

The argument has two halves:

1. **An analytic doubling criterion.** Write the partition function of the ring as a trace of
   powers of a spatial transfer operator X_β of an exact continuous-time thermal MPO. Bond
   reflection positivity makes X_β self-adjoint, so space can be squared like imaginary time.
   Two purities, Pₛ(n, β) = Z₂ₙ(β)/Zₙ(β)² and P_τ(n, β) = Zₙ(2β)/Zₙ(β)², then improve
   quadratically under doubling of length and inverse temperature, and two purity bounds at a
   single scale give a gap bound for every larger ring.
2. **Certified finite inputs.** The purities at the starting scale are bounded with exact
   integer and rational arithmetic, reproducing the certificates of OpenAI's preprint
   [*The periodic spin-one Haldane gap*](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026).

## Results

| statement | status | where |
|---|---|---|
| Δ_L > (log 20)/784 ≈ 0.00382 for every even L ≥ 2304 | proof in the main note, given the purity bounds at (2304, 784) | [`outputs/thermal-mpo-gap-note.pdf`](outputs/thermal-mpo-gap-note.pdf) |
| the purity bounds Pₛ(2304, 784), P_τ(4608, 784) ≥ 0.99 | reproduced in exact arithmetic (status `passed`) | [`spin1/`](spin1/) |
| the doubling criterion, abstractly | formalised in Lean 4; the two purity bounds are its only unproved hypotheses | [`lean/`](lean/) |
| Δ_L > 0.0729 for every even L ≥ 36 | exact-arithmetic computation with a new temperature-lifting argument; not yet independently reviewed | [`spin1/experiments/`](spin1/experiments/) |
| open chains with endpoint fields h = 3/5: Δ_n > (log 10)/392 for odd n ≥ 1921 | proof in appendices A–C, given four finite estimates of the boundary-field paper | main note |
| which Hamiltonians the method covers | appendices D–H | [`outputs/thermal-mpo-gap-scope.pdf`](outputs/thermal-mpo-gap-scope.pdf) |
| how the starting estimate is certified, in diagrams | appendices I–M | [`outputs/thermal-mpo-gap-certificate.pdf`](outputs/thermal-mpo-gap-certificate.pdf) |

## Layout

| directory | content |
|---|---|
| [`outputs/`](outputs/) | the LaTeX notes and their PDFs, the shared diagram package `tnnotes.sty`, notation and diagram conventions |
| [`julia/`](julia/) | `ThermalMPOGap`, the shared Julia core: the doubling recurrence with outward-rounded and exact-rational bounds, exact arithmetic, reflection-positivity tests, exact diagonalisation, certified MPO arithmetic (Nemo/Arb), thermal purification (ITensors) |
| [`spin1/`](spin1/) | `Spin1Certificates`: exact reproduction of the spin-1 starting estimates (integer trace recursions on twisted rings, polynomial filters, MPS trial contractions) and the improved bounds |
| [`spin2/`](spin2/) | the spin-2 Heisenberg chain: transfer algebra, DMRG baselines, certified MPO bounds; no uniform gap yet |
| [`aklt/`](aklt/) | `AKLTPilot`: numerical test of the 2+1-dimensional recursion on the hexagonal AKLT model |
| [`lean/`](lean/) | Lean 4 + Mathlib formalisation of the abstract criterion |
| [`tools/`](tools/) | `fetch_upstream.jl`: downloads the preprint's evidence files at a pinned commit and checks their hashes |

[`AGENTS.md`](AGENTS.md) states what counts as a gap certificate here and how new ones are added.

## Reproduce

Requirements: Julia 1.12, elan (Lean), a LaTeX installation with REVTeX 4.2 and `latexmk`.

```sh
# Julia core and its tests
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# spin-1 certificates: fetch the preprint's evidence (pinned commit, hash-checked), then test
julia tools/fetch_upstream.jl
julia --project=spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=spin1 spin1/test/runtests.jl
julia --project=spin1 --threads=4 spin1/scripts/reproduce.jl all /tmp/spin1-fresh   # ≈20–45 min

# Lean
cd lean && lake exe cache get && lake build

# notes
latexmk -pdf -outdir=work/latex-build outputs/thermal-mpo-gap-note.tex
```

Each directory's README gives its own commands, run times, and what is and is not rigorous.

## Status of the claims

* The proof of the main theorem is complete on paper given the purity bounds at (2304, 784);
  those bounds are reproduced here in exact finite arithmetic, and the abstract doubling
  argument is checked in Lean. The transfer-operator representation and the finite
  computations are not formalised.
* Numbers labelled numerical (DMRG, Lanczos, floating-point exact diagonalisation, thermal
  purification) are diagnostics, never certificates.
