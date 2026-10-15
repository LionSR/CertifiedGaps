# A diagrammatic proof of the Haldane gap

The spin-1 antiferromagnetic Heisenberg chain H_L = Σⱼ **S**ⱼ·**S**ⱼ₊₁ has a spectral gap that
stays open as the ring grows (Haldane's conjecture; numerically Δ ≈ 0.4105). This repository
proves a uniform lower bound on that gap with spacetime tensor networks, drawn as diagrams,
checks every finite input by computer, and applies the same method to other spin chains.

The method, the **spacetime bootstrap**, has two halves:

1. **An analytic criterion.** Write the partition function of the ring as a trace of powers
   of a spatial transfer operator X_β of an exact continuous-time thermal MPO. Bond reflection
   positivity makes X_β self-adjoint, so space can be squared like imaginary time. Two
   purities, Pₛ(n, β) = Z₂ₙ(β)/Zₙ(β)² and P_τ(n, β) = Zₙ(2β)/Zₙ(β)², then improve
   quadratically under doubling of length and inverse temperature, and two purity bounds at a
   single scale (the *seed*) give a gap bound for every larger even ring.
2. **Certified finite inputs.** The seed purities are bounded by computer: for the spin-1
   chain with exact integer and rational arithmetic, reproducing the certificates of OpenAI's
   preprint [*The periodic spin-one Haldane gap*](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-periodic-spin-one-Haldane-gap-September-24-2026);
   for the other models with verified eigenvalue enclosures of complete seed-ring spectra.

## Spin-1 Heisenberg chain

| statement | status | where |
|---|---|---|
| Δ_L > (log 20)/784 ≈ 0.00382 for every even L ≥ 2304 | proof in the main note, given the purity bounds at (2304, 784) | [`notes/thermal-mpo-gap-note.pdf`](notes/thermal-mpo-gap-note.pdf) |
| the purity bounds Pₛ(2304, 784), P_τ(4608, 784) ≥ 0.99 | reproduced in exact arithmetic (status `passed`) | [`models/spin1/`](models/spin1/) |
| the criterion, and the theorem above given the two purity bounds | formalised in Lean 4 | [`lean/`](lean/) |
| Δ_L > 0.0729 for every even L ≥ 36 | exact-arithmetic computation with a temperature-lifting argument; not yet independently reviewed | [`models/spin1/`](models/spin1/) (`experiments/`) |
| open chains with endpoint fields h = 3/5: Δ_n > (log 10)/392 for odd n ≥ 1921 | proof in appendices A–C, given four finite estimates of the boundary-field paper | main note |
| which Hamiltonians the method covers | appendices D–H | [`notes/thermal-mpo-gap-scope.pdf`](notes/thermal-mpo-gap-scope.pdf) |
| how the starting estimate is certified, in diagrams | appendices I–M | [`notes/thermal-mpo-gap-certificate.pdf`](notes/thermal-mpo-gap-certificate.pdf) |
| smaller seed inputs for a Lean-checked seed: gap > 0.039 for every even L ≥ 60 from traces at even n ≤ 12 with real twists only | exact rational chain, conditional on the finite trace tables and the analytic seed lemmas | [`models/spin1_lean_seed/`](models/spin1_lean_seed/) (`seed_design/`) |
| proving the seed inside the Lean kernel (packed-integer arithmetic, counting anchors) | feasibility studies; no seed proved yet | [`models/spin1_lean_seed/`](models/spin1_lean_seed/), [`models/lean_seed_small/`](models/lean_seed_small/) |

## Other models

Each bound holds, with a unique ground state, for every even ring length at or above the seed
length n₀ (lengths count unit cells). Bounds are rounded down here; the model's README gives
the exact values, the RP basis, the exact seeds, run times, and what is and is not rigorous.
[`lean/`](lean/) states each certified model's gap theorem with its two seed inequalities as
the only open hypotheses.

| model | seed (n₀, β₀) | certified uniform gap bound | directory |
|---|---|---|---|
| spin-½ Heisenberg + staggered field, h = 1, 3/4, 1/2, 7/20, 1/4 | (8, 21/4 … 25/4) | 1.148, 0.890, 0.594, 0.383, 0.214 | [`heisenberg_staggered_field`](models/heisenberg_staggered_field/) |
| spin-½ dimerized Heisenberg, δ = 1/2, 3/10, 1/4 (2-site cells; δ = 1/4 is the smallest certifiable grid point with rings ≤ 16 spins) | (4, 6 … 29/4) | 0.534, 0.173, 0.058 | [`dimerized_heisenberg`](models/dimerized_heisenberg/) |
| spin-1 dimerized Heisenberg, δ = 9/10, 7/10 (2-site cells) | (2, 7/2), (3, 17/4) | 0.528, 0.218 | [`spin1_dimerized`](models/spin1_dimerized/) |
| spin-1 Heisenberg + D Σ(Sᶻ)², D = 4 and D ∈ [3.2059998, 3.206] | (12, 29/10), (6, 381/125) | 0.899, 0.0010 (every even L ≥ 6) | [`spin1_largeD`](models/spin1_largeD/) |
| two-leg Heisenberg ladder, J⊥ = 4, 7/2, 3, 5/2 (rung cells) | (4, 2 … 3) | 1.102, 0.778, 0.454, 0.042 | [`heisenberg_ladder`](models/heisenberg_ladder/) |
| two-leg transverse-field Ising ladder, g = 3 (rung cells) | (4, 3/2) | 0.636 | [`two_leg_tfim`](models/two_leg_tfim/) |
| Ising chain −ZZ − gX − hZ, 30 points g ∈ [1/4, 2], h ∈ [1/20, 1] | (6, per point) | ≥ 0.0763 at every point (smallest at g = 1, h = 1/20) | [`tfim_longitudinal`](models/tfim_longitudinal/) |
| spin-1 bilinear–biquadratic, θ ∈ [−7/10, 0] | no admissible even seed with n₀ ≤ 6 | negative result; seed scale estimated | [`blbq_haldane`](models/blbq_haldane/) |
| spin-2 Heisenberg chain | — | investigation, no certificate | [`spin2`](models/spin2/) |
| hexagonal spin-3/2 AKLT (2+1 dimensions) | — | numerical pilot | [`aklt`](models/aklt/) |

These seeds come from complete spectra of the seed rings, every eigenvalue enclosed by a
residual/Weyl bound whose floating-point rounding is bounded a priori under the IEEE binary64
model, including the BLAS products; that arithmetic model is the trusted base.

## Layout

| directory | content |
|---|---|
| [`notes/`](notes/) | the LaTeX notes and their PDFs, the shared diagram package `tnnotes.sty`, notation and diagram conventions |
| [`julia/`](julia/) | `CertifiedGaps`, the shared Julia core: the doubling recurrence with outward-rounded and exact-rational bounds, exact arithmetic, results I/O, reflection-positivity tests and exact RP certificates, exact diagonalisation, verified ring spectra, certified MPO arithmetic (Nemo/Arb), thermal purification (ITensors) |
| [`models/`](models/) | one directory per model in a common format (AGENTS.md, "Directory format"): certificate, `certify.jl`, `check.jl`, tests and README; `spin1/` is the exact-arithmetic reproduction of the spin-1 certificate; `spin1_lean_seed/` and `lean_seed_small/` study kernel-checked seeds |
| [`lean/`](lean/) | Lean 4 + Mathlib formalisation of the criterion and of each certified model's gap theorem, with the seeds as hypotheses |
| [`tools/`](tools/) | `fetch_upstream.jl`: downloads the preprint's evidence files at a pinned commit and checks their hashes |

[`AGENTS.md`](AGENTS.md) states what counts as a gap certificate here and how new ones are added.

## Reproduce

Requirements: Julia 1.12, elan (Lean), a LaTeX installation with REVTeX 4.2 and `latexmk`.

```sh
# Julia core and its tests
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# spin-1 certificates: fetch the preprint's evidence (pinned commit, hash-checked), then test
julia tools/fetch_upstream.jl
julia --project=models/spin1 -e 'using Pkg; Pkg.instantiate()'
julia --project=models/spin1 models/spin1/test/runtests.jl
julia --project=models/spin1 --threads=4 models/spin1/scripts/reproduce.jl all /tmp/spin1-fresh   # ≈20–45 min

# any other model: fast checker (add --full for the slow replay), and its tests
OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder models/heisenberg_ladder/check.jl
OPENBLAS_NUM_THREADS=1 julia --project=models/heisenberg_ladder models/heisenberg_ladder/test/runtests.jl

# Lean
cd lean && lake exe cache get && lake build

# notes
latexmk -pdf -outdir=work/latex-build notes/thermal-mpo-gap-note.tex
```

Each directory's README gives its own commands, run times, and what is and is not rigorous.

## Status of the claims

* The proof of the spin-1 theorem is complete on paper given the purity bounds at
  (2304, 784); those bounds are reproduced here in exact finite arithmetic, and the criterion
  with its application to the spin-1 ring is checked in Lean. The finite computations are not
  formalised: in Lean each model's two seeds are named `sorry` lemmas. Kernel-checked seed
  arithmetic exists so far only as a study (`models/lean_seed_small/`, one Ising point,
  conditional on four log Z enclosures).
* The other models' seeds rest on verified eigenvalue enclosures whose rounding analysis
  assumes IEEE binary64 arithmetic and conventional BLAS; each was independently replayed or
  reviewed (see the model's README and `METHOD.md`).
* Numbers labelled numerical (DMRG, Lanczos, floating-point exact diagonalisation, thermal
  purification, extrapolations) are diagnostics, never certificates.
