# Toward a Lean-checked seed

The Lean gap proofs of `lean/` are complete except for their seed inequalities (`sorry` in the
model files). This directory collects five studies of how to prove a seed inside the Lean
kernel with `decide +kernel` only (no `native_decide`, no new axioms): four for the spin-1
Heisenberg seed and one carrying the pipeline furthest on the cheapest model, the longitudinal
Ising chain. Julia is an untrusted generator; the Lean side is the opt-in library
`SeedKernel` of the repository's Lake package (`lake build SeedKernel`).

**One unconditional gap.** The Ising study closes: `SeedKernel.Ising.ising_gap` proves, with no
hypothesis and no `sorry`, that the Ising ring at g = 1/4, h = 1/20 has gap > log(13480)/(97/4)
≈ 0.3921 for every even L ≥ 6. The spin-1 studies remain exact but conditional:

| study | question | best result (status) |
|---|---|---|
| `seed_design` | which finite inputs does the spin-1 seed need? | design C: gap > 0.039036761 for every even L ≥ 60 from even n ≤ 12, real twists, table error 10⁻³ (exact replay of `ChainCertificate.Valid`; conditional on trace tables and analytic lemmas) |
| `spectrum` | how many low states, how precise, how expensive? | 26 (I) / 34 (P) states below the tail cutoff at n = 12, β = 49/4; Q48 filtering meets every residual target (Float64, non-certifying) |
| `anchor` | can a cheap comparison anchor the eigenvalue count? | exact local PSD comparison; 121 → 49 comparison levels at n = 12, Sᶻ = 0 against 8 physical (exact counts; anchor does not close) |
| `spin1_kernel` | does packed arithmetic fit the kernel? | yes: n = 12 matvec + dot in 20–26 s; Ritz witnesses for H₆, H₈, H₁₀ kernel-checked (arithmetic only) |
| `ising` | the whole pipeline on the cheapest model | **proved**: seed of `IsingLongitudinal` at g = 1/4, h = 1/20, (n₀, β₀) = (6, 97/4), and the gap > log(13480)/(97/4) ≈ 0.3921 for every even L ≥ 6 (`ising_gap`; standard axioms only) |

## Model

- **Spin-1 Heisenberg** (`seed_design`, `spectrum`, `anchor`, `spin1_kernel`):
  `H = Σⱼ Sⱼ·Sⱼ₊₁` on a ring of n sites (lengths count sites), shifted by
  `a = 700741/500000` per bond as in `models/spin1`, with the closing bond twisted by a rotation
  about z of 0 (I), π (P) or 2π/3 (C). The numerical and kernel studies use the Sᶻ product
  basis, index `Σⱼ (mⱼ + 1) 3ʲ`, integer bond `h = SᶻSᶻ + (S⁺S⁻ + S⁻S⁺)/2` (eigenvalues
  −2, −1, 1). Source 1 is β = 49/4 on n = 4–12 with I, P; source 2 is β = 21/2 on n = 6, 8,
  10 with I, P, C (`models/spin1/results/source{1,2}.json`).
- **Longitudinal Ising** (`ising`): `H_L = −Σⱼ ZⱼZⱼ₊₁ − g Σⱼ Xⱼ − h Σⱼ Zⱼ` at g = 1/4,
  h = 1/20, exactly `CertifiedGaps.IsingLongitudinal.model (1/4) (1/20)`. The checkers use
  the integer matrix `A = 20H` in the computational basis; `test/ising.jl` pins its
  periodic-bond convention at L = 2 to that of `RPRing.hamiltonian`.

## Reflection positivity

Nothing new is proved here. The spin-1 route uses the Cartesian basis of `models/spin1`
(`−h = Σ_α S^α ⊗ conj(S^α)` up to the shift; `lean/CertifiedGaps/Heisenberg/Main.lean`); the
unitary bridge from it to the Sᶻ-basis integer matrices used here is open. The Ising model is
the repository's `RPRing` with `C = Z` (Pauli-Z basis, `−h_bond = Z ⊗ Z`, real on-site term),
whose RP structure Lean checks in `lean/`.

## Result

### seed_design: smaller finite seed inputs

Three seed designs replayed in exact rationals down to `ChainCertificate.Valid`
(`results/seed_design/certificate.json`; 9-step chains on the grid 10⁻⁵¹²). **Conditional** on
the trace tables to the stated absolute error and on the analytic sector, filter and lifting
lemmas (METHOD.md §1), none checked in Lean.

| design | thermal inputs (lengths; twists; abs. error) | b → β₀ | k; trial | εₛ (p) | εₜ (q) | gap > | lengths |
|---|---|---|---|---:|---:|---:|---|
| A | 6, 8, 10; 1, P, C; 3/10000 | 21/2 → 105/4 | 30; 120 | 0.124048525 | 0.125580476 | 0.005966152 | even L ≥ 60 |
| B | 6, 8, 10; 1, P, C; 4/10000 | 21/2 → 525/16 | 30; 144 | 0.117346685 | 0.128925715 | 0.006309519 | even L ≥ 72 |
| C | 4–12 even; 1, P; 1/1000 | 49/4 → 49/2 | 30; 120 | 0.047914768 | 0.112578548 | 0.039036761 | even L ≥ 60 |

All three exceed the original target `log 20 / 784`. B's εₜ exceeds the diagonal threshold
0.12731, but its unequal-defect chain is valid. Negative results (`diagnostics.json`, Float64):
length 8 fails even at zero error (best max defect ≈ 0.5288 over 132,075 lift cases); length
10 without C fails (≈ 0.2956); length 10 with C at error 10⁻³ narrowly fails (≈ 0.1338). The
numerical Pareto front (`pareto.csv`) has 39 rows from 94 grid representatives.

### spectrum: seed spectra, filtering precision, kernel costs (Float64, non-certifying)

| source | β | rings | states kept at n_max (I / P / C) | adaptive Q48 Hv | Lean vector + Gram forecast |
|---|---:|---|---|---:|---|
| 1 | 49/4 | n = 4–12, I/P | 26 / 34 / – at n = 12 | 34,388 | 38–83 min |
| 2 | 21/2 | n = 6, 8, 10, I/P/C | 22 / 30 / 26 at n = 10 | 17,800 | 5–11 min |

- 48-bit fixed-point Chebyshev filtering of degree 16 meets every residual target; degree-32
  Q48 plateaus at 2.99·10⁻⁶ at n = 12/I/Sᶻ = 0, Q112 reaches 2.34·10⁻⁹.
- Compressed Sᶻ indexing needs thousands of shift masks (3,816 at n = 12, Sᶻ = 0; 13,620 at
  Sᶻ = 2) instead of about 50, so keep full 3ⁿ vectors with 2n masks.
- The naive two-cut anchor leaves 242 Sᶻ = 0 levels below ρ at n = 12/I where the ring keeps
  10.

### anchor: exact local comparison anchors

`H ≥ H0 = A_l ⊗ I + I ⊗ A_l + 2c` with `A_l = H_open,l + B₁₂ + B_(l−1,l)`, `B = u h + v h²`,
holds because the four-site remainder `h₂₃ − B₁₂ − B₃₄ − c` is PSD, checked by exact LDL in
every Sᶻ block (`results/anchor/certificate.json`; half-chain spectra enclosed exactly, radius
≤ 2.8·10⁻⁷).

| n | count-opt (u, v, c) | Sᶻ = 0 count: bare / ground-opt / count-opt | binned-tail Sᶻ = 0 (total) |
|---:|---|---|---|
| 8 | (−3/10, −1/10, −5076841/2500000) | 44 / 63 / 34 | 28 (98) |
| 10 | (−1/5, 0, −5538721/2500000) | 86 / 115 / 62 | 38 (140) |
| 12 | (−1/5, 1/20, −48533/20000) | 191 / 212 / 121 | 49 (183) |

The physical ring has 2, 5, 8 Sᶻ = 0 levels below the cutoff (n = 8, 10, 12). A one-shot count
does not give ordered lower bounds (`H = diag(−2, −1, 10)` counterexample); a Float64
Lehmann–Goerisch continuation completed in all 26 nonzero sectors (658 stages, ≈ 40 min of
Julia), estimated at 80–97 kernel hours naively at n = 12.

### spin1_kernel: packed `Nat` arithmetic in the kernel

Kernel-checked (`lean/SeedKernel/Spin1/H*.lean`, Std only): integer inequalities `ritzCheck` on
the moments G = vᵀv, A = vᵀFv, C = (Fv)ᵀFv of the trial v = (nI − H_M)^q e_Néel.

| n | degree q | lane bits | Ritz witness θ ≤ | residual witness ≤ |
|---|---:|---:|---:|---:|
| 6 | 24 | 262 | −8495447/10⁶ | 274967/10⁶ |
| 8 | 32 | 345 | −11179504/10⁶ | 279904/10⁶ |
| 10 | 16 | 222 | −13548710/10⁶ | 632355/10⁶ |

Dense 128-bit matvec + dot: 0.5 s (n = 6, 8), 1.3 s (n = 10), 26 s (n = 12) per fresh Lean
process on an Apple M5 (`results/spin1_kernel/benchmarks.tsv`). The untuned cut anchor fails
the separation test at n = 6, 8 by more than one energy unit. Mathlib (pinned c55e6e78) lacks
Courant–Fischer, the Weyl inequalities and Cauchy interlacing.

### ising: the pipeline on the cheapest model

Kernel-checked (`lean/SeedKernel/Ising/`): with `S = 2⁵²`, generated integer vectors v₁, v₂
and integers λ₁, λ₂ on the 6- and 12-site rings, `|vᵢᵀvᵢ − S²|, |v₁ᵀv₂| ≤ 10⁻¹² S²`,
`Σᵢ ‖S·A vᵢ − λᵢ vᵢ‖² ≤ 10⁻¹⁸ S⁴`, `λ₁ + S < λ₂`, `λ₂ + S < (82 − 26n) S` (packed and direct,
the direct 12-site check in 64 blocks with a coverage proof), exactly two classical
configurations below the anchor, and rejection of corrupted proposals.

| ring | dimension | candidate E₀ | candidate E₁ | anchor ρ | normalized residual |
|---|---:|---:|---:|---:|---:|
| 6 sites | 64 | −6.391774806847206 | −5.796599647006008 | −3.7 | 5.97×10⁻¹⁵ |
| 12 sites | 4096 | −12.78354960854822 | −11.59319932209873 | −11.5 | 1.49×10⁻¹³ |

`Ising/Seed.lean` turns these checks into `first_point_seed`: exactly the first of the thirty
conjunctions of `IsingLongitudinal.seeds` (`1 − 9064·10⁻⁹ ≤ spatialPurity 6 (97/4)`,
`1 − 3141·10⁻⁸ ≤ thermalPurity 12 (97/4)`), with no hypothesis, and then `ising_gap`. The
steps: a residual bound gives an eigenvalue within 10⁻¹⁰ of each candidate (`Ising/Spectral`);
`Re ⟪x, Hx⟫ ≥ ∑ (E(σ) − L/4) |x(σ)|²` and the classical count give at most two eigenvalues below
the anchor (counting by comparison); the checker's integer sums are identified with the matrix
`IsingLongitudinal.model (1/4) (1/20)` (`Ising/Hamiltonian`, `Ising/Bridge`); then
`e^{-β(θ₀+r)} ≤ Z ≤ e^{-β(θ₀−r)}(1 + e^{-β(θ₁−θ₀)} + 2^L e^{-β(ρ−θ₀+r)})` with `e > 2.7182818283`
bounds both purities. It does not use `seeds`; the earlier conditional bridge (`LogEnclosures`) is
retired. A 43-target inventory (`results/ising/costs.json`, Float64) shows the full k×k
Gram route reaching hours to ~1,000 hours for the 65,536- and 531,441-dimensional seeds.

### Open obligations

1. **Anchor.** A kernel-checkable counting anchor ρ ≤ λ_{k+1} for the spin-1 rings (better
   comparisons, fewer continuation stages, Lehmann–Goerisch with a priori bounds); the
   bottleneck of the whole route.
2. **Spectral linear algebra in Lean.** Lehmann bounds and ordered-eigenvalue statements
   (absent from Mathlib). The residual bound and counting by comparison are now proved
   (`Ising/Spectral.lean`) and apply to any Hermitian matrix.
3. **Packing correctness.** The `Prop`s in `lean/SeedKernel/Obligations.lean` and
   `Spin1/Obligations.lean` (`DotCorrect`, `HApplyCorrect`, `FilterCorrect`, ...) and the lane
   decoding of `Ising/Packed.lean`; nothing assumes them.
4. **Hamiltonian bridges.** Identify the checked integer matrices with the Lean Hamiltonians
   (Cartesian ↔ Sᶻ for spin 1, including twisted rings); done for Ising (`Ising/Bridge.lean`).
5. **Scalar enclosures.** The trace tables with their sector/filter lemmas (spin 1); the Ising
   case is done in `Ising/Seed.lean`.

Estimates: 900–1,950 Lean lines for the Ising point, 8–15k for the spin-1 route (METHOD.md).
Target design C first (real arithmetic only, largest ring 12).

## Reproduce

Julia 1.12, `export OPENBLAS_NUM_THREADS=1`, and `lake` (elan) for the pinned toolchain.
From this directory; times measured on an Apple M5 laptop.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl         # 30 s: 1,582 tests of all five studies (no Lean)
julia --project=. check.jl                 # 1 min: exact replays + data checks, then
                                           #   `lake build SeedKernel` and the axiom audit
julia --project=. check.jl --full          # 4.7 min: also regenerate the spin-1 witnesses,
                                           #   rebuild SeedKernel from clean, run every benchmark
julia --project=. certify.jl [study ...]   # 75 s: rewrite the exact records and Lean witnesses
```

`certify.jl` takes `seed_design` (5 s), `anchor` (3 s), `ising` (9 s) and `spin1_kernel`
(55 s; run it last, it hashes the shared drivers). `check.jl` needs the Lean build cache of
`lean/` (`lake exe cache get` there, or copy an existing `.lake`); with Mathlib cached,
`lake build SeedKernel` from clean takes about 80 s, the audit (`lake env lean
SeedKernel/Audit.lean`) 5–30 s, and the `SeedKernelBench` files about 1.7 min.

Non-certifying drivers, by study (outputs under `results/<study>/`; caches and raw rows under
the gitignored `cache/`, `upstream/`, `.cache/`, `runtime/` and the paths in `.gitignore`):

```sh
julia scripts/fetch_paper.jl         # paper sources -> upstream/ (network; tools/fetch_upstream.jl)
# seed_design: spectra (70 s for n <= 10), searches, diagnostics, trial contraction, summary
julia --project=. scripts/seed_design/spectra.jl 4 5 6 7 8 9 10
julia --project=. scripts/seed_design/search_existing.jl   # results/seed_design/raw/
julia --project=. scripts/seed_design/search_new.jl 8      # also 10
julia --project=. scripts/seed_design/refine.jl
julia --project=. scripts/seed_design/diagnostics.jl
julia --project=. scripts/seed_design/shift_search.jl      # 8 s
julia --project=. scripts/seed_design/trial_scan.jl        # network
julia --project=. scripts/seed_design/summarize.jl         # needs raw/, spectra, upstream/
# spectrum: 21 rings (2.6 min), analysis, filters, packing, costs, provenance (~10 s besides)
julia --project=. scripts/spectrum/spectra.jl
julia --project=. scripts/spectrum/analyze.jl
julia --project=. scripts/spectrum/filters.jl float        # also fixed, sweep
julia --project=. scripts/spectrum/packing.jl
julia --project=. scripts/spectrum/anchor_probe.jl
julia --project=. scripts/spectrum/summarize.jl
julia --project=. scripts/spectrum/runtime_estimates.jl
julia --project=. scripts/spectrum/provenance.jl
# anchor: scans, ring spectra (minutes), analysis, bins, continuation (~40 min), tables
julia --project=. scripts/anchor/scan.jl
julia --project=. scripts/anchor/count_scan.jl
julia --project=. scripts/anchor/spectra.jl 8 10 12
julia --project=. scripts/anchor/analyze.jl
julia --project=. scripts/anchor/tail_bins.jl
ANCHOR_FRESH=1 julia --project=. scripts/anchor/homotopy.jl 12 I binned   # each n, g
julia --project=. scripts/anchor/homotopy_sectors.jl      # n = 8, 10; then `I`, `P` for n = 12
julia --project=. scripts/anchor/report_tables.jl         # METHOD.md §3 tables, homotopy_summary
# spin1_kernel: Float64 structure/anchor tables (5 s), kernel benchmarks, Mathlib audit
julia --project=. scripts/spin1_kernel/experiments.jl
julia --startup-file=no scripts/spin1_kernel/benchmark.jl 300   # ~2.5 min
julia scripts/audit_mathlib.jl       # Mathlib audit; reads lean/.lake/packages/mathlib, needs rg
# ising: two-level feasibility, cost model, 43-row table, inventory
julia --project=. scripts/ising/feasibility.jl
julia --project=. scripts/ising/costs.jl --reuse           # without --reuse: re-run benchmarks
julia --project=. scripts/ising/report_tables.jl
julia --project=. scripts/ising/inventory.jl               # see "Layout" for its inputs
```

## What is rigorous

- **Lean kernel, standard axioms only.** Every declaration in `lean/SeedKernel/Audit.lean`
  (40, checked by `check.jl`) depends on at most `propext`, `Classical.choice`, `Quot.sound`;
  none uses `sorryAx` or `Lean.ofReduceBool`. The computational checks (spin-1 witnesses and
  tests, Ising Gram/residual/anchor checks) use only `propext` or nothing. `Ising/Seed` imports
  `IsingLongitudinal`, whose `seeds` is a `sorry`, but no audited declaration depends on it.
- **Unconditional (Ising).** `seed_spatial`, `seed_thermal`, `first_point_seed` and `ising_gap`
  have no hypotheses. **Conditional (spin 1).** The spin-1 witnesses are arithmetic facts about
  packed integers until the obligations above are proved.
- **Exact in Julia (not Lean).** `seed_design` (BigInt rationals, integer-root and series
  enclosures; a line-by-line transcription of `ChainCertificate.Valid`) and `anchor` (rational
  LDL, half-chain enclosures from Float64 eigenvector proposals by bounded Int128 polar
  decomposition and Weyl, comparison counts, Taylor-degree-256 tail sums). The length-144 trial
  of design B is a stored exact Nemo contraction (`results/seed_design/trial_scan.json`), not
  recomputed by the checker.
- **Untrusted generators.** Float64 LAPACK proposes the Ising vectors and the anchor
  eigenvectors; wrong proposals fail the exact checks. The Julia replays of kernel-checked
  facts (`check_spin1_kernel`, `check_ising`) are diagnostics.
- **Floating point only.** All spectra, counts, cutoffs, filter residuals, searches, Pareto
  fronts, continuations, inventories and cost/time forecasts.
- **Platform.** LAPACK proposals differ across platforms; the stored Ising and anchor records
  were regenerated on the M5 machine and reproduce bit for bit there. Other platforms may give
  different but equally valid proposals.

## Layout

| path | content |
|---|---|
| `certify.jl`, `check.jl` | exact records and Lean witnesses by study; the checker (fast / `--full`) |
| `src/LeanSeeds.jl` | module `LeanSeeds`, one submodule per study |
| `src/Common.jl` | shared paths, source hashes, Sᶻ basis, Lake driver and axiom-audit parser |
| `src/seed_design.jl` | `SeedDesign`: Float64 search, exact filters/seeds/chains, certify and check |
| `src/spectrum.jl` | `Spectrum`: twisted Sᶻ Hamiltonians, low spectra, fixed-point filters, data checks |
| `src/anchor.jl` | `AnchorStudy`: integer bonds, comparison matrices, exact LDL and enclosures, certify and check |
| `src/spin1_kernel.jl` | `Spin1Kernel`: witness generator, Lean text, witness inventory and check |
| `src/ising.jl` | `IsingSeed`: proposal generator, `Generated*.lean` text, exact replay and check |
| `scripts/<study>/` | the drivers under "Reproduce" |
| `scripts/fetch_paper.jl` | shared: runs `tools/fetch_upstream.jl --paper` (paper sources into `upstream/`) |
| `scripts/audit_mathlib.jl` | spin1_kernel Mathlib audit; kept outside `scripts/spin1_kernel/`, whose `.jl` files are hash-pinned by `results/spin1_kernel/witnesses.json` |
| `test/runtests.jl` | runs `test/<study>.jl` (one module each) |
| `results/seed_design/` | `certificate.json`, `compact_witnesses.json` (exact); `trial_scan.json`, `pareto.csv`, `summary.json`, `diagnostics.json`, `shift_search.json` (search records) |
| `results/spectrum/` | `spectrum_<n>_<twist>.json`, `analysis.json`, `summary.json`, CSV tables, kernel logs, `provenance.json` |
| `results/anchor/` | `certificate.json` (exact); `optimized.json`, `count_optimized.json`, `binned_tails.json` (Float64 inputs of `certify.jl`, the continuation); `sectors.csv`, `homotopy_summary.csv` |
| `results/spin1_kernel/` | `witnesses.json`, `exact_witnesses.tsv`; `axioms.txt` (audit of all of `SeedKernel`), benchmarks, Mathlib audit, structure and anchor scans, paper hashes |
| `results/ising/` | `certificate.json` (integer proposals); `feasibility.json`, `inventory.{json,tsv}`, `costs.json` |
| `lean/SeedKernel/` (repository `lean/`) | the Lean side; file map in `lean/README.md` |

Hash pinning: `check.jl` compares `results/spin1_kernel/witnesses.json` with the current
SeedKernel spin-1 Lean files, Lake configuration and this study's Julia sources (including
`certify.jl`, `check.jl`, `src/Common.jl`), and `results/ising/certificate.json` with
`src/ising.jl`; after editing any of these, rerun `certify.jl spin1_kernel` (and `ising`). The
`metadata` inside the seed_design search records (`summary.json`, `trial_scan.json`,
`diagnostics.json`, `shift_search.json`), the spectrum data files (`provenance.json` hashes the
current layout) and the anchor inputs `optimized.json`, `count_optimized.json` name the source
layout at their generation time. They were kept: some need hours, the network or gitignored
caches, a rerun of `shift_search.jl` here changes Float64 values in the last digits only, and
rerunning the anchor scans could move the exact parameters (u, v, c). `binned_tails.json` was
regenerated (rows identical). `scripts/ising/inventory.jl` reads other models' stored spectra
as of commit `463c4a5`; two of them are now gitignored regenerable data
(`heisenberg_ladder/results/verified_L8.json.gz`, `spin1_dimerized/results/spectrum_d9_c4.json`).

## History

- PR #24 (`seed_design`), PR #22 (`spectrum`), PR #27 (`anchor`), PR #21 (`spin1_kernel`, then
  `lean_prototype`) and PR #25 (`ising`, then `models/lean_seed_small`), each first in its own
  directory and Julia project, with standalone Lake projects for the two Lean prototypes.
- Consolidation: one Julia project here (union of the environments, same package versions);
  one module per study with shared paths, Sᶻ basis, Lake driver and source hashing; the paper
  fetchers merged into `scripts/fetch_paper.jl`. The two Lake projects became the opt-in
  library `SeedKernel` of `lean/`, with one shared packed-lane module (`SeedKernel/Packed.lean`,
  now also under the Ising `Packed` checks) and `SeedBridge` importing `CertifiedGaps`
  directly. All exact records were regenerated: every exact number is identical (only
  provenance hashes, timings and the seed_design trial-file paths changed), and the Lean
  witness files are byte-identical apart from module names.
- Layout tidy: `scripts/fetch_paper.jl` now forwards to the shared
  `tools/fetch_upstream.jl --paper` (one place for the pinned commit; it also skips files that
  are already present with the pinned hash), and `scripts/spin1_kernel/audit_mathlib.sh` was
  ported to `scripts/audit_mathlib.jl` (byte-identical output). The port sits outside
  `scripts/spin1_kernel/` so that the spin1_kernel source-hash set is unchanged.
  `scripts/spectrum/KernelBench.lean` stays where it is: `scripts/spectrum/provenance.jl` hashes
  it into `results/spectrum/provenance.json`, and the Lean benchmarks under `lean/` are not
  edited from model directories.
- The Ising seed closed: `Ising/Spectral`, `Ising/Hamiltonian`, `Ising/Bridge` and `Ising/Seed`
  prove `first_point_seed` and the unconditional `ising_gap` from the existing kernel checks
  (no new witnesses); `SeedBridge`, `Arithmetic` and `ResidualLemma` were retired.
