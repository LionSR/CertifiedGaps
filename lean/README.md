# CertifiedGaps: Lean 4 formalization of spectral gaps via the spacetime bootstrap

This Lean 4 + Mathlib package formalizes the spacetime-bootstrap proof that reflection-positive quantum
spin rings have a uniform spectral gap for all sufficiently large even ring lengths, following the main
note ([`../notes/thermal-mpo-gap-note.tex`](../notes/thermal-mpo-gap-note.tex)) and scope note
([`../notes/thermal-mpo-gap-scope.tex`](../notes/thermal-mpo-gap-scope.tex)). **Every analytic and
algebraic step is proved from standard axioms except the finite numerical seed inputs** (two purity
inequalities per model at the seed scale $(n_0, \beta_0)$). Those finite seed bounds appear as explicit
named `sorry` lemmas whose values come from the computer certificates in [`../models/`](../models/).
One seed is now proved outright in the kernel: `SeedKernel.Ising.ising_gap` is an unconditional gap
theorem (below). The open-chain argument of the note (Appendices A–C) is formalized as well: it
gives a uniform gap for odd open chains with boundary fields from four finite inputs.

## Overview

### The headline theorem

```lean
-- CertifiedGaps/Heisenberg/Main.lean
theorem SpinOneHeisenberg.spin_one_heisenberg_gap {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784)
```

For every even ring length `L ≥ 2304`, the periodic spin-1 Heisenberg Hamiltonian
$H_L = \sum_{j=1}^L \mathbf{S}_j \cdot \mathbf{S}_{j+1}$ acting on `(Fin L → Fin 3)` has a unique ground
state, and every excited eigenvalue lies strictly more than `log 20 / 784 ≈ 0.00382` above the ground
energy (`GapGT e c` asserts that some basis index `m₀` has `e m₀ + c < e m` for all `m ≠ m₀`). The proof
follows OpenAI's Proposition 5.5. The seeds `seed_spatial : P_s(2304, 784) ≥ 0.99` and
`seed_thermal : P_τ(4608, 784) ≥ 0.99` are proved from two finite seed hypotheses at the ring of
36 sites by six kernel-checked doubling steps (`chain36`):
- `seed36_spatial : 1 - P_s(36, 49/4) ≤ 0.047916`, and
- `seed36_thermal : 1 - P_τ(72, 49/4) ≤ 0.181521`.

The first bound of OpenAI's Theorem 1.1 is formalized too:

```lean
theorem SpinOneHeisenberg.spin_one_heisenberg_gap_sixty {L : ℕ} (hL : 60 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (4 / 105 * log (80 / 79))
```

It uses OpenAI's general criterion (their Proposition 4.3, `gapGT_of_coefficient`: seed defect
`u₀ < 1/6`, coefficient `C ≥ G(u₀)` with `C u₀ < 1` and `C (1 - 3u₀) > 3`), here with `u₀ = 1/8` and
`C = 79/10`, and the seeds `seed60_spatial : P_s(60, 105/4) ≥ 7/8` and
`seed60_thermal : P_τ(120, 105/4) ≥ 7/8` (their Proposition 5.4). All four periodic seeds are
reproduced in exact integer and rational arithmetic in [`../models/spin1/`](../models/spin1/).

### The open chain with boundary fields

```lean
-- CertifiedGaps/Heisenberg/OpenChain.lean
theorem SpinOneHeisenberg.spin_one_boundary_field_gap {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum (3 / 5) n) (Real.log 10 / 392)
```

For every odd `n ≥ 1921`, the open spin-1 Heisenberg chain with equal end fields
$`H_n = \sum_{i=0}^{n-2} \mathbf{S}_i \cdot \mathbf{S}_{i+1} - \tfrac35 (S_0^z + S_{n-1}^z)`$ has a
unique ground state and a gap `> log 10 / 392 ≈ 0.005874` (Theorem 2 of the main note, Appendices
A–C; the theorem and its inputs are from OpenAI's boundary-field preprint). The proof is complete
modulo four finite inputs at `β = 49` (the note's eq. (C1)):
- `seed_ring_spatial : P_s(120, 49) ≥ 0.948` and `seed_ring_thermal : P_τ(240, 49) ≥ 0.948` for the
  periodic ring;
- `seed_boundary`: the boundary defect `e₁₄₂ ≤ 0.14` of the opposite-field chain, a statement about
  traces through `boundaryDefect` (the leading eigenvalue is `inf_N Z_{2N}^{1/(2N)}` and the leading
  boundary weight is `inf_N G_{2N} / μ^{2N}`);
- `seed_open_thermal`: the opposite-field thermal purity `P_τ(240, 49) ≥ 0.79`.

`gapGT_of_boundaryFieldInputs` proves the same gap from these four inputs as hypotheses, with no
`sorry`.

The criterion itself is general. `gapGT_of_openChain` takes any seed scale `β₀`, ring seeds at
lengths `2N₀` and `4N₀`, a boundary seed at length `2K₀`, an open-chain seed at length `4N₀`, and
`m` base steps, and gives `Δ_n > log(1/(10x)) / (β₀ 2^m)` for odd `n ≥ 4 N₀ 2^m`. Its finite data
form an `OpenCertificate`, whose rational conditions `OpenCertificate.Valid` are checked by
`decide +kernel`; `gapGT_of_openCertificate` applies one. The script
[`scripts/open_certificate.jl`](scripts/open_certificate.jl) searches for the data from the seeds,
and Lean re-checks it. The note's numbers are `boundaryFieldCertificate`. From the same four
inputs, the script finds `sharpBoundaryFieldCertificate`, which hands over at the defect actually
reached (`x = 2.241 · 10⁻⁴` in place of `1/100`). The resulting
`spin_one_boundary_field_gap_sharp` gives `Δ_n > log(10⁶/2241)/392 ≈ 0.01556` for every odd
`n ≥ 1921`, 2.65 times the note's bound.

The open-chain layer follows the ring layer. At every finite Trotter depth, the ring
traces and the open-chain traces are moments of one real symmetric transfer matrix `X = F J` (`J`
reverses imaginary time): `Z_{2N} = Tr X^{2N}` and `G_n = ⟨v, Xⁿ v⟩`. The opposite-field chain at
even `n` and the equal-field chain at odd `n` are both of this form. The closed inequalities
`OpenMomentBounds` pass to the limit by compactness. The sign of the leading eigenvalue at odd
lengths is fixed by positivity of the trace, so the note's `U_P` symmetry argument is not needed.

### The first unconditional gap

```lean
-- SeedKernel/Ising/Seed.lean (lake build SeedKernel)
theorem SeedKernel.Ising.ising_gap {L : ℕ} (hL : 6 ≤ L) (hLe : Even L) :
    GapGT ((IsingLongitudinal.model (1 / 4) (1 / 20)).spectrum L) (log 13480 / (97 / 4))
```

For every even ring length `L ≥ 6`, the longitudinal-field Ising ring
$-\sum_j Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac{1}{20} \sum_j Z_j$ has a unique ground state and
a spectral gap `> 0.392`. The theorem has no `sorry`, no hypothesis, and no `native_decide`; it
depends only on the axioms `propext`, `Classical.choice`, and `Quot.sound`. Its seed is proved from
kernel-checked integer residual and counting checks on the 6-site and 12-site rings (see
[SeedKernel](#seedkernel-kernel-checked-seed-arithmetic) below).

### Architecture of the formalization

The library is structured in five layers, each depending only on the preceding ones:

1. **From purities to a spectral gap** (`PurityCriterion.lean`). Write the partition function as
   `Z_L(β) = Re Tr e^{-βH_L}` and define the spatial purity `P_s(n, β) = Z_{2n}(β)/Z_n(β)²` and thermal
   purity `P_τ(n, β) = Z_n(2β)/Z_n(β)²`. Assuming the *spatial moment inequalities* (which capture the
   consequences of `Z_n(β)` being the trace of the $n$-th power of a self-adjoint transfer operator),
   both purities improve quadratically when ring length and inverse temperature double together. Two
   purity bounds at a single seed scale $(n_0, \beta_0)$ therefore imply a uniform spectral gap for all
   even $L \ge n_0$: `gapGT_of_purities` (for equal seed purities $\ge 0.99$) and `gapGT_of_seeds`.
2. **Unequal seed defects and rational chains** (`Chain.lean`, `ChainCheck.lean`). Certified seeds at
   small $n_0$ typically have unequal spatial and thermal defects $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau)$
   larger than $0.01$. The theorem `gapGT_of_chain` tracks the exact asymmetric doubling recurrence for
   $J$ steps before handing off to the equal-defect tail bound. A `ChainCertificate` encodes all step
   conditions as rational inequalities that Lean's kernel verifies by `decide +kernel`, yielding a gap
   lower bound of the form `log r / β₀`.
3. **Proving the spatial moment inequalities** (`Transfer/*`, `Ring/*`). This layer contains the
   operator-theoretic core of the formalization and involves no numerics. For any reflection-positive
   ring (`RPRing`), where the nearest-neighbor bond interaction has the form
   $-h = \sum_\alpha C_\alpha \otimes F C_\alpha F$ with real symmetric/antisymmetric $C_\alpha$, a real
   symmetric involution $F$, and a real symmetric $F$-invariant on-site term $o$:
   * Trotterize $e^{-\beta H_L}$ into $M$ slices; each Trotter trace equals $\operatorname{Tr}(B^T B)^{L/2}$
     for a real matrix $B$ (`Ring/Transfer.lean`, `Transfer/PiTensor.lean`).
   * Any positive semidefinite matrix $B^T B$ satisfies the spatial moment form at even powers
     (`Transfer/MomentForm.lean`).
   * As $M \to \infty$, the Trotter traces converge to $Z_L(\beta)$ by a matrix product formula proved
     from scratch (`Transfer/ProductLimit.lean`).
   * Because the spatial moment inequalities are closed (non-strict) relations among finitely many
     traces, they pass to the $M \to \infty$ limit (`RPRing.spatialMomentBounds`). Combining this with
     Layer 2 yields `RPRing.gapGT_of_chainCertificate`.
4. **Concrete spin models** (`Heisenberg/*`, `Models/*`). Each model is instantiated as an `RPRing` in a
   real basis, paired with its seed hypotheses and kernel-checked `ChainCertificate`s to obtain its
   uniform gap theorem. Models that are reflection positive only after a sublattice unitary rotation are
   linked to their standard physical Hamiltonians via `Ring/Conjugation.lean`, which proves that
   sitewise unitarily conjugate Hamiltonians have identical spectra (`physical_gap`).
5. **Kernel-checked seed arithmetic** (`SeedKernel/*`). The seed inequalities carry the model-specific
   finite-size physics, which the bootstrap cannot produce. The opt-in library `SeedKernel` verifies
   seed bounds directly inside the Lean kernel using `decide +kernel` (backed by GMP integer
   arithmetic, with no `native_decide`): Julia generates witnesses and the kernel checks them.

### Summary of proved vs. open statements

| component | status |
|---|---|
| two-purity criterion (Layers 1–2), spatial moment inequalities for every `RPRing` (Layer 3), and unitary frame conjugations | completely proved; depends only on the standard Lean axioms `propext`, `Classical.choice`, and `Quot.sound` |
| uniform gap theorems for the spin-1 Heisenberg ring and 8 additional model families (table below) | proved **conditional on their seed lemmas** (`sorryAx` enters exclusively through the named seed declarations) |
| the finite seed lemmas | one proved: the longitudinal-field Ising chain at `g = 1/4`, `h = 1/20` (`SeedKernel.Ising.first_point_seed`, giving the unconditional `ising_gap`). All others remain open; [`../models/lean_seeds/`](../models/lean_seeds/) has the plan and the open obligations of the spin-1 seed |

### Build instructions

```sh
cd lean
lake exe cache get          # download prebuilt Mathlib oleans (never build Mathlib from source)
lake build                  # build CertifiedGaps (≈4 min from scratch; only warnings are the seed sorrys)
lake build SeedKernel       # opt-in: kernel-checked seed arithmetic (≈80 s)
lake build SeedKernelBench  # opt-in: kernel benchmarks up to n = 12 (≈2 min, up to 4 GB RAM per file)
```

Lean is pinned at `v4.36.0-rc1` and Mathlib at commit
`807d3836b165983210ee2ecf3cd28f239436949c` (October 19, 2026). `CertifiedGaps` comprises
roughly 6,800 lines of Lean
and `SeedKernel` roughly 1,200 lines (excluding generated witness files). All declarations live in
the `CertifiedGaps` namespace; declaration names below omit that prefix.

### Navigation

* [Status](#status): table of all model gap theorems, their open seed hypotheses, and the
  `sorry`-free theorems.
* [SeedKernel](#seedkernel-kernel-checked-seed-arithmetic): scope and file map of the kernel-checked
  seed library.
* [Adding a model](#adding-a-model): six-step recipe for formalizing a new model certificate, and
  the naming scheme.
* [Files](#files): file-by-file guide to `CertifiedGaps`.
* [Why each piece is needed in Lean](#why-each-piece-is-needed-in-lean): how the formal proof adapts
  the continuous-time MPO argument of the LaTeX notes.

## Status

`GapGT e c` states that the minimum of the eigenvalue family `e` is nondegenerate and every other
eigenvalue lies strictly more than `c` above it. Every model theorem below holds for all even ring
lengths `L ≥ n₀` (where `L` counts unit cells).

### Model theorems (conditional on two seed lemmas)

| Theorem | Model | Certified bound | Open seed hypotheses (the only `sorry`s) |
|---|---|---|---|
| `SpinOneHeisenberg.spin_one_heisenberg_gap` | spin-1 Heisenberg ring | `Δ_L > log 20 / 784` for `L ≥ 2304` | `seed36_spatial`: `1 − P_s(36, 49/4) ≤ 0.047916`; `seed36_thermal`: `1 − P_τ(72, 49/4) ≤ 0.181521` (then six doubling steps) |
| `SpinOneHeisenberg.spin_one_heisenberg_gap_sixty` | spin-1 Heisenberg ring | `Δ_L > (4/105) log (80/79) ≈ 0.000479` for `L ≥ 60` | `seed60_spatial`: `P_s(60, 105/4) ≥ 7/8`; `seed60_thermal`: `P_τ(120, 105/4) ≥ 7/8` |
| `SpinOneLargeD.gap` | spin-1 Heisenberg `+ D ∑ (Sᶻ)²`, `D = 4, 1603/500` | `log 4.755 / (9/4) ≈ 0.693`, `log 1.0065429 / (381/125) ≈ 0.00214`; `L ≥ 6` | `seeds` at `(6, β₀)`: `P_s ≥ 1 − 0.06411`, `P_τ ≥ 1 − 0.04053` at `D = 4` |
| `SpinOneDimerized.gap` | spin-1 dimerized chain, `δ = 9/10` (cells = weak pairs) | `log 6.339 / (7/2) ≈ 0.528`, `L ≥ 2` | `seeds`: `P_s(2, 7/2) ≥ 1 − 0.04972`, `P_τ(4, 7/2) ≥ 1 − 0.03206` |
| `HeisenbergLadder.physical_gap` (rotated: `gap`) | spin-½ ladder `∑_{legs} S·S + 4 ∑_{rungs} S·S` | `log 9.033 / 2 ≈ 1.100`, `L ≥ 4` rungs | `seeds`: `P_s(4, 2) ≥ 1 − 0.03509`, `P_τ(8, 2) ≥ 1 − 0.02463` |
| `StaggeredHeisenberg.physical_gap` (rotated: `gap`) | `∑ S·S − h ∑ (−1)ʲ Sᶻⱼ`, `h = 1, 3/4, 1/2, 7/20, 1/4` | `1.147, 0.890, 0.594, 0.383, 0.214`, `L ≥ 8` | `seeds` at `(8, β₀)` |
| `DimerizedHeisenberg.physical_gap` (rotated: `gap`) | spin-½ dimerized chain, `δ = 1/2, 3/10, 1/4` (cells = strong pairs) | `0.532, 0.172, 0.057`, `L ≥ 4` | `seeds` at `(4, β₀)` |
| `TwoLegIsing.gap` | two-leg TFIM ladder, `g = 3` | `log 2.585 / (3/2) ≈ 0.633`, `L ≥ 4` rungs | `seeds`: `P_s(4, 3/2) ≥ 1 − 0.1125`, `P_τ(8, 3/2) ≥ 1 − 0.05162` |
| `IsingLongitudinal.gap` | Ising chain `−ZZ − gX − hZ`, 30 points `g ∈ {1/4, …, 2}`, `h ∈ {1/20, …, 1}` | `log r / β₀` per point (`0.075` at `g = 1, h = 1/20` up to `3.70`), `L ≥ 6` | `seeds` (single lemma covering all 30 points) |
| `SpinOneHeisenberg.spin_one_boundary_field_gap` | open spin-1 Heisenberg chain, equal end fields `h = 3/5`, odd `n` sites | `Δ_n > log 10 / 392 ≈ 0.005874` for odd `n ≥ 1921` | four inputs at `β = 49`: `P_s(120) ≥ 0.948`, `P_τ(240) ≥ 0.948` (ring), `e₁₄₂ ≤ 0.14`, `P_τ^{opp}(240) ≥ 0.79` (open chain) |
| `SpinOneHeisenberg.spin_one_boundary_field_gap_sharp` | the same chain | `Δ_n > log(10⁶/2241) / 392 ≈ 0.01556` for odd `n ≥ 1921` | the same four inputs (`sharpBoundaryFieldCertificate`) |

Every model except the spin-1 Heisenberg ring lists its parameter points in `points`; `gap hp` (for
`hp : p ∈ points`) gives `Δ_L > log r / β₀` with `r = p.cert.r`, and the model's single `sorry`,
`seeds`, states both seed inequalities at every point. Running `#print axioms` on any of these
theorems reports only `propext`, `Classical.choice`, `Quot.sound`, and `sorryAx` (entering solely
through the seed hypotheses). The `physical_gap` theorems apply to the Hamiltonian in the standard
physical spin basis, while the rotated `gap` theorems apply to the unitarily rotated
reflection-positive frame, which has the same spectrum.

Note that the Lean model files record the seed bounds that were current when each file was formalized.
Where a model certificate in `../models/` was subsequently tightened (for example, `models/spin1_largeD`
now also certifies a sharper bound `0.899` at `D = 4` from an `n₀ = 12` seed, and
`models/dimerized_heisenberg` has slightly tighter enclosures than the rounded decimals in
`DimerizedHeisenberg.lean`), the weaker seed inequalities assumed in Lean remain a fortiori valid.

### Unconditional theorems (no `sorry`)

Every theorem in the table below is completely proved using only `propext`, `Classical.choice`, and
`Quot.sound`:

| Theorem | Statement |
|---|---|
| `gapGT_of_purities` | Proposition 1 of the main note: spatial moment form and purities `P_s(n₀, β₀), P_τ(2n₀, β₀) ≥ 0.99` imply `Δ_L > log 20 / β₀` for all even `L ≥ n₀` |
| `gapGT_of_coefficient` | OpenAI's Proposition 4.3: seed defects `≤ u₀ < 1/6`, `G(u₀) ≤ C`, `r = C u₀ < 1`, `C (1 - 3u₀) > 3` imply `Δ_L > log (1/r) / β₀`; `gapGT_of_seeds` is the case `C = 5` |
| `gapGT_of_seeds` | general equal-defect criterion under `SpatialMomentBounds`: seed defects `≤ ε₀ ≤ 1/100` imply `Δ_L > log (1/(5ε₀)) / β₀` |
| `gapGT_of_chain`, `gapGT_of_chainCertificate` | unequal-defect criterion via the exact asymmetric doubling recurrence; `ChainCertificate` validity is decided in the kernel |
| `RPRing.spatialMomentBounds` | every reflection-positive ring (`RPRing`) satisfies `SpatialMomentBounds` (proved via Trotter transfer matrices and the product formula) |
| `RPRing.gapGT_of_chainCertificate` | end-to-end criterion for an `RPRing`: seed bounds plus a valid `ChainCertificate` imply `Δ_L > log r / β₀` |
| `gapGT_of_openChain`, `gapGT_of_openCertificate` | the boundary-field criterion of Appendices A–C for abstract ring and open-chain spectra under `OpenMomentBounds`, with general seeds at `(2N₀, 4N₀, 2K₀, β₀)` and `m` base steps: `Δ_n > log(1/(10x)) / (β₀ 2^m)` for odd `n ≥ 4 N₀ 2^m`; the certificate form checks the rational data by `decide +kernel` |
| `gapGT_of_openSeeds` | the case of the note's eq. (C1) (`boundaryFieldCertificate`): four seeds imply `Δ_n > log 10 / 392` for odd `n ≥ 1921` |
| `RPRing.openMomentBounds`, `RPRing.openGapGT_of_certificate`, `RPRing.openGapGT_of_seeds` | every reflection-positive ring without internal reflection and every boundary term `b = -∑ w_α C_α` satisfy `OpenMomentBounds` (Trotter transfer matrix `X = F J` and the product formula); the end-to-end open-chain criterion |
| `SpinOneHeisenberg.gapGT_of_boundaryFieldInputs_certificate` | the spin-1 boundary-field gap from the four packaged inputs `BoundaryFieldInputs` and any valid certificate for those seeds; `gapGT_of_boundaryFieldInputs` (note's certificate) and `gapGT_of_boundaryFieldInputs_sharp` (`sharpBoundaryFieldCertificate`) are its cases |
| `RPRing.eigenvalues_eq_spectrum` | any Hamiltonian conjugate to `R.hamiltonian L` by a sitewise unitary has the exact same eigenvalue spectrum |
| `SpinOneHeisenberg.gapGT_log_twenty` | Theorem 1 of the main note deduced from packaged `ThermalMPOInputs`: `Δ_L > log 20 / 784` |
| `SpinOneHeisenberg.gapGT_sixty_of_seeds`, `SpinOneHeisenberg.seeds2304_of_seeds36` | the `L ≥ 60` bound from its two seeds; the seeds at `(2304, 784)` from those at `(36, 49/4)` via `chain36` (`ChainCertificate.StepsValid`, decided in the kernel) |
| `*.physicalSpectrum_eq` | equality of spectra between each physical Hamiltonian and its rotated reflection-positive counterpart |
| `*.points_valid` | kernel verification (`decide +kernel`) of every model's rational `ChainCertificate` |

The model seeds come from the certificates in [`../models/`](../models/) (complete spectra of the
seed rings with rigorous eigenvalue enclosures; each model file names its directory). The spin-1
Heisenberg seeds are Proposition 5.5 of the OpenAI preprint, reproduced in
[`../models/spin1/`](../models/spin1/).

## SeedKernel: kernel-checked seed arithmetic

`SeedKernel` is an opt-in library (`lake build SeedKernel`; not built by `lake build`) for proving
finite seed inequalities inside the Lean kernel. Every computational step is discharged by
`decide +kernel` on `Nat` or `Int` (evaluated by GMP inside the kernel), without `native_decide`,
`sorry`, or additional axioms. `SeedKernel/Audit.lean` (a root of the library) verifies via
`#print axioms` that all 40 main declarations depend at most on `propext`, `Classical.choice`, and
`Quot.sound`. The corresponding Julia generators, witness inventory, and verification driver
(`check.jl`: Julia replay, `lake build SeedKernel`, and the axiom audit) live in
[`../models/lean_seeds/`](../models/lean_seeds/), whose README states the results and the open
obligations.

**Scope of what is currently proved in `SeedKernel`** (one seed is closed, for the Ising chain; the
spin-1 work is arithmetic only):

* **Shared packed-integer lanes (`SeedKernel/Packed.lean`, `PackingProofs.lean`).** Packs a vector
  of `d` natural numbers into a single `Nat` in radix `2^w`, providing lane extraction, shifts,
  butterfly bit-reversal, and a single-multiplication dot product (`dot` / `gram`), together with
  radix encoding theorems (`radixPack_extract`, `pack_eq_radixPack`) and kernel-checked test cases.
  The general semantic equivalence between packed operations and unpacked sums is formulated as
  explicit target propositions (`DotCorrect`, etc.) in `SeedKernel/Obligations.lean`, which nothing
  assumes.
* **Spin-1 Heisenberg (`SeedKernel/Spin1/*`, Std only).** Implements the packed Hamiltonian action
  `hApply` and shifted power filter $nI - H_M$ on full $3^n$-lane vectors, checks agreement with an
  independent row-wise specification on small rings, and verifies exact integer Ritz/residual
  witnesses (`h6_ritz_witness`, `h8_ritz_witness`, `h10_ritz_witness`) for the filtered Néel trial
  vector (integer facts, not eigenvalue theorems), along with a conditional one-Ritz-value
  trace-counting lemma (`Counting.lean`).
* **Longitudinal Ising at $g = 1/4, h = 1/20$ (`SeedKernel/Ising/*`).** Verifies exact integer Gram,
  residual, and anchor inequalities for two generated approximate eigenvectors on the 6-site and
  12-site rings (both in packed form and, for 12 sites, in 64 direct chunks with a chunk-coverage
  theorem and negative rejection tests), and turns them into a proof of the seed:
  - `Spectral`: for a Hermitian matrix, a residual bound `∃ i, (λᵢ − θ)² ‖x‖² ≤ ‖Ax − θx‖²` and
    counting by comparison (`Re ⟪x, Ax⟫ ≥ ρ ‖x‖²` off a set `B` ⇒ at most `|B|` eigenvalues below
    `ρ`), derived from the unitary diagonalization;
  - `Hamiltonian`: `IsingLongitudinal.model g h` acts as `E(σ) x(σ) − g ∑ⱼ x(flipⱼ σ)`, and
    `Re ⟪x, Hx⟫ ≥ ∑ (E(σ) − gL) |x(σ)|²`;
  - `Bridge`: the checker's bit-indexed sums are exactly these quadratic forms (`xor 2ʲ` is a spin
    flip);
  - `Seed`: one eigenvalue lies within `10⁻¹⁰` of each candidate and at most two lie below the
    anchor, hence two-sided bounds on `Z₆(β)`, `Z₁₂(β)`, `Z₁₂(2β)` and the theorems `seed_spatial`,
    `seed_thermal`, `first_point_seed`, and `ising_gap`.

| File | Role |
|---|---|
| `SeedKernel.lean` | Library root importing all `SeedKernel` modules |
| `SeedKernel/Packed.lean` | Shared packed-`Nat` lane operations: `lane`, `repunit`, `pack`, `move`, `Signed`, floors, `reverseLanes`, `dot`, `gram`, and fast variants |
| `SeedKernel/PackingProofs.lean` | Radix encoding proofs: lane extraction and `pack = radixPack` |
| `SeedKernel/Obligations.lean`, `SeedKernel/PackedTests.lean` | Target correctness propositions for packed operations and kernel-checked examples |
| `SeedKernel/Spin1/Packed.lean`, `Spin1/Optimized.lean` | Packed spin-1 `hApply`, `filterApply`, trial moments, `ritzCheck`, and prepared-mask variants |
| `SeedKernel/Spin1/Obligations.lean`, `Spin1/Tests.lean` | Row-wise specification, open semantic propositions (`HApplyCorrect`, `FilterCorrect`, ...), and small-ring checks |
| `SeedKernel/Spin1/Counting.lean` | Conditional one-Ritz-value trace-counting lemma and formal counterexamples to naive two-value variants |
| `SeedKernel/Spin1/H6.lean`, `H8.lean`, `H10.lean` | Generated Ritz/residual witnesses (`models/lean_seeds/certify.jl spin1_kernel`) |
| `SeedKernel/Ising/Checker.lean`, `Ising/Packed.lean`, `Ising/Chunked.lean` | Integer $20H$ action and direct, packed (on the shared lanes), and 64-chunk certificate checkers |
| `SeedKernel/Ising/Generated.lean`, `Ising/GeneratedChunks.lean` | Generated eigenvector proposals and kernel checks (`models/lean_seeds/certify.jl ising`) |
| `SeedKernel/Ising/ChunkCoverage.lean`, `Ising/NegativeChecks.lean` | Proof that the 64 chunks cover `{0, ..., 4095}` and rejection tests for corrupted witnesses |
| `SeedKernel/Ising/Spectral.lean` | Residual bound and counting by comparison for Hermitian matrices |
| `SeedKernel/Ising/Hamiltonian.lean`, `Ising/Bridge.lean` | The Ising Hamiltonian acting on vectors and its identification with the integer checker |
| `SeedKernel/Ising/Seed.lean` | The unconditional seed and `ising_gap` |
| `SeedKernel/Audit.lean` | `#print axioms` audit across all main declarations |
| `SeedKernel/Bench/*.lean` | Kernel performance benchmarks up to $n = 12$ (`lake build SeedKernelBench`) |

## Adding a model

A gap certificate enters `CertifiedGaps` in exactly this form:

1. **Instantiate the model as an `RPRing d ι`:** Specify real matrices `C_α` (each symmetric or
   antisymmetric), a real symmetric on-site matrix `o`, and optionally an internal cell reflection
   `refl` (a real symmetric involution commuting with `o`), so that the bond interaction is
   $h = -\sum_\alpha C_\alpha \otimes (\text{refl} \cdot C_\alpha \cdot \text{refl})$. A multi-site
   unit cell is represented as a single larger local space `d` with all intra-cell terms placed in
   `o`. The matrix symmetry conditions are finite decidable identities; the helper lemmas
   `transpose_kronecker_eq_or`, `transpose_sum_kronecker_self`, and
   `cellSwap_conj_smul_sum_kronecker_self` handle the standard cases.
2. **List the parameter points:** Define a `structure Point` holding the model parameters (as `ℚ`),
   the seed inverse temperature `β₀ : ℚ`, and a `cert : ChainCertificate`, and
   `def points : List Point`. Generate each certificate with
   `julia lean/scripts/chain_certificate.jl εₛ εₜ β₀ [J] [digits]` (the arguments are exact literals
   such as `a//b` or `0.064`; nothing is evaluated) and paste it in. Its first entries `s[0], t[0]`
   are the seed defects. `points_valid` and `points_β₀_pos` are checked by `decide +kernel`; the gap
   bound is `log r / β₀`.
3. **State the seeds:** Declare one `sorry` lemma
   `seeds : ∀ p ∈ points, 1 − s[0] ≤ (model p.x).spatialPurity n₀ p.β₀ ∧ 1 − t[0] ≤ (model p.x).thermalPurity (2 n₀) p.β₀`,
   with even `n₀ ≥ 2` written as a literal. These are the only statements left to numerics.
4. **Apply the criterion:** `gap hp hL hLe` is
   `(model p.x).gapGT_of_chainCertificate (n₀ := n₀) (by norm_num) (by decide) (by exact_mod_cast points_β₀_pos p hp) (points_valid p hp) ⟨rfl, rfl, rfl⟩ (seeds p hp).1 (seeds p hp).2 hL hLe`.
   See `Models/TwoLegIsing.lean` for a complete example (roughly 80 lines).
5. **Transfer across a unitary frame rotation (if needed):** When reflection positivity requires
   one-site unitaries `a`, define the physical Hamiltonian in `Models/<Name>Physical.lean`, prove
   that it equals `piTensor a * R.hamiltonian L * piTensor aᴴ` using `RPRing.conj_hamiltonian`, and
   transfer the spectral gap via `RPRing.eigenvalues_eq_spectrum`.
6. **Check the seeds against the certificate:** Add the model to `lean/scripts/check_seeds.jl` and
   run `julia --project=julia lean/scripts/check_seeds.jl` from the repository root. It checks that
   every Lean `s[0], t[0]` is at least the certified defect at the same `(n₀, β₀)`.

A certificate should therefore report its seeds as exact rationals together with `(n₀, β₀)` and the
exact `C_α, o`. Building blocks shared by several models (one-site matrices, the cell swap) go in
`Local/`, not in `Models/`.

**Naming scheme.** Every model uses the same names, so that each model file reads the same way:

| Name | Meaning |
|---|---|
| `Models/<Name>.lean`, namespace `CertifiedGaps.<Name>` | `<Name>` is the Julia directory `models/<slug>` in UpperCamelCase (table below) |
| `model` | the model as an `RPRing`, a function of its real parameters |
| `Point`, `points` | the certified parameter points with `β₀` and `cert` |
| `points_valid`, `points_β₀_pos`, `points_<field>_…` | decidable facts about all points (`decide +kernel`) |
| `seeds` | the seed inequalities at every point: the model's single `sorry` |
| `gap` | the gap at every point, in the frame of `model` |
| `physicalHamiltonian`, `physicalSpectrum`, `physicalHamiltonian_eq`, `physicalSpectrum_eq`, `physical_gap` | `Models/<Name>Physical.lean`: the physical frame, when `model` is rotated |

| Lean namespace | Julia directory |
|---|---|
| `SpinOneHeisenberg` (`Heisenberg/`) | `models/spin1` (the headline theorem keeps its own names: `seed_spatial`, `seed_thermal`, `spin_one_heisenberg_gap`) |
| `SpinOneLargeD` | `models/spin1_largeD` |
| `SpinOneDimerized` | `models/spin1_dimerized` |
| `HeisenbergLadder` | `models/heisenberg_ladder` |
| `StaggeredHeisenberg` | `models/heisenberg_staggered_field` |
| `DimerizedHeisenberg` | `models/dimerized_heisenberg` |
| `TwoLegIsing` | `models/two_leg_tfim` |
| `IsingLongitudinal` | `models/tfim_longitudinal` |

## Files

| File | Role |
|---|---|
| `PurityCriterion.lean` | Abstract two-purity criterion: squaring inequality, thermal-spatial coupling, one exact doubling step for arbitrary defects, $`\beta`$-interpolation, the purity-to-gap bound, the criterion with a coefficient `C ≥ G(u₀)` (`gapGT_of_coefficient`) and its equal-defect case (seed defects $\le 1/100$, `C = 5`) |
| `Chain.lean` | Exact asymmetric doubling recurrence for unequal seed defects and `gapGT_of_chain` with bound `log r / β₀` |
| `ChainCheck.lean` | Decidable rational `ChainCertificate` structure discharging the hypotheses of `gapGT_of_chain`, and `StepsValid` for the purities along a chain alone |
| `Transfer/PiTensor.lean` | Sitewise tensor products (`piTensor`), expansion of $`\mathrm{Tr}\, T^N`$ as a sum over closed label paths and of matrix elements of ordered products as sums over open paths |
| `Transfer/MomentForm.lean` | Proof that traces of even powers of a real matrix factorization $B^T B \succeq 0$ satisfy the spatial moment form |
| `Transfer/ContinuousKernel/Basic.lean` | Actual bounded measurable kernel integral operator on `L²`, operator-norm bound, and Hermitian-kernel self-adjointness |
| `Transfer/ContinuousKernel/HilbertSchmidt.lean` | Basis-free Hilbert–Schmidt criterion; operator-norm convergence of finite-rank approximations and compactness |
| `Transfer/ContinuousKernel/Compact.lean`, `Audit.lean` | Integrated Bessel bound, Hilbert–Schmidt and compact kernel operators, axiom audit and degenerate/complex regression examples |
| `Transfer/ProductLimit.lean` | Lie–Trotter product formula: $\lim_{M\to\infty} (\prod_k (1 + (c/M) h_k))^M = e^{c \sum_k h_k}$ |
| `Ring/Model.lean` | Periodic ring geometry, `bondOp`/`siteOp`, `RPRing` structure, Hermiticity of $H_L$, and $\operatorname{Re} \operatorname{Tr} e^{-\beta H_L} = \sum_m e^{-\beta E_m}$ |
| `Ring/Transfer.lean` | Even-bond, on-site, and odd-bond Trotter slices and the reflection-positivity identity $\operatorname{Tr} S^M = \operatorname{Tr} (B^T B)^{L/2}$ with real $B$ |
| `Ring/Gap.lean` | Proof of `RPRing.spatialMomentBounds` and `RPRing.gapGT_of_chainCertificate` |
| `Ring/Conjugation.lean` | Sitewise unitary conjugation of ring Hamiltonians and spectral invariance |
| `Heisenberg/Hamiltonian.lean` | The spin-1 Heisenberg ring $H_L = \sum_j \mathbf{S}_j \cdot \mathbf{S}_{j+1}$ in the real Cartesian basis |
| `Heisenberg/Gap.lean`, `Heisenberg/Main.lean` | Both bounds of OpenAI's Theorem 1.1: `spin_one_heisenberg_gap` (`L ≥ 2304`, from the seeds at 36 sites and `chain36`) and `spin_one_heisenberg_gap_sixty` (`L ≥ 60`) |
| `Heisenberg/OpenChain.lean` | Theorem 2 of the main note (`spin_one_boundary_field_gap`) and its sharper form (`spin_one_boundary_field_gap_sharp`): the equal- and opposite-field open chains as open chains of the Heisenberg ring, and the four inputs |
| `Open/MomentBounds.lean` | `OpenMomentBounds`: closed inequalities of the boundary moment form, from a finite Hermitian transfer matrix and through limits; `leadingRoot`, `leadingOverlap`, `boundaryDefect` |
| `Open/Criterion.lean` | The boundary-field criterion `gapGT_of_openChain` for general seeds: the ring at the open-chain scales, the three moves of Appendix A, the base steps (`OpenStep`) and the uniform phase |
| `Open/Certificate.lean` | Rational certificates `OpenCertificate` checked by `decide +kernel`, `gapGT_of_openCertificate`, the note's `boundaryFieldCertificate` with `gapGT_of_openSeeds`, and `sharpBoundaryFieldCertificate` |
| `Open/Model.lean` | Open chains of an `RPRing` with boundary terms (`Boundary`, `openHamiltonian`, `openSpectrum`) |
| `Open/Transfer.lean` | Open-chain Trotter slices, time reversal `J`, and `Tr (slice)^M = ⟨v, Xⁿ v⟩` with `X = F J` |
| `Open/Bounds.lean`, `Open/Gap.lean` | `RPRing.openMomentBounds` and the end-to-end `RPRing.openGapGT_of_certificate` and `RPRing.openGapGT_of_seeds` |
| `Local/Pauli.lean`, `Local/SpinHalf.lean` | Real Pauli matrices $X, Z$ and real spin-½ factors $S^x, iS^y, S^z$ (`SpinHalf.P`) |
| `Local/SpinHalfRotation.lean` | Sublattice rotation $e^{i\pi S^y}$ and complex spin-½ operators |
| `Local/SpinOne.lean` | Real antisymmetric spin-1 generators `SpinOne.G α` $= G_\alpha = i S^\alpha$ |
| `Local/CellSwap.lean` | Two-site cell swap involution $F$ satisfying $F (A \otimes B) F = B \otimes A$ |
| `Models/SpinOneLargeD.lean` | Spin-1 chain with single-ion anisotropy $\sum \mathbf{S}_j \cdot \mathbf{S}_{j+1} + D \sum (S_j^z)^2$ as an `RPRing` ($D = 0$ is the Heisenberg ring), at $D = 4$ and $D = 1603/500$ |
| `Models/SpinOneDimerized.lean` | Dimerized spin-1 chain with internal cell reflection $F$ at $\delta = 9/10$ |
| `Models/TwoLegIsing.lean` | Two-leg transverse-field Ising ladder at $g = 3$ |
| `Models/IsingLongitudinal.lean` | Transverse-field Ising chain with longitudinal field across 30 parameter points |
| `Models/HeisenbergLadder.lean`, `Models/StaggeredHeisenberg.lean`, `Models/DimerizedHeisenberg.lean` | Reflection-positive rotated Hamiltonians, seeds, and gap theorems for the three spin-½ models |
| `Models/*Physical.lean` | Physical spin-½ Hamiltonians, unitary sublattice rotations, and transferred gap theorems |
| `scripts/chain_certificate.jl` | Julia generator that constructs an exact rational `ChainCertificate` from $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$ (exact literal parser; Lean re-checks) |
| `scripts/open_certificate.jl` | Julia generator that searches an exact rational `OpenCertificate` from the four open-chain seeds (reuses the parser of `chain_certificate.jl`; Lean re-checks) |
| `scripts/check_seeds.jl` | Checks every Lean seed defect against the certificate JSON in `../models/*/results` |

## Why each piece is needed in Lean

On paper, the main note obtains the spatial moment representation in a single conceptual step: it
constructs an exact continuous-time thermal MPO on Poisson bond histories whose spatial transfer operator
$X_\beta$ is a compact self-adjoint integral operator on $L^2(\mu_\beta)$, giving
$Z_n(\beta) = \operatorname{Tr} X_\beta^n = \sum_i |\lambda_i|^n$. Formalizing that infinite-dimensional
path-space construction directly in Lean still needs Poisson history measures, the concrete
chronologically merged history kernel, and the later cycle-trace identification. The generic
bounded-kernel operator bridge is now proved in `Transfer/ContinuousKernel/`: from joint strong
measurability, a pointwise bound $`|k(x,y)| \le B`$, and a finite measure, it constructs the actual
integral operator $`Kf(x)=\int k(x,y)f(y)\,d\mu(y)`$ on real or complex $L^2$. It proves
$`\|K\| \le B\mu(\Omega)`$, the squared Hilbert–Schmidt bound
$`B^2\mu(\Omega)^2`$, and compactness; Hermitian symmetry of the kernel implies self-adjointness.
The Hilbert–Schmidt property is expressed by a uniform bound on squared image norms of finite
orthonormal families, and compactness is derived from operator-norm convergence of rank-one
truncations. No separability, positive total mass, or pre-existing compactness assumption is needed.
For the paper's bound $B=3$, this gives $`9\mu(\Omega)^2`$; identifying the physical history kernel
and applying its scalar normalization remain separate steps. This layer asserts no cycle-trace,
partition-function, moment, twist, or sector identity.

The current end-to-end gap formalization continues to use the finite-dimensional Trotterization
route, which needs only finite matrix algebra and the matrix exponential limit. Each module below
addresses a specific distinction between the paper proof and the formal proof:

**Spectral gap as a predicate on unordered levels (`GapGT`).** The paper writes $\Delta_L = E_1 - E_0$
with eigenvalues sorted in nondecreasing order. In Mathlib, the eigenvalues of a Hermitian matrix are
indexed by the finite basis type rather than returned as a sorted list. Formulating the conclusion as
`GapGT e c` ("there exists a ground level `m₀` such that every other level `m ≠ m₀` satisfies
`e m₀ + c < e m`") avoids sorting permutations and simultaneously captures both ground-state uniqueness
and the spectral gap bound.

**Connecting matrix traces to eigenvalue sums (`partitionFunction_eq_sum`).** The transfer-operator
identities and seed purities are statements about matrix traces $\operatorname{Tr} e^{-\beta H_L}$,
whereas `GapGT` is a statement about the eigenvalues of $H_L$. Diagonalizing the Hermitian matrix $H_L$
and applying the spectral mapping theorem for the matrix exponential bridges these two representations.

**Spatial moment inequalities instead of an infinite-dimensional transfer operator (`SpatialMomentBounds`).**
Rather than constructing a single limiting operator $X_\beta$, `PurityCriterion.lean` axiomatizes the
three trace inequalities that the doubling recurrence actually uses:
* $Z_{2n}(\beta) \le Z_n(\beta)^2$;
* the spatial squaring inequality bounding $1 - P_{\mathrm{s}}(2n, \beta)$ by $f(1 - P_{\mathrm{s}}(n, \beta))$; and
* $\ell^p$ moment monotonicity across even lengths.

Because these are non-strict inequalities among a finite number of partition functions $Z_n(\beta)$,
they are preserved under pointwise limits of traces.

**Finite-dimensional Trotter transfer matrices (`Ring/Transfer.lean`, `Transfer/PiTensor.lean`).**
Replacing $e^{-\beta H_L}$ by $M$ first-order Trotter slices (even-bond factors $1 - \tau h_{2k}$,
on-site factors $1 - \tau o_j$, and odd-bond factors $1 - \tau h_{2k+1}$, with $\tau = \beta/M$) produces
a finite-dimensional tensor network. Expanding each bond factor as
$1 \otimes 1 + \tau \sum_\alpha C_\alpha \otimes F C_\alpha F$ and regrouping the bond indices around
even and odd unit cells rewrites the $M$-slice Trotter trace on an even ring of $L = 2N$ cells as
$\operatorname{Tr}(B^T B)^N$. Because the factors $C_\alpha$, the on-site term $o$, and the cell
reflection $F$ have real entries, the transfer block $B$ is a real matrix, so $B^T B = B^\dagger B$ is
Hermitian positive semidefinite. This is why every `RPRing` is formulated in a real local basis (such as
the Cartesian basis $G_\alpha = iS^\alpha$ for spin-1 or the $\{S^x, iS^y, S^z\}$ basis for rotated
spin-½ models).

**From positive semidefiniteness to the moment form (`Transfer/MomentForm.lean`).** Because $T = B^T B$
is positive semidefinite, its eigenvalues $\mu_i$ are real and nonnegative, so
$\operatorname{Tr} T^N = \sum_i (\sqrt{\mu_i})^{2N}$ satisfies the spatial moment inequalities at every
finite Trotter number $M$.

**Passing to the continuous-time limit (`Transfer/ProductLimit.lean`, `SpatialMomentBounds.of_tendsto_hasMomentForm`).**
`ProductLimit.lean` proves that the $M$-slice first-order Trotter product converges to $e^{-\beta H_L}$
as $M \to \infty$ for each fixed ring length $L$. Because the limit is taken in the trace inequalities
 *before* running the doubling bootstrap, no uniform-in-$L$ Trotter error bound or eigenvalue
perturbation argument is required.

**Why the seed inequalities cannot be eliminated (`seed_spatial`, `seed_thermal`).** The spatial moment
inequalities hold for *every* reflection-positive ring, including gapless chains (such as the spin-½
Heisenberg ring) and chains with spontaneous symmetry breaking (such as the classical ferromagnetic Ising
ring, which has $P_{\mathrm{s}} \to 1$ but $P_\tau \to 1/2$ and a twofold degenerate ground state). The
two finite-scale purity bounds $P_{\mathrm{s}}(n_0, \beta_0) \ge 1 - \varepsilon_{\mathrm{s}}$ and
$P_\tau(2n_0, \beta_0) \ge 1 - \varepsilon_\tau$ inside the admissible basin are the physical inputs
that distinguish a gapped, unique-ground-state phase from a gapless or symmetry-broken one.

**Decidable rational certificates (`Chain.lean`, `ChainCheck.lean`).** Writing the target gap bound as
$\log r / \beta_0$ for a rational base $r > 1$ converts the transcendental inequality
$\Delta_L > \log r / \beta_0$ into purely rational inequalities along the $J$ doubling steps
($r^{2^j} z_j < 1 - z_j$ and $r^{2^J} \cdot 5\varepsilon_J \le 1$). Lean's kernel verifies those rational
inequalities in seconds via `decide +kernel` without needing interval arithmetic for transcendental
functions inside Lean.
