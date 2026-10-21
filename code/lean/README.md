# CertifiedGaps: Lean 4 formalization of spectral gaps via the spacetime bootstrap

This Lean 4 + Mathlib package formalizes the spacetime-bootstrap proof that reflection-positive quantum
spin rings have a uniform spectral gap for all sufficiently large even ring lengths, following the main
note ([`../../paper/thermal-mpo-gap-note.tex`](../../paper/thermal-mpo-gap-note.tex)) and scope note
([`../../paper/thermal-mpo-gap-scope.tex`](../../paper/thermal-mpo-gap-scope.tex)). **Every analytic and
algebraic step is proved from standard axioms except the finite numerical seed inputs** (two purity
inequalities per model at the seed scale $(n_0, \beta_0)$). Those finite seed bounds appear as explicit
named `sorry` lemmas whose values come from the computer certificates in [`../models/`](../models/).
Two gaps are now proved outright in the kernel: `SeedKernel.Ising.ising_gap` (Ising rings) and
`SeedKernel.IsingOpen.ising_free_gap` (Ising chains with free ends) are unconditional (below). The
open-chain argument of the note (Appendices A–C) is formalized as well: it gives a uniform gap for
open chains, with boundary fields or free ends, from four finite inputs.

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
`sorry`. `SeedKernel/Heisenberg/Integers.lean` proves the four inputs themselves from the
paper's printed integers (see "The spin-1 gaps from the printed integers" below):
`spin_one_boundary_field_gap_of_integers` has no `sorry` and two hypotheses: the eighteen
totals of the periodic table and the six boundary integers with `n = 5, 6`.

The criterion itself is general. `gapGT_of_openChain` takes any seed scale `β₀`, ring seeds at
lengths `2N₀` and `4N₀`, a boundary seed at length `2K₀`, an open-chain seed at length `4N₀`, and
`m` base steps, and gives `Δ_n > log(1/(10x)) / (β₀ 2^m)` for every `n ≥ 4 N₀ 2^m`, odd or even
(`oddChain_purity`, `evenChain_purity`). Its finite data
form an `OpenCertificate`, whose rational conditions `OpenCertificate.Valid` are checked by
`decide +kernel`; `gapGT_of_openCertificate` applies one. The script
[`scripts/open_certificate.jl`](scripts/open_certificate.jl) searches for the data from the seeds,
and Lean re-checks it. The note's numbers are `boundaryFieldCertificate`. From the same four
inputs, the script finds `sharpBoundaryFieldCertificate`, which hands over at the defect actually
reached (`x = 2.241 · 10⁻⁴` in place of `1/100`). The resulting
`spin_one_boundary_field_gap_sharp` gives `Δ_n > log(10⁶/2241)/392 ≈ 0.01556` for every odd
`n ≥ 1921`, 2.65 times the note's bound. At even lengths the same open chain is the opposite-field
chain, and `spin_one_opposite_field_gap(_sharp)` gives the same two bounds for every even
`n ≥ 1920`.

All four seeds can be stated for finite partition functions. The boundary defect involves the
leading eigenvalue and overlap of the transfer operator, but `boundaryDefect_le_of_traces` bounds
it by three finite traces: `G(2K - 2N₀ r) (d(u₀) Z(2N₀))^r ≤ E G(2K)` implies `e_{2K} ≤ E`.
`gapGT_of_openCertificate_traces` and `RPRing.openGapGT_of_traces` use this form.

Open boundary conditions are the case `b = 0`. Every bond factor `C_α` of an `RPRing` is symmetric
or antisymmetric (`C_symm`), so `RPRing.freeBoundary` is always an admissible end term, its open chain is the chain with free ends at every length (`openHamiltonian_free`), and
`RPRing.freeGapGT_of_traces` gives its gap at every `n ≥ 4 N₀ 2^m` from a certificate and four
finite seeds. No model instance with free ends is formalized yet: the free spin-1 Heisenberg chain
has spin-½ edge states whose splitting closes exponentially in `n`, so it has no uniform gap; that
is why the boundary fields are there.

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

### An unconditional gap with free ends

```lean
-- SeedKernel/IsingOpen/Seed.lean (lake build SeedKernel)
theorem SeedKernel.IsingOpen.ising_free_gap {n : ℕ} (hn : 24 ≤ n) :
    GapGT (isingModel.openSpectrum isingModel.freeBoundary n) (log (10 ^ 12 / 3801) / (97 / 2))
```

For every length `n ≥ 24`, odd or even, the same Ising chain with open boundary conditions,
$`-\sum_{j=1}^{n-1} Z_j Z_{j+1} - \tfrac14 \sum_j X_j - \tfrac{1}{20} \sum_j Z_j`$, has a unique ground state
and a spectral gap `> 0.3998` (and `> 0.4723` for `n ≥ 96`, `ising_free_gap_long`), again with no
`sorry` and only the standard axioms. It is `RPRing.freeGapGT_of_traces` with one base step
(`isingFreeCertificate`). Its four seeds come from low spectra of the rings of 6 and 12 sites and
the free chains of 2, 8 and 12 sites (kernel-checked residual, counting and comparison-weight
checks), which do not depend on `β`; the general seed theorems of `SeedKernel.Seeds` turn them into
the seeds at any rational `β₀` through one decidable rational check (`free_gap_of_check`).

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
| the finite seed lemmas | one proved: the longitudinal-field Ising chain at `g = 1/4`, `h = 1/20` (`SeedKernel.Ising.first_point_seed`, giving the unconditional `ising_gap`), together with the open-chain seeds of the same chain with free ends (`SeedKernel.IsingOpen`, giving `ising_free_gap`). All others remain open. The spin-1 seeds of both OpenAI papers reduce to the printed integer totals of their trace tables (`SeedKernel/Heisenberg/Integers.lean`); [`../models/lean_seeds/`](../models/lean_seeds/) has the plan and the open obligations of the spin-1 seed |

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
| `SpinOneHeisenberg.spin_one_opposite_field_gap`, `spin_one_opposite_field_gap_sharp` | open spin-1 Heisenberg chain, opposite end fields `±3/5`, even `n` sites | `Δ_n > log 10 / 392` and `Δ_n > log(10⁶/2241) / 392` for even `n ≥ 1920` | the same four inputs |

Every model except the spin-1 Heisenberg ring lists its parameter points in `points`; `gap hp` (for
`hp : p ∈ points`) gives `Δ_L > log r / β₀` with `r = p.cert.r`, and the model's single `sorry`,
`seeds`, states both seed inequalities at every point. Running `#print axioms` on any of these
theorems reports only `propext`, `Classical.choice`, `Quot.sound`, and `sorryAx` (entering solely
through the seed hypotheses). The `physical_gap` theorems apply to the Hamiltonian in the standard
physical spin basis, while the rotated `gap` theorems apply to the unitarily rotated
reflection-positive frame, which has the same spectrum.

Note that the Lean model files record the seed bounds that were current when each file was formalized.
Where a model certificate in `../models/` was subsequently tightened (for example, `code/models/spin1_largeD`
now also certifies a sharper bound `0.899` at `D = 4` from an `n₀ = 12` seed, and
`code/models/dimerized_heisenberg` has slightly tighter enclosures than the rounded decimals in
`DimerizedHeisenberg.lean`), the weaker seed inequalities assumed in Lean remain a fortiori valid.

### The spin-1 gap from its trace tables (no `sorry`)

The four spin-1 seed `sorry`s of `Heisenberg/Main.lean` reduce to two tables of short-ring
partition functions. Everything else in OpenAI's two initializations is proved: the transfer
spectrum in the `D₂` and `I/O/N` sectors at every ring length (odd included), the polynomial
filters, capped masses and dominant eigenvalue, the temperature lift of Proposition 5.4, and the
integer matrix-product trial state on 72 and 120 sites, checked in the kernel
(`SeedKernel/Heisenberg/Trial*.lean`). The theorems below take the table as their only hypothesis and
use only `propext`, `Classical.choice`, `Quot.sound` (`SeedKernel/Audit.lean`).

| Theorem (`SpinOneHeisenberg.`) | Hypothesis | Bound |
|---|---|---|
| `gap_of_table` | `TraceTable`: `Z_n(49/4)`, `Z_n(49/4, P)`, `n = 4, …, 12`, within `3 × 10⁻⁵` | `Δ_L > log(1.05882)/(49/4) ≈ 0.00467` for every even `L ≥ 36` |
| `gap_log_twenty_of_table` | the same | OpenAI's eq. (1.4): `Δ_L > log 20 / 784` for `L ≥ 2304` |
| `gap_sixty_chain_of_table` | `TraceTable60`: `Z_n(21/2, g)`, `g = 1, P, C`, `n = 6, 8, 10`, within `10⁻⁵` | `Δ_L > log(1.16321)/(105/4) ≈ 0.00576` for every even `L ≥ 60` |
| `gap_sixty_of_table` | the same | OpenAI's eq. (1.3): `Δ_L > (4/105) log(80/79)` for `L ≥ 60` |

The first and third bounds come from the exact unequal-defect recurrence (`chain36Gap`, `chain60`)
applied to OpenAI's own seed values; their eq. (1.3) and (1.4) use a single coefficient and are
weaker (`Main.lean`).

### The spin-1 gaps from the printed integers

Both trace tables follow from the integers the papers print. The papers compute each trace with
the same algorithm, and Lean formalizes it:
1. Shift the Hamiltonian to a contraction `R` (energy bounds).
2. Truncate the Poisson expansion of `e^{m(R - 1)}`.
3. Run a floored integer Horner recursion on each representative column.
4. Sum the integer totals over orbit representatives with multiplicities.

Lean proves that the computable integer recurrence, if it returns the printed totals, implies the
rational enclosures of the table (`ThermalTable.lean`, `ThermalTable60.lean`, the `Trace/` layer).
The kernel cannot evaluate these recurrences at the paper's sizes, so the totals remain
hypotheses; every other step, including all trial states, is proved
(`SeedKernel/Heisenberg/Integers.lean`, standard axioms only). For the `b = 49/4` table,
`ThermalPacked.lean` evaluates the recurrence with packed big integers on the representative
columns and proves the result equal to `hornerTotal` (`hornerTotal_eq_of_packed`). With it the
kernel checks the two totals with `n = 4` (`SeedKernel/Heisenberg/ThermalTraces.lean`), and
compiled evaluation (`SeedKernel/Eval/ThermalTotals.lean`, not a proof) reproduces the twelve
with `n ≤ 9`. For the `b = 21/2` table, `ThermalPacked60.lean` reduces all three twists to the
same packed run (the twist `C` through the real form of its `ℤ[ω]` matrix) and proves
`hornerTotals60_of_packed`, and compiled evaluation (`SeedKernel/Eval/ThermalTotals60.lean`)
reproduces the six totals with `n = 6, 8`.

| Theorem (`SpinOneHeisenberg.`) | Hypotheses | Bound |
|---|---|---|
| `gap_of_hornerTotals`, `gap_log_twenty_of_hornerTotals` | `HornerTotals`: the eighteen totals `W_{n,g}` of the `b = 49/4` table | as `gap_of_table`, `gap_log_twenty_of_table` |
| `gap_of_hornerTotalsLarge`, `gap_log_twenty_of_hornerTotalsLarge` | `HornerTotalsLarge`: the sixteen totals with `n ≥ 5` (the two with `n = 4` are kernel-checked) | as `gap_of_table`, `gap_log_twenty_of_table` |
| `gap_sixty_of_hornerTotals60` | `HornerTotals60`: the nine totals of the `b = 21/2` table (`ℤ[ω]` arithmetic for the `C` twist) | as `gap_sixty_chain_of_table` |
| `spin_one_boundary_field_gap_of_integers`, `spin_one_opposite_field_gap_of_integers` | `HornerTotals`, `BoundaryTotalsLarge` (the six boundary integers `k_{n,g}` with `n = 5, 6`) | `Δ_n > log(10⁶/2241)/392` for odd `n ≥ 1921` (equal fields) and even `n ≥ 1920` (opposite fields) |
| `spin_one_boundary_field_gap_of_hornerTotalsLarge`, `spin_one_opposite_field_gap_of_hornerTotalsLarge` | `HornerTotalsLarge`, `BoundaryTotalsLarge` | as the previous row |

For the boundary-field paper, the chain from integers to inputs is the following:
1. Reseeding to `(120, 49)` (`Reseed120.lean`) gives `1 - P_s(120, 49) ≤ 0.047` and
   `1 - P_τ(240, 49) ≤ 0.0179`; the paper has `0.0502` and `0.0198`.
2. The initialization at `b₀ = 49/4` (`BoundaryInit.lean`, `BoundaryInitTable.lean`) yields
   `e₁₄₂ ≤ 0.1306` and the open thermal purity.
3. The fifteen short-chain traces `Z_{n,g}` (`BoundaryTraces.lean`, `BoundaryThermal.lean`) give
   `Y₆(b₀) < 0.027` and the leading-weight polynomial `< 0.003745`. The nine entries with `n ≤ 4`
   are checked by `decide +kernel`.
4. The four open-chain trial vectors on 34, 78, 142 and 240 sites are checked in the kernel
   (`SeedKernel/Heisenberg/OpenTrial*.lean`).
5. The spectral witness at `b₀` (`BoundaryLeading.lean`, `thermalCore_of_tables`) comes from both
   tables: on the single-site Trotter matrices, the `D₂` sector data of the seed at `(36, 49/4)`
   give a leading line with `λ₀ ≥ 0.999` and `ρ₀ ≤ 0.874`, fixed by the `A₄` average `Π`. The
   two-boundary transfer identifies the averaged pairings with `I(n)`, so that
   `c₀ P_*(λ₀)² ≤ 0.003745`, and closed inequalities pass all of it to the limit.

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
| `gapGT_of_openChain`, `gapGT_of_openCertificate` | the boundary-field criterion of Appendices A–C for abstract ring and open-chain spectra under `OpenMomentBounds`, with general seeds at `(2N₀, 4N₀, 2K₀, β₀)` and `m` base steps: `Δ_n > log(1/(10x)) / (β₀ 2^m)` for every `n ≥ 4 N₀ 2^m`, odd or even; the certificate form checks the rational data by `decide +kernel` |
| `boundaryDefect_le_of_traces`, `gapGT_of_openCertificate_traces` | the boundary seed from three finite traces, `G(2K - 2N₀ r) (d(u₀) Z(2N₀))^r ≤ E G(2K)`, so that all four seeds concern finite partition functions |
| `gapGT_of_openSeeds` | the case of the note's eq. (C1) (`boundaryFieldCertificate`): four seeds imply `Δ_n > log 10 / 392` for every `n ≥ 1920` |
| `RPRing.openMomentBounds`, `RPRing.openGapGT_of_certificate`, `RPRing.openGapGT_of_traces`, `RPRing.freeGapGT_of_traces`, `RPRing.openGapGT_of_seeds` | every reflection-positive ring without internal reflection and every boundary term `b = -∑ w_α C_α` satisfy `OpenMomentBounds` (Trotter transfer matrix `X = F J` and the product formula); the end-to-end open-chain criterion |
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
`#print axioms` that all 108 main declarations depend at most on `propext`, `Classical.choice`, and
`Quot.sound`. The corresponding Julia generators, witness inventory, and verification driver
(`check.jl`: Julia replay, `lake build SeedKernel`, and the axiom audit) live in
[`../models/lean_seeds/`](../models/lean_seeds/), whose README states the results and the open
obligations.

**Scope of what is currently proved in `SeedKernel`** (the Ising seeds are closed, for rings and for
free chains; the spin-1 work is arithmetic only):

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
  - the spectral input is the shared layer `CertifiedGaps/Spectral`: the residual bound
    `∃ i, (λᵢ − θ)² ‖x‖² ≤ ‖Ax − θx‖²` (`Spectral.exists_eigenvalue_near`) and counting by
    comparison (`Re ⟪x, Ax⟫ ≥ ρ ‖x‖²` off a set `B` ⇒ at most `|B|` eigenvalues below `ρ`,
    `Spectral.eigCount_le_card_of_vanishing`);
  - `Hamiltonian`: `IsingLongitudinal.model g h` acts as `E(σ) x(σ) − g ∑ⱼ x(flipⱼ σ)` (from the
    generic `siteOp_mulVec_apply` and `bondOp_diag_apply` of `Ring/Hamiltonian.lean`), and
    `Re ⟪x, Hx⟫ ≥ ∑ (E(σ) − gL) |x(σ)|²`;
  - `Bridge`: the checker's bit-indexed sums are exactly these quadratic forms (`xor 2ʲ` is a spin
    flip);
  - `Seed`: one eigenvalue lies within `10⁻¹⁰` of each candidate and at most two lie below the
    anchor, hence two-sided bounds on `Z₆(β)`, `Z₁₂(β)`, `Z₁₂(2β)` and the theorems `seed_spatial`,
    `seed_thermal`, `first_point_seed`, and `ising_gap`.
* **General seeds from low spectra (`SeedKernel/Spectral`, `Comparison`, `LowSpectrum`, `ExpBounds`,
  `Seeds`).** Model-independent: the residual bound and counting by comparison for any Hermitian
  matrix; the weighted comparison `le_re_quadratic_of_weights` (positive weights `w` off a set `B`
  with `∑_{τ ∉ B, τ ≠ σ} |A σ τ| w(τ) ≤ w(σ) (Re A σ σ − ρ)` give `Re ⟪x, Ax⟫ ≥ ρ ‖x‖²` on vectors
  vanishing on `B`); the low spectrum `LowSpectrum` and its Boltzmann-sum bounds; rational upper
  bounds `expUpQ` on `e^y`; and, for any `RPRing` and any boundary, the spatial, thermal,
  open-thermal and boundary-trace seeds from low spectra under decidable rational conditions
  (`SpatialCheck`, `ThermalCheck`, `TraceCheck`).
* **The Ising ring at any `β` (`SeedKernel/Ising/Checks`).** `ring_gap_of_check`: the gap from any
  chain certificate whose seeds pass `RingCheck` at `β₀`.
* **The free Ising chain (`SeedKernel/IsingOpen/*`).** The same pipeline on the open chains of 2, 8
  and 12 sites with free ends:
  - `Hamiltonian`: the free chain acts as `E(σ) x(σ) − g ∑ⱼ x(flipⱼ σ)` with the open-chain classical
    energy (`freeH_apply` gives its entries), and `freeH_weights` reduces the weighted comparison to
    `g ∑ⱼ w(flipⱼ σ) ≤ w(σ) (E(σ) − ρ)`. The ring bound `E(σ) − gL` fails here: flipping an end spin
    breaks one bond, not two;
  - `Bridge`: the checker sums, the classical count and the integer weight checks (`weightOK`) are
    these quadratic forms and comparison hypotheses;
  - `Seed`: the low spectra `lowSpectrum2`, `near8`, `lowSpectrum12`; `free_gap_of_check` at any `β`
    and certificate passing `FreeCheck`; and at `β₀ = 97/4` `ising_free_gap` (`n ≥ 24`) and
    `ising_free_gap_long` (`n ≥ 96`).

| File | Role |
|---|---|
| `SeedKernel.lean` | Library root importing all `SeedKernel` modules |
| `SeedKernel/Packed.lean` | Shared packed-`Nat` lane operations: `lane`, `repunit`, `pack`, `move`, `Signed`, floors, `reverseLanes`, `dot`, `gram`, and fast variants |
| `SeedKernel/PackingProofs.lean` | Radix encoding proofs: lane extraction and `pack = radixPack` |
| `SeedKernel/Obligations.lean`, `SeedKernel/PackedTests.lean` | Target correctness propositions for packed operations and kernel-checked examples |
| `SeedKernel/Spin1/Packed.lean`, `Spin1/Optimized.lean` | Packed spin-1 `hApply`, `filterApply`, trial moments, `ritzCheck`, and prepared-mask variants |
| `SeedKernel/Spin1/Obligations.lean`, `Spin1/Tests.lean` | Row-wise specification, open semantic propositions (`HApplyCorrect`, `FilterCorrect`, ...), and small-ring checks |
| `SeedKernel/Spin1/Counting.lean` | Conditional one-Ritz-value trace-counting lemma and formal counterexamples to naive two-value variants |
| `SeedKernel/Spin1/H6.lean`, `H8.lean`, `H10.lean` | Generated Ritz/residual witnesses (`code/models/lean_seeds/certify.jl spin1_kernel`) |
| `SeedKernel/Heisenberg/Trial*.lean` | The D = 26 integer trial state checked in the kernel (`trial_check_72`, `trial_check_120`): SU(2)-and-swap block reduction by an integer intertwiner, packed integer matrix powers with a correctness proof; data from `code/models/lean_seeds/scripts/trial/generate.jl` |
| `SeedKernel/Heisenberg/Seed36.lean`, `Heisenberg/Seed60.lean` | The spin-1 seeds and gaps with the trace tables as the only hypotheses (`gap_of_table`, `gap_sixty_chain_of_table`) |
| `SeedKernel/Heisenberg/Reseed120.lean` | The ring inputs of the boundary-field gap at `(120, 49)` from the trace table (`seed_ring_spatial_of_table`, `seed_ring_thermal_of_table`) |
| `SeedKernel/Heisenberg/OpenTrial*.lean` | The four open-chain trial vectors of the boundary paper on 34, 78, 142 and 240 sites (`oppositeField_trial_*`); data from `code/models/lean_seeds/scripts/open_trial/generate.jl` |
| `SeedKernel/Heisenberg/BoundaryTraces.lean` | The nine boundary integers `k_{n,g}` with `n ≤ 4` by `decide +kernel` (about 45 s, 6 GB) |
| `SeedKernel/Heisenberg/ThermalTraces.lean` | The two periodic totals `W_{4,g}` by `decide +kernel` on the packed evaluator; `HornerTotalsLarge` (the sixteen with `n ≥ 5`) and `hornerTotals_of_large` |
| `SeedKernel/Heisenberg/Integers.lean` | Both papers' gap theorems from the printed integers (`gap_of_hornerTotals`, `gap_sixty_of_hornerTotals60`, `spin_one_boundary_field_gap_of_integers`) |
| `SeedKernel/Ising/Checker.lean`, `Ising/Packed.lean`, `Ising/Chunked.lean` | Integer $20H$ action and direct, packed (on the shared lanes), and 64-chunk certificate checkers |
| `SeedKernel/Ising/Generated.lean`, `Ising/GeneratedChunks.lean` | Generated eigenvector proposals and kernel checks (`code/models/lean_seeds/certify.jl ising`) |
| `SeedKernel/Ising/ChunkCoverage.lean`, `Ising/NegativeChecks.lean` | Proof that the 64 chunks cover `{0, ..., 4095}` and rejection tests for corrupted witnesses |
| `SeedKernel/Spectral.lean`, `SeedKernel/Comparison.lean` | Adapters to the shared spectral residual and counting bounds, and weighted comparison for any Hermitian matrix |
| `SeedKernel/LowSpectrum.lean`, `SeedKernel/ExpBounds.lean` | Low spectra, Boltzmann-sum and purity-ratio bounds, and rational bounds on the exponential |
| `SeedKernel/Seeds.lean` | The seeds of any `RPRing` (ring and open chain) from low spectra and decidable rational checks |
| `SeedKernel/Ising/Hamiltonian.lean`, `Ising/Bridge.lean` | The Ising Hamiltonian acting on vectors and its identification with the integer checker (spectral input from `CertifiedGaps/Spectral`) |
| `SeedKernel/Ising/Seed.lean` | The unconditional seed and `ising_gap` |
| `SeedKernel/Ising/Checks.lean` | The ring seeds and gap at any `β` (`RingCheck`, `ring_gap_of_check`) |
| `SeedKernel/IsingOpen/Checker.lean`, `IsingOpen/Generated.lean`, `IsingOpen/GeneratedChunks.lean` | Integer checks on the free chains and their generated witnesses (`code/models/lean_seeds/certify.jl ising_open`) |
| `SeedKernel/IsingOpen/Hamiltonian.lean`, `IsingOpen/Bridge.lean` | The free chain acting on vectors, its entries and weight condition, and the identification with the checker |
| `SeedKernel/IsingOpen/Seed.lean` | The free-chain low spectra, `free_gap_of_check`, the certificates, `ising_free_gap` and `ising_free_gap_long` |
| `SeedKernel/Audit.lean` | `#print axioms` audit across all main declarations |
| `SeedKernel/Eval/BoundaryTotals.lean`, `Eval/ThermalTotals.lean`, `Eval/ThermalTotals60.lean` | Compiled evaluation (`#guard`, not a proof; opt-in `lake build SeedKernelEval`) of the fifteen boundary integers, of the `b = 49/4` totals with `n ≤ 9`, and of the six `b = 21/2` totals with `n = 6, 8` |
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
   `julia code/lean/scripts/chain_certificate.jl εₛ εₜ β₀ [J] [digits]` (the arguments are exact literals
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
6. **Check the seeds against the certificate:** Add the model to `code/lean/scripts/check_seeds.jl` and
   run `julia --project=code/julia code/lean/scripts/check_seeds.jl` from the repository root. It checks that
   every Lean `s[0], t[0]` is at least the certified defect at the same `(n₀, β₀)`.
7. **Discharge the seeds from a trace table (optional):** The `seeds` lemma can be proved instead
   of assumed, from a table of short-ring partition functions and a trial state. The generic
   layer is in `Seed/`, `Transfer/SingleSite*.lean` and `Transfer/MatrixProduct.lean`; a model
   supplies only data and model-specific identities:
   - **The scheme** (`Seed/Scheme.lean`): which twisted traces are tabulated and how they split
     into symmetry sectors, $`Z_g = \sum_s c_{gs} M_s`$ with a rational left inverse. The
     schemes `SectorScheme.plain` (untwisted traces, one sector), `SectorScheme.d2` (`D₂`
     sectors, $`Z = J + 3N`$, $`Z(P) = J - N`$) and `SectorScheme.ion` (twists $`1, P, C`$) are
     provided.
   - **The sector representation** (`SectorApprox`): free for an `RPRing` with `refl = 1`, from
     `RPRing.sectorApprox_plain` or, given the symmetry data `SingleSite.D2Symmetry` (two sign
     symmetries and a label permutation of the bond factors, `Transfer/SingleSiteD2.lean`), from
     `RPRing.sectorApprox_d2`.
   - **The certificate** (`TableCertificate`, `Seed/Certificate.lean`): the table, one filter
     per capped sector, the dominant-sector data `SectorSeed` and `(β, κ)`; `Valid` is checked by
     `decide +kernel`. `RPRing.seeds_of_certificate` (`Seed/Ring.lean`) then gives the
     conclusion of `seeds`. The data come from `julia code/lean/scripts/seed_certificate.jl
     config.toml` (scheme, table, trace error, `(q, a, ℓ)`, `(β, κ)`; see the header of the
     script), which searches the filters, checks every condition in exact rationals and prints
     the Lean certificate; changing a parameter point means rerunning it.
   - **The remaining hypotheses** are numerical: the table bounds (`TraceBounds.HoldsStrict` for
     the shifted traces $`e^{-\beta\kappa n} Z_g(n, \beta)`$; a missing lower or upper bound is
     `none`, and only the side that a filter coefficient uses is needed) and the trial-state
     bounds at $\beta$ and $2\beta$. The lemmas for certifying tables by eigenvalue counting
     are in `Spectral/`; for a trial state, `MPS.exp_lt_gibbsTrace_int`
     (`Transfer/MatrixProduct.lean`) turns one integer inequality between two transfer-matrix
     traces of an integer MPS into $`e^{\beta\kappa N} < Z_N(\beta)`$.
   The spin-1 seed `Heisenberg/Seed36.lean` is the worked example.

A certificate should therefore report its seeds as exact rationals together with `(n₀, β₀)` and the
exact `C_α, o`. Building blocks shared by several models (one-site matrices, the cell swap) go in
`Local/`, not in `Models/`.

**Naming scheme.** Every model uses the same names, so that each model file reads the same way:

| Name | Meaning |
|---|---|
| `Models/<Name>.lean`, namespace `CertifiedGaps.<Name>` | `<Name>` is the Julia directory `code/models/<slug>` in UpperCamelCase (table below) |
| `model` | the model as an `RPRing`, a function of its real parameters |
| `Point`, `points` | the certified parameter points with `β₀` and `cert` |
| `points_valid`, `points_β₀_pos`, `points_<field>_…` | decidable facts about all points (`decide +kernel`) |
| `seeds` | the seed inequalities at every point: the model's single `sorry` |
| `gap` | the gap at every point, in the frame of `model` |
| `physicalHamiltonian`, `physicalSpectrum`, `physicalHamiltonian_eq`, `physicalSpectrum_eq`, `physical_gap` | `Models/<Name>Physical.lean`: the physical frame, when `model` is rotated |

| Lean namespace | Julia directory |
|---|---|
| `SpinOneHeisenberg` (`Heisenberg/`) | `code/models/spin1` (the headline theorem keeps its own names: `seed_spatial`, `seed_thermal`, `spin_one_heisenberg_gap`) |
| `SpinOneLargeD` | `code/models/spin1_largeD` |
| `SpinOneDimerized` | `code/models/spin1_dimerized` |
| `HeisenbergLadder` | `code/models/heisenberg_ladder` |
| `StaggeredHeisenberg` | `code/models/heisenberg_staggered_field` |
| `DimerizedHeisenberg` | `code/models/dimerized_heisenberg` |
| `TwoLegIsing` | `code/models/two_leg_tfim` |
| `IsingLongitudinal` | `code/models/tfim_longitudinal` |

### Replacing a seed by a kernel proof

The seeds of every gap theorem above are inequalities between partition functions of a few small
rings or chains. `SeedKernel` proves them in three model-independent steps, so a new model or
parameter point only supplies the model-specific first step:

1. **Low spectra.** For each length the seeds use, prove a `LowSpectrum` (one eigenvalue within
   `r` of a Ritz value `θ₀`, a second at least `θ₁ - r`, all others at least `ρ`) of
   `R.spectrum L` or `R.openSpectrum b n`. The tools are model-independent: `exists_near_of_residual`
   for a trial vector, counting by comparison (`card_eigenvalues_lt_le` with
   `le_re_quadratic_of_weights`), and `lowSpectrum_of_near`. What is model-specific is an integer
   checker for the Hamiltonian in a product basis, its identification with `R.hamiltonian L` or
   `R.openHamiltonian b n` (as in `Ising/Bridge.lean` and `IsingOpen/Bridge.lean`), and the
   generated witnesses (`code/models/lean_seeds/certify.jl`). The low spectra do not depend on `β`.
2. **Seeds at any `β₀`.** `spatialPurity_ge_of_check`, `thermalPurity_ge_of_check`,
   `openThermalPurity_ge_of_check` and `openTrace_le_of_check` (`SeedKernel/Seeds.lean`) turn the
   low spectra into the seeds when a rational inequality (`SpatialCheck`, `ThermalCheck`,
   `TraceCheck`) holds at the chosen rational `β₀`; `decide +kernel` checks it.
3. **The gap.** Feed the seeds to `RPRing.gapGT_of_chainCertificate` (rings) or
   `RPRing.openGapGT_of_traces` and `RPRing.freeGapGT_of_traces` (open chains).

With the low spectra in place, a new inverse temperature or a new certificate is one theorem and
one `decide +kernel`: `SeedKernel.Ising.ring_gap_of_check` and
`SeedKernel.IsingOpen.free_gap_of_check` take `β₀` and the certificate as arguments, and
`ising_free_gap` and `ising_free_gap_long` are two such instances.

## Files

| File | Role |
|---|---|
| `PurityCriterion.lean` | Abstract two-purity criterion: squaring inequality, thermal-spatial coupling, one exact doubling step for arbitrary defects, $`\beta`$-interpolation, the purity-to-gap bound, the criterion with a coefficient `C ≥ G(u₀)` (`gapGT_of_coefficient`) and its equal-defect case (seed defects $\le 1/100$, `C = 5`) |
| `Chain.lean` | Exact asymmetric doubling recurrence for unequal seed defects and `gapGT_of_chain` with bound `log r / β₀` |
| `ChainCheck.lean` | Decidable rational `ChainCertificate` structure discharging the hypotheses of `gapGT_of_chain`, and `StepsValid` for the purities along a chain alone |
| `Transfer/PiTensor.lean` | Sitewise tensor products (`piTensor`), expansion of $`\mathrm{Tr}\, T^N`$ as a sum over closed label paths and of matrix elements of ordered products as sums over open paths |
| `Transfer/MomentForm.lean` | Proof that traces of even powers of a real matrix factorization $B^T B \succeq 0$ satisfy the spatial moment form |
| `Transfer/ProductLimit.lean` | Lie–Trotter product formula: $\lim_{M\to\infty} (\prod_k (1 + (c/M) h_k))^M = e^{c \sum_k h_k}$ |
| `Transfer/PolynomialFilter.lean` | Polynomial filters: eigenvalue caps from finitely many moments (`CapCheck`, `momentBound`, one-sided `pairingBound`), checked by `decide` |
| `Transfer/MatrixLemmas.lean` | Traces of powers, block matrices `blockOf` by a key, and reindexing invariance |
| `Transfer/SingleSite.lean`, `Transfer/SingleSiteLimit.lean` | Single-site Trotter transfer matrix of any ring with $`-h = \sum_\alpha C_\alpha \otimes C_\alpha`$ and on-site `o`: exact trace identity at every length (odd included), Trotter limit, and `RPRing.exists_transfer_approx` |
| `Transfer/SingleSiteSectors.lean`, `Transfer/SingleSiteD2.lean` | Sign and label symmetries of the transfer matrix, `D2Symmetry`, and the sector representation `D2Approx` (`RPRing.d2Approx`) |
| `Transfer/SpectralForm.lean` | Partition functions and both purities as limits of power sums of transfer eigenvalues, at every length (`RPRing.tendsto_sum_transferEigenvalues_pow`, `RPRing.tendsto_spatialPurity`, `RPRing.tendsto_thermalPurity`) |
| `Transfer/CyclicSectors.lean`, `Transfer/SingleSiteA4.lean` | ℤ₃ sectors of a symmetric matrix by projectors; `A4Symmetry` and the `I/O/N` sector representation `RPRing.sectorApprox_ion` |
| `Local/SpinOneRotation.lean` | Rotated spin-1 exchanges `∑ R_{αβ} S^α ⊗ S^β` and their conjugation by rotations |
| `Seed/Lift.lean`, `Seed/LiftedSeed.lean` | Temperature lifting of Boltzmann sums; the seed with a temperature lift and `I/O/N` multiplicities (OpenAI's Proposition 5.4 in generic form) |
| `Heisenberg/SectorTransferC.lean`, `Heisenberg/Seed60.lean` | The `2π/3`-twisted spin-1 ring and `exists_ion_approx`; the seed at `(60, 105/4)` from its inputs |
| `Transfer/MatrixProduct.lean` | Integer matrix-product states for any ring Hamiltonian: norm and energy as transfer-matrix traces, and the trial-state bound `MPS.exp_lt_gibbsTrace` |
| `Ring/Variational.lean` | Rayleigh-quotient bound on the Gibbs trace `gibbsTrace H β` for Hermitian `H` |
| `Spectral/*.lean` | Eigenvalue counting (`eigCount`; counting by a subspace, a kernel, or vanishing coordinates), the residual bound `exists_eigenvalue_near`, rational `LDLᵀ` inertia certificates, the Schur compression bound, Weyl monotonicity, Ritz bounds and binned trace bounds |
| `Spectral/OperatorBounds.lean` | Operator inequalities `a ≤ A ≤ b` through quadratic forms: sums, unitary conjugates, the weighted Schur test `abs_qform_le_of_schur`, bounds from an annihilating polynomial, and the Gram factorization of a nonnegative Hermitian matrix |
| `Spectral/FunctionalCalculus.lean` | Real functions `hA.cfc f = U diag(f(λᵢ)) Uᴴ` of a Hermitian matrix (Mathlib's `Matrix.IsHermitian.cfc`) without continuity side conditions: algebra, `exp`, and the norm bound `‖f(A) v‖ ≤ (max_i ‖f(λᵢ)‖) ‖v‖` |
| `Seed/Scheme.lean` | Sector schemes (`plain`, `d2`, `ion`), the sector representation `SectorApprox`, trace tables `TraceBounds` and the sector moment bounds |
| `Seed/CappedMass.lean`, `Seed/Dominant.lean` | Capped-mass bounds and the dominant-sector seed `SectorSeed` |
| `Seed/Limits.lean` | The limit `δ → 0` of a seed, shifted traces, and shift invariance of the purities |
| `Seed/Certificate.lean`, `Seed/Ring.lean` | `TableCertificate` and its purities; `RPRing.seeds_of_certificate` and the plain and `D₂` sector representations of rings |
| `Heisenberg/SzBridge.lean` | Spin-1 twisted ring and `S^z`-basis definitions, their identification with the generic `twistedRingHamiltonian`, and the Cartesian ↔ `S^z` bridge for partition functions |
| `Heisenberg/SectorTransfer.lean`, `Heisenberg/TrialBound.lean` | Spin-1 instances: the `D₂` symmetry of the Heisenberg ring and the trial-state bound |
| `Heisenberg/TraceTable.lean`, `Heisenberg/Seed36.lean` | The spin-1 trace table and the seed at `(36, 49/4)` (`seed36_of_table_mps`: from the table and one integer MPS inequality) |
| `Heisenberg/EnergyBound.lean` | The Lemma "Energy bounds" of the spin-1 paper: `-(3L/2) ≤ H_L(R) ≤ L` for every rotation `R` of the closing bond, and the untwisted, `P`- and `C`-twisted cases |
| `Ring/Hamiltonian.lean` | The one Hamiltonian API: ring geometry, `bondOp`/`siteOp` (`bondOp j k` a `*`-algebra embedding for `j ≠ k`; entries `siteOp_apply`, `siteOp_mulVec_apply`, `bondOp_diag_apply`), `twistedRingHamiltonian B Bc L` (closing bond `Bc`), `ringHamiltonian B L`, `siteSum A L`, and the Gibbs trace `gibbsTrace H β` $= \mathrm{Re}\,\mathrm{Tr}\, e^{-\beta H}$ (eigenvalue sum, positivity, unitary invariance) |
| `Ring/Model.lean` | `RPRing` structure, `H_L = ringHamiltonian h L + siteSum o L`, Hermiticity of $H_L$, and $Z_L = \sum_m e^{-\beta E_m}$ |
| `Ring/Transfer.lean` | Even-bond, on-site, and odd-bond Trotter slices and the reflection-positivity identity $\operatorname{Tr} S^M = \operatorname{Tr} (B^T B)^{L/2}$ with real $B$ |
| `Ring/Gap.lean` | Proof of `RPRing.spatialMomentBounds` and `RPRing.gapGT_of_chainCertificate` |
| `Ring/Conjugation.lean` | Sitewise unitary conjugation of ring Hamiltonians (twisted or not, `conj_twistedRingHamiltonian`), invariance of Gibbs traces and spectra |
| `Ring/Symmetry.lean` | Signed-permutation symmetries stated entrywise over any commutative semiring (`Intertwines`, `signedPerm`): site permutations and on-site relabelings of `bondOp`, sign gauges, the rotation and reflection of a ring with a seam (`seamHamiltonian`), and the signed rotation of a twisted ring (`intertwines_twisted_rotate`) |
| `Ring/EnergyBound.lean` | Uniform energy bounds `a L ≤ H_L ≤ b L` of twisted and untwisted rings from bounds on one bond or on pairs of consecutive bonds, and row sums of `bondOp` for local Schur tests |
| `Trace/Poisson.lean` | Conditioned Poisson weights, the Fox–Glynn truncation bound for `e^{m(R - 1)}` of a Hermitian contraction (`nsq_exp_sub_poisson_mulVec_le`, instances `(98, 280)` and `(84, 256)`), and the Taylor enclosure of `e^x` |
| `Trace/Horner.lean` | Integer Horner columns `X ← ⌊D X / s⌋ + cⱼ e_r` (`hornerFloor`) and their rounding error against `(∑ aⱼ Rʲ) e_r` for a contraction `R = D / s` |
| `Trace/Frobenius.lean` | Gibbs traces from representative columns: `Z = e^{bc} ∑_r ‖F e_r‖²` with `F = e^{-(b/2)(H + c)}`, symmetric columns, weighted Minkowski, and the rational enclosure `gibbsTrace_enclosure` |
| `Trace/Orbit.lean` | Column symmetries (`IsColumnSymmetry`) as a submonoid, and the representative hypothesis of `gibbsTrace_enclosure_of_rep` from generating symmetries (`nsq_gibbsHalf_rep`) |
| `Heisenberg/ThermalTable.lean` | The `b = 49/4` table from the paper's integer algorithm: `D = 32R`, coefficients `⌊2⁵⁵ pⱼ⌋`, representatives and `hornerTotal`; `traceTable_of_hornerTotals : HornerTotals → TraceTable` (the eighteen rational endpoints checked by `decide +kernel`) |
| `Trace/Eisenstein.lean`, `Trace/HornerEisenstein.lean`, `Trace/OrbitConj.lean`, `Ring/SymmetryPhase.lean` | Eisenstein integers `ℤ[ω]` with floored division and the Horner error over them; antiunitary column symmetries; symmetries up to a phase, for the `2π/3`-twisted ring |
| `Heisenberg/ThermalTable60.lean` | The `b = 21/2` table from the paper's integer algorithm (`g = 1, P, C`, `n = 6, 8, 10`): `traceTable60_of_hornerTotals60 : HornerTotals60 → TraceTable60`, and the seed at `(60, 105/4)` |
| `Trace/WeightedHorner.lean`, `Trace/WeightedFrobenius.lean` | The Horner rounding error and the Gibbs-trace enclosure in a weighted norm, for integer matrices `D = s J A J⁻¹` that are contractions only after a diagonal similarity |
| `Trace/PackedLanes.lean`, `Trace/PackedHorner.lean` | A packed big-integer evaluator of `hornerFloor` for banded integer matrices with bands of either sign (one `Nat` per row, one lane per computed column, any list of columns, division by a constant via multiply and shift) proved equal to `hornerFloor` |
| `Trace/PackedRealify.lean` | Four banded matrices as the blocks of one (`realBand`, `realBand_entry`), for the real form of a matrix over `ℤ[ω]`, and the packed numerator of the norms `u² - uv + v²` (`runCols_numeratorE`) |
| `Heisenberg/WordBand.lean` | Spin-one words as base-3 numbers, the spin-flip hops of an exchange bond as bands (`hopBand_term`), shared by the boundary and periodic recurrences |
| `Heisenberg/ThermalPacked.lean` | The `b = 49/4` recurrence as a packed run on the representative columns: `D = (32 - 3n) I - 2H` as a banded matrix, representatives on word numbers, and `hornerTotal_eq_of_packed` |
| `Heisenberg/ThermalPacked60.lean` | The `b = 21/2` recurrence as the same packed run (`2B` with divisor `32`: the band of `ThermalPacked.lean` for the twists `1`, `P`, and the real form `ringBandC` for `C`): `hornerTotal60_eq_of_packedTotal60`, `hornerTotals60_of_packed` |
| `Seed/Reseed.lean`, `Heisenberg/Reseed120.lean` | Reseeding from one ring seed to a longer one at four times the inverse temperature; the ring inputs `P_s(120, 49) ≥ 1 − 0.047`, `P_τ(240, 49) ≥ 1 − 0.0179` from the trace table |
| `Transfer/RatioMonotone.lean`, `Open/Initialization.lean` | Monotonicity of `x^m / P(x)²` on an interval from a decidable `RatioCheck`; the scalar initialization of the boundary seeds (`move_one` steps, defect bounds from moment witnesses) |
| `Heisenberg/BoundaryInit.lean`, `Heisenberg/BoundaryInitTable.lean` | The boundary paper's initialization at `b₀ = 49/4`: `seed_boundary` and `seed_open_thermal` from `ThermalCore`, the trial bounds `TrialData` and the trace table |
| `Transfer/OpenMatrixProduct.lean`, `Heisenberg/OpenTrialBound.lean` | Open matrix-product trial vectors `bᵀ A ⋯ A c` for open chains with end terms, and the variational eigenvalue bound for the opposite-field chain |
| `Open/TwoBoundary.lean`, `Open/LeadingLimit.lean` | A Hermitian transfer matrix with an isolated leading line (`IsLeadingLine`): the boundary moment witness from even and odd pairings, the leading-weight polynomial bound `c f(e₀) ≤ ∑ aₙ ⟨v, Xⁿ P v⟩` for an invariant projection `P`, and their passage to the limit (`LeadingWitness.of_tendsto`) |
| `Transfer/OpenSingleSite.lean`, `Transfer/SingleSiteLeading.lean` | Open chains with two arbitrary end terms in the single-site Trotter matrix (`⟨v_L, Xⁿ v_R⟩` tends to the Gibbs trace); sign and relabeling actions, the `A₄` average `Π`, and the leading line from the `D₂` block eigenvalues |
| `Heisenberg/BoundaryLeading.lean` | The spectral witness at `b₀ = 49/4` (`λ₀ ≥ 0.999`, `ρ₀ ≤ 0.874`, `c₀ P_*(λ₀)² ≤ 0.003745`) and `thermalCore_of_tables` from the trace table and the boundary trace table |
| `Heisenberg/BoundaryData.lean`, `Heisenberg/BoundaryTraces.lean`, `Heisenberg/BoundaryThermal.lean` | The fifteen boundary traces `Z_{n,g}` (`n = 2, …, 6`, right field rotated by `g`): the paper's integer recurrence, the table `BoundaryTraceTable` from `BoundaryTotals`, `Y₆(b₀) < 0.027` and the leading-weight polynomial `< 0.003745` |
| `Heisenberg/Hamiltonian.lean` | The spin-1 Heisenberg ring $H_L = \sum_j \mathbf{S}_j \cdot \mathbf{S}_{j+1}$ (`ringHamiltonian SpinOne.exchange L`) in the Cartesian basis |
| `Heisenberg/Gap.lean`, `Heisenberg/Main.lean` | Both bounds of OpenAI's Theorem 1.1: `spin_one_heisenberg_gap` (`L ≥ 2304`, from the seeds at 36 sites and `chain36`) and `spin_one_heisenberg_gap_sixty` (`L ≥ 60`) |
| `Heisenberg/OpenChain.lean` | Theorem 2 of the main note (`spin_one_boundary_field_gap`) and its sharper form (`spin_one_boundary_field_gap_sharp`): the equal- and opposite-field open chains as open chains of the Heisenberg ring, the four inputs, and the opposite-field gap at even lengths (`spin_one_opposite_field_gap`) |
| `Open/MomentBounds.lean` | `OpenMomentBounds`: closed inequalities of the boundary moment form, from a finite Hermitian transfer matrix and through limits; `leadingRoot`, `leadingOverlap`, `boundaryDefect` |
| `Open/Criterion.lean` | The boundary-field criterion `gapGT_of_openChain` for general seeds: the ring at the open-chain scales, the three moves of Appendix A with the even-length purity, the boundary seed from finite traces (`boundaryDefect_le_of_traces`), the base steps (`OpenStep`) and the uniform phase |
| `Open/Certificate.lean` | Rational certificates `OpenCertificate` checked by `decide +kernel`, `gapGT_of_openCertificate`, the note's `boundaryFieldCertificate` with `gapGT_of_openSeeds`, and `sharpBoundaryFieldCertificate` |
| `Open/Model.lean` | Open chains of an `RPRing` with boundary terms (`Boundary`, `openHamiltonian`, `openSpectrum`), and free ends (`freeBoundary`, `openHamiltonian_free`) |
| `Open/Transfer.lean` | Open-chain Trotter slices, time reversal `J`, and `Tr (slice)^M = ⟨v, Xⁿ v⟩` with `X = F J` |
| `Open/Bounds.lean`, `Open/Gap.lean` | `RPRing.openMomentBounds` and the end-to-end `RPRing.openGapGT_of_certificate`, `RPRing.openGapGT_of_traces`, `RPRing.freeGapGT_of_traces` and `RPRing.openGapGT_of_seeds` |
| `Local/Pauli.lean`, `Local/SpinHalf.lean` | Real Pauli matrices $X, Z$ and real spin-½ factors $S^x, iS^y, S^z$ (`SpinHalf.P`) |
| `Local/SpinHalfRotation.lean` | Sublattice rotation $e^{i\pi S^y}$ and complex spin-½ operators |
| `Local/SpinOne.lean` | Spin-1 building blocks: Levi-Civita symbol, Cartesian `SpinOne.S α`, real antisymmetric generators `SpinOne.G α` $= i S^\alpha$, the exchange `SpinOne.exchange`, the $\pi$-rotation signs `rotX`, `rotY`, `rotZ`, and the $S^z$ basis `szBasis`, `spinZ` |
| `Local/CellSwap.lean` | Two-site cell swap involution $F$ satisfying $F (A \otimes B) F = B \otimes A$ |
| `Models/SpinOneLargeD.lean` | Spin-1 chain with single-ion anisotropy $\sum \mathbf{S}_j \cdot \mathbf{S}_{j+1} + D \sum (S_j^z)^2$ as an `RPRing` ($D = 0$ is the Heisenberg ring), at $D = 4$ and $D = 1603/500$ |
| `Models/SpinOneDimerized.lean` | Dimerized spin-1 chain with internal cell reflection $F$ at $\delta = 9/10$ |
| `Models/TwoLegIsing.lean` | Two-leg transverse-field Ising ladder at $g = 3$ |
| `Models/IsingLongitudinal.lean` | Transverse-field Ising chain with longitudinal field across 30 parameter points |
| `Models/HeisenbergLadder.lean`, `Models/StaggeredHeisenberg.lean`, `Models/DimerizedHeisenberg.lean` | Reflection-positive rotated Hamiltonians, seeds, and gap theorems for the three spin-½ models |
| `Models/*Physical.lean` | Physical spin-½ Hamiltonians, unitary sublattice rotations, and transferred gap theorems |
| `scripts/chain_certificate.jl` | Julia generator that constructs an exact rational `ChainCertificate` from $(\varepsilon_{\mathrm{s}}, \varepsilon_\tau, \beta_0)$ (exact literal parser; Lean re-checks) |
| `scripts/open_certificate.jl` | Julia generator that searches an exact rational `OpenCertificate` from the four open-chain seeds (reuses the parser of `chain_certificate.jl`; Lean re-checks) |
| `scripts/seed_certificate.jl` | Julia generator that searches the filters and prints an exact `TableCertificate` from a trace table, its sector scheme and the seed parameters (Lean re-checks) |
| `scripts/check_seeds.jl` | Checks every Lean seed defect against the certificate JSON in `../models/*/results` |

## Why each piece is needed in Lean

On paper, the main note obtains the spatial moment representation in a single conceptual step: it
constructs an exact continuous-time thermal MPO on Poisson bond histories whose spatial transfer operator
$X_\beta$ is a compact self-adjoint integral operator on $L^2(\mu_\beta)$, giving
$Z_n(\beta) = \operatorname{Tr} X_\beta^n = \sum_i |\lambda_i|^n$. Formalizing that infinite-dimensional
path-space construction directly in Lean still needs Poisson history measures and the bridge from
integral kernels to Hilbert–Schmidt operators and their cycle-trace identities. Mathlib already has
$L^2$ integration tools and the compact self-adjoint spectral theorem, but not this ready-to-use
path-space and trace infrastructure. Instead, the Lean formalization uses a finite-dimensional
Trotterization route that needs only finite matrix algebra and the matrix exponential limit. Each module below addresses a specific distinction between the paper proof and the formal proof:

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
