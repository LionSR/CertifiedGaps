# CertifiedGaps: Lean proofs of spectral gaps by the spacetime bootstrap

This Lean 4 + Mathlib package proves that quantum spin rings have a spectral gap that stays
open as the ring grows. It follows the spacetime-bootstrap argument of the main note
(`../notes/thermal-mpo-gap-note.tex`). **Everything is proved except the finite numerical
inputs** (two "seed" inequalities per model). They are named `sorry` lemmas, and their values
come from the certified computations in `../models/`.

## The rundown

### The headline theorem

```lean
-- CertifiedGaps/Heisenberg/Main.lean
theorem SpinOneHeisenberg.spin_one_heisenberg_gap {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784)
```

For every even `L ≥ 2304`, the periodic spin-1 Heisenberg ring `H_L = ∑ⱼ Sⱼ·Sⱼ₊₁` has a unique
ground state. Every other level lies more than `log 20 / 784 ≈ 0.0038` above it. The statement
uses the actual Hamiltonian matrix on `(Fin L → Fin 3)`. `GapGT e c` means: some level `m₀` has
every other level more than `c` above it. The proof is complete except for two `sorry`s,
`seed_spatial : P_s(2304, 784) ≥ 0.99` and `seed_thermal : P_τ(4608, 784) ≥ 0.99`. These are
OpenAI's Proposition 5.5, reproduced in exact arithmetic in `../models/spin1/`.

### How the proof is put together

Read the files in this order. Each step uses only the ones before it.

1. **From purities to a gap** (`PurityCriterion.lean`). Write the partition function as
   `Z_L(β) = Re Tr e^{-βH_L}` and define two purities, `P_s(n,β) = Z_{2n}/Z_n²` and
   `P_τ(n,β) = Z_n(2β)/Z_n²`. Assume the *spatial moment inequalities*: traces behave like
   moments `∑ᵢ λᵢⁿ` of a positive transfer operator. Then both purities improve quadratically
   when length and inverse temperature double. Two purity bounds at one scale therefore give a
   gap for every larger even ring: `gapGT_of_purities` (seeds ≥ 0.99) and `gapGT_of_seeds`.
2. **Realistic seeds** (`Chain.lean`, `ChainCheck.lean`). Certified seeds are usually worse
   than 0.99. `gapGT_of_chain` follows the exact doubling recurrence for unequal defects. A
   `ChainCertificate` turns the whole chain into rational inequalities, which `decide` checks.
   The bound comes out as `log r / β₀`.
3. **The moment inequalities are true** (`Transfer/*`, `Ring/*`). This is the core of the
   formalisation and contains no numerics. It covers every reflection-positive ring
   (`RPRing`): bond `−h = ∑_α C_α ⊗ C_α` with real `C_α`, plus a real on-site term.
   * Trotterise `e^{-βH}`; each Trotter trace equals `Tr (FᵀF)^{L/2}` with real `F`
     (`Ring/Transfer`, `Transfer/PiTensor`).
   * A PSD matrix has the moment form (`Transfer/MomentForm`).
   * The Trotter traces converge to `Z_L(β)` by a product formula proved here
     (`Transfer/ProductLimit`).
   * The inequalities are closed, so they pass to the limit. The result is
     `RPRing.spatialMomentBounds`, and with step 2, `RPRing.gapGT_of_chainCertificate`.
4. **Models** (`Heisenberg/*`, `Models/*`). Each model is written as an `RPRing` in a real
   basis, plus its two seeds and a chain certificate, which together give its gap theorem.
   Some models are reflection positive only after a sublattice rotation. For those,
   `Ring/Conjugation` transfers the gap to the physical Hamiltonian (`physical_gap`).
5. **Seeds** (still `sorry`). These carry the physics, and the bootstrap cannot produce them.
   The opt-in library `SeedKernel` works toward proving them in the kernel. Julia generates
   witnesses, `decide +kernel` checks them, and no `native_decide` is used.

### What is proved, and what is open

| | status |
|---|---|
| the criterion (steps 1–2), the moment inequalities for every `RPRing` (step 3), the frame changes | proved; axioms `propext`, `Classical.choice`, `Quot.sound` only |
| gap theorems for the spin-1 Heisenberg ring and 8 other model families (table below) | proved **given their seeds** (`sorryAx` appears only through the seed lemmas) |
| the seeds | open. `SeedKernel` checks seed arithmetic in the kernel. Its furthest result is one longitudinal-Ising point, proved from four explicit `log Z` enclosure hypotheses. See `../models/lean_seeds/` for the plan and the open obligations |

### Build

```sh
cd lean
lake exe cache get          # prebuilt Mathlib; never build Mathlib from source
lake build                  # CertifiedGaps: about 4 min from scratch, warning-free but the seeds
lake build SeedKernel       # opt-in: kernel-checked seed arithmetic, about 80 s more
lake build SeedKernelBench  # opt-in: kernel benchmarks, about 2 min, up to 4 GB per file
```

Lean and Mathlib are pinned at `v4.35.0-rc3`. `CertifiedGaps` has about 4300 lines, `SeedKernel`
about 1200 (plus generated witness files). All declarations live in the namespace
`CertifiedGaps`; names below are relative to it.

### Where to read further

* [Status](#status): every model theorem with its bound and its open seeds; the sorry-free theorems.
* [SeedKernel](#seedkernel-kernel-checked-seed-arithmetic): what the kernel-checked seed work proves.
* [Adding a model](#adding-a-model): the five-step recipe.
* [Files](#files): one line per file.
* [Why each piece is needed in Lean](#why-each-piece-is-needed-in-lean): where the formal
  proof differs from the note, and why.

## Status

`GapGT e c` says that the lowest level of the spectrum `e` is nondegenerate and every other level
lies more than `c` above it. Every gap theorem below holds for every even ring length
`L ≥ n₀` (counted in cells).

### Model theorems (seeds open)

| Theorem | Model | Bound | Open seeds (the only `sorry`s) |
|---|---|---|---|
| `SpinOneHeisenberg.spin_one_heisenberg_gap` | spin-1 Heisenberg ring | `Δ_L > log 20 / 784` for `L ≥ 2304` | `seed_spatial`: `P_s(2304, 784) ≥ 0.99`; `seed_thermal`: `P_τ(4608, 784) ≥ 0.99` |
| `SpinOneLargeD.gap` | spin-1, `+ D ∑ (Sᶻ)²`, `D = 4, 1603/500` | `log 4.755 / (9/4) ≈ 0.693`, `log 1.0065429 / (381/125) ≈ 0.00214`; `L ≥ 6` | `seeds` at `(6, β₀)`: `P_s ≥ 1 − 0.06411`, `P_τ ≥ 1 − 0.04053` at `D = 4` |
| `SpinOneDimerized.gap` | spin-1 dimerized chain, `δ = 9/10` (cells = weak pairs) | `log 6.339 / (7/2) ≈ 0.528`, `L ≥ 2` | `seeds`: `P_s(2, 7/2) ≥ 1 − 0.04972`, `P_τ(4, 7/2) ≥ 1 − 0.03206` |
| `HeisenbergLadder.physical_gap` (rotated: `gap`) | spin-½ ladder `∑_{legs} S·S + 4 ∑_{rungs} S·S` | `log 9.033 / 2 ≈ 1.100`, `L ≥ 4` rungs | `seeds`: `P_s(4, 2) ≥ 1 − 0.03509`, `P_τ(8, 2) ≥ 1 − 0.02463` |
| `StaggeredHeisenberg.physical_gap` (rotated: `gap`) | `∑ S·S − h ∑ (−1)ʲ Sᶻⱼ`, `h = 1, 3/4, 1/2, 7/20, 1/4` | `1.147, 0.890, 0.594, 0.383, 0.214`, `L ≥ 8` | `seeds` at `(8, β₀)` |
| `DimerizedHeisenberg.physical_gap` (rotated: `gap`) | spin-½ dimerized chain, `δ = 1/2, 3/10, 1/4` (cells = strong pairs) | `0.532, 0.172, 0.057`, `L ≥ 4` | `seeds` at `(4, β₀)` |
| `TwoLegIsing.gap` | two-leg TFIM, `g = 3` | `log 2.585 / (3/2) ≈ 0.633`, `L ≥ 4` rungs | `seeds`: `P_s(4, 3/2) ≥ 1 − 0.1125`, `P_τ(8, 3/2) ≥ 1 − 0.05162` |
| `IsingLongitudinal.gap` | Ising chain `−ZZ − gX − hZ`, 30 points `g ∈ {1/4, …, 2}`, `h ∈ {1/20, …, 1}` | `log r / β₀` per point (`0.075` at `g = 1, h = 1/20` to `3.70`), `L ≥ 6` | `seeds` (one lemma for all 30 points) |

Every model except the spin-1 Heisenberg ring lists its parameter points in `points`; `gap hp`
(for `hp : p ∈ points`) gives `Δ_L > log r / β₀` with `r = p.cert.r`, and the one `sorry` of
the model, `seeds`, states both seed inequalities at every point. `#print axioms` on each of
these lists `propext`, `Classical.choice`, `Quot.sound` and `sorryAx` (from its seeds). The
`physical_gap` theorems are about the physical Hamiltonian; the rotated `gap` theorems about the
reflection-positive frame, which has the same spectrum.

The Lean theorems use the seeds that were current when they were written. Some model
certificates have since improved (for example `models/spin1_largeD` now certifies 0.899 at
`D = 4` from an `n₀ = 12` seed, and `models/dimerized_heisenberg` has slightly stronger seeds
than the ones used here, which therefore remain valid); the model READMEs and the root README give the current bounds.
Kernel-checked seed arithmetic is in the opt-in library `SeedKernel` (below), whose generators
and checker are in `models/lean_seeds/`.

### Proved theorems (no `sorry`)

These use only `propext`, `Classical.choice` and `Quot.sound`.

| Theorem | Statement |
|---|---|
| `gapGT_of_purities` | Proposition 1 of the note: moment form and both purities `≥ 0.99` at `(n₀, β₀)` give `Δ_L > log 20 / β₀` |
| `gapGT_of_seeds` | the same under the moment inequalities, with seed defect `ε₀ ≤ 1/100` and bound `log (1/(5ε₀)) / β₀` |
| `gapGT_of_chain`, `gapGT_of_chainCertificate` | unequal seed defects, via the exact recurrence; a `ChainCertificate` is checked by `decide` |
| `RPRing.spatialMomentBounds` | the moment inequalities for every reflection-positive ring (Trotter transfer matrices + product formula) |
| `RPRing.gapGT_of_chainCertificate` | the criterion for a reflection-positive ring: seeds + valid certificate give `Δ_L > log r / β₀` |
| `RPRing.eigenvalues_eq_spectrum` | a Hamiltonian conjugate to the ring by one-site unitaries has the ring's spectrum |
| `SpinOneHeisenberg.gapGT_log_twenty` | Theorem 1 of the note from its inputs (`ThermalMPOInputs`): `Δ_L > log 20 / 784` |
| `*.physicalSpectrum_eq` | each physical Hamiltonian has the spectrum of its rotated reflection-positive ring |
| `*.points_valid` | every chain certificate, by `decide +kernel` |

The seeds of the models come from the certificates in `../models/` (complete spectra of the seed
rings with rigorous eigenvalue enclosures; each model file names its directory). The spin-1
Heisenberg seeds are OpenAI's Proposition 5.5, reproduced in `../models/spin1/`.

## SeedKernel: kernel-checked seed arithmetic

An opt-in library (`lake build SeedKernel`; not built by `lake build`) toward proving seeds in
the kernel. Every computational statement is checked by `decide +kernel` on `Nat`/`Int`, which
the kernel evaluates with GMP; there is no `native_decide`, no `sorry` and no new axiom.
`SeedKernel/Audit.lean` (a root of the library) prints the axioms of 36 main declarations: all
use at most `propext`, `Classical.choice`, `Quot.sound`. The Julia generators, the witness
inventory and the checker (`check.jl`: Julia replay, `lake build SeedKernel`, audit) are in
`models/lean_seeds/`, whose README states the results and the open obligations.

**What it proves.** No seed is closed.

* *Shared packed lanes.* A vector of `d` naturals packed into one `Nat` in radix `2^w`, with
  lane extraction, shifts, butterfly reversal and a one-multiplication dot product; the radix
  encoding lemmas (`radixPack_extract`, `pack_eq_radixPack`) and kernel-checked examples. The
  general correctness of the packed operations is stated as `Prop`s (`DotCorrect`, ...) that
  nothing assumes.
* *Spin-1 Heisenberg (Std only).* Packed `H` and power filter `nI − H_M` on full `3ⁿ`-lane
  vectors; comparison with an independent row specification on small rings; arithmetic
  Ritz/residual witnesses `h6_ritz_witness`, `h8_ritz_witness`, `h10_ritz_witness` (integer
  facts about the filtered Néel trial, not eigenvalue theorems); a conditional one-Ritz-value
  counting step.
* *Longitudinal Ising, g = 1/4, h = 1/20.* Exact integer Gram/residual/anchor inequalities for
  two generated vectors on the 6- and 12-site rings (packed and direct, the 12-site direct check
  in 64 blocks with a coverage proof; rejection of corrupted proposals), and
  `first_point_seed_instance`: the first of the thirty seed conjunctions of
  `IsingLongitudinal.seeds`, proved from four explicit log Z enclosures (`LogEnclosures`), which
  remain hypotheses.

| File | Role |
|---|---|
| `SeedKernel.lean` | Root: imports the library |
| `SeedKernel/Packed.lean` | Shared packed-`Nat` lanes: `lane`, `repunit`, `pack`, `move`, `Signed`, floors, `reverseLanes`, `dot`, `gram`, fast variants |
| `SeedKernel/PackingProofs.lean` | Radix encoding: extraction and `pack = radixPack` |
| `SeedKernel/Obligations.lean`, `SeedKernel/PackedTests.lean` | Open correctness `Prop`s of the lanes; kernel-checked examples |
| `SeedKernel/Spin1/Packed.lean`, `Spin1/Optimized.lean` | Packed spin-1 `hApply`, `filterApply`, trial moments, `ritzCheck`; prepared-mask variants |
| `SeedKernel/Spin1/Obligations.lean`, `Spin1/Tests.lean` | Row specification and open `Prop`s (`HApplyCorrect`, `FilterCorrect`, ...); small-ring comparisons |
| `SeedKernel/Spin1/Counting.lean` | Conditional one-Ritz-value counting step and two obstructions |
| `SeedKernel/Spin1/H6.lean`, `H8.lean`, `H10.lean` | Generated Ritz/residual witnesses (`models/lean_seeds/certify.jl spin1_kernel`) |
| `SeedKernel/Ising/Checker.lean`, `Ising/Packed.lean`, `Ising/Chunked.lean` | Integer `20H`, direct, packed (on the shared lanes) and blockwise checkers |
| `SeedKernel/Ising/Generated.lean`, `Ising/GeneratedChunks.lean` | Generated proposals and their checks (`models/lean_seeds/certify.jl ising`) |
| `SeedKernel/Ising/ChunkCoverage.lean`, `Ising/NegativeChecks.lean` | Coverage of the 64 blocks; rejection of wrong proposals |
| `SeedKernel/Ising/Arithmetic.lean`, `Ising/ResidualLemma.lean` | Endpoint arithmetic of the log enclosures; the scalar weighted-residual lemma |
| `SeedKernel/Ising/SeedBridge.lean` | `LogEnclosures` and the conditional seed of `IsingLongitudinal` at its first point |
| `SeedKernel/Audit.lean` | `#print axioms` for every main declaration |
| `SeedKernel/Bench/*.lean` | Kernel benchmarks up to n = 12 (library `SeedKernelBench`) |

## Adding a model

A gap certificate enters Lean in exactly this form:

1. **The model as an `RPRing d ι`**: real matrices `C_α` (each symmetric or antisymmetric), a
   real symmetric on-site term `o` and optionally an internal reflection `refl` of the cell (a
   real symmetric involution commuting with `o`), with bond `h = −∑_α C_α ⊗ (refl C_α refl)`.
   A multi-site unit cell is a bigger local space `d` with the intra-cell terms in `o`. The
   symmetry conditions are finite matrix computations; `transpose_kronecker_eq_or`,
   `transpose_sum_kronecker_self` and `cellSwap_conj_smul_sum_kronecker_self` cover the usual
   cases.
2. **The points**: a `structure Point` with the model parameters (as `ℚ`), the seed inverse
   temperature `β₀ : ℚ` and a `cert : ChainCertificate`, and `def points : List Point`. Generate
   each certificate with `julia lean/scripts/chain_certificate.jl εₛ εₜ β₀ [J] [digits]` (the
   arguments are exact literals such as `a//b` or `0.064`; nothing is evaluated) and paste it in.
   Its first entries `s[0], t[0]` are the seed defects. `points_valid` and `points_β₀_pos` are
   checked by `decide +kernel`; the bound is `log r / β₀`.
3. **The seeds**: one `sorry` lemma `seeds : ∀ p ∈ points, 1 − s[0] ≤ (model p.x).spatialPurity n₀
   p.β₀ ∧ 1 − t[0] ≤ (model p.x).thermalPurity (2 n₀) p.β₀`, `n₀ ≥ 2` even and written as a
   literal. These are the only statements left to numerics.
4. **The gap**: `gap hp hL hLe` is `(model p.x).gapGT_of_chainCertificate (n₀ := n₀) (by norm_num)
   (by decide) (by exact_mod_cast points_β₀_pos p hp) (points_valid p hp) ⟨rfl, rfl, rfl⟩
   (seeds p hp).1 (seeds p hp).2 hL hLe`. See `Models/TwoLegIsing.lean` for a complete example
   (about 80 lines).
5. **A rotated frame** (if the model is reflection positive only after one-site unitaries `a`):
   in `Models/<Name>Physical.lean`, define the physical Hamiltonian, prove it is
   `piTensor a * R.hamiltonian L * piTensor aᴴ` with `RPRing.conj_hamiltonian`, and transfer the
   gap with `RPRing.eigenvalues_eq_spectrum`.
6. **Check the seeds against the certificate**: add the model to `lean/scripts/check_seeds.jl`
   and run `julia --project=julia lean/scripts/check_seeds.jl` from the repository root. It
   checks that every Lean `s[0], t[0]` is at least the certified defect at the same `(n₀, β₀)`.

The seeds of a certificate should therefore be reported as exact rationals together with
`(n₀, β₀)` and the exact `C_α, o`. Building blocks shared by several models (one-site matrices,
the cell swap) go in `Local/`, not in `Models/`.

**Naming scheme.** Every model uses the same names, so that a model is read the same way in
every file:

| Name | Meaning |
|---|---|
| `Models/<Name>.lean`, namespace `CertifiedGaps.<Name>` | `<Name>` is the Julia directory `models/<slug>` in UpperCamelCase (table below) |
| `model` | the model as an `RPRing`, a function of its real parameters |
| `Point`, `points` | the certified parameter points with `β₀` and `cert` |
| `points_valid`, `points_β₀_pos`, `points_<field>_…` | decidable facts about all points (`decide +kernel`) |
| `seeds` | the seed inequalities at every point: the one `sorry` of the model |
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
| `PurityCriterion.lean` | The abstract two-purity criterion: squaring, coupling, doubling, interpolation, purity to gap, seeds `≤ 1/100` |
| `Chain.lean` | The exact doubling recurrence with unequal defects; `gapGT_of_chain` with bound `log r / β₀` |
| `ChainCheck.lean` | `ChainCertificate`: the rational conditions of `gapGT_of_chain`, decidable |
| `Transfer/PiTensor.lean` | Site-wise tensor products; trace of `T^N` as a sum over closed paths |
| `Transfer/MomentForm.lean` | Traces of powers of a PSD matrix have the moment form |
| `Transfer/ProductLimit.lean` | Product formula: `(∏ₖ (1 + (c/M) hₖ))^M → e^{c ∑ hₖ}` |
| `Ring/Model.lean` | Ring geometry, `bondOp`/`siteOp`, reflection-positive rings `RPRing`: Hamiltonian, Hermiticity, `Re Tr e^{-βH} = ∑ e^{-βEₘ}` |
| `Ring/Transfer.lean` | Trotter slices (even bonds, on-site terms, odd bonds) and `Tr S^M = Tr (FᵀF)^{L/2}` with real `F` |
| `Ring/Gap.lean` | Moment inequalities for every `RPRing`; `RPRing.gapGT_of_chainCertificate` |
| `Ring/Conjugation.lean` | Sitewise unitary conjugation of the ring Hamiltonian; conjugate Hamiltonians have the same spectrum |
| `Heisenberg/Hamiltonian.lean` | The spin-1 ring `H_L = ∑ⱼ Sⱼ·Sⱼ₊₁` in the Cartesian basis |
| `Heisenberg/Gap.lean`, `Heisenberg/Main.lean` | Theorem 1 of the note, with the two seed `sorry`s |
| `Local/Pauli.lean`, `Local/SpinHalf.lean` | The real factors `X, Z` and `Sˣ, iSʸ, Sᶻ` (`SpinHalf.P`) |
| `Local/SpinHalfRotation.lean` | The rotation `e^{iπSʸ}` and the complex spin-½ matrices |
| `Local/SpinOne.lean` | The real antisymmetric spin-1 factors `SpinOne.G α = i S^α` |
| `Local/CellSwap.lean` | The swap `F` of a two-site cell, `F (A ⊗ B) F = B ⊗ A` |
| `Models/SpinOneLargeD.lean` | `∑ Sⱼ·Sⱼ₊₁ + D ∑ (Sᶻⱼ)²` as an `RPRing` (`D = 0` is the Heisenberg ring); `D = 4`, `D = 1603/500` |
| `Models/SpinOneDimerized.lean` | The dimerized spin-1 chain, internal reflection `F`; `δ = 9/10` |
| `Models/TwoLegIsing.lean` | The two-leg transverse-field Ising ladder; `g = 3` |
| `Models/IsingLongitudinal.lean` | The Ising chain in transverse and longitudinal fields; 30 grid points |
| `Models/HeisenbergLadder.lean`, `Models/StaggeredHeisenberg.lean`, `Models/DimerizedHeisenberg.lean` | The three rotated spin-½ models; seeds and gaps |
| `Models/*Physical.lean` | Their physical Hamiltonians, the rotations, and the physical gap theorems |
| `scripts/chain_certificate.jl` | Searches a `ChainCertificate` from the seeds (exact rationals, literal parser; Lean re-checks) |
| `scripts/check_seeds.jl` | Checks every Lean seed defect against the certificate JSON in `../models/*/results` |

## Why each piece is needed in Lean

On paper, the note gets its spatial input in one line: an exact thermal MPO whose transfer
operator `X_β` is self-adjoint, so `Z_n(β) = Tr X_β^n = ∑ᵢ |λᵢ|^n`. Lean needs every hidden step
written out. Each piece below closes one such gap.

**Gap as a statement about levels (`GapGT`).** The note writes `Δ_L = E₁ − E₀`. Mathlib returns
the eigenvalues of a Hermitian matrix as an unordered function on the basis, not as a sorted
list. So the conclusion is stated as "there is a level `m₀` with every other level more than
`c` above it". This also includes the uniqueness of the ground state, which the note proves
separately.

**Trace = spectral sum (`partitionFunction_eq_sum`).** The transfer operator and the seeds
concern matrix traces `Tr e^{-βH}`. The gap concerns eigenvalues. Diagonalizing `H` and
exponentiating connects the two. The note uses this silently.

**Moment inequalities instead of the moment form (`SpatialMomentBounds`).** The note's
`X_β` is an integral operator on bond histories in `L²(μ_β)`. Using it needs:
* a measure on histories;
* the spectral theorem for compact self-adjoint operators;
* traces of trace-class products.

Mathlib has little of this ready to use. The doubling argument uses `X_β` only through three
inequalities between traces:
* `Z_{2n} ≤ Z_n²`;
* the squaring inequality;
* `ℓ^p` monotonicity.

Stating exactly these makes the input smaller. They are also non-strict inequalities between
finitely many numbers, so they survive limits. That is what allows the finite route below.

**Trotter slices (`Ring/Transfer`).** These are the finite substitute for the continuous-time
MPO. Replacing `e^{-βH}` by `M` slices of first-order factors (even bonds `1 − τ hⱼ`, then the
on-site factors `1 − τ oᵢ`, then odd bonds) gives an exact tensor network of finite bond
dimension for an approximate Gibbs operator.

**The factorization (`Transfer/PiTensor`, `Ring/Transfer`).** The note draws the ring
of transfer vertices and reads off `Tr X^n`. In Lean this picture is a theorem:
* Every bond factor is expanded as `1 ⊗ 1 + τ ∑_α C_α ⊗ C_α`; every cell carries the word
  `∏ₜ C_{p t} O C_{q t}` with `O = 1 − τ o` between its even-bond and odd-bond factors.
* Products of site-wise tensors multiply site by site, and their trace factorizes over sites.
* The labels are regrouped by bond.
* "Sum over closed label paths = `Tr T^N`" turns the sum into a matrix power.

The proof also needs positivity, `T = FᵀF ⪰ 0`. Summing the even-bond labels pairs the even and
odd cells into `FᵀF`, and `Fᵀ = Fᴴ` because `C_α` and `o` are real, so `F` is. That is why the
model must be written in a real basis (the Cartesian basis `G_α = iS^α` for spin one, `iSʸ` for
spin one-half). This is the reflection positivity the note appeals to when it says the transfer
operator is self-adjoint. With an internal reflection `R` of the cell (dimerized
chains, where a bond joins the second site of a cell to the first site of the next), the word
on an odd cell is `R W R` for the word `W` on an even cell, so it has the same trace and the
transfer matrix is still `FᵀF`.

**Moment form from positivity (`Transfer/MomentForm`).** A PSD `T` has eigenvalues `λᵢ ≥ 0`, so
`Tr T^N = ∑ (√λᵢ)^{2N}`, which is the moment form at even lengths `2N`. The note uses
`|λᵢ|^n` for even `n`. Here the eigenvalues are nonnegative outright.

**Product formula (`Transfer/ProductLimit`).** The Trotter traces must converge to the true
`Tr e^{-βH_L}` at each fixed `L`. Mathlib has no Lie–Trotter or Chernoff formula, so it is
proved here for any ordered product of first-order factors. QICLean's
`Channel/Semigroup/ProductFormula.lean` has only the error estimate under an assumed step bound,
stated for superoperators.

**Why the limit is harmless (`SpatialMomentBounds.of_tendsto_hasMomentForm`).** The limit
`M → ∞` is taken in the input, not in the conclusion. So nothing about convergence of
eigenvalues or gaps of the Trotter operators is needed. This avoids Weyl's inequality and
Trotter error bounds. Only convergence of traces at each fixed `L` is used.

**The seeds (still `sorry`).** The doubling argument produces a gap from two purities near one
at a single scale. It cannot produce those purities itself. They carry the physics:
* The moment form alone holds for gapless models too.
* The classical Ising ring satisfies it, has `P_s → 1`, and still has a twofold degenerate
  ground state.
* Neither seed can be dropped (Theorem 2.24 of the local survey notes, not in this repository).

The seeds `P_s(2304, 784), P_τ(4608, 784) ≥ 0.99` are Proposition 5.5 of the OpenAI preprint.
They come from certified numerics, which `../models/spin1/` reproduces in Julia. Proving them in Lean
is the remaining open problem, and the hard one.

**Unequal seed defects (`Chain`, `ChainCheck`).** Certified seeds are rarely below `1/100`; for
the large-D chain at `D = 4` they are `(0.064, 0.041)`. `gapGT_of_chain` follows the exact
recurrence `s' ≥ f(z(s, t))`, `t' ≥ f(z(t, s'))` for `J` steps and then hands over to the
`≤ 1/100` theorem. Writing the bound as `log r / β₀` turns every remaining condition into an
inequality between rationals (`r^{2^j} z_j < 1 − z_j`, `r^{2^J} · 5ε ≤ 1`), so a certificate is a
short list of decimals that `decide` checks in seconds, with no logarithms or interval
arithmetic in the kernel. The loss against the floating-point infimum of the rates is about
`0.2%` at `D = 4` and for the ladder, and `3%` at the borderline `D = 1603/500`.

**Rotated frames (`Ring/Conjugation`, `Models/*Physical`).** The ladder, the staggered-field
chain and the dimerized spin-½ chain are reflection positive only after the rotation `e^{iπSʸ}`
on a sublattice. Each `*Physical.lean` file defines the physical Hamiltonian, proves it equals
the rotated one conjugated by that product of one-site unitaries (on even rings where the
sublattice requires it), and transfers the gap: unitarily equivalent Hermitian matrices have
the same characteristic polynomial, hence the same eigenvalues with multiplicity.
