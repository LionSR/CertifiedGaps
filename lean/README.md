# Lean checks of two-purity gap certificates

A Lean 4 project, with Mathlib as its only dependency, that proves uniform spectral gaps of
quantum spin rings with the two-purity doubling argument of the thermal-MPO note
(`../notes/thermal-mpo-gap-note.tex`). Everything is proved except the seed inequalities: two
numbers per model, stated as named `sorry` lemmas. All declarations live in the namespace
`CertifiedGaps`; the names below are relative to it.

```sh
cd lean
lake exe cache get   # prebuilt Mathlib; never build Mathlib from source
lake build           # about 4 minutes from scratch, warning-free except the seed `sorry`s
```

Lean and Mathlib are pinned at `v4.35.0-rc3`. The package has about 4400 lines.

## Status

`GapGT e c` says that the lowest level of the spectrum `e` is nondegenerate and every other level
lies more than `c` above it. Every gap theorem below holds for every even ring length
`L ≥ n₀` (counted in cells).

### Model theorems (seeds open)

| Theorem | Model | Bound | Open seeds (the only `sorry`s) |
|---|---|---|---|
| `SpinOneHeisenberg.spin_one_heisenberg_gap` | spin-1 Heisenberg ring | `Δ_L > log 20 / 784` for `L ≥ 2304` | `seed_spatial`: `P_s(2304, 784) ≥ 0.99`; `seed_thermal`: `P_τ(4608, 784) ≥ 0.99` |
| `SpinOneLargeD.gap_D4` | spin-1, `+ 4 ∑ (Sᶻ)²` | `log 4.755 / (9/4) ≈ 0.693`, `L ≥ 6` | `seed_spatial_D4`, `seed_thermal_D4`: `P_s(6, 9/4) ≥ 1 − 0.06411`, `P_τ(12, 9/4) ≥ 1 − 0.04053` |
| `SpinOneLargeD.gap_D3206` | spin-1, `+ (1603/500) ∑ (Sᶻ)²` | `log 1.0065429 / (381/125) ≈ 0.00214`, `L ≥ 6` | `seed_spatial_D3206`, `seed_thermal_D3206` at `(6, 381/125)` |
| `SpinOneDimerized.gap_d910` | spin-1 dimerized chain, `δ = 9/10` (cells = weak pairs) | `log 6.339 / (7/2) ≈ 0.528`, `L ≥ 2` | `seed_spatial_d910`, `seed_thermal_d910`: `P_s(2, 7/2) ≥ 1 − 0.04972`, `P_τ(4, 7/2) ≥ 1 − 0.03206` |
| `TwoLegLadder.physical_gap_J4` (rotated: `gap_J4`) | spin-½ ladder `∑_{legs} S·S + 4 ∑_{rungs} S·S` | `log 9.033 / 2 ≈ 1.100`, `L ≥ 4` rungs | `seed_spatial_J4`, `seed_thermal_J4`: `P_s(4, 2) ≥ 1 − 0.03509`, `P_τ(8, 2) ≥ 1 − 0.02463` |
| `StaggeredHeisenberg.physical_gap_h1`, `_h34`, `_h12`, `_h720`, `_h14` (rotated: `gap_h*`) | `∑ S·S − h ∑ (−1)ʲ Sᶻⱼ`, `h = 1, 3/4, 1/2, 7/20, 1/4` | `1.147, 0.890, 0.594, 0.383, 0.214`, `L ≥ 8` | `seed_spatial_h*`, `seed_thermal_h*` at `(8, β₀)` |
| `DimerizedHeisenberg.physical_gap_d12`, `_d310`, `_d14` (rotated: `gap_d*`) | spin-½ dimerized chain, `δ = 1/2, 3/10, 1/4` (cells = strong pairs) | `0.532, 0.172, 0.057`, `L ≥ 4` | `seed_spatial_d*`, `seed_thermal_d*` at `(4, β₀)` |
| `IsingLadder.gap_G3` | two-leg TFIM, `g = 3` | `log 2.585 / (3/2) ≈ 0.633`, `L ≥ 4` rungs | `seed_spatial_G3`, `seed_thermal_G3`: `P_s(4, 3/2) ≥ 1 − 0.1125`, `P_τ(8, 3/2) ≥ 1 − 0.05162` |
| `IsingLongitudinal.gap` | Ising chain `−ZZ − gX − hZ`, 30 points `g ∈ {1/4, …, 2}`, `h ∈ {1/20, …, 1}` | `log r / β₀` per point (`0.075` at `g = 1, h = 1/20` to `3.70`), `L ≥ 6` | `seeds` (one lemma for all 30 points) |

`#print axioms` on each of these lists `propext`, `Classical.choice`, `Quot.sound` and `sorryAx`
(from its seeds). The `physical_*` theorems are about the physical Hamiltonian; the rotated
`gap_*` theorems about the reflection-positive frame, which has the same spectrum.

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
| `*_valid`, `IsingLongitudinal.points_valid` | every chain certificate, by `decide +kernel` |

The seeds of the models come from the certificates in `../models/` (complete spectra of the seed
rings with rigorous eigenvalue enclosures; each model file names its directory). The spin-1
Heisenberg seeds are OpenAI's Proposition 5.5, reproduced in `../models/spin1/`.

## Adding a model

A gap certificate enters Lean in exactly this form:

1. **The model as an `RPRing d ι`**: real matrices `C_α` (each symmetric or antisymmetric), a
   real symmetric on-site term `o` and optionally an internal reflection `refl` of the cell (a
   real symmetric involution commuting with `o`), with bond `h = −∑_α C_α ⊗ (refl C_α refl)`.
   A multi-site unit cell is a bigger local space `d` with the intra-cell terms in `o`. The
   symmetry conditions are finite matrix computations; `transpose_kronecker_eq_or`,
   `transpose_sum_kronecker_self` and `cellSwap_conj_smul_sum_kronecker_self` cover the usual
   cases.
2. **The seeds**: `1 − εₛ ≤ R.spatialPurity n₀ β₀` and `1 − εₜ ≤ R.thermalPurity (2 n₀) β₀` with
   decimal `εₛ, εₜ`, `n₀ ≥ 2` even. These are the only statements left to numerics.
3. **A `ChainCertificate`**: generate it with
   `julia lean/scripts/chain_certificate.jl εₛ εₜ β₀ [J] [digits]` and paste it in. Lean checks it
   with `decide +kernel`. The resulting bound is `log r / β₀`.
4. **The gap**: `R.gapGT_of_chainCertificate (by norm_num) (by decide) (by norm_num) cert_valid
   (by norm_num [cert, ChainCertificate.sAt, ChainCertificate.tAt]) seed_spatial seed_thermal hL hLe`.
   See `Models/TwoLegLadder.lean` for a complete example (about 80 lines).
5. **A rotated frame** (if the model is reflection positive only after one-site unitaries `a`):
   define the physical Hamiltonian, prove it is `piTensor a * R.hamiltonian L * piTensor aᴴ` with
   `RPRing.conj_hamiltonian`, and transfer the gap with `RPRing.eigenvalues_eq_spectrum`.

The seeds of a certificate should therefore be reported as exact rationals together with
`(n₀, β₀)` and the exact `C_α, o`.

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
| `Models/Pauli.lean`, `Models/SpinHalf.lean` | The real factors `X, Z` and `Sˣ, iSʸ, Sᶻ` |
| `Models/SpinHalfRotation.lean` | The rotation `e^{iπSʸ}` and the complex spin-½ matrices |
| `Models/CellSwap.lean` | The swap `F` of a two-site cell, `F (A ⊗ B) F = B ⊗ A` |
| `Models/SpinOneLargeD.lean` | `∑ Sⱼ·Sⱼ₊₁ + D ∑ (Sᶻⱼ)²` as an `RPRing` (`D = 0` is the Heisenberg ring); `D = 4`, `D = 1603/500` |
| `Models/SpinOneDimerized.lean` | The dimerized spin-1 chain, internal reflection `F`; `δ = 9/10` |
| `Models/IsingLadder.lean` | The two-leg transverse-field Ising ladder; `g = 3` |
| `Models/IsingLongitudinal.lean` | The Ising chain in transverse and longitudinal fields; 30 grid points |
| `Models/TwoLegLadder.lean`, `Models/StaggeredHeisenberg.lean`, `Models/DimerizedHeisenberg.lean` | The three rotated spin-½ models; seeds and gaps |
| `Models/*Physical.lean` | Their physical Hamiltonians, the rotations, and the physical gap theorems |
| `scripts/chain_certificate.jl` | Searches a `ChainCertificate` from the seeds (exact rationals; Lean re-checks) |

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
