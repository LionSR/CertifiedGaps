# Lean check of the spin-1 Heisenberg gap

A minimal Lean 4 project, with Mathlib as its only dependency, that checks Theorem 1 of the
thermal-MPO note (`../outputs/thermal-mpo-gap-note.tex`). The periodic spin-1 Heisenberg ring
has a unique ground state and gap `Δ_L > (log 20)/784` for every even `L ≥ 2304`. The proof is
the note's two-purity doubling argument. It needs none of the twisted rings, symmetry sectors,
polynomial filters or trial states of OpenAI's *The periodic spin-one Haldane gap*.

## Status

```lean
theorem SpinOneHeisenberg.spin_one_heisenberg_gap {L : ℕ} (hL : 2304 ≤ L) (hLe : Even L) :
    GapGT (spectrum L) (log 20 / 784)
```

Everything is proved except the two numerical seeds, `seed_spatial` and `seed_thermal` in
`Heisenberg/Main.lean`. Until those are proved, `#print axioms` on this theorem shows `sorryAx`
alongside `propext`, `Classical.choice` and `Quot.sound`. The transfer-operator input
`spatialMomentBounds_partitionFunction` is proved and uses only the three standard axioms.

```sh
cd lean
lake exe cache get   # prebuilt Mathlib; never build Mathlib from source
lake build
```

Lean and Mathlib are pinned at `v4.35.0-rc3`. The project has about 1900 lines.

## Files

| File | Role |
|---|---|
| `PurityCriterion.lean` | The abstract two-purity criterion: squaring, coupling, doubling, interpolation, purity to gap |
| `Heisenberg/Hamiltonian.lean` | The ring `H_L = ∑ⱼ Sⱼ·Sⱼ₊₁`, its Hermiticity, and `Re Tr e^{-βH} = ∑ e^{-βEₘ}` |
| `Heisenberg/Gap.lean` | `ThermalMPOInputs` and the conditional theorem `gapGT_log_twenty` |
| `Heisenberg/Main.lean` | The unconditional theorem, with the two seed `sorry`s |
| `Transfer/Defs.lean` | Trotter slices `∏ⱼ (1 − τ hⱼ)` (even bonds, then odd) and their traces |
| `Transfer/PiTensor.lean` | Site-wise tensor products of `3 × 3` matrices; trace of `T^N` as a sum over closed paths |
| `Transfer/Factorization.lean` | `Tr S^M = Tr (FᴴF)^{L/2}` for a finite matrix `F` |
| `Transfer/MomentForm.lean` | Traces of powers of a PSD matrix have the moment form |
| `Transfer/ProductLimit.lean` | Product formula: `(∏ₖ (1 + (c/M) hₖ))^M → e^{c ∑ hₖ}` |
| `Transfer/Moments.lean` | Assembly: the moment inequalities for `Re Tr e^{-βH_L}` |

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

**Trotter slices (`Transfer/Defs`).** These are the finite substitute for the continuous-time
MPO. Replacing `e^{-βH}` by `M` slices of first-order bond factors `1 − τ hⱼ`, even bonds then
odd bonds, gives an exact tensor network of finite bond dimension `4^M` for an approximate Gibbs
operator.

**The factorization (`Transfer/PiTensor`, `Transfer/Factorization`).** The note draws the ring
of transfer vertices and reads off `Tr X^n`. In Lean this picture is a theorem:
* Every bond factor is expanded as `1 ⊗ 1 + τ ∑_α G_α ⊗ G_α`, with `G_α = iS^α`.
* Products of site-wise tensors multiply site by site, and their trace factorizes over sites.
* The labels are regrouped by bond.
* "Sum over closed label paths = `Tr T^N`" turns the sum into a matrix power.

The proof also needs positivity, `T = FᴴF ⪰ 0`. In the even/odd ordering, odd sites carry `F`
and even sites carry `Fᵀ`, and `Fᵀ = Fᴴ` because the entries of `F` are real. That is why the
Cartesian basis is used: `G_α = iS^α` is real. This is the reflection positivity the note
appeals to when it says the transfer operator is self-adjoint.

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
