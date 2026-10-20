/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Conjugation
import CertifiedGaps.Ring.Variational

/-!
# Matrix-product trial states on a ring

A periodic matrix-product state `ψ(σ) = Tr (A σ₀ ⋯ A σ_{N-1})` with integer tensors bounds the
partition function of any translation-invariant ring Hamiltonian `H_N = ∑ⱼ hⱼ,ⱼ₊₁` with a
symmetric integer bond `h` from below: if its energy density lies below `-κ`, then some
eigenvalue of `H_N` lies below `-κN` and `Re Tr e^{-βH_N} > e^{βκN}`. Both the norm and the
energy of `ψ` are traces of integer transfer matrices, so the hypothesis is a finite integer
inequality.

The computation runs in three steps.

1. **Sandwich of a product operator.** For `⨂ⱼ Mⱼ`, `⟨ψ, (⨂ⱼ Mⱼ) ψ⟩` is the trace of the
   ordered product of `Eⱼ = ∑ Mⱼ(s, t) Aₛ ⊗ Aₜ`, because `ψ(σ) ψ(τ)` is the trace of the
   Kronecker product of the two ordered products.
2. **One bond.** Expanding `h` in matrix units gives, for the bond `j → j + 1`, the transfer
   product with `Aₛ ⊗ Aₛ'` at `j`, `Aₜ ⊗ Aₜ'` at `j + 1` and `T` elsewhere. Rotating the ring
   (cyclicity of the trace) moves every bond to `0 → 1`, where the product is
   `(Aₛ Aₜ ⊗ Aₛ' Aₜ') Tⁿ`; summing against `h` gives `Tr (Tⁿ O)`.
3. **Variational principle** (`exp_lt_gibbsTrace_of_re_dotProduct_lt`).

## Main definitions

* `MPS.state A σ`, `MPS.transfer A = ∑ₛ Aₛ ⊗ Aₛ`, and the two-site energy insertion
  `MPS.insertion h A = ∑ h(st; s't') (Aₛ Aₜ) ⊗ (Aₛ' Aₜ')`.

## Main results

* `MPS.sum_state_sq`: `∑_σ ψ(σ)² = Tr Tᴺ`.
* `MPS.dotProduct_ringHamiltonian_mulVec`: `⟨ψ, H_N ψ⟩ = N Tr (T^(N-2) O)` for `N ≥ 2`.
* `MPS.exp_lt_gibbsTrace`: `Tr (T^(N-2) O) < -κ Tr Tᴺ → e^{βκN} < Z_N(β)`, with the
  integer form `MPS.exp_lt_gibbsTrace_int`.
-/

open Matrix NormedSpace
open scoped Kronecker

namespace CertifiedGaps

/-! ### Ordered products -/

section Kronecker

variable {R D : Type*} [CommSemiring R] [Fintype D] [DecidableEq D]

/-- An ordered product of Kronecker products is the Kronecker product of the ordered products. -/
theorem list_prod_ofFn_kronecker :
    ∀ (m : ℕ) (f g : Fin m → Matrix D D R),
      (List.ofFn fun j ↦ f j ⊗ₖ g j).prod = (List.ofFn f).prod ⊗ₖ (List.ofFn g).prod
  | 0, f, g => by simp
  | m + 1, f, g => by
    simp only [List.ofFn_succ, List.prod_cons, list_prod_ofFn_kronecker m, mul_kronecker_mul]

/-- The trace of an ordered product is invariant under a cyclic shift of the factors. -/
theorem trace_list_prod_ofFn_comp_finRotate {N : ℕ} (f : Fin N → Matrix D D R) :
    trace (List.ofFn (f ∘ finRotate N)).prod = trace (List.ofFn f).prod := by
  cases N with
  | zero => rfl
  | succ n =>
    rw [List.ofFn_succ', List.prod_concat, List.ofFn_succ, List.prod_cons, trace_mul_comm]
    simp

/-- The factors of a ring transfer product with one bond `j → k`: `X` at `j`, `Y` at `k`, `T`
elsewhere. -/
def bondFactors {M : Type*} {N : ℕ} (j k : Fin N) (X Y T : M) (i : Fin N) : M :=
  if i = j then X else if i = k then Y else T

/-- Rotating the bond and the factors together leaves the factors unchanged. -/
theorem bondFactors_finRotate {M : Type*} {N : ℕ} (j : Fin N) (X Y T : M) :
    bondFactors (finRotate N j) (finRotate N (finRotate N j)) X Y T ∘ finRotate N =
      bondFactors j (finRotate N j) X Y T := by
  ext1 i
  simp [bondFactors]

/-- By cyclicity of the trace, every bond of the ring gives the same transfer-matrix trace. -/
theorem trace_bondFactors {n : ℕ} (j : Fin (n + 1)) (X Y T : Matrix D D R) :
    trace (List.ofFn (bondFactors j (finRotate _ j) X Y T)).prod =
      trace (List.ofFn (bondFactors (0 : Fin (n + 1)) (finRotate (n + 1) 0) X Y T)).prod := by
  obtain ⟨m, hm⟩ := j
  induction m with
  | zero => rfl
  | succ m ih =>
    have hj : (⟨m + 1, hm⟩ : Fin (n + 1)) = finRotate _ ⟨m, by omega⟩ := by
      rw [finRotate_of_lt (by omega)]
    rw [hj, ← ih (by omega), ← trace_list_prod_ofFn_comp_finRotate, bondFactors_finRotate]

/-- The bond `0 → 1` on a ring of `n + 2` sites gives the product `X Y Tⁿ`. -/
theorem list_prod_bondFactors_zero (n : ℕ) (X Y T : Matrix D D R) :
    (List.ofFn (bondFactors (0 : Fin (n + 2)) (finRotate _ 0) X Y T)).prod = X * (Y * T ^ n) := by
  simp [List.ofFn_succ, bondFactors, Fin.succ_ne_zero,
    Fin.succ_succ_ne_one, List.ofFn_const, List.prod_replicate]

end Kronecker

/-! ### Integer matrices as complex matrices -/

section IntCast

variable {D : Type*} [Fintype D] [DecidableEq D]

/-- Casting an integer matrix to complex entries commutes with ordered products. -/
theorem map_list_prod_ofFn_intCast {m : ℕ} (f : Fin m → Matrix D D ℤ) :
    (List.ofFn f).prod.map (Int.cast : ℤ → ℂ) = (List.ofFn fun j ↦ (f j).map (↑)).prod := by
  change (Int.castRingHom ℂ).mapMatrix (List.ofFn f).prod = _
  rw [map_list_prod, List.map_ofFn]
  rfl

omit [DecidableEq D] in
/-- The trace of a cast integer matrix is the cast of the trace. -/
theorem trace_map_intCast (X : Matrix D D ℤ) :
    trace (X.map (Int.cast : ℤ → ℂ)) = ((trace X : ℤ) : ℂ) :=
  (AddMonoidHom.map_trace (Int.castRingHom ℂ) X).symm

omit [Fintype D] [DecidableEq D] in
/-- A symmetric integer matrix is Hermitian as a complex matrix. -/
theorem isHermitian_map_intCast {X : Matrix D D ℤ} (hX : Xᵀ = X) :
    (X.map (Int.cast : ℤ → ℂ)).IsHermitian := by
  ext p q
  simp only [conjTranspose_apply, map_apply, star_intCast]
  rw [← transpose_apply X, hX]

omit [Fintype D] in
/-- A two-site matrix as a combination of tensor products of matrix units. -/
theorem map_intCast_eq_sum_single_kronecker {d : Type*} [Fintype d] [DecidableEq d]
    (h : Matrix (d × d) (d × d) ℤ) :
    h.map (↑) = ∑ p : d × d, ∑ q : d × d,
      (h p q : ℂ) • (single p.1 q.1 (1 : ℂ) ⊗ₖ single p.2 q.2 (1 : ℂ)) := by
  conv_lhs => rw [matrix_eq_sum_single (h.map (Int.cast : ℤ → ℂ))]
  refine Finset.sum_congr rfl fun p _ ↦ Finset.sum_congr rfl fun q _ ↦ ?_
  rw [single_kronecker_single, smul_single, mul_one, smul_eq_mul, mul_one]
  rfl

end IntCast

/-! ### Matrix-product states -/

namespace MPS

variable {d : Type} [Fintype d] [DecidableEq d] {D : Type*} [Fintype D] [DecidableEq D]

/-- The periodic matrix-product state `ψ(σ) = Tr (A σ₀ ⋯ A σ_{N-1})`. -/
def state {N : ℕ} (A : d → Matrix D D ℤ) (σ : Fin N → d) : ℤ :=
  trace (List.ofFn fun j ↦ A (σ j)).prod

/-- The transfer matrix `T = ∑ₛ Aₛ ⊗ Aₛ` of a real matrix-product state. -/
def transfer (A : d → Matrix D D ℤ) : Matrix (D × D) (D × D) ℤ :=
  ∑ s, A s ⊗ₖ A s

/-- The two-site energy insertion `O = ∑ h(st; s't') (Aₛ Aₜ) ⊗ (Aₛ' Aₜ')` of the bond `h`, with
rows of `h` the output pair and columns the input pair. -/
def insertion (h : Matrix (d × d) (d × d) ℤ) (A : d → Matrix D D ℤ) :
    Matrix (D × D) (D × D) ℤ :=
  ∑ p : d × d, ∑ q : d × d, h p q • ((A p.1 * A p.2) ⊗ₖ (A q.1 * A q.2))

variable (A : d → Matrix D D ℤ)

/-- The matrix-product tensors as complex matrices. -/
def tensorC (s : d) : Matrix D D ℂ := (A s).map (↑)

/-- The matrix-product state as a complex vector on configurations. -/
def vec (N : ℕ) : (Fin N → d) → ℂ := fun σ ↦ (state A σ : ℂ)

omit [Fintype d] [DecidableEq d] in
/-- The state is the trace of the ordered product of the complex tensors. -/
theorem vec_apply {N : ℕ} (σ : Fin N → d) :
    vec A N σ = trace (List.ofFn fun j ↦ tensorC A (σ j)).prod := by
  rw [vec, state, ← trace_map_intCast, map_list_prod_ofFn_intCast]
  rfl

omit [Fintype d] [DecidableEq d] in
/-- The state is real. -/
theorem star_vec (N : ℕ) : star (vec A N) = vec A N := by
  ext σ
  simp [vec]

omit [DecidableEq d] in
/-- **Sandwich of a product operator.** For a product operator `⨂ⱼ Mⱼ`, the matrix element in
the state is the trace of the ordered product of the transfer matrices
`Eⱼ = ∑ M_j(s, t) Aₛ ⊗ Aₜ`. -/
theorem dotProduct_piTensor_mulVec {N : ℕ} (M : Fin N → Matrix d d ℂ) :
    star (vec A N) ⬝ᵥ (piTensor M *ᵥ vec A N) =
      trace (List.ofFn fun j ↦ ∑ s, ∑ t, M j s t • (tensorC A s ⊗ₖ tensorC A t)).prod := by
  have hsum : ∀ j, ∑ s, ∑ t, M j s t • (tensorC A s ⊗ₖ tensorC A t) =
      ∑ x : d × d, M j x.1 x.2 • (tensorC A x.1 ⊗ₖ tensorC A x.2) := fun j ↦
    (Fintype.sum_prod_type' (fun s t ↦ M j s t • (tensorC A s ⊗ₖ tensorC A t))).symm
  simp_rw [hsum]
  rw [list_prod_ofFn_sum, trace_sum,
    ← (Equiv.arrowProdEquivProdArrow (Fin N) (fun _ ↦ d) (fun _ ↦ d)).symm.sum_comp,
    Fintype.sum_prod_type, star_vec]
  simp only [dotProduct, mulVec, piTensor_apply, Finset.mul_sum]
  refine Finset.sum_congr rfl fun σ _ ↦ Finset.sum_congr rfl fun τ _ ↦ ?_
  rw [list_prod_ofFn_smul, trace_smul, list_prod_ofFn_kronecker, trace_kronecker,
    ← vec_apply, ← vec_apply]
  simp [smul_eq_mul, Equiv.arrowProdEquivProdArrow]
  ring

/-! ### One bond -/

/-- The transfer matrix `T = ∑ₛ Aₛ ⊗ Aₛ` over the complex numbers. -/
def transferC : Matrix (D × D) (D × D) ℂ := ∑ s, tensorC A s ⊗ₖ tensorC A s

/-- The energy insertion `O` over the complex numbers. -/
def insertionC (h : Matrix (d × d) (d × d) ℤ) : Matrix (D × D) (D × D) ℂ :=
  ∑ p : d × d, ∑ q : d × d,
    (h p q : ℂ) • ((tensorC A p.1 * tensorC A p.2) ⊗ₖ (tensorC A q.1 * tensorC A q.2))

omit [DecidableEq d] [Fintype D] [DecidableEq D] in
/-- The complex transfer matrix is the cast of the integer one. -/
theorem transferC_eq : transferC A = (transfer A).map (↑) := by
  ext x y
  simp [transferC, transfer, tensorC, Matrix.sum_apply, kroneckerMap_apply]

omit [DecidableEq d] [DecidableEq D] in
/-- The complex insertion is the cast of the integer one. -/
theorem insertionC_eq (h : Matrix (d × d) (d × d) ℤ) :
    insertionC A h = (insertion h A).map (↑) := by
  ext x y
  simp only [insertionC, insertion, tensorC, Matrix.map_apply, Matrix.sum_apply,
    Matrix.smul_apply, smul_eq_mul, kroneckerMap_apply, mul_apply]
  push_cast
  rfl

omit [Fintype D] [DecidableEq D] in
/-- The local transfer matrix `∑ 1(s, t) Aₛ ⊗ Aₜ` of the identity is `T`. -/
theorem sum_one_smul_kronecker :
    ∑ s, ∑ t, (1 : Matrix d d ℂ) s t • (tensorC A s ⊗ₖ tensorC A t) = transferC A := by
  simp [one_apply, ite_smul, transferC]

omit [Fintype D] [DecidableEq D] in
/-- The local transfer matrix of the matrix unit `|a⟩⟨b|` is `Aₐ ⊗ A_b`. -/
theorem sum_single_smul_kronecker (a b : d) :
    ∑ s, ∑ t, (single a b (1 : ℂ)) s t • (tensorC A s ⊗ₖ tensorC A t) =
      tensorC A a ⊗ₖ tensorC A b := by
  simp [single_apply, ite_smul, ite_and]

/-- Matrix element of one bond in the state, as transfer-matrix traces. -/
theorem dotProduct_bondOp_mulVec (h : Matrix (d × d) (d × d) ℤ) {N : ℕ} {j k : Fin N}
    (hjk : j ≠ k) :
    star (vec A N) ⬝ᵥ (bondOp j k (h.map (↑)) *ᵥ vec A N) =
      ∑ p : d × d, ∑ q : d × d, (h p q : ℂ) *
        trace (List.ofFn (bondFactors j k (tensorC A p.1 ⊗ₖ tensorC A q.1)
          (tensorC A p.2 ⊗ₖ tensorC A q.2) (transferC A))).prod := by
  simp only [map_intCast_eq_sum_single_kronecker h, bondOp_sum, bondOp_smul,
    bondOp_kronecker hjk, sum_mulVec, smul_mulVec, dotProduct_sum, dotProduct_smul,
    dotProduct_piTensor_mulVec, smul_eq_mul]
  refine Finset.sum_congr rfl fun p _ ↦ Finset.sum_congr rfl fun q _ ↦ ?_
  congr 4
  ext1 i
  unfold bondFactors
  split_ifs
  · exact sum_single_smul_kronecker A _ _
  · exact sum_single_smul_kronecker A _ _
  · exact sum_one_smul_kronecker A

/-- **One bond of the energy.** On a ring of `n + 2` sites, each bond contributes
`Tr (Tⁿ O)`. -/
theorem dotProduct_bond_mulVec (h : Matrix (d × d) (d × d) ℤ) (n : ℕ) (j : Fin (n + 2)) :
    star (vec A (n + 2)) ⬝ᵥ (bondOp j (finRotate _ j) (h.map (↑)) *ᵥ vec A (n + 2)) =
      trace (transferC A ^ n * insertionC A h) := by
  rw [dotProduct_bondOp_mulVec A h (ne_finRotate (by omega) j), insertionC]
  simp only [trace_bondFactors j, list_prod_bondFactors_zero, Finset.mul_sum, trace_sum,
    mul_smul_comm, trace_smul, smul_eq_mul]
  refine Finset.sum_congr rfl fun p _ ↦ Finset.sum_congr rfl fun q _ ↦ ?_
  rw [← mul_assoc, trace_mul_comm, mul_kronecker_mul]

/-! ### The norm and the energy -/

omit [DecidableEq d] in
/-- **Norm of the state.** `∑_σ ψ(σ)² = Tr Tᴺ`. -/
theorem sum_state_sq (N : ℕ) :
    ∑ σ : Fin N → d, state A σ ^ 2 = trace (transfer A ^ N) := by
  have h : ∀ σ : Fin N → d, state A σ ^ 2 =
      trace (List.ofFn fun j ↦ A (σ j) ⊗ₖ A (σ j)).prod := fun σ ↦ by
    rw [list_prod_ofFn_kronecker, trace_kronecker, state, sq]
  simp_rw [h, ← trace_sum, ← list_prod_ofFn_sum N (fun _ s ↦ A s ⊗ₖ A s), List.ofFn_const,
    List.prod_replicate, transfer]

omit [DecidableEq d] in
/-- The norm `⟨ψ, ψ⟩ = Tr Tᴺ` of the state as a complex vector. -/
theorem dotProduct_vec (N : ℕ) :
    star (vec A N) ⬝ᵥ vec A N = ((trace (transfer A ^ N) : ℤ) : ℂ) := by
  rw [star_vec, ← sum_state_sq]
  simp [dotProduct, vec, sq]

omit [DecidableEq d] in
/-- `Tr (Tⁿ O)` over the complex numbers is the cast of the integer trace. -/
theorem trace_transferC_pow_mul (h : Matrix (d × d) (d × d) ℤ) (n : ℕ) :
    trace (transferC A ^ n * insertionC A h) =
      ((trace (transfer A ^ n * insertion h A) : ℤ) : ℂ) := by
  rw [transferC_eq, insertionC_eq, ← trace_map_intCast]
  congr 1
  change (Int.castRingHom ℂ).mapMatrix _ ^ n * (Int.castRingHom ℂ).mapMatrix _ =
    (Int.castRingHom ℂ).mapMatrix _
  rw [map_mul, map_pow]

/-- **Energy of the state.** On a ring of `N ≥ 2` sites, `⟨ψ, H ψ⟩ = N Tr (T^(N-2) O)`: every
bond contributes the same by cyclicity of the trace. -/
theorem dotProduct_ringHamiltonian_mulVec (h : Matrix (d × d) (d × d) ℤ) {N : ℕ} (hN : 2 ≤ N) :
    star (vec A N) ⬝ᵥ (ringHamiltonian (h.map (↑)) N *ᵥ vec A N) =
      (N : ℂ) * ((trace (transfer A ^ (N - 2) * insertion h A) : ℤ) : ℂ) := by
  obtain ⟨n, rfl⟩ : ∃ n, N = n + 2 := ⟨N - 2, by omega⟩
  rw [ringHamiltonian, sum_mulVec, dotProduct_sum]
  simp only [dotProduct_bond_mulVec, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul, Nat.add_sub_cancel, trace_transferC_pow_mul]

/-! ### The trial-state bound on the partition function -/

/-- **Trial-state bound.** For a symmetric integer bond `h` and integer tensors `A` on a ring of
`N ≥ 2` sites: if the state has energy density below `-κ`, i.e. `Tr (T^(N-2) O) < -κ Tr Tᴺ`,
then `Z_N(β) > e^{βκN}` for every `β > 0`. -/
theorem exp_lt_gibbsTrace (h : Matrix (d × d) (d × d) ℤ) (hh : hᵀ = h) {N : ℕ}
    (hN : 2 ≤ N) (κ : ℚ)
    (hU : ((trace (transfer A ^ (N - 2) * insertion h A) : ℤ) : ℚ) <
      -κ * ((trace (transfer A ^ N) : ℤ) : ℚ))
    {β : ℝ} (hβ : 0 < β) :
    Real.exp (β * κ * N) < gibbsTrace (ringHamiltonian (h.map (↑)) N) β := by
  set U := trace (transfer A ^ (N - 2) * insertion h A)
  set V := trace (transfer A ^ N)
  have hUR : (U : ℝ) < -(κ : ℝ) * V := by exact_mod_cast hU
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hgibbs := exp_lt_gibbsTrace_of_re_dotProduct_lt
    (ringHamiltonian_isHermitian (isHermitian_map_intCast hh) N) (vec A N)
    (E := -(κ : ℝ) * N) (by
      rw [dotProduct_ringHamiltonian_mulVec A h hN, dotProduct_vec]
      simp only [← Complex.ofReal_natCast, ← Complex.ofReal_intCast, ← Complex.ofReal_mul,
        Complex.ofReal_re]
      nlinarith) hβ
  rw [show -(β * (-(κ : ℝ) * N)) = β * κ * N by ring] at hgibbs
  exact hgibbs

/-- The trial-state bound with an integer certificate `a Tr (T^(N-2) O) < -b Tr Tᴺ`, `a > 0`:
then `Z_N(β) > e^{β (b/a) N}`. -/
theorem exp_lt_gibbsTrace_int (h : Matrix (d × d) (d × d) ℤ) (hh : hᵀ = h) {N : ℕ}
    (hN : 2 ≤ N) {a b : ℤ} (ha : 0 < a)
    (hU : a * trace (transfer A ^ (N - 2) * insertion h A) < -b * trace (transfer A ^ N))
    {β : ℝ} (hβ : 0 < β) :
    Real.exp (β * (b / a) * N) < gibbsTrace (ringHamiltonian (h.map (↑)) N) β := by
  set U := trace (transfer A ^ (N - 2) * insertion h A)
  set V := trace (transfer A ^ N)
  have haQ : (0 : ℚ) < a := by exact_mod_cast ha
  have hUQ : (a : ℚ) * U < -(b : ℚ) * V := by exact_mod_cast hU
  have hZ := exp_lt_gibbsTrace A h hh hN ((b : ℚ) / a) (by
    rw [show -((b : ℚ) / a) * V = (-(b : ℚ) * V) / a by ring, lt_div_iff₀ haQ]
    linarith) hβ
  simpa using hZ

/-- **Trial-state bound on the ground-state energy.** If `Tr (T^(N-2) O) < -κ Tr Tᴺ` on a ring of
`N ≥ 2` sites, some eigenvalue of the ring Hamiltonian lies below `-κN` (the variational
principle; the bound on `Z_N` follows from it). -/
theorem exists_eigenvalue_lt (h : Matrix (d × d) (d × d) ℤ) (hh : hᵀ = h) {N : ℕ}
    (hN : 2 ≤ N) (κ : ℚ)
    (hU : ((trace (transfer A ^ (N - 2) * insertion h A) : ℤ) : ℚ) <
      -κ * ((trace (transfer A ^ N) : ℤ) : ℚ)) :
    ∃ i, (ringHamiltonian_isHermitian (isHermitian_map_intCast hh) N).eigenvalues i <
      -(κ : ℝ) * N := by
  set U := trace (transfer A ^ (N - 2) * insertion h A)
  set V := trace (transfer A ^ N)
  have hUR : (U : ℝ) < -(κ : ℝ) * V := by exact_mod_cast hU
  have hNpos : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  refine exists_eigenvalue_lt_of_re_dotProduct_lt _ (vec A N) ?_
  rw [dotProduct_ringHamiltonian_mulVec A h hN, dotProduct_vec]
  simp only [← Complex.ofReal_natCast, ← Complex.ofReal_intCast, ← Complex.ofReal_mul,
    Complex.ofReal_re]
  nlinarith

/-- The eigenvalue bound with an integer certificate `a Tr (T^(N-2) O) < -b Tr Tᴺ`, `a > 0`. -/
theorem exists_eigenvalue_lt_int (h : Matrix (d × d) (d × d) ℤ) (hh : hᵀ = h) {N : ℕ}
    (hN : 2 ≤ N) {a b : ℤ} (ha : 0 < a)
    (hU : a * trace (transfer A ^ (N - 2) * insertion h A) < -b * trace (transfer A ^ N)) :
    ∃ i, (ringHamiltonian_isHermitian (isHermitian_map_intCast hh) N).eigenvalues i <
      -((b : ℝ) / a) * N := by
  have haQ : (0 : ℚ) < a := by exact_mod_cast ha
  have hUQ : (a : ℚ) * trace (transfer A ^ (N - 2) * insertion h A) <
      -(b : ℚ) * trace (transfer A ^ N) := by exact_mod_cast hU
  obtain ⟨i, hi⟩ := exists_eigenvalue_lt A h hh hN ((b : ℚ) / a) (by
    set V : ℚ := ((trace (transfer A ^ N) : ℤ) : ℚ)
    rw [show -((b : ℚ) / a) * V = (-(b : ℚ) * V) / a by ring, lt_div_iff₀ haQ]
    linarith)
  exact ⟨i, by push_cast at hi; exact hi⟩

end MPS

end CertifiedGaps
