/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.MatrixProduct

/-!
# Matrix-product trial states on an open chain

An open matrix-product state `ψ(σ) = bᵀ A_{σ₀} ⋯ A_{σ_{N-1}} c` with integer tensors and
integer end vectors bounds the ground-state energy of an open chain

  `H = ∑ⱼ h_{j,j+1} + f₀ ⊗ 1 ⊗ ⋯ ⊗ 1 + 1 ⊗ ⋯ ⊗ 1 ⊗ f₁`

(`openChainHamiltonian`) from above: if its Rayleigh quotient lies below `E`, some eigenvalue of
`H` lies below `E`. Its norm and its energy are matrix elements `L₀ ⋯ u` of integer transfer
matrices between `L₀ = b ⊗ b` and `u = c ⊗ c`, so the hypothesis is a finite integer inequality.
This is the open-chain analogue of `CertifiedGaps/Transfer/MatrixProduct.lean` and reuses its
transfer matrix `MPS.transfer A = ∑ₛ Aₛ ⊗ Aₛ` and bond insertion `MPS.insertion h A`.

The computation runs in three steps.

1. **Sandwich of a product operator.** For `⨂ⱼ Mⱼ`, `⟨ψ, (⨂ⱼ Mⱼ) ψ⟩ = L₀ (∏ⱼ Eⱼ) u` with
   `Eⱼ = ∑ Mⱼ(s, t) Aₛ ⊗ Aₜ`, because `ψ(σ) ψ(τ) = (b ⊗ b)ᵀ (∏ A_{σⱼ} ⊗ A_{τⱼ}) (c ⊗ c)`.
2. **Bonds and end sites.** The bond `j, j + 1` gives the factors `Tʲ O T^(n-1-j)` and an end
   site the factors `G Tⁿ` or `Tⁿ G`, with `G = ∑ f(s, t) Aₛ ⊗ Aₜ` (`siteInsertion`).
3. **Variational principle** (`exists_eigenvalue_lt_of_re_dotProduct_lt`).

## Main definitions

* `openChainHamiltonian h f₀ f₁ N`: the open chain of `N` sites.
* `OpenMPS.state b c A σ`, `OpenMPS.pair b = b ⊗ b`, `OpenMPS.siteInsertion f A`.
* The integer contractions `OpenMPS.norm`, `OpenMPS.bondEnergy`, `OpenMPS.leftField`,
  `OpenMPS.rightField`.

## Main results

* `conj_openChainHamiltonian`: a sitewise change of basis of the open chain.
* `OpenMPS.sum_state_sq`: `∑_σ ψ(σ)² = L₀ Tᴺ u`.
* `OpenMPS.dotProduct_openChainHamiltonian_mulVec`: `⟨ψ, H ψ⟩` on `n + 1` sites, with end
  fields `μ₀ f₀` and `μ₁ f₁`, is `∑ⱼ L₀ Tʲ O T^(n-1-j) u + μ₀ L₀ G₀ Tⁿ u + μ₁ L₀ Tⁿ G₁ u`.
* `OpenMPS.exists_eigenvalue_lt`: the variational bound with these integers.
-/

open Matrix
open scoped Kronecker

namespace CertifiedGaps

/-! ### The open chain -/

section OpenChain

variable {d : Type} [Fintype d] [DecidableEq d]

/-- The open chain `∑ⱼ h_{j,j+1} + f₀ ⊗ 1 ⊗ ⋯ ⊗ 1 + 1 ⊗ ⋯ ⊗ 1 ⊗ f₁` on `N` sites, with end
terms `f₀` on the first and `f₁` on the last site (`0` for `N = 0`; both end terms act on the
same site for `N = 1`). -/
def openChainHamiltonian (h : Matrix (d × d) (d × d) ℂ) (f₀ f₁ : Matrix d d ℂ) :
    (N : ℕ) → Matrix (Fin N → d) (Fin N → d) ℂ
  | 0 => 0
  | n + 1 => ∑ j : Fin n, bondOp j.castSucc j.succ h + siteOp 0 f₀ + siteOp (Fin.last n) f₁

omit [Fintype d] in
/-- The open chain of Hermitian terms is Hermitian. -/
theorem openChainHamiltonian_isHermitian {h : Matrix (d × d) (d × d) ℂ} {f₀ f₁ : Matrix d d ℂ}
    (hh : h.IsHermitian) (hf₀ : f₀.IsHermitian) (hf₁ : f₁.IsHermitian) (N : ℕ) :
    (openChainHamiltonian h f₀ f₁ N).IsHermitian := by
  cases N with
  | zero => exact isHermitian_zero
  | succ n =>
    exact ((isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian _ _ hh).add
      (siteOp_isHermitian _ hf₀)).add (siteOp_isHermitian _ hf₁)

/-- **Sitewise change of basis of the open chain**: the same unitary `V` on every site turns
each bond `h` into `(V ⊗ V)ᴴ h (V ⊗ V)` and each end term `f` into `Vᴴ f V`. -/
theorem conj_openChainHamiltonian (V : Matrix d d ℂ) (hV : Vᴴ * V = 1)
    (h : Matrix (d × d) (d × d) ℂ) (f₀ f₁ : Matrix d d ℂ) (N : ℕ) :
    (piTensor fun _ : Fin N ↦ V)ᴴ * openChainHamiltonian h f₀ f₁ N *
        piTensor (fun _ : Fin N ↦ V) =
      openChainHamiltonian ((V ⊗ₖ V)ᴴ * h * (V ⊗ₖ V)) (Vᴴ * f₀ * V) (Vᴴ * f₁ * V) N := by
  cases N with
  | zero => simp [openChainHamiltonian]
  | succ n =>
    have hsite (j : Fin (n + 1)) (f : Matrix d d ℂ) :
        (piTensor fun _ : Fin (n + 1) ↦ V)ᴴ * siteOp j f * piTensor (fun _ ↦ V) =
          siteOp j (Vᴴ * f * V) := by
      have h := conj_siteOp (fun _ : Fin (n + 1) ↦ Vᴴ)
        (fun _ ↦ by rw [conjTranspose_conjTranspose, hV]) j f
      simp only [conjTranspose_conjTranspose] at h
      rw [conjTranspose_piTensor, h]
    simp only [openChainHamiltonian, Matrix.mul_add, Matrix.add_mul, Finset.mul_sum,
      Finset.sum_mul, hsite]
    congr 2
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    exact conj_bondOp_const V hV (show j.castSucc < j.succ from Fin.castSucc_lt_succ).ne _

end OpenChain

/-! ### Ordered products with one or two marked factors -/

section MarkedProducts

variable {M : Type*} [Monoid M]

/-- An ordered product of `T`s with `X` at position `j`: `Tʲ X T^(N-1-j)`. -/
theorem list_prod_ofFn_site : ∀ (N : ℕ) (j : Fin N) (X T : M),
    (List.ofFn fun i : Fin N ↦ if i = j then X else T).prod = T ^ (j : ℕ) * X * T ^ (N - 1 - j)
  | 0, j, _, _ => j.elim0
  | N + 1, j, X, T => by
    rw [List.ofFn_succ, List.prod_cons]
    refine Fin.cases ?_ (fun j ↦ ?_) j
    · have h : (List.ofFn fun i : Fin N ↦ if i.succ = (0 : Fin (N + 1)) then X else T) =
          List.ofFn fun _ : Fin N ↦ T := by
        congr 1
      rw [h, List.ofFn_const, List.prod_replicate]
      simp
    · have h : (List.ofFn fun i : Fin N ↦ if i.succ = j.succ then X else T) =
          List.ofFn fun i ↦ if i = j then X else T := by
        congr 1; ext1 i; simp [Fin.succ_inj]
      have e : N + 1 - 1 - (j.succ : ℕ) = N - 1 - j := by simp only [Fin.val_succ]; omega
      rw [h, list_prod_ofFn_site N j X T, ite_eq_right_iff.2 (fun h ↦ absurd h
        (Fin.succ_ne_zero j).symm), e, Fin.val_succ, pow_succ', mul_assoc, mul_assoc, mul_assoc]

/-- An ordered product of `T`s with `X Y` at the bond `j, j + 1`: `Tʲ X Y T^(n-1-j)`. -/
theorem list_prod_ofFn_bond : ∀ (n : ℕ) (j : Fin n) (X Y T : M),
    (List.ofFn (bondFactors j.castSucc j.succ X Y T)).prod =
      T ^ (j : ℕ) * (X * (Y * T ^ (n - 1 - j)))
  | 0, j, _, _, _ => j.elim0
  | n + 1, j, X, Y, T => by
    refine Fin.cases ?_ (fun j ↦ ?_) j
    · rw [List.ofFn_succ, List.prod_cons, List.ofFn_succ, List.prod_cons]
      have h : (List.ofFn fun i : Fin n ↦ bondFactors (0 : Fin (n + 1)).castSucc
          (0 : Fin (n + 1)).succ X Y T i.succ.succ) = List.ofFn fun _ : Fin n ↦ T := by
        congr 1
      rw [h, List.ofFn_const, List.prod_replicate]
      simp [bondFactors]
    · rw [List.ofFn_succ, List.prod_cons]
      have h0 : bondFactors j.succ.castSucc j.succ.succ X Y T 0 = T := by
        simp [bondFactors, Fin.ext_iff]
      have hs : (List.ofFn fun i : Fin (n + 1) ↦
          bondFactors j.succ.castSucc j.succ.succ X Y T i.succ) =
          List.ofFn (bondFactors j.castSucc j.succ X Y T) := by
        congr 1
        ext1 i
        simp only [bondFactors, ← Fin.succ_castSucc, Fin.succ_inj]
      have e : n + 1 - 1 - (j.succ : ℕ) = n - 1 - j := by simp only [Fin.val_succ]; omega
      rw [h0, hs, list_prod_ofFn_bond n j X Y T, e, Fin.val_succ, pow_succ', mul_assoc]

end MarkedProducts

/-! ### Matrix-product states on an open chain -/

namespace OpenMPS

variable {d : Type} [Fintype d] [DecidableEq d] {D : Type*} [Fintype D] [DecidableEq D]

/-- The open matrix-product state `ψ(σ) = bᵀ A_{σ₀} ⋯ A_{σ_{N-1}} c`. -/
def state (b c : D → ℤ) (A : d → Matrix D D ℤ) {N : ℕ} (σ : Fin N → d) : ℤ :=
  b ⬝ᵥ ((List.ofFn fun j ↦ A (σ j)).prod *ᵥ c)

/-- The doubled vector `b ⊗ b`. -/
def pair {R : Type*} [Mul R] (b : D → R) : D × D → R := fun p ↦ b p.1 * b p.2

/-- The one-site insertion `G = ∑ f(s, t) Aₛ ⊗ Aₜ` of a one-site operator `f`. -/
def siteInsertion (f : Matrix d d ℤ) (A : d → Matrix D D ℤ) : Matrix (D × D) (D × D) ℤ :=
  ∑ s, ∑ t, f s t • (A s ⊗ₖ A t)

/-- The norm `L₀ Tᴺ u` of the state on `N` sites. -/
def norm (b c : D → ℤ) (A : d → Matrix D D ℤ) (N : ℕ) : ℤ :=
  pair b ⬝ᵥ (MPS.transfer A ^ N *ᵥ pair c)

/-- The bond energy `∑ⱼ L₀ Tʲ O T^(n-1-j) u` on `n + 1` sites (`n` bonds). -/
def bondEnergy (h : Matrix (d × d) (d × d) ℤ) (b c : D → ℤ) (A : d → Matrix D D ℤ) (n : ℕ) :
    ℤ :=
  ∑ j : Fin n, pair b ⬝ᵥ
    ((MPS.transfer A ^ (j : ℕ) * MPS.insertion h A * MPS.transfer A ^ (n - 1 - j)) *ᵥ pair c)

/-- The first-site term `L₀ G Tⁿ u` on `n + 1` sites. -/
def leftField (f : Matrix d d ℤ) (b c : D → ℤ) (A : d → Matrix D D ℤ) (n : ℕ) : ℤ :=
  pair b ⬝ᵥ ((siteInsertion f A * MPS.transfer A ^ n) *ᵥ pair c)

/-- The last-site term `L₀ Tⁿ G u` on `n + 1` sites. -/
def rightField (f : Matrix d d ℤ) (b c : D → ℤ) (A : d → Matrix D D ℤ) (n : ℕ) : ℤ :=
  pair b ⬝ᵥ ((MPS.transfer A ^ n * siteInsertion f A) *ᵥ pair c)

/-! ### Doubled vectors -/

section Pair

variable {R : Type*} [CommSemiring R]

omit [Fintype d] [DecidableEq d] [DecidableEq D] in
/-- `(b ⊗ b)ᵀ (X ⊗ Y) (c ⊗ c) = (bᵀ X c) (bᵀ Y c)`. -/
theorem dotProduct_kronecker_mulVec (b c : D → R) (X Y : Matrix D D R) :
    pair b ⬝ᵥ ((X ⊗ₖ Y) *ᵥ pair c) = (b ⬝ᵥ (X *ᵥ c)) * (b ⬝ᵥ (Y *ᵥ c)) := by
  simp only [pair, dotProduct, mulVec, kroneckerMap_apply, Fintype.sum_prod_type]
  rw [Fintype.sum_mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ Finset.sum_congr rfl fun i' _ ↦ ?_
  have e : ∀ (u v : R) (f g : D → R),
      (u * ∑ x, f x) * (v * ∑ y, g y) = u * v * ∑ x, ∑ y, f x * g y := by
    intro u v f g
    rw [← Finset.sum_mul_sum]
    ring
  rw [e]
  refine congrArg _ (Finset.sum_congr rfl fun j _ ↦ Finset.sum_congr rfl fun j' _ ↦ ?_)
  ring

omit [Fintype d] [DecidableEq d] [DecidableEq D] in
/-- `(b ⊗ b)ᵀ (X ⊗ Y) = (bᵀ X) ⊗ (bᵀ Y)`. -/
theorem vecMul_pair_kronecker (b : D → R) (X Y : Matrix D D R) :
    pair b ᵥ* (X ⊗ₖ Y) = fun p ↦ (b ᵥ* X) p.1 * (b ᵥ* Y) p.2 := by
  ext ⟨i, j⟩
  simp only [pair, vecMul, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type]
  rw [Fintype.sum_mul_sum]
  refine Finset.sum_congr rfl fun k _ ↦ Finset.sum_congr rfl fun l _ ↦ ?_
  ring

omit [Fintype d] [DecidableEq d] [DecidableEq D] in
/-- `(X ⊗ Y) (c ⊗ c) = (X c) ⊗ (Y c)`. -/
theorem kronecker_mulVec_pair (c : D → R) (X Y : Matrix D D R) :
    (X ⊗ₖ Y) *ᵥ pair c = fun p ↦ (X *ᵥ c) p.1 * (Y *ᵥ c) p.2 := by
  ext ⟨i, j⟩
  simp only [pair, mulVec, dotProduct, kroneckerMap_apply, Fintype.sum_prod_type]
  rw [Fintype.sum_mul_sum]
  refine Finset.sum_congr rfl fun k _ ↦ Finset.sum_congr rfl fun l _ ↦ ?_
  ring

omit [DecidableEq D] in
/-- Casting to `ℂ` commutes with the matrix element `(b ⊗ b)ᵀ M (c ⊗ c)`. -/
theorem pair_dotProduct_map_intCast (b c : D → ℤ) (M : Matrix (D × D) (D × D) ℤ) :
    pair (fun i ↦ (b i : ℂ)) ⬝ᵥ (M.map (↑) *ᵥ pair fun i ↦ (c i : ℂ)) =
      ((pair b ⬝ᵥ (M *ᵥ pair c) : ℤ) : ℂ) := by
  simp only [pair, dotProduct, mulVec, map_apply]
  push_cast
  rfl

end Pair

variable (b c : D → ℤ) (A : d → Matrix D D ℤ)

/-- The state as a complex vector on configurations. -/
def vec (N : ℕ) : (Fin N → d) → ℂ := fun σ ↦ (state b c A σ : ℂ)

omit [Fintype d] [DecidableEq d] in
/-- The state is the matrix element of the ordered product of the complex tensors. -/
theorem vec_apply {N : ℕ} (σ : Fin N → d) :
    vec b c A N σ = (fun i ↦ (b i : ℂ)) ⬝ᵥ
      ((List.ofFn fun j ↦ MPS.tensorC A (σ j)).prod *ᵥ fun i ↦ (c i : ℂ)) := by
  have h := map_list_prod_ofFn_intCast (fun j ↦ A (σ j))
  simp only [MPS.tensorC]
  rw [← h, vec, state]
  simp only [dotProduct, mulVec, map_apply]
  push_cast
  rfl

omit [Fintype d] [DecidableEq d] in
/-- The state is real. -/
theorem star_vec (N : ℕ) : star (vec b c A N) = vec b c A N := by
  ext σ
  simp [vec]

omit [DecidableEq d] in
/-- **Sandwich of a product operator.** For a product operator `⨂ⱼ Mⱼ`, the matrix element in
the state is `L₀ (∏ⱼ Eⱼ) u` with the transfer matrices `Eⱼ = ∑ Mⱼ(s, t) Aₛ ⊗ Aₜ`. -/
theorem dotProduct_piTensor_mulVec {N : ℕ} (M : Fin N → Matrix d d ℂ) :
    star (vec b c A N) ⬝ᵥ (piTensor M *ᵥ vec b c A N) =
      pair (fun i ↦ (b i : ℂ)) ⬝ᵥ ((List.ofFn fun j ↦
        ∑ s, ∑ t, M j s t • (MPS.tensorC A s ⊗ₖ MPS.tensorC A t)).prod *ᵥ
          pair fun i ↦ (c i : ℂ)) := by
  have hsum : ∀ j, ∑ s, ∑ t, M j s t • (MPS.tensorC A s ⊗ₖ MPS.tensorC A t) =
      ∑ x : d × d, M j x.1 x.2 • (MPS.tensorC A x.1 ⊗ₖ MPS.tensorC A x.2) := fun j ↦
    (Fintype.sum_prod_type' (fun s t ↦ M j s t • (MPS.tensorC A s ⊗ₖ MPS.tensorC A t))).symm
  have hL : star (vec b c A N) ⬝ᵥ (piTensor M *ᵥ vec b c A N) =
      ∑ σ : Fin N → d, ∑ τ : Fin N → d,
        (∏ i, M i (σ i) (τ i)) * (vec b c A N σ * vec b c A N τ) := by
    rw [star_vec]
    simp only [dotProduct, mulVec, piTensor_apply, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ ↦ Finset.sum_congr rfl fun τ _ ↦ ?_
    ring
  simp_rw [hsum]
  rw [hL, list_prod_ofFn_sum, sum_mulVec, dotProduct_sum,
    ← (Equiv.arrowProdEquivProdArrow (Fin N) (fun _ ↦ d) (fun _ ↦ d)).symm.sum_comp,
    Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun σ _ ↦ Finset.sum_congr rfl fun τ _ ↦ ?_
  simp only [Equiv.arrowProdEquivProdArrow, Equiv.coe_fn_symm_mk]
  rw [list_prod_ofFn_smul, list_prod_ofFn_kronecker, smul_mulVec, dotProduct_smul, smul_eq_mul,
    dotProduct_kronecker_mulVec, ← vec_apply, ← vec_apply]

/-! ### The norm -/

omit [Fintype d] [DecidableEq d] in
/-- `ψ(σ)²` is a matrix element of the doubled product. -/
theorem state_sq {N : ℕ} (σ : Fin N → d) :
    state b c A σ ^ 2 =
      pair b ⬝ᵥ ((List.ofFn fun j ↦ A (σ j) ⊗ₖ A (σ j)).prod *ᵥ pair c) := by
  rw [list_prod_ofFn_kronecker, dotProduct_kronecker_mulVec, state, sq]

omit [DecidableEq d] in
/-- **Norm of the state.** `∑_σ ψ(σ)² = L₀ Tᴺ u`. -/
theorem sum_state_sq (N : ℕ) : ∑ σ : Fin N → d, state b c A σ ^ 2 = norm b c A N := by
  simp_rw [state_sq, ← dotProduct_sum, ← sum_mulVec,
    ← list_prod_ofFn_sum N (fun _ s ↦ A s ⊗ₖ A s), List.ofFn_const, List.prod_replicate]
  rfl

omit [DecidableEq d] in
/-- The norm `⟨ψ, ψ⟩ = L₀ Tᴺ u` of the state as a complex vector. -/
theorem dotProduct_vec (N : ℕ) : star (vec b c A N) ⬝ᵥ vec b c A N = ((norm b c A N : ℤ) : ℂ) := by
  rw [star_vec, ← sum_state_sq]
  simp [dotProduct, vec, sq]

/-! ### Bonds and end sites -/

/-- The complex one-site insertion. -/
def siteInsertionC (f : Matrix d d ℂ) : Matrix (D × D) (D × D) ℂ :=
  ∑ s, ∑ t, f s t • (MPS.tensorC A s ⊗ₖ MPS.tensorC A t)

omit [DecidableEq d] [Fintype D] [DecidableEq D] in
/-- The insertion of `μ f` for an integer `f` is `μ` times the cast of the integer insertion. -/
theorem siteInsertionC_smul (μ : ℂ) (f : Matrix d d ℤ) :
    siteInsertionC A (μ • f.map (Int.cast : ℤ → ℂ)) =
      μ • (siteInsertion f A).map (Int.cast : ℤ → ℂ) := by
  ext x y
  simp only [siteInsertionC, siteInsertion, MPS.tensorC, Matrix.sum_apply, Matrix.smul_apply,
    smul_eq_mul, map_apply, kroneckerMap_apply]
  push_cast
  simp only [Finset.mul_sum]
  refine Finset.sum_congr rfl fun s _ ↦ Finset.sum_congr rfl fun t _ ↦ ?_
  ring

/-- **One end site.** The term `f` at site `j` of `N` sites contributes `L₀ Tʲ G T^(N-1-j) u`. -/
theorem dotProduct_siteOp_mulVec {N : ℕ} (j : Fin N) (f : Matrix d d ℂ) :
    star (vec b c A N) ⬝ᵥ (siteOp j f *ᵥ vec b c A N) =
      pair (fun i ↦ (b i : ℂ)) ⬝ᵥ ((MPS.transferC A ^ (j : ℕ) * siteInsertionC A f *
        MPS.transferC A ^ (N - 1 - j)) *ᵥ pair fun i ↦ (c i : ℂ)) := by
  have h : (fun i ↦ ∑ s, ∑ t, (if i = j then f else 1) s t •
      (MPS.tensorC A s ⊗ₖ MPS.tensorC A t)) =
      fun i ↦ if i = j then siteInsertionC A f else MPS.transferC A := by
    funext i
    split_ifs
    · rfl
    · exact MPS.sum_one_smul_kronecker A
  rw [siteOp, dotProduct_piTensor_mulVec, h, list_prod_ofFn_site]

/-- **One bond.** The bond `j, j + 1` of `n + 1` sites contributes `L₀ Tʲ O T^(n-1-j) u`. -/
theorem dotProduct_bondOp_mulVec (h : Matrix (d × d) (d × d) ℤ) {n : ℕ} (j : Fin n) :
    star (vec b c A (n + 1)) ⬝ᵥ (bondOp j.castSucc j.succ (h.map (↑)) *ᵥ vec b c A (n + 1)) =
      pair (fun i ↦ (b i : ℂ)) ⬝ᵥ ((MPS.transferC A ^ (j : ℕ) * MPS.insertionC A h *
        MPS.transferC A ^ (n - 1 - j)) *ᵥ pair fun i ↦ (c i : ℂ)) := by
  have hjk : j.castSucc ≠ j.succ := (show j.castSucc < j.succ from Fin.castSucc_lt_succ).ne
  simp only [map_intCast_eq_sum_single_kronecker h, bondOp_sum, bondOp_smul,
    bondOp_kronecker hjk, sum_mulVec, smul_mulVec, dotProduct_sum, dotProduct_smul,
    dotProduct_piTensor_mulVec, smul_eq_mul]
  have hfac : ∀ p q : d × d, (List.ofFn fun i ↦ ∑ s, ∑ t,
      (if i = j.castSucc then single p.1 q.1 (1 : ℂ) else if i = j.succ then single p.2 q.2 1
        else 1) s t • (MPS.tensorC A s ⊗ₖ MPS.tensorC A t)) =
      List.ofFn (bondFactors j.castSucc j.succ (MPS.tensorC A p.1 ⊗ₖ MPS.tensorC A q.1)
        (MPS.tensorC A p.2 ⊗ₖ MPS.tensorC A q.2) (MPS.transferC A)) := by
    intro p q
    congr 1
    ext1 i
    unfold bondFactors
    split_ifs
    · exact MPS.sum_single_smul_kronecker A _ _
    · exact MPS.sum_single_smul_kronecker A _ _
    · exact MPS.sum_one_smul_kronecker A
  simp only [hfac, list_prod_ofFn_bond, MPS.insertionC, Finset.mul_sum, Finset.sum_mul,
    sum_mulVec, dotProduct_sum, mul_smul_comm, smul_mul_assoc, smul_mulVec, dotProduct_smul,
    smul_eq_mul]
  refine Finset.sum_congr rfl fun p _ ↦ Finset.sum_congr rfl fun q _ ↦ ?_
  rw [mul_kronecker_mul, Matrix.mul_assoc, Matrix.mul_assoc]

/-! ### The energy and the variational bound -/

omit [DecidableEq d] in
/-- Powers of the complex transfer matrix are casts of the integer powers. -/
theorem transferC_pow (k : ℕ) : MPS.transferC A ^ k = (MPS.transfer A ^ k).map (↑) := by
  rw [MPS.transferC_eq]
  change (Int.castRingHom ℂ).mapMatrix _ ^ k = (Int.castRingHom ℂ).mapMatrix _
  rw [map_pow]

/-- **Energy of the state.** On `n + 1` sites, with end terms `μ₀ f₀` and `μ₁ f₁`,
`⟨ψ, H ψ⟩ = ∑ⱼ L₀ Tʲ O T^(n-1-j) u + μ₀ L₀ G₀ Tⁿ u + μ₁ L₀ Tⁿ G₁ u`. -/
theorem dotProduct_openChainHamiltonian_mulVec (h : Matrix (d × d) (d × d) ℤ)
    (f₀ f₁ : Matrix d d ℤ) (μ₀ μ₁ : ℂ) (n : ℕ) :
    star (vec b c A (n + 1)) ⬝ᵥ (openChainHamiltonian (h.map (↑)) (μ₀ • f₀.map (↑))
        (μ₁ • f₁.map (↑)) (n + 1) *ᵥ vec b c A (n + 1)) =
      ((bondEnergy h b c A n : ℤ) : ℂ) + μ₀ * ((leftField f₀ b c A n : ℤ) : ℂ) +
        μ₁ * ((rightField f₁ b c A n : ℤ) : ℂ) := by
  have hcast (M : Matrix (D × D) (D × D) ℤ) :=
    pair_dotProduct_map_intCast b c M
  rw [openChainHamiltonian, add_mulVec, add_mulVec, dotProduct_add, dotProduct_add, sum_mulVec,
    dotProduct_sum]
  refine congrArg₂ (· + ·) (congrArg₂ (· + ·) ?_ ?_) ?_
  · rw [bondEnergy]
    push_cast
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [dotProduct_bondOp_mulVec, transferC_pow, transferC_pow, MPS.insertionC_eq,
      ← hcast]
    congr 2
    change (Int.castRingHom ℂ).mapMatrix _ * (Int.castRingHom ℂ).mapMatrix _ *
      (Int.castRingHom ℂ).mapMatrix _ = (Int.castRingHom ℂ).mapMatrix _
    rw [map_mul, map_mul]
  · rw [dotProduct_siteOp_mulVec, siteInsertionC_smul, leftField, ← hcast]
    simp only [Fin.val_zero, pow_zero, Matrix.one_mul, Nat.add_sub_cancel, Nat.sub_zero,
      transferC_pow, Matrix.smul_mul, smul_mulVec, dotProduct_smul, smul_eq_mul]
    congr 3
    change (Int.castRingHom ℂ).mapMatrix _ * (Int.castRingHom ℂ).mapMatrix _ =
      (Int.castRingHom ℂ).mapMatrix _
    rw [map_mul]
  · rw [dotProduct_siteOp_mulVec, siteInsertionC_smul, rightField, ← hcast]
    simp only [Fin.val_last, Nat.add_sub_cancel, Nat.sub_self, pow_zero, Matrix.mul_one,
      transferC_pow, Matrix.mul_smul, smul_mulVec, dotProduct_smul, smul_eq_mul]
    congr 3
    change (Int.castRingHom ℂ).mapMatrix _ * (Int.castRingHom ℂ).mapMatrix _ =
      (Int.castRingHom ℂ).mapMatrix _
    rw [map_mul]

omit [Fintype d] [DecidableEq d] in
/-- A real multiple of a symmetric integer matrix is Hermitian over `ℂ`. -/
theorem isHermitian_smul_map_intCast {f : Matrix d d ℤ} (hf : fᵀ = f) (μ : ℝ) :
    ((μ : ℂ) • f.map (Int.cast : ℤ → ℂ)).IsHermitian := by
  change ((μ : ℂ) • f.map (Int.cast : ℤ → ℂ))ᴴ = _
  rw [conjTranspose_smul, (isHermitian_map_intCast hf).eq, Complex.star_def,
    Complex.conj_ofReal]

/-- **Trial-state bound on the open chain.** For integer data `b, c, A`, a symmetric integer
bond `h` and symmetric integer end terms `f₀, f₁` with real coefficients `μ₀, μ₁`: if the state
on `n + 1` sites has positive norm and Rayleigh quotient below `E`, i.e.
`∑ⱼ L₀ Tʲ O T^(n-1-j) u + μ₀ L₀ G₀ Tⁿ u + μ₁ L₀ Tⁿ G₁ u < E L₀ T^(n+1) u`, then some
eigenvalue of the open chain lies below `E`. -/
theorem exists_eigenvalue_lt (h : Matrix (d × d) (d × d) ℤ) (hh : hᵀ = h)
    {f₀ f₁ : Matrix d d ℤ} (hf₀ : f₀ᵀ = f₀) (hf₁ : f₁ᵀ = f₁) (μ₀ μ₁ : ℝ) (n : ℕ) {E : ℝ}
    (hE : (bondEnergy h b c A n : ℝ) + μ₀ * leftField f₀ b c A n + μ₁ * rightField f₁ b c A n <
      E * norm b c A (n + 1)) :
    ∃ i, (openChainHamiltonian_isHermitian (isHermitian_map_intCast hh)
      (isHermitian_smul_map_intCast hf₀ μ₀) (isHermitian_smul_map_intCast hf₁ μ₁)
        (n + 1)).eigenvalues i < E := by
  refine exists_eigenvalue_lt_of_re_dotProduct_lt _ (vec b c A (n + 1)) ?_
  rw [dotProduct_openChainHamiltonian_mulVec, dotProduct_vec]
  simp only [← Complex.ofReal_intCast, ← Complex.ofReal_mul, ← Complex.ofReal_add,
    Complex.ofReal_re]
  exact hE

end OpenMPS

end CertifiedGaps
