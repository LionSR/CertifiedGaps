/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.SingleSite

/-!
# On-site symmetries of the single-site transfer matrix

Let `Φ` be a trace-preserving algebra endomorphism of the site matrices, for example the
conjugation `A ↦ U A Uᵀ` by an orthogonal matrix `U` (`trace_conjStarAlgAut`) or a relabeling
of the site basis (`Matrix.reindexAlgEquiv`). If `Φ` fixes the on-site term `o`, it acts on the
transfer matrix `X` of `Transfer/SingleSite.lean` through its action on the bond factors:

* **Sign symmetries.** If `Φ (C_α) = s_α C_α` with signs `s_α = ±1`, then `Φ` multiplies the
  site factor `M(a, b) (1 - τ o)` by `s(a) s(b)`, so the site trace of two time-label sequences
  picks up `D_s(γ) D_s(γ')`, with `D_s(γ) = ∏ₜ s(γ t)`. Since the trace is preserved,
  `X(γ, γ') = 0` whenever `D_s(γ) ≠ D_s(γ')`: `X` is block diagonal by the character `D_s`.
* **Label permutations.** If `Φ (C_α) = C_{π α}`, then relabeling every time label by `π`
  preserves `X`.

Combined with `CertifiedGaps.trace_diagonal_mul_pow_eq_sum_blockOf` and
`CertifiedGaps.trace_blockOf_pow_of_reindex`, these split the ring traces into sector traces and
identify the sectors that a label permutation exchanges.

## Main results

* `transfer_eq_char_mul`, `transfer_eq_zero_of_labelChar_ne`: sign symmetries.
* `transfer_map_labels`, `labelChar_map_labels`: label permutations.
* `trace_conjStarAlgAut`, `trace_conjAlgHom`: conjugation by an orthogonal matrix preserves the
  trace; `conjAlgHom_signUnitary_apply`: conjugation by a diagonal sign matrix, entrywise.
* `labelChar_mul`: characters are multiplicative in the sign vector.
* `SignSymmetry`, `LabelSymmetry`: the data a model supplies for each kind of symmetry, with
  `SignSymmetry.transfer_eq_zero` and `LabelSymmetry.transfer_map`.
-/

open Matrix

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d] (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

/-- Conjugation `A ↦ U A Uᵀ` by an orthogonal matrix preserves the trace. -/
theorem trace_conjStarAlgAut (U : unitary (Matrix d d ℝ)) (A : Matrix d d ℝ) :
    trace (Unitary.conjStarAlgAut ℝ (Matrix d d ℝ) U A) = trace A := by
  rw [Unitary.conjStarAlgAut_apply, trace_mul_cycle, Unitary.coe_star_mul_self, Matrix.one_mul]

/-- Conjugation `A ↦ U A Uᵀ` by an orthogonal matrix, as an algebra endomorphism (the form the
symmetry lemmas below take). -/
noncomputable def conjAlgHom (U : unitary (Matrix d d ℝ)) : Matrix d d ℝ →ₐ[ℝ] Matrix d d ℝ :=
  AlgHom.ofClass (Unitary.conjStarAlgAut ℝ _ U)

/-- Conjugation by an orthogonal matrix preserves the trace. -/
theorem trace_conjAlgHom (U : unitary (Matrix d d ℝ)) (A : Matrix d d ℝ) :
    trace (conjAlgHom U A) = trace A :=
  trace_conjStarAlgAut U A

/-- The diagonal matrix of signs `r a = ±1`, as an orthogonal matrix. In a basis where a
symmetry acts diagonally (for example the `π` rotations in a Cartesian basis) this is the
symmetry itself. -/
def signUnitary (r : d → ℝ) (hr : ∀ a, r a * r a = 1) : unitary (Matrix d d ℝ) :=
  ⟨diagonal r, by
    rw [Unitary.mem_iff, star_eq_conjTranspose, diagonal_conjTranspose, star_trivial,
      diagonal_mul_diagonal, ← diagonal_one]
    exact ⟨congrArg diagonal (funext hr), congrArg diagonal (funext hr)⟩⟩

/-- Conjugation by a diagonal sign matrix multiplies the entry `(a, b)` by `r a r b`, so it
reduces a sign symmetry `Φ (C_α) = s_α C_α` to an entrywise check. -/
theorem conjAlgHom_signUnitary_apply (r : d → ℝ) (hr : ∀ a, r a * r a = 1) (A : Matrix d d ℝ) :
    conjAlgHom (signUnitary r hr) A = of fun a b ↦ r a * r b * A a b := by
  change diagonal r * A * star (diagonal r) = _
  rw [star_eq_conjTranspose, diagonal_conjTranspose, star_trivial]
  ext a b
  simp only [mul_diagonal, diagonal_mul, of_apply]
  ring

variable (Φ : Matrix d d ℝ →ₐ[ℝ] Matrix d d ℝ)

/-- A trace-preserving algebra endomorphism can be applied to every factor of a site trace. -/
theorem siteTrace_eq_map (htr : ∀ A, trace (Φ A) = trace A) (τ : ℝ) {m : ℕ}
    (γ γ' : Fin m → Option ι) :
    siteTrace C o τ γ γ' = trace (List.ofFn fun t ↦ Φ (siteFactor C o τ (γ t) (γ' t))).prod := by
  rw [siteTrace, ← htr, map_list_prod, List.map_ofFn]
  rfl

/-- `Φ` fixes `1 - τ o` when it fixes `o`. -/
theorem map_one_sub_smul (ho : Φ o = o) (τ : ℝ) : Φ (1 - τ • o) = 1 - τ • o := by
  rw [map_sub, map_one, map_smul, ho]

/-! ### Sign symmetries -/

section Sign

variable {C o Φ} {s : ι → ℝ}

/-- Under a sign symmetry the factor of a label picks up its sign. -/
theorem map_labelGen (hC : ∀ α, Φ (C α) = s α • C α) (a : Option ι) :
    Φ (labelGen C a) = labelSign s a • labelGen C a := by
  rcases a with _ | α
  · simp [labelGen, labelSign]
  · exact hC α

/-- Under a sign symmetry the site factor picks up the signs of its two labels. -/
theorem map_siteFactor_of_sign (hC : ∀ α, Φ (C α) = s α • C α) (ho : Φ o = o) (τ : ℝ)
    (a b : Option ι) :
    Φ (siteFactor C o τ a b) = (labelSign s a * labelSign s b) • siteFactor C o τ a b := by
  rw [siteFactor, map_mul, map_one_sub_smul _ Φ ho, symFactor, map_smul, map_add, map_mul,
    map_mul, map_labelGen hC, map_labelGen hC, smul_mul_smul_comm, smul_mul_smul_comm,
    mul_comm (labelSign s b), ← smul_add, smul_comm, smul_mul_assoc]

/-- Under a sign symmetry the site trace picks up the characters of both sequences. -/
theorem siteTrace_eq_char_mul (htr : ∀ A, trace (Φ A) = trace A)
    (hC : ∀ α, Φ (C α) = s α • C α) (ho : Φ o = o) (τ : ℝ) {m : ℕ}
    (γ γ' : Fin m → Option ι) :
    siteTrace C o τ γ γ' = (labelChar s γ * labelChar s γ') * siteTrace C o τ γ γ' := by
  conv_lhs => rw [siteTrace_eq_map C o Φ htr]
  simp only [map_siteFactor_of_sign hC ho]
  rw [list_prod_ofFn_smul, trace_smul, smul_eq_mul, siteTrace, labelChar, labelChar,
    ← Finset.prod_mul_distrib]

/-- Under a sign symmetry the transfer matrix picks up the characters of both sequences. -/
theorem transfer_eq_char_mul (htr : ∀ A, trace (Φ A) = trace A)
    (hC : ∀ α, Φ (C α) = s α • C α) (ho : Φ o = o) (m : ℕ) (τ : ℝ)
    (γ γ' : Fin m → Option ι) :
    transfer C o m τ γ γ' = (labelChar s γ * labelChar s γ') * transfer C o m τ γ γ' := by
  simp only [transfer, of_apply]
  conv_lhs => rw [siteTrace_eq_char_mul htr hC ho]
  ring

omit [Fintype d] [DecidableEq d] in
/-- With signs `±1`, a character is `±1`. -/
theorem labelChar_eq_one_or_neg_one (hs : ∀ α, s α * s α = 1) {m : ℕ}
    (γ : Fin m → Option ι) : labelChar s γ = 1 ∨ labelChar s γ = -1 := by
  refine mul_self_eq_one_iff.1 ?_
  rw [labelChar, ← Finset.prod_mul_distrib]
  refine Finset.prod_eq_one fun t _ ↦ ?_
  rcases γ t with _ | α
  · simp [labelSign]
  · exact hs α

/-- **Block diagonality.** Under a sign symmetry with signs `±1`, the transfer matrix vanishes
between sequences of different characters. -/
theorem transfer_eq_zero_of_labelChar_ne (htr : ∀ A, trace (Φ A) = trace A)
    (hs : ∀ α, s α * s α = 1) (hC : ∀ α, Φ (C α) = s α • C α) (ho : Φ o = o) (m : ℕ) (τ : ℝ)
    {γ γ' : Fin m → Option ι} (h : (labelChar s γ = 1) ≠ (labelChar s γ' = 1)) :
    transfer C o m τ γ γ' = 0 := by
  have hm : labelChar s γ * labelChar s γ' = -1 := by
    rcases labelChar_eq_one_or_neg_one hs γ with e1 | e1 <;>
      rcases labelChar_eq_one_or_neg_one hs γ' with e2 | e2 <;> simp_all
  have := transfer_eq_char_mul htr hC ho m τ γ γ'
  rw [hm] at this
  linarith

end Sign

/-! ### Label permutations -/

section Permutation

variable {C o Φ} {π : ι → ι}

/-- Under a label permutation the factor of a label is the factor of the relabeled label. -/
theorem map_labelGen_of_perm (hC : ∀ α, Φ (C α) = C (π α)) (a : Option ι) :
    Φ (labelGen C a) = labelGen C (a.map π) := by
  rcases a with _ | α
  · simp [labelGen]
  · exact hC α

/-- Under a label permutation the site factor is the site factor of the relabeled labels. -/
theorem map_siteFactor_of_perm (hC : ∀ α, Φ (C α) = C (π α)) (ho : Φ o = o) (τ : ℝ)
    (a b : Option ι) :
    Φ (siteFactor C o τ a b) = siteFactor C o τ (a.map π) (b.map π) := by
  rw [siteFactor, map_mul, map_one_sub_smul _ Φ ho, symFactor, map_smul, map_add, map_mul,
    map_mul, map_labelGen_of_perm hC, map_labelGen_of_perm hC, siteFactor, symFactor]

/-- **Relabeling invariance.** Relabeling every time label by `π` preserves the transfer
matrix. -/
theorem transfer_map_labels (htr : ∀ A, trace (Φ A) = trace A) (hC : ∀ α, Φ (C α) = C (π α))
    (ho : Φ o = o) (m : ℕ) (τ : ℝ) (γ γ' : Fin m → Option ι) :
    transfer C o m τ (fun t ↦ (γ t).map π) (fun t ↦ (γ' t).map π) = transfer C o m τ γ γ' := by
  have hw : ∀ γ : Fin m → Option ι, halfWeight τ (fun t ↦ (γ t).map π) = halfWeight τ γ :=
    fun γ ↦ Finset.prod_congr rfl fun t _ ↦ by
      change ((γ t).map π).elim _ _ = (γ t).elim _ _
      rcases γ t with _ | α <;> rfl
  have hs : siteTrace C o τ (fun t ↦ (γ t).map π) (fun t ↦ (γ' t).map π) =
      siteTrace C o τ γ γ' := by
    rw [siteTrace_eq_map C o Φ htr τ γ γ']
    simp only [map_siteFactor_of_perm hC ho, siteTrace]
  simp only [transfer, of_apply, hw, hs]

omit [Fintype d] [DecidableEq d] in
/-- The character of a relabeled sequence: `D_s(π γ) = D_{s ∘ π}(γ)`. -/
theorem labelChar_map_labels (s : ι → ℝ) (π : ι → ι) {m : ℕ} (γ : Fin m → Option ι) :
    labelChar s (fun t ↦ (γ t).map π) = labelChar (s ∘ π) γ :=
  Finset.prod_congr rfl fun t _ ↦ by
    change labelSign s ((γ t).map π) = labelSign (s ∘ π) (γ t)
    rcases γ t with _ | α <;> rfl

end Permutation

omit [Fintype d] [DecidableEq d] in
/-- The character of a product of sign vectors is the product of the characters. -/
theorem labelChar_mul (s s' : ι → ℝ) {m : ℕ} (γ : Fin m → Option ι) :
    labelChar (fun α ↦ s α * s' α) γ = labelChar s γ * labelChar s' γ := by
  rw [labelChar, labelChar, labelChar, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  rcases γ t with _ | α
  · simp [labelSign]
  · rfl

/-! ### Bundled symmetries

A model supplies its on-site symmetries in these two forms; the sector theorems
(`Transfer/SingleSiteD2.lean`) only use the bundled fields. -/

/-- A **sign symmetry** of the bond factors `C` and the on-site term `o`: a trace-preserving
algebra endomorphism `Φ` (typically `conjAlgHom U` for an orthogonal `U`) with
`Φ (C_α) = s_α C_α`, `s_α = ±1`, and `Φ o = o`. -/
structure SignSymmetry (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ) where
  /-- The sign `s_α` that the symmetry puts on the bond factor `C_α`. -/
  sign : ι → ℝ
  sign_sq : ∀ α, sign α * sign α = 1
  /-- The action of the symmetry on the site matrices. -/
  map : Matrix d d ℝ →ₐ[ℝ] Matrix d d ℝ
  trace_map : ∀ A, trace (map A) = trace A
  map_C : ∀ α, map (C α) = sign α • C α
  map_o : map o = o

/-- A **label symmetry** of `C` and `o`: a trace-preserving algebra endomorphism `Φ` (an
orthogonal conjugation or a relabeling of the site basis) with `Φ (C_α) = C_{π α}` for a
permutation `π` of the labels, and `Φ o = o`. -/
structure LabelSymmetry (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ) where
  /-- The permutation of the bond labels. -/
  perm : Equiv.Perm ι
  /-- The action of the symmetry on the site matrices. -/
  map : Matrix d d ℝ →ₐ[ℝ] Matrix d d ℝ
  trace_map : ∀ A, trace (map A) = trace A
  map_C : ∀ α, map (C α) = C (perm α)
  map_o : map o = o

variable {C o}

/-- The transfer matrix is block diagonal by the character of a sign symmetry. -/
theorem SignSymmetry.transfer_eq_zero (S : SignSymmetry C o)
    (m : ℕ) (τ : ℝ) {γ γ' : Fin m → Option ι}
    (h : (labelChar S.sign γ = 1) ≠ (labelChar S.sign γ' = 1)) : transfer C o m τ γ γ' = 0 :=
  transfer_eq_zero_of_labelChar_ne S.trace_map S.sign_sq S.map_C S.map_o m τ h

/-- The transfer matrix is invariant under relabeling every time label by a label symmetry. -/
theorem LabelSymmetry.transfer_map (S : LabelSymmetry C o) (m : ℕ)
    (τ : ℝ) (γ γ' : Fin m → Option ι) :
    transfer C o m τ (fun t ↦ (γ t).map S.perm) (fun t ↦ (γ' t).map S.perm) =
      transfer C o m τ γ γ' :=
  transfer_map_labels S.trace_map S.map_C S.map_o m τ γ γ'

end CertifiedGaps.SingleSite
