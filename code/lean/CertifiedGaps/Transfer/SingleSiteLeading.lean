/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.TwoBoundary
import CertifiedGaps.Spectral.OperatorBounds
import CertifiedGaps.Transfer.OpenSingleSite
import CertifiedGaps.Transfer.SingleSiteA4

/-!
# The rotation group on the single-site transfer matrix and its leading line

OpenAI's *A boundary-field gap for the spin-one Heisenberg chain* (2026) lets the rotation group
`𝒢 = ⟨D₂, C⟩ ≅ A₄` act on the spatial transfer operator by relabeling the bond events with signs
(eq. (tr:rotations)), averages over it (`Π = |𝒢|⁻¹ ∑_g U_g`, thermal appendix) and uses that `𝒢`
fixes the leading eigenvector (Proposition per:concentration). This file does the same for the
single-site transfer matrix `X` of `Transfer/SingleSite.lean` and an `A4Symmetry`.

**The group.** On time-label sequences a sign symmetry acts by the diagonal matrix of its
character, `charMat r`, and the cyclic relabeling by the permutation matrix `relabelMat π`. Both
commute with `X` (`charMat_mul_comm`, `relabelMat_mul_comm`), are unitary, and move the boundary
vectors of `Transfer/OpenSingleSite.lean` to the boundary vectors of the rotated end fields
(`charMat_mulVec_capVec`, `relabelMat_mulVec_capVec`).

**The average.** `avgProj = Π_D Π_C` with `Π_D = ¼ (1 + D_x)(1 + D_y)` (the projection onto the
trivial `D₂` character `J`, diagonal) and `Π_C = ⅓ (1 + P + P²)`. It is the average of the twelve
group elements, an orthogonal projection, and commutes with `X` (`avgProj_isHermitian`,
`avgProj_mul_self`, `avgProj_mul_comm`).

**The leading line.** Suppose the trivial block `Y_J` of `t X` has an eigenvalue `λ₀` with every
other eigenvalue of `t Y_J`, and every eigenvalue of the block `N`, at most `r < |λ₀|` in modulus.
The three nontrivial `D₂` blocks are conjugate under the cyclic relabeling, so all of `t X` off the
eigenvector `ψ` of `λ₀` has quadratic form at most `r` (`isLeadingLine_of_blocks`). The eigenvector
lies in `J`, so `D₂` fixes it; the relabeling maps it to an eigenvector of the same eigenvalue,
hence to `a ψ` with `a` real (everything is real) and `a³ = 1`, so `a = 1`. Thus `Π ψ = ψ`
(`avgProj_mulVec_of_fixed`): the paper's leading-line invariance, at the level of `X`.

## Main results

* `charMat`, `relabelMat`, `avgProj` and their algebra.
* `isLeadingLine_of_blocks`: the leading line of `t X` from the eigenvalues of two blocks, with
  `ψ` in the trivial `D₂` character and fixed by the relabeling.
* `avgProj_mulVec_of_fixed`: such a `ψ` is fixed by `Π`.
-/

open Matrix Finset

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d]

/-! ### Label actions on sequences -/

section Actions

variable {m : ℕ} [Fintype ι] [DecidableEq ι]

/-- The diagonal matrix of the character `D_r(γ) = ∏ₜ r(γ t)`. -/
noncomputable def charMat (r : ι → ℝ) : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  diagonal fun γ ↦ (labelChar r γ : ℂ)

/-- The relabeling of a sequence: every time label is relabeled by `π`. -/
def relabelSeq (π : ι → ι) (γ : Fin m → Option ι) : Fin m → Option ι := fun t ↦ (γ t).map π

/-- The permutation matrix of a relabeling, `(P f)(γ) = f(π γ)`. -/
noncomputable def relabelMat (π : ι → ι) :
    Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  of fun γ γ' ↦ if γ' = relabelSeq π γ then 1 else 0

omit [DecidableEq d] [Fintype d] in
/-- A character acts on a vector by multiplying each entry by its sign. -/
theorem charMat_mulVec (r : ι → ℝ) (f : (Fin m → Option ι) → ℂ) :
    charMat r *ᵥ f = fun γ ↦ (labelChar r γ : ℂ) * f γ := by
  ext γ; simp [charMat, mulVec_diagonal]

omit [DecidableEq d] [Fintype d] in
/-- The relabeling matrix acts on a vector by relabeling its argument. -/
theorem relabelMat_mulVec (π : ι → ι) (f : (Fin m → Option ι) → ℂ) :
    relabelMat π *ᵥ f = fun γ ↦ f (relabelSeq π γ) := by
  ext γ
  simp [relabelMat, mulVec, dotProduct]

omit [DecidableEq d] [Fintype d] in
/-- A character moves a boundary vector to the boundary vector of the signed weights. -/
theorem charMat_mulVec_capVec (r : ι → ℝ) (τ : ℝ) (u : ι → ℂ) :
    charMat r *ᵥ capVec (m := m) τ u = capVec τ fun α ↦ (r α : ℂ) * u α := by
  rw [charMat_mulVec]
  ext γ
  simp only [labelChar, capVec, Complex.ofReal_prod, ← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  rcases γ t with _ | α
  · simp [labelSign, capWeight]
  · simp only [labelSign, capWeight, Option.elim]; ring

omit [DecidableEq d] [Fintype d] in
/-- A relabeling moves a boundary vector to the boundary vector of the relabeled weights. -/
theorem relabelMat_mulVec_capVec (π : ι → ι) (τ : ℝ) (u : ι → ℂ) :
    relabelMat π *ᵥ capVec (m := m) τ u = capVec τ (u ∘ π) := by
  rw [relabelMat_mulVec]
  ext γ
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  simp only [relabelSeq]
  rcases γ t with _ | α <;> rfl

omit [DecidableEq d] [Fintype d] in
/-- A character with signs `±1` squares to one. -/
theorem charMat_mul_self {r : ι → ℝ} (hr : ∀ α, r α * r α = 1) :
    charMat (m := m) r * charMat r = 1 := by
  rw [charMat, diagonal_mul_diagonal, ← diagonal_one]
  congr 1; ext γ
  have h := mul_self_eq_one_iff.2 (labelChar_eq_one_or_neg_one hr γ)
  rw [← Complex.ofReal_mul, h, Complex.ofReal_one]

omit [DecidableEq d] [Fintype d] [Fintype ι] in
/-- A character matrix is Hermitian. -/
theorem charMat_conjTranspose (r : ι → ℝ) : (charMat (m := m) r)ᴴ = charMat r := by
  rw [charMat, diagonal_conjTranspose]
  congr 1; ext γ; simp [Complex.conj_ofReal]

omit [DecidableEq d] [Fintype d] in
/-- A real matrix whose entries pick up the character of both indices commutes with the
character matrix. -/
theorem charMat_mul_comm {r : ι → ℝ} (hr : ∀ α, r α * r α = 1)
    {T : Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ}
    (hT : ∀ γ γ', T γ γ' = (labelChar r γ * labelChar r γ') * T γ γ') :
    charMat r * T.map (↑) = T.map (↑) * charMat r := by
  ext γ γ'
  simp only [charMat, diagonal_mul, mul_diagonal, map_apply]
  have h := mul_self_eq_one_iff.2 (labelChar_eq_one_or_neg_one hr γ)
  have e : labelChar r γ * T γ γ' = T γ γ' * labelChar r γ' := by
    conv_lhs => rw [hT]
    linear_combination (labelChar r γ' * T γ γ') * h
  exact_mod_cast e

omit [DecidableEq d] [Fintype d] [Fintype ι] [DecidableEq ι] in
/-- Relabeling by `π` and by `π⁻¹` are inverse. -/
theorem relabelSeq_symm (π : Equiv.Perm ι) (γ : Fin m → Option ι) :
    relabelSeq π (relabelSeq (⇑π.symm) γ) = γ ∧ relabelSeq (⇑π.symm) (relabelSeq π γ) = γ := by
  constructor <;> (funext t; simp only [relabelSeq, Option.map_map]; rcases γ t with _ | α <;> simp)

omit [DecidableEq d] [Fintype d] in
/-- The relabeling matrix of a permutation is unitary. -/
theorem relabelMat_conjTranspose_mul_self (π : Equiv.Perm ι) :
    (relabelMat (m := m) π)ᴴ * relabelMat π = 1 := by
  ext γ γ'
  simp only [mul_apply, conjTranspose_apply, relabelMat, of_apply, one_apply]
  rw [Finset.sum_eq_single (relabelSeq (⇑π.symm) γ)]
  · rw [(relabelSeq_symm π γ).1]
    by_cases h : γ = γ' <;> simp [h, eq_comm]
  · intro δ _ hδ
    have : γ ≠ relabelSeq π δ := fun h ↦ hδ (by rw [h, (relabelSeq_symm π δ).2])
    simp [this]
  · simp

omit [DecidableEq d] [Fintype d] in
/-- A real matrix invariant under relabeling both indices by a permutation commutes with the
relabeling matrix. -/
theorem relabelMat_mul_comm (π : Equiv.Perm ι)
    {T : Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ}
    (hT : ∀ γ γ', T (relabelSeq π γ) (relabelSeq π γ') = T γ γ') :
    relabelMat π * T.map (↑) = T.map (↑) * relabelMat π := by
  ext γ γ'
  simp only [mul_apply, relabelMat, of_apply, map_apply, ite_mul, one_mul, zero_mul,
    mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rw [Finset.sum_eq_single (relabelSeq (⇑π.symm) γ')]
  · rw [ite_eq_left (relabelSeq_symm π γ').1.symm, ← hT γ (relabelSeq (⇑π.symm) γ'),
      (relabelSeq_symm π γ').1]
  · intro δ _ hδ
    have : γ' ≠ relabelSeq π δ := fun h ↦ hδ (by rw [h, (relabelSeq_symm π δ).2])
    simp [this]
  · simp

omit [DecidableEq d] [Fintype d] in
/-- Relabeling twice is relabeling by the composite. -/
theorem relabelMat_mul (π π' : ι → ι) :
    relabelMat (m := m) π * relabelMat π' = relabelMat (π' ∘ π) := by
  ext γ γ'
  simp only [mul_apply, relabelMat, of_apply, ite_mul, one_mul, zero_mul, Finset.sum_ite_eq',
    Finset.mem_univ, ite_true]
  congr 1
  apply propext
  have : relabelSeq π' (relabelSeq π γ) = relabelSeq (π' ∘ π) γ := by
    funext t; simp only [relabelSeq, Option.map_map]
  rw [this]

omit [DecidableEq d] [Fintype d] [Fintype ι] in
/-- Relabeling by the identity is the identity. -/
theorem relabelMat_id : relabelMat (m := m) (id : ι → ι) = 1 := by
  ext γ γ'
  simp only [relabelMat, of_apply, one_apply]
  have : relabelSeq (id : ι → ι) γ = γ := by funext t; simp [relabelSeq]
  rw [this]
  by_cases h : γ = γ' <;> simp [h, eq_comm]

end Actions

/-! ### The average over the rotation group -/

section Average

variable [Fintype ι] [DecidableEq ι] {C : ι → Matrix d d ℝ} {o : Matrix d d ℝ}
  (S : A4Symmetry C o) (m : ℕ)

/-- The weight `¼ (1 + D_x)(1 + D_y)` of the trivial `D₂` character `J`: `1` on `J`, `0`
elsewhere. -/
noncomputable def trivWeight (γ : Fin m → Option ι) : ℝ :=
  (1 + labelChar S.x.sign γ) * (1 + labelChar S.y.sign γ) / 4

/-- The projection `Π_D = ¼ (1 + D_x)(1 + D_y)` onto the trivial `D₂` character. -/
noncomputable def projD : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  diagonal fun γ ↦ (trivWeight S m γ : ℂ)

/-- The cyclic relabeling as a matrix. -/
noncomputable abbrev cycleMat : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  relabelMat S.cycle.perm

/-- The average `Π_C = ⅓ (1 + P + P²)` over the cyclic relabeling. -/
noncomputable def projC : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  (1 / 3 : ℂ) • (1 + cycleMat S m + cycleMat S m * cycleMat S m)

/-- **The average over the rotation group** `Π = Π_D Π_C`. -/
noncomputable def avgProj : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ :=
  projD S m * projC S m

omit [Fintype ι] [DecidableEq ι] in
/-- The weight of the trivial character is the indicator of `J`. -/
theorem trivWeight_eq (γ : Fin m → Option ι) :
    trivWeight S m γ = if S.key γ = D2Symmetry.sectorJ then 1 else 0 := by
  simp only [trivWeight, D2Symmetry.key, D2Symmetry.sectorJ]
  rcases labelChar_eq_one_or_neg_one S.x.sign_sq γ with h1 | h1 <;>
    rcases labelChar_eq_one_or_neg_one S.y.sign_sq γ with h2 | h2 <;> norm_num [h1, h2]

omit [Fintype ι] [DecidableEq ι] in
/-- The cyclic relabeling preserves the trivial character. -/
theorem trivWeight_relabel (γ : Fin m → Option ι) :
    trivWeight S m (relabelSeq S.cycle.perm γ) = trivWeight S m γ := by
  have hJ : S.key (S.relabel m γ) = D2Symmetry.sectorJ ↔ S.key γ = D2Symmetry.sectorJ := by
    rw [S.key_relabel]
    constructor
    · intro h; exact d2KeyPerm_injective (h.trans ((d2KeyPerm_eq_self_iff _).2 rfl).symm)
    · intro h; rw [h]; rfl
  rw [trivWeight_eq, trivWeight_eq, show relabelSeq S.cycle.perm γ = S.relabel m γ from rfl]
  simp only [hJ]

/-- The cyclic relabeling matrix has order three. -/
theorem cycleMat_cube : cycleMat S m * cycleMat S m * cycleMat S m = 1 := by
  rw [cycleMat, relabelMat_mul, relabelMat_mul]
  have : ⇑S.cycle.perm ∘ ⇑S.cycle.perm ∘ ⇑S.cycle.perm = id := funext S.cycle_cube
  rw [this, relabelMat_id]

/-- `Pᴴ = P²`. -/
theorem cycleMat_conjTranspose : (cycleMat S m)ᴴ = cycleMat S m * cycleMat S m := by
  have h1 := relabelMat_conjTranspose_mul_self (m := m) S.cycle.perm
  calc (cycleMat S m)ᴴ = (cycleMat S m)ᴴ * (cycleMat S m * cycleMat S m * cycleMat S m) := by
        rw [cycleMat_cube, Matrix.mul_one]
    _ = ((cycleMat S m)ᴴ * cycleMat S m) * (cycleMat S m * cycleMat S m) := by
        simp only [Matrix.mul_assoc]
    _ = _ := by rw [h1, Matrix.one_mul]

omit [Fintype ι] in
/-- `Π_D` is Hermitian. -/
theorem projD_conjTranspose : (projD S m)ᴴ = projD S m := by
  rw [projD, diagonal_conjTranspose]
  congr 1; ext γ; simp [Complex.conj_ofReal]

/-- `Π_D` is idempotent. -/
theorem projD_mul_self : projD S m * projD S m = projD S m := by
  rw [projD, diagonal_mul_diagonal]
  congr 1; ext γ
  rw [trivWeight_eq]; split_ifs <;> simp

/-- `Π_D` commutes with the cyclic relabeling. -/
theorem projD_mul_cycleMat : projD S m * cycleMat S m = cycleMat S m * projD S m := by
  ext γ γ'
  simp only [projD, cycleMat, diagonal_mul, mul_diagonal, relabelMat, of_apply]
  split_ifs with h
  · rw [h, trivWeight_relabel]; ring
  · simp

/-- `Π_C` is Hermitian. -/
theorem projC_conjTranspose : (projC S m)ᴴ = projC S m := by
  have h3 := cycleMat_cube S m
  rw [projC, conjTranspose_smul, conjTranspose_add, conjTranspose_add, conjTranspose_one,
    conjTranspose_mul, cycleMat_conjTranspose]
  have : cycleMat S m * cycleMat S m * (cycleMat S m * cycleMat S m) = cycleMat S m := by
    rw [← Matrix.mul_assoc, h3, Matrix.one_mul]
  rw [this]
  congr 1
  · simp
  · abel

/-- `Π_C` is idempotent. -/
theorem projC_mul_self : projC S m * projC S m = projC S m := by
  have h3 := cycleMat_cube S m
  have key : ∀ P : Matrix (Fin m → Option ι) (Fin m → Option ι) ℂ, P * P * P = 1 →
      (1 + P + P * P) * (1 + P + P * P) = (3 : ℂ) • (1 + P + P * P) := fun P hP ↦ by
    simp only [Matrix.add_mul, Matrix.mul_add, Matrix.one_mul, Matrix.mul_one,
      ← Matrix.mul_assoc, hP]
    rw [show (3 : ℂ) = 1 + 1 + 1 by norm_num, add_smul, add_smul, one_smul]
    abel
  rw [projC, smul_mul_smul_comm, key _ h3, smul_smul]
  congr 1
  norm_num

/-- `Π_D` commutes with `Π_C`. -/
theorem projD_mul_projC : projD S m * projC S m = projC S m * projD S m := by
  have h := projD_mul_cycleMat S m
  have h2 : projD S m * (cycleMat S m * cycleMat S m) =
      cycleMat S m * cycleMat S m * projD S m := by
    rw [← Matrix.mul_assoc, h, Matrix.mul_assoc, h, Matrix.mul_assoc]
  rw [projC, mul_smul_comm, smul_mul_assoc]
  congr 1
  simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, Matrix.one_mul, h, h2]

/-- **`Π` is Hermitian.** -/
theorem avgProj_conjTranspose : (avgProj S m)ᴴ = avgProj S m := by
  rw [avgProj, conjTranspose_mul, projD_conjTranspose, projC_conjTranspose, projD_mul_projC]

/-- **`Π` is idempotent.** -/
theorem avgProj_mul_self : avgProj S m * avgProj S m = avgProj S m := by
  rw [avgProj, Matrix.mul_assoc, ← Matrix.mul_assoc (projC S m), ← projD_mul_projC,
    Matrix.mul_assoc, projC_mul_self, ← Matrix.mul_assoc, projD_mul_self]

/-- **`Π` commutes with every real matrix that is block diagonal in the `D₂` characters and
invariant under the cyclic relabeling**, such as the transfer matrix. -/
theorem avgProj_mul_comm {T : Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ}
    (hkey : ∀ γ γ', S.key γ ≠ S.key γ' → T γ γ' = 0)
    (hrel : ∀ γ γ', T (relabelSeq S.cycle.perm γ) (relabelSeq S.cycle.perm γ') = T γ γ') :
    avgProj S m * T.map (↑) = T.map (↑) * avgProj S m := by
  have hD : projD S m * T.map (↑) = T.map (↑) * projD S m := by
    ext γ γ'
    simp only [projD, diagonal_mul, mul_diagonal, map_apply]
    by_cases h : S.key γ = S.key γ'
    · rw [trivWeight_eq, trivWeight_eq, h]; ring
    · rw [hkey γ γ' h]; simp
  have hP : cycleMat S m * T.map (↑) = T.map (↑) * cycleMat S m := relabelMat_mul_comm _ hrel
  have hC : projC S m * T.map (↑) = T.map (↑) * projC S m := by
    rw [projC, smul_mul_assoc, mul_smul_comm, Matrix.add_mul, Matrix.add_mul, Matrix.mul_add,
      Matrix.mul_add, Matrix.one_mul, Matrix.mul_one, hP, Matrix.mul_assoc, hP,
      ← Matrix.mul_assoc, hP, Matrix.mul_assoc]
  rw [avgProj, Matrix.mul_assoc, hC, ← Matrix.mul_assoc, hD, Matrix.mul_assoc]

/-- A vector in the trivial character fixed by the relabeling is fixed by `Π`. -/
theorem avgProj_mulVec_of_fixed {ψ : (Fin m → Option ι) → ℂ}
    (hJ : ∀ γ, S.key γ ≠ D2Symmetry.sectorJ → ψ γ = 0) (hP : cycleMat S m *ᵥ ψ = ψ) :
    avgProj S m *ᵥ ψ = ψ := by
  have hD : projD S m *ᵥ ψ = ψ := by
    ext γ
    rw [projD, mulVec_diagonal, trivWeight_eq]
    split_ifs with h
    · simp
    · simp [hJ γ h]
  have hC : projC S m *ᵥ ψ = ψ := by
    rw [projC, smul_mulVec, add_mulVec, add_mulVec, one_mulVec, ← mulVec_mulVec, hP, hP]
    ext γ; simp only [Pi.smul_apply, Pi.add_apply, smul_eq_mul]; ring
  rw [avgProj, ← mulVec_mulVec, hC, hD]

/-- `Π_D` on a vector: the average of the four `D₂` characters. -/
theorem projD_mulVec (f : (Fin m → Option ι) → ℂ) :
    projD S m *ᵥ f = (1 / 4 : ℂ) • (f + charMat S.x.sign *ᵥ f + charMat S.y.sign *ᵥ f +
      charMat (fun α ↦ S.x.sign α * S.y.sign α) *ᵥ f) := by
  ext γ
  simp only [projD, mulVec_diagonal, charMat_mulVec, trivWeight, Pi.smul_apply, Pi.add_apply,
    smul_eq_mul, labelChar_mul]
  push_cast
  ring

/-- `Π_C` on a vector: the average of the three powers of the relabeling. -/
theorem projC_mulVec (f : (Fin m → Option ι) → ℂ) :
    projC S m *ᵥ f =
      (1 / 3 : ℂ) • (f + cycleMat S m *ᵥ f + cycleMat S m *ᵥ (cycleMat S m *ᵥ f)) := by
  rw [projC, smul_mulVec, add_mulVec, add_mulVec, one_mulVec, ← mulVec_mulVec]

end Average

/-! ### Quadratic forms of block matrices -/

section BlockForms

open Spectral

variable {D κ : Type*} [Fintype D] [DecidableEq κ]

/-- The restriction of a real vector to the indices of key `c`. -/
def restrictKey (key : D → κ) (c : κ) (y : D → ℝ) : D → ℝ :=
  fun x ↦ if key x = c then y x else 0

/-- The form of a matrix supported on the indices of key `c` sees only the restriction. -/
theorem qform_restrictKey {A : Matrix D D ℝ} {key : D → κ} {c : κ}
    (hA : ∀ x z, A x z ≠ 0 → key x = c ∧ key z = c) (y : D → ℝ) :
    qform A y = qform A (restrictKey key c y) := by
  simp only [qform, dotProduct, mulVec, Finset.mul_sum, RCLike.re_to_real, star_trivial]
  refine Finset.sum_congr rfl fun x _ ↦ Finset.sum_congr rfl fun z _ ↦ ?_
  by_cases hxz : A x z = 0
  · simp [hxz]
  · obtain ⟨hx, hz⟩ := hA x z hxz
    simp [restrictKey, hx, hz]

/-- The squared norm splits over the keys. -/
theorem nsq_eq_sum_restrictKey [Fintype κ] (key : D → κ) (y : D → ℝ) :
    nsq y = ∑ c, nsq (restrictKey key c y) := by
  simp only [nsq_real, restrictKey]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  simp [ite_pow]

variable [DecidableEq D]

/-- Eigenvalues at most `b` in modulus bound the form by `b`. -/
theorem abs_qform_le_of_eigenvalues {A : Matrix D D ℝ} (hA : A.IsHermitian) {b : ℝ}
    (h : ∀ i, |hA.eigenvalues i| ≤ b) (y : D → ℝ) : |qform A y| ≤ b * nsq y := by
  refine abs_le.2 ⟨?_, qform_le_of_eigenvalues_le hA (fun i ↦ (le_abs_self _).trans (h i)) y⟩
  have hlo : ∀ i, -b ≤ hA.eigenvalues i := fun i ↦ by
    have := neg_abs_le (hA.eigenvalues i)
    linarith [h i]
  have := le_qform_of_le_eigenvalues hA hlo y
  linarith

omit [DecidableEq D] in
/-- Reindexing a matrix reindexes its form. -/
theorem qform_reindex (A : Matrix D D ℝ) (σ : D ≃ D) (y : D → ℝ) :
    qform (reindex σ σ A) y = qform A (y ∘ σ) := by
  simp only [qform, dotProduct, mulVec, reindex_apply, submatrix_apply, RCLike.re_to_real,
    star_trivial, Function.comp_apply]
  rw [← σ.sum_comp]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [← σ.sum_comp]
  simp

omit [DecidableEq D] in
/-- Reindexing preserves the squared norm. -/
theorem nsq_comp_equiv (σ : D ≃ D) (y : D → ℝ) : nsq (y ∘ σ) = nsq y := by
  simp only [nsq_real, Function.comp_apply]
  exact σ.sum_comp (fun x ↦ y x ^ 2)

/-- **A Hermitian matrix without one eigenvector.** With `U` the eigenvector unitary of `A` and
`ψ = U e_{j₀}`, the form of `t A - t λ_{j₀} |ψ⟩⟨ψ|` is at most `b` in modulus as soon as every
other `|t λⱼ|` is. -/
theorem abs_qform_sub_proj_le {A : Matrix D D ℝ} (hA : A.IsHermitian) (j₀ : D) (t : ℝ)
    {b : ℝ} (hb : 0 ≤ b) (h : ∀ j, j ≠ j₀ → |t * hA.eigenvalues j| ≤ b) (y : D → ℝ) :
    |qform (t • A - (t * hA.eigenvalues j₀) • vecMulVec (eigU hA *ᵥ Pi.single j₀ 1)
      (eigU hA *ᵥ Pi.single j₀ 1)) y| ≤ b * nsq y := by
  set U := eigU hA
  set ψ := U *ᵥ Pi.single j₀ 1
  set z := Uᴴ *ᵥ y
  have hz : z j₀ = ψ ⬝ᵥ y := by
    simp only [z, ψ, mulVec, dotProduct, conjTranspose_apply, star_trivial,
      Pi.single_apply, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ,
      ite_true]
  have hA' : qform A y = ∑ j, hA.eigenvalues j * z j ^ 2 := by
    rw [← eigU_mulVec_star hA y, qform_eigU]
    simp [z, U, Real.norm_eq_abs, sq_abs]
  have hP : qform (vecMulVec ψ ψ) y = z j₀ ^ 2 := by
    simp only [qform, vecMulVec_mulVec, RCLike.re_to_real, star_trivial, hz]
    rw [dotProduct_smul, MulOpposite.smul_eq_mul_unop]
    simp [dotProduct_comm, sq]
  have hsm : ∀ (c : ℝ) (B : Matrix D D ℝ), qform (c • B) y = c * qform B y := fun c B ↦ by
    simpa using qform_smul_left (𝕜 := ℝ) c B y
  rw [qform_sub_left, hsm, hsm, hA', hP, Finset.mul_sum,
    ← Finset.add_sum_erase _ _ (Finset.mem_univ j₀)]
  have hn : nsq y = ∑ j, z j ^ 2 := by
    rw [← nsq_conjTranspose_mulVec (eigU_mul_star hA) y, nsq_real]
  rw [hn, Finset.mul_sum]
  have e : t * (hA.eigenvalues j₀ * z j₀ ^ 2) - t * hA.eigenvalues j₀ * z j₀ ^ 2 = 0 := by ring
  rw [add_sub_right_comm, e, zero_add]
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  refine (Finset.sum_le_sum fun j hj ↦ ?_).trans
    (Finset.sum_le_sum_of_subset_of_nonneg (Finset.erase_subset j₀ Finset.univ)
      fun j _ _ ↦ mul_nonneg hb (sq_nonneg _))
  rw [← mul_assoc, abs_mul, abs_of_nonneg (sq_nonneg (z j))]
  exact mul_le_mul_of_nonneg_right (h j (Finset.ne_of_mem_erase hj)) (sq_nonneg _)

omit [DecidableEq D] in
/-- The complex form of a real matrix is the sum of the forms of the real and imaginary
parts. -/
theorem qform_map_ofReal (A : Matrix D D ℝ) (y : D → ℂ) :
    qform (A.map (↑) : Matrix D D ℂ) y =
      qform A (fun x ↦ (y x).re) + qform A (fun x ↦ (y x).im) := by
  simp only [qform, dotProduct, mulVec, map_apply, RCLike.re_to_real, star_trivial,
    Pi.star_apply, Finset.mul_sum, RCLike.re_eq_complex_re, Complex.re_sum,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun x _ ↦ Finset.sum_congr rfl fun z _ ↦ ?_
  simp only [Complex.star_def, Complex.mul_re, Complex.conj_re, Complex.conj_im,
    Complex.ofReal_re, Complex.ofReal_im, Complex.mul_im]
  ring

omit [DecidableEq D] in
/-- The squared norm of a complex vector is that of its real and imaginary parts. -/
theorem nsq_complex (y : D → ℂ) :
    nsq y = nsq (fun x ↦ (y x).re) + nsq (fun x ↦ (y x).im) := by
  simp only [nsq_eq_sum, Real.norm_eq_abs, sq_abs, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun x _ ↦ ?_
  rw [Complex.sq_norm, Complex.normSq_apply]
  ring

end BlockForms

/-! ### The leading line of the transfer matrix -/

section Leading

open Spectral D2Symmetry

variable [Fintype ι] [DecidableEq ι] {C : ι → Matrix d d ℝ} {o : Matrix d d ℝ}
  (S : A4Symmetry C o) (m : ℕ) (τ : ℝ)

omit [Fintype ι] [DecidableEq ι] in
/-- The transfer matrix is the sum of its four `D₂` blocks. -/
theorem transfer_eq_sum_blocks : transfer C o m τ = ∑ c, S.block m τ c := by
  ext γ γ'
  simp only [Matrix.sum_apply, blockOf, of_apply]
  rw [Finset.sum_eq_single (S.key γ)]
  · by_cases h : S.key γ = S.key γ'
    · simp [h]
    · simp [S.transfer_eq_zero_of_key_ne τ h]
  · intro c _ hc; simp [Ne.symm hc]
  · simp

omit [DecidableEq ι] in
/-- A bound on the form of a block passes to the block of the relabeled character. -/
theorem abs_qform_block_relabel (c : Bool × Bool) {b : ℝ}
    (h : ∀ y, |qform (S.block m τ c) y| ≤ b * nsq y) (y : (Fin m → Option ι) → ℝ) :
    |qform (S.block m τ (d2KeyPerm c)) y| ≤ b * nsq y := by
  change |qform (blockOf (transfer C o m τ) S.key (d2KeyPerm c)) y| ≤ b * nsq y
  rw [blockOf_reindex (transfer C o m τ) S.key (S.relabel m) d2KeyPerm_injective
    (S.cycle.transfer_map m τ) S.key_relabel c, qform_reindex, ← nsq_comp_equiv (S.relabel m) y]
  exact h _

omit [Fintype ι] [DecidableEq ι] in
/-- The block of a character is supported on that character. -/
theorem block_support (c : Bool × Bool) (x z : Fin m → Option ι)
    (h : S.block m τ c x z ≠ 0) : S.key x = c ∧ S.key z = c := by
  by_contra hc
  exact h (by simp [blockOf, hc])

/-- **The leading line of the transfer matrix from two blocks.** Let `t > 0`. If the trivial
block of `X` has an eigenvalue `λ_{j₀}` with `|t λ_{j₀}| > r ≥ 0`, and `|t λ| ≤ r` for all its
other eigenvalues and all eigenvalues of the block `N`, then `t X` has an isolated leading line
`ψ` with eigenvalue `t λ_{j₀}` and bound `r`; `ψ` lies in the trivial `D₂` character and is fixed
by the cyclic relabeling. -/
theorem isLeadingLine_of_blocks {t : ℝ} (ht : 0 < t) (j₀ : Fin m → Option ι) {r : ℝ}
    (hr0 : 0 ≤ r)
    (hJ : ∀ j, j ≠ j₀ → |t * (S.block_isHermitian m τ sectorJ).eigenvalues j| ≤ r)
    (hN : ∀ j, |t * (S.block_isHermitian m τ sectorN).eigenvalues j| ≤ r)
    (hr : r < |t * (S.block_isHermitian m τ sectorJ).eigenvalues j₀|) :
    ∃ ψ : (Fin m → Option ι) → ℂ, IsLeadingLine ((t • transfer C o m τ).map (↑)) ψ
        (t * (S.block_isHermitian m τ sectorJ).eigenvalues j₀) r ∧
      (∀ γ, S.key γ ≠ sectorJ → ψ γ = 0) ∧ cycleMat S m *ᵥ ψ = ψ := by
  set hJh := S.block_isHermitian m τ sectorJ
  set hNh := S.block_isHermitian m τ sectorN
  set e₀ := hJh.eigenvalues j₀
  set ψr : (Fin m → Option ι) → ℝ := eigU hJh *ᵥ Pi.single j₀ 1
  have hBψ : S.block m τ sectorJ *ᵥ ψr = e₀ • ψr := by
    have := mulVec_eigU_single hJh j₀
    convert this using 2
  have he0 : e₀ ≠ 0 := by
    intro h; rw [h, mul_zero, abs_zero] at hr; linarith
  -- The eigenvector lies in the trivial character.
  have hsupp : ∀ γ, S.key γ ≠ sectorJ → ψr γ = 0 := by
    intro γ hγ
    have h1 := congrFun hBψ γ
    simp only [mulVec, dotProduct, blockOf, of_apply, hγ, false_and, ite_false, zero_mul,
      Finset.sum_const_zero, Pi.smul_apply, smul_eq_mul] at h1
    exact (mul_eq_zero.1 h1.symm).resolve_left he0
  have hTψ : transfer C o m τ *ᵥ ψr = e₀ • ψr := by
    rw [← hBψ]; ext γ; simp only [mulVec, dotProduct]
    refine Finset.sum_congr rfl fun z _ ↦ ?_
    by_cases hz : S.key z = sectorJ
    · by_cases hγ : S.key γ = sectorJ
      · simp [blockOf, hz, hγ]
      · rw [S.transfer_eq_zero_of_key_ne (γ := γ) (γ' := z) τ (by rw [hz]; exact hγ)]
        simp [blockOf, hγ]
    · simp [hsupp z hz]
  have hunit : nsq ψr = 1 := by
    rw [nsq_eigU, nsq_real]
    simp [Pi.single_apply]
  -- The real rest, block by block.
  set Pr := vecMulVec ψr ψr
  have hsm : ∀ (c : ℝ) (B : Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ) y,
      qform (c • B) y = c * qform B y := fun c B y ↦ by
    simpa using qform_smul_left (𝕜 := ℝ) c B y
  have hNb : ∀ y, |qform (S.block m τ sectorN) y| ≤ r / t * nsq y :=
    abs_qform_le_of_eigenvalues hNh fun j ↦ by
      rw [le_div_iff₀ ht, mul_comm]; simpa [abs_mul, abs_of_pos ht] using hN j
  have hN1 : ∀ y, |qform (S.block m τ (true, false)) y| ≤ r / t * nsq y :=
    abs_qform_block_relabel S m τ sectorN hNb
  have hN2 : ∀ y, |qform (S.block m τ (false, true)) y| ≤ r / t * nsq y :=
    abs_qform_block_relabel S m τ (true, false) hN1
  have hpiece : ∀ (c : Bool × Bool) y, (∀ y, |qform (S.block m τ c) y| ≤ r / t * nsq y) →
      |qform (t • S.block m τ c) y| ≤ r * nsq (restrictKey S.key c y) := fun c y hc ↦ by
    rw [hsm, qform_restrictKey (block_support S m τ c) y, abs_mul, abs_of_pos ht]
    calc t * |qform (S.block m τ c) (restrictKey S.key c y)|
        ≤ t * (r / t * nsq (restrictKey S.key c y)) := mul_le_mul_of_nonneg_left (hc _) ht.le
      _ = r * nsq (restrictKey S.key c y) := by field_simp
  have hJpiece : ∀ y, |qform (t • S.block m τ sectorJ - (t * e₀) • Pr) y| ≤
      r * nsq (restrictKey S.key sectorJ y) := fun y ↦ by
    have hsuppJ : ∀ x z, (t • S.block m τ sectorJ - (t * e₀) • Pr) x z ≠ 0 →
        S.key x = sectorJ ∧ S.key z = sectorJ := fun x z hxz ↦ by
      by_contra hc
      apply hxz
      rcases not_and_or.1 hc with hx | hz
      · simp [Pr, vecMulVec_apply, blockOf, hx, hsupp x hx]
      · simp [Pr, vecMulVec_apply, blockOf, hz, hsupp z hz]
    rw [qform_restrictKey hsuppJ y]
    exact abs_qform_sub_proj_le hJh j₀ t hr0 hJ _
  have hrest : ∀ y, |qform (t • transfer C o m τ - (t * e₀) • Pr) y| ≤ r * nsq y := fun y ↦ by
    have hdecomp : t • transfer C o m τ - (t * e₀) • Pr =
        (t • S.block m τ sectorJ - (t * e₀) • Pr) + t • S.block m τ (true, false) +
          t • S.block m τ (false, true) + t • S.block m τ sectorN := by
      rw [transfer_eq_sum_blocks S m τ]
      simp only [Fintype.sum_prod_type, Fintype.sum_bool, smul_add, sectorJ, sectorN]
      abel
    rw [hdecomp, qform_add_left, qform_add_left, qform_add_left,
      nsq_eq_sum_restrictKey S.key y]
    simp only [Fintype.sum_prod_type, Fintype.sum_bool]
    have h1 := hJpiece y
    have h2 := hpiece (true, false) y hN1
    have h3 := hpiece (false, true) y hN2
    have h4 := hpiece sectorN y hNb
    simp only [sectorJ, sectorN] at h1 h4 ⊢
    set A := qform (t • S.block m τ (true, true) - (t * e₀) • Pr) y
    set B := qform (t • S.block m τ (true, false)) y
    set C' := qform (t • S.block m τ (false, true)) y
    set D' := qform (t • S.block m τ (false, false)) y
    have t1 := abs_add_le (A + B + C') D'
    have t2 := abs_add_le (A + B) C'
    have t3 := abs_add_le A B
    linarith
  -- Complexify.
  set ψ : (Fin m → Option ι) → ℂ := fun γ ↦ (ψr γ : ℂ)
  have hmap : ((t • transfer C o m τ).map (↑) : Matrix _ _ ℂ) -
      ((t * e₀ : ℝ) : ℂ) • lineProj ψ = (t • transfer C o m τ - (t * e₀) • Pr).map (↑) := by
    ext x z
    simp [lineProj, Pr, ψ, vecMulVec_apply, Complex.conj_ofReal]
  have hL : IsLeadingLine ((t • transfer C o m τ).map (↑)) ψ (t * e₀) r := by
    refine ⟨?_, ?_, ?_, fun y ↦ ?_⟩
    · ext x z
      simp [conjTranspose_apply, transfer_symm C o m τ z x, Complex.conj_ofReal]
    · have : nsq ψr = ∑ γ, ψr γ ^ 2 := nsq_real ψr
      simp only [dotProduct, Pi.star_apply, ψ, Complex.star_def, Complex.conj_ofReal,
        ← Complex.ofReal_mul, ← Complex.ofReal_sum]
      have h1 : ∑ γ, ψr γ * ψr γ = 1 := by rw [← hunit, this]; simp [sq]
      rw [h1, Complex.ofReal_one]
    · ext γ
      have := congrFun hTψ γ
      simp only [mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at this
      simp only [mulVec, dotProduct, map_apply, Matrix.smul_apply, smul_eq_mul, ψ,
        Pi.smul_apply, mul_assoc, ← Finset.mul_sum, ← Complex.ofReal_sum,
        ← Complex.ofReal_mul, this]
    · rw [hmap, qform_map_ofReal, nsq_complex]
      have h1 := hrest (fun x ↦ (y x).re)
      have h2 := hrest (fun x ↦ (y x).im)
      calc _ ≤ |qform (t • transfer C o m τ - (t * e₀) • Pr) fun x ↦ (y x).re| +
            |qform (t • transfer C o m τ - (t * e₀) • Pr) fun x ↦ (y x).im| := abs_add_le _ _
        _ ≤ _ := by linarith
  refine ⟨ψ, hL, fun γ hγ ↦ by simp [ψ, hsupp γ hγ], ?_⟩
  -- The relabeling maps `ψ` to `a ψ` with `a` real and `a³ = 1`.
  set P := cycleMat S m
  have hPX : P * (t • transfer C o m τ).map (↑) = (t • transfer C o m τ).map (↑) * P :=
    relabelMat_mul_comm S.cycle.perm
    fun γ γ' ↦ by
      rw [Matrix.smul_apply, Matrix.smul_apply]
      exact congrArg (t • ·) (S.cycle.transfer_map m τ γ γ')
  have hφ : (t • transfer C o m τ).map (↑) *ᵥ (P *ᵥ ψ) = ((t * e₀ : ℝ) : ℂ) • (P *ᵥ ψ) := by
    rw [mulVec_mulVec, ← hPX, ← mulVec_mulVec, hL.eigen, mulVec_smul]
  have hlt : r < |t * e₀| := hr
  have heq := hL.eq_smul_of_eigen hlt hφ
  set a := star ψ ⬝ᵥ (P *ᵥ ψ)
  have hreal : a = ((∑ γ, ψr γ * ψr (relabelSeq S.cycle.perm γ) : ℝ) : ℂ) := by
    simp only [a, P, cycleMat, relabelMat_mulVec, dotProduct, Pi.star_apply, ψ,
      Complex.star_def, Complex.conj_ofReal, Complex.ofReal_sum, Complex.ofReal_mul]
  have hcube : P *ᵥ (P *ᵥ (P *ᵥ ψ)) = ψ := by
    rw [mulVec_mulVec, mulVec_mulVec, cycleMat_cube, one_mulVec]
  have h2 : P *ᵥ (P *ᵥ ψ) = (a * a) • ψ := by rw [heq, mulVec_smul, heq, smul_smul]
  have h3 : P *ᵥ (P *ᵥ (P *ᵥ ψ)) = (a * (a * a)) • ψ := by
    rw [heq, mulVec_smul, mulVec_smul, h2, smul_smul]
  rw [h3] at hcube
  have ha3 : a * (a * a) = 1 := by
    have := congrArg (fun f ↦ star ψ ⬝ᵥ f) hcube
    simp only [dotProduct_smul, hL.unit, smul_eq_mul, mul_one] at this
    exact this
  set x := ∑ γ, ψr γ * ψr (relabelSeq S.cycle.perm γ)
  have hx3 : x * (x * x) = 1 := by
    have := ha3
    rw [hreal] at this
    exact_mod_cast this
  have hx : x = 1 := by
    have h1 : (x - 1) * (x ^ 2 + x + 1) = 0 := by linear_combination hx3
    rcases mul_eq_zero.1 h1 with h | h
    · linarith
    · nlinarith [sq_nonneg (x + 1 / 2)]
  rw [heq, hreal, hx, Complex.ofReal_one, one_smul]

end Leading

end CertifiedGaps.SingleSite
