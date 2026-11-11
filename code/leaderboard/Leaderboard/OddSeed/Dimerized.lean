/-
Copyright (c) 2026 CertifiedGaps contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: CertifiedGaps contributors
-/
import Leaderboard.OddSeed.Moment
import CertifiedGaps.Models.SpinOneDimerized
import CertifiedGaps.Transfer.ProductLimit
import Mathlib.Analysis.Matrix.Normed
import Mathlib.Analysis.Normed.Algebra.Exponential

/-!
# The positive cell-transfer operator of the dimerized spin-1 chain (odd seeds in Lean)

For the dimerized spin-1 chain `H = ∑ⱼ (1 + (-1)ʲ δ) Sⱼ · Sⱼ₊₁`, every bond is reflection positive
with real factors and no internal reflection at the single-site level:
`-J S·S = J ∑_α G_α ⊗ G_α` with `G_α = i S^α` real antisymmetric. A ring of `L` two-spin cells is a
ring of `2L` sites, and `2L` is always even, so the standard transfer factorization applies at
**every** cell count `L`: cutting the cell across its internal (weak) bond factors the cell
transfer operator as `X_β = K K* ⪰ 0` (`paper/thermal-mpo-gap-scope.tex`, "Positive cell-transfer
operators and odd seed lengths"). Concretely, the Trotter slice of the `2L`-site ring factorizes as
`Tr (slice^M) = Tr ((Fᵀ F)^L)` for a real transfer matrix `F` and every `L ≥ 1`
(`trace_siteSlice_pow_eq_trace_transfer`), so all transfer eigenvalues are nonnegative and the
moment form holds at all lengths. Transporting along the cell/site relabeling (`siteHam_reindex`)
identifies the site Hamiltonian with the cell model's Hamiltonian, giving

* `hasSpatialMomentBoundsAll : HasSpatialMomentBoundsAll (model δ).spectrum 2` for `-1 ≤ δ ≤ 1`.

## Main results

* `trace_siteSlice_pow_eq_trace_transfer`: `Tr (siteSlice δ L τ ^ M) = Tr ((Fᵀ F) ^ L)`.
* `hasSpatialMomentBoundsAll`: the all-lengths moment inequalities of `model δ`.
-/

open Matrix Real Filter NormedSpace
open scoped Kronecker ComplexOrder Matrix.Norms.Operator

namespace CertifiedGaps.SpinOneDimerized

open SpinOne

/-! ### The site-level ring -/

/-- The complex spin generator `G_α = i S^α`. -/
noncomputable def Gc (α : Fin 3) : Matrix (Fin 3) (Fin 3) ℂ := (G α).map (↑)

/-- The site-level bond matrix of coupling `J`: `h = J S·S = -J ∑_α G_α ⊗ G_α`. -/
noncomputable def siteBondM (J : ℝ) : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ :=
  -(J : ℂ) • ∑ α : Fin 3, Gc α ⊗ₖ Gc α

/-- The real local factor `√(τJ) G_α`. -/
noncomputable def realSFactor (τ J : ℝ) : Option (Fin 3) → Matrix (Fin 3) (Fin 3) ℝ
  | none => 1
  | some α => Real.sqrt (τ * J) • G α

/-- The complex local factor `√(τJ) G_α`. -/
noncomputable def sFactor (τ J : ℝ) : Option (Fin 3) → Matrix (Fin 3) (Fin 3) ℂ
  | none => 1
  | some α => (Real.sqrt (τ * J) : ℂ) • Gc α

theorem sFactor_eq_map (τ J : ℝ) (γ : Option (Fin 3)) :
    sFactor τ J γ = (realSFactor τ J γ).map (↑) := by
  cases γ with
  | none => ext a b; simp only [sFactor, realSFactor, map_apply, one_apply]; split_ifs <;> simp
  | some α => ext a b; simp [sFactor, realSFactor, Gc, map_apply, Matrix.smul_apply,
      Complex.ofReal_mul]

/-- The site-level ring Hamiltonian of `2L` spins: weak bonds `(2m, 2m+1)` of coupling `1 - δ`
and strong bonds `(2m+1, 2m+2)` of coupling `1 + δ`. -/
noncomputable def siteHam (δ : ℝ) (L : ℕ) :
    Matrix (Fin (2 * L) → Fin 3) (Fin (2 * L) → Fin 3) ℂ :=
  (List.ofFn fun m ↦ bondOp (evenSite L m) (finRotate (2 * L) (evenSite L m))
      (siteBondM (1 - δ))).sum +
    (List.ofFn fun m ↦ bondOp (oddSite L m) (finRotate (2 * L) (oddSite L m))
      (siteBondM (1 + δ))).sum

/-- One first-order Trotter slice of the site ring: even (weak) bonds, then odd (strong) bonds. -/
noncomputable def siteSlice (δ : ℝ) (L : ℕ) (τ : ℝ) :
    Matrix (Fin (2 * L) → Fin 3) (Fin (2 * L) → Fin 3) ℂ :=
  (List.ofFn fun m ↦ 1 - (τ : ℂ) • bondOp (evenSite L m) (finRotate (2 * L) (evenSite L m))
      (siteBondM (1 - δ))).prod *
    (List.ofFn fun m ↦ 1 - (τ : ℂ) • bondOp (oddSite L m) (finRotate (2 * L) (oddSite L m))
      (siteBondM (1 + δ))).prod

/-- The Trotter approximant `Re Tr (siteSlice δ L (β/M))^M`. -/
noncomputable def siteTrotterTrace (δ : ℝ) (M : ℕ) (β : ℝ) (L : ℕ) : ℝ :=
  (trace (siteSlice δ L (β / M) ^ M)).re

/-! ### Bond factorization at the site level (no internal reflection) -/

/-- `√(τJ) G_α ⊗ √(τJ) G_α = τJ (G_α ⊗ G_α)`. -/
theorem kronecker_sFactor_some {τ J : ℝ} (hτJ : 0 ≤ τ * J) (α : Fin 3) :
    sFactor τ J (some α) ⊗ₖ sFactor τ J (some α) = ((τ : ℂ) * (J : ℂ)) • (Gc α ⊗ₖ Gc α) := by
  have hs : sFactor τ J (some α) = (Real.sqrt (τ * J) : ℂ) • Gc α := by
    ext a b; simp [sFactor, Gc]
  have h : (Real.sqrt (τ * J) : ℂ) * Real.sqrt (τ * J) = (τ : ℂ) * (J : ℂ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτJ]
    push_cast; ring
  rw [hs, Matrix.smul_kronecker, Matrix.kronecker_smul, smul_smul, h]

/-- The tensor-product factor of label `γ` of the site bond `j` of coupling `J` (no internal
reflection: both ends carry the same factor). -/
noncomputable def siteBondSite (L : ℕ) (τ J : ℝ) (j : Fin (2 * L)) (γ : Option (Fin 3)) :
    Fin (2 * L) → Matrix (Fin 3) (Fin 3) ℂ :=
  fun i ↦ if i = j then sFactor τ J γ else if i = finRotate (2 * L) j then sFactor τ J γ else 1

/-- A first-order site bond factor as a sum over the bond label: `1 - τ h = ∑_γ C_γ ⊗ C_γ`. -/
theorem one_sub_smul_siteBond {L : ℕ} (hL : 2 ≤ 2 * L) {τ J : ℝ} (hτJ : 0 ≤ τ * J)
    (j : Fin (2 * L)) :
    1 - (τ : ℂ) • bondOp j (finRotate (2 * L) j) (siteBondM J) =
      ∑ γ : Option (Fin 3), piTensor (siteBondSite L τ J j γ) := by
  have hjk := ne_finRotate hL j
  have hb : bondOp j (finRotate (2 * L) j) (siteBondM J) =
      -∑ α : Fin 3, (J : ℂ) • bondOp j (finRotate (2 * L) j) (Gc α ⊗ₖ Gc α) := by
    rw [siteBondM, neg_smul, bondOp_neg, bondOp_smul, bondOp_sum, Finset.smul_sum]
  rw [hb, Fintype.sum_option]
  have h1 : siteBondSite L τ J j none = fun _ ↦ 1 := funext fun i ↦ by simp [siteBondSite, sFactor]
  rw [h1, piTensor_one, smul_neg, sub_neg_eq_add, Finset.smul_sum]
  have hterm : ∀ α : Fin 3, piTensor (siteBondSite L τ J j (some α)) =
      ((τ : ℂ) * (J : ℂ)) • bondOp j (finRotate (2 * L) j) (Gc α ⊗ₖ Gc α) := by
    intro α
    change piTensor (fun i ↦ if i = j then sFactor τ J (some α) else
      if i = finRotate (2 * L) j then sFactor τ J (some α) else 1) = _
    rw [← bondOp_kronecker hjk, kronecker_sFactor_some hτJ, bondOp_smul]
  congr 1
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [smul_smul, hterm α]

/-! ### The site slice as a sum of tensor products -/

/-- The local word at site `i` of one slice with even-bond labels `bc.1` and odd-bond labels
`bc.2`. -/
noncomputable def ssliceSite (δ : ℝ) (L : ℕ) (τ : ℝ)
    (bc : (Fin L → Option (Fin 3)) × (Fin L → Option (Fin 3))) (i : Fin (2 * L)) :
    Matrix (Fin 3) (Fin 3) ℂ :=
  (List.ofFn fun k ↦ siteBondSite L τ (1 - δ) (evenSite L k) (bc.1 k) i).prod *
    (List.ofFn fun k ↦ siteBondSite L τ (1 + δ) (oddSite L k) (bc.2 k) i).prod

/-- One site-level Trotter slice on a ring of `2L` sites is a sum of tensor products. -/
theorem siteSlice_eq_sum (δ : ℝ) (L : ℕ) (hL : 1 ≤ L) {τ : ℝ} (hτ : 0 ≤ τ) (hmδ : -1 ≤ δ)
    (hδ : δ ≤ 1) :
    siteSlice δ L τ = ∑ bc, piTensor (ssliceSite δ L τ bc) := by
  have h2L : 2 ≤ 2 * L := by omega
  have hτw : 0 ≤ τ * (1 - δ) := by nlinarith
  have hτs : 0 ≤ τ * (1 + δ) := by nlinarith
  rw [siteSlice]
  simp only [one_sub_smul_siteBond h2L hτw, one_sub_smul_siteBond h2L hτs]
  rw [list_prod_ofFn_sum, list_prod_ofFn_sum, Finset.sum_mul_sum, ← Finset.sum_product']
  refine Finset.sum_congr rfl fun bc _ ↦ ?_
  rw [list_prod_ofFn_piTensor, list_prod_ofFn_piTensor, piTensor_mul]
  rfl

/-- The word at the even site `2m`: the weak factor of its even bond and the strong factor of the
preceding odd bond. -/
theorem ssliceSite_evenSite (δ : ℝ) (L : ℕ) (τ : ℝ)
    (bc : (Fin L → Option (Fin 3)) × (Fin L → Option (Fin 3))) (m : Fin L) :
    ssliceSite δ L τ bc (evenSite L m) =
      sFactor τ (1 - δ) (bc.1 m) * sFactor τ (1 + δ) (bc.2 ((finRotate L).symm m)) := by
  simp only [ssliceSite, siteBondSite, evenSite_ne_rot_evenSite, evenSite_ne_oddSite, ite_false]
  rw [list_prod_ofFn_ite _ _ m (evenSite_eq_evenSite_iff m),
    list_prod_ofFn_ite _ _ _ (evenSite_eq_rot_oddSite_iff m)]

/-- The word at the odd site `2m+1`: the weak factor of its even bond and the strong factor of its
odd bond. -/
theorem ssliceSite_oddSite (δ : ℝ) (L : ℕ) (τ : ℝ)
    (bc : (Fin L → Option (Fin 3)) × (Fin L → Option (Fin 3))) (m : Fin L) :
    ssliceSite δ L τ bc (oddSite L m) =
      sFactor τ (1 - δ) (bc.1 m) * sFactor τ (1 + δ) (bc.2 m) := by
  simp only [ssliceSite, siteBondSite, oddSite_ne_evenSite, oddSite_ne_rot_oddSite, ite_false]
  rw [list_prod_ofFn_ite _ _ m (oddSite_eq_rot_evenSite_iff m),
    list_prod_ofFn_ite _ _ m (oddSite_eq_oddSite_iff m)]

/-! ### The site transfer matrix and its factorization -/

/-- The time-ordered word over `M` slices with even-bond labels `p` and odd-bond labels `q`. -/
noncomputable def siteTimeWord (δ : ℝ) (M : ℕ) (τ : ℝ) (p q : Fin M → Option (Fin 3)) :
    Matrix (Fin 3) (Fin 3) ℂ :=
  (List.ofFn fun t ↦ sFactor τ (1 - δ) (p t) * sFactor τ (1 + δ) (q t)).prod

/-- The one-site transfer matrix `F p q = Tr (siteTimeWord p q)` on time labels. -/
noncomputable def siteTransfer (δ : ℝ) (M : ℕ) (τ : ℝ) :
    Matrix (Fin M → Option (Fin 3)) (Fin M → Option (Fin 3)) ℂ :=
  of fun p q ↦ trace (siteTimeWord δ M τ p q)

/-- The trace of `M` site slices as a sum over space-time bond labels. -/
theorem trace_siteSlice_pow_eq_sum (δ : ℝ) (L : ℕ) (hL : 1 ≤ L) {τ : ℝ} (hτ : 0 ≤ τ) (hmδ : -1 ≤ δ)
    (hδ : δ ≤ 1) (M : ℕ) :
    trace (siteSlice δ L τ ^ M) =
      ∑ a : Fin M → (Fin L → Option (Fin 3)) × (Fin L → Option (Fin 3)),
        ∏ m, siteTransfer δ M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 ((finRotate L).symm m)) *
          siteTransfer δ M τ (fun t ↦ (a t).1 m) (fun t ↦ (a t).2 m) := by
  rw [← List.prod_replicate, ← List.ofFn_const, siteSlice_eq_sum δ L hL hτ hmδ hδ,
    list_prod_ofFn_sum, trace_sum]
  refine Finset.sum_congr rfl fun a _ ↦ ?_
  rw [list_prod_ofFn_piTensor, trace_piTensor, prod_fin_two_mul]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  simp only [ssliceSite_evenSite, ssliceSite_oddSite, siteTransfer, of_apply, siteTimeWord,
    List.ofFn_eq_map]

/-- **Site transfer factorization (the `X_β = K K* ⪰ 0` of the scope note).** `Tr (siteSlice^M) =
`Tr ((Fᵀ F)^L)` on the ring of `2L` sites, for every `L ≥ 1` — the cell count `L` of the dimerized
chain may be odd. -/
theorem trace_siteSlice_pow_eq_trace_transfer (δ : ℝ) (L : ℕ) (hL : 1 ≤ L) {τ : ℝ} (hτ : 0 ≤ τ)
    (hmδ : -1 ≤ δ) (hδ : δ ≤ 1) (M : ℕ) :
    trace (siteSlice δ L τ ^ M) =
      trace (((siteTransfer δ M τ)ᵀ * siteTransfer δ M τ) ^ L) := by
  classical
  rw [trace_siteSlice_pow_eq_sum δ L hL hτ hmδ hδ M, trace_pow_eq_sum_closedPath _ hL]
  let F := siteTransfer δ M τ
  let e := (Equiv.arrowProdEquivProdArrow (Fin M) (fun _ ↦ Fin L → Option (Fin 3))
    (fun _ ↦ Fin L → Option (Fin 3))).trans
    (Equiv.prodCongr (Equiv.piComm _) (Equiv.piComm _))
  refine (Fintype.sum_equiv e _
    (fun pq ↦ ∏ m, F (pq.1 m) (pq.2 ((finRotate L).symm m)) * F (pq.1 m) (pq.2 m))
    fun a ↦ rfl).trans ?_
  rw [Fintype.sum_prod_type, Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ ↦ ?_
  rw [← Fintype.prod_sum fun m x ↦ F x (q ((finRotate L).symm m)) * F x (q m),
    ← Equiv.prod_comp (finRotate L)]
  refine Finset.prod_congr rfl fun m _ ↦ ?_
  simp only [Equiv.symm_apply_apply, mul_apply, transpose_apply]
  rfl

/-- The site transfer matrix has real entries. -/
theorem star_siteTransfer_apply (δ : ℝ) (M : ℕ) (τ : ℝ) (p q : Fin M → Option (Fin 3)) :
    star (siteTransfer δ M τ p q) = siteTransfer δ M τ p q := by
  have h : siteTimeWord δ M τ p q = (List.ofFn fun t ↦ realSFactor τ (1 - δ) (p t) *
      realSFactor τ (1 + δ) (q t)).prod.map (↑) := by
    change _ = (Complex.ofRealHom.mapMatrix : Matrix (Fin 3) (Fin 3) ℝ →+* _) _
    rw [map_list_prod, List.map_ofFn]
    simp only [siteTimeWord, sFactor_eq_map, Function.comp_def, RingHom.mapMatrix_apply,
      Matrix.map_mul]
    rfl
  simp only [siteTransfer, of_apply, h, trace, diag_apply, map_apply, ← Complex.ofReal_sum]
  exact Complex.conj_ofReal _

/-- The conjugate transpose of the site transfer matrix is its transpose. -/
theorem conjTranspose_siteTransfer (δ : ℝ) (M : ℕ) (τ : ℝ) :
    (siteTransfer δ M τ)ᴴ = (siteTransfer δ M τ)ᵀ := by
  ext p q
  exact star_siteTransfer_apply δ M τ q p

/-- **Positive site transfer operator.** For every number of slices `M`, every `τ ≥ 0` and
`-1 ≤ δ ≤ 1` there is a positive semidefinite `T`, independent of the ring size, with
`Tr (siteSlice δ L τ ^ M) = Tr (T ^ L)` for every ring of `2L ≥ 2` sites — that is, every cell
count `L ≥ 1` of the dimerized chain, odd included. -/
theorem exists_posSemidef_siteTransfer (δ : ℝ) (hmδ : -1 ≤ δ) (hδ : δ ≤ 1) (M : ℕ) (τ : ℝ)
    (hτ : 0 ≤ τ) :
    ∃ (D : Type) (_ : Fintype D) (_ : DecidableEq D) (T : Matrix D D ℂ), T.PosSemidef ∧
      ∀ L : ℕ, 1 ≤ L → trace (siteSlice δ L τ ^ M) = trace (T ^ L) := by
  classical
  refine ⟨Fin M → Option (Fin 3), inferInstance, inferInstance,
    (siteTransfer δ M τ)ᴴ * siteTransfer δ M τ, posSemidef_conjTranspose_mul_self _,
    fun L hL ↦ ?_⟩
  rw [conjTranspose_siteTransfer, trace_siteSlice_pow_eq_trace_transfer δ L hL hτ hmδ hδ M]

/-! ### The all-lengths moment form of the site Trotter traces, and convergence -/

/-- Each site Trotter approximant has the all-lengths moment form from `L = 1` on. -/
theorem hasMomentFormAll_siteTrotter (δ : ℝ) (hmδ : -1 ≤ δ) (hδ : δ ≤ 1) (M : ℕ) {β : ℝ}
    (hβ : 0 < β) :
    HasMomentFormAll (fun L ↦ siteTrotterTrace δ M β L) 1 := by
  obtain ⟨D, _, _, T, hT, htr⟩ := exists_posSemidef_siteTransfer δ hmδ hδ M (β / M) (by positivity)
  exact hasMomentFormAll_of_trace_pow hT fun L hL ↦ by rw [siteTrotterTrace, htr L hL]

/-- The site Trotter traces converge to the Gibbs trace of the site Hamiltonian. -/
theorem tendsto_siteTrotterTrace (δ : ℝ) (β : ℝ) (L : ℕ) :
    Tendsto (fun M : ℕ ↦ siteTrotterTrace δ M β L) atTop (nhds (gibbsTrace (siteHam δ L) β)) := by
  let hs : List (Matrix (Fin (2 * L) → Fin 3) (Fin (2 * L) → Fin 3) ℂ) :=
    (List.ofFn fun m ↦ bondOp (evenSite L m) (finRotate (2 * L) (evenSite L m))
        (siteBondM (1 - δ))) ++
      (List.ofFn fun m ↦ bondOp (oddSite L m) (finRotate (2 * L) (oddSite L m))
        (siteBondM (1 + δ)))
  have hslice : ∀ τ : ℝ, siteSlice δ L τ = (hs.map fun h ↦ 1 - (τ : ℂ) • h).prod := by
    intro τ
    simp [siteSlice, hs, List.map_append, List.prod_append, List.map_ofFn, Function.comp_def]
  have hsum : hs.sum = siteHam δ L := by
    simp [hs, List.sum_append, List.sum_ofFn, siteHam]
  have h := tendsto_list_prod_one_add_smul_pow hs (-β)
  rw [hsum] at h
  have hc : Continuous fun A : Matrix (Fin (2 * L) → Fin 3) (Fin (2 * L) → Fin 3) ℂ ↦
      (trace A).re :=
    Complex.continuous_re.comp (continuous_id.matrix_trace)
  have hcomp := (hc.tendsto _).comp h
  simpa [siteTrotterTrace, hslice, gibbsTrace, Function.comp_def, neg_div, sub_eq_add_neg]
    using hcomp

/-- The site Hamiltonian's partition functions have the all-lengths moment inequalities. -/
theorem siteHam_spatialMomentBoundsAll (δ : ℝ) (hmδ : -1 ≤ δ) (hδ : δ ≤ 1) (β : ℝ) (hβ : 0 < β) :
    SpatialMomentBoundsAll (fun L ↦ gibbsTrace (siteHam δ L) β) 1 :=
  SpatialMomentBoundsAll.of_tendsto_hasMomentFormAll one_pos
    (Eventually.of_forall fun M ↦ hasMomentFormAll_siteTrotter δ hmδ hδ M hβ)
    fun L _ ↦ tendsto_siteTrotterTrace δ β L

/-! ### Transport to the cell model -/

/-- The `a` spin of cell `m`, as a site: `2m`. -/
def siteA (L : ℕ) (m : Fin L) : Fin (2 * L) := ⟨2 * m.val, by omega⟩

/-- The `b` spin of cell `m`, as a site: `2m + 1`. -/
def siteB (L : ℕ) (m : Fin L) : Fin (2 * L) := ⟨2 * m.val + 1, by omega⟩

/-- The cell containing a site. -/
def cellOfSite (L : ℕ) (i : Fin (2 * L)) : Fin L := ⟨i.val / 2, by omega⟩

/-- The rotation of the `a` spin of cell `m` is the `b` spin of cell `m`. -/
theorem finRotate_siteA (L : ℕ) (m : Fin L) : finRotate (2 * L) (siteA L m) = siteB L m := by
  have hm : (m : ℕ) < L := m.isLt
  have hlt : 2 * (m : ℕ) + 1 < 2 * L := by omega
  apply Fin.ext
  rw [val_finRotate_eq_mod, show (siteA L m : ℕ) = 2 * (m : ℕ) from rfl, Nat.mod_eq_of_lt hlt,
    show (siteB L m : ℕ) = 2 * (m : ℕ) + 1 from rfl]

theorem finRotate_siteB (L : ℕ) (m : Fin L) :
    finRotate (2 * L) (siteB L m) = siteA L (finRotate L m) := by
  have hm : (m : ℕ) < L := m.isLt
  apply Fin.ext
  rw [val_finRotate_eq_mod, show (siteB L m : ℕ) = 2 * (m : ℕ) + 1 from rfl,
    show (siteA L (finRotate L m) : ℕ) = 2 * ((finRotate L m : ℕ)) from rfl, val_finRotate_eq_mod]
  by_cases h : (m : ℕ) + 1 < L
  · rw [Nat.mod_eq_of_lt h, Nat.mod_eq_of_lt (by omega)]
    ring
  · have h' : (m : ℕ) + 1 = L := by omega
    have h2 : 2 * (m : ℕ) + 2 = 2 * L := by omega
    rw [h2, h']
    simp

/-- The `siteA`/`siteB` maps are injective and have disjoint images. -/
theorem siteA_inj (L : ℕ) : Function.Injective (siteA L) := by
  intro a b h; ext; have := congrArg Fin.val h; simp only [siteA, Fin.val_mk] at this; omega

theorem siteB_inj (L : ℕ) : Function.Injective (siteB L) := by
  intro a b h; ext; have := congrArg Fin.val h; simp only [siteB, Fin.val_mk] at this; omega

theorem siteA_ne_siteB (L : ℕ) (i j : Fin L) : siteA L i ≠ siteB L j := by
  intro h; have := congrArg Fin.val h; simp only [siteA, siteB, Fin.val_mk] at this; omega

theorem siteB_ne_siteA (L : ℕ) (i j : Fin L) : siteB L i ≠ siteA L j := by
  intro h; have := congrArg Fin.val h; simp only [siteA, siteB, Fin.val_mk] at this; omega

/-- Every site is the `a` or the `b` spin of its cell. -/
theorem site_eq_siteA_or_siteB (L : ℕ) (s : Fin (2 * L)) :
    s = siteA L (cellOfSite L s) ∨ s = siteB L (cellOfSite L s) := by
  have hs : s.val < 2 * L := s.isLt
  by_cases h : s.val % 2 = 0
  · left; apply Fin.ext; simp only [siteA, cellOfSite, Fin.val_mk]; omega
  · right; apply Fin.ext; simp only [siteB, cellOfSite, Fin.val_mk]; omega

/-- The cell/site relabeling: a ring of `L` two-spin cells is a ring of `2L` spins. -/
def cellSiteEquiv (L : ℕ) : (Fin L → (Fin 3 × Fin 3)) ≃ (Fin (2 * L) → Fin 3) :=
  (Equiv.piCongrRight fun _ : Fin L ↦ (piFinTwoEquiv (fun _ : Fin 2 ↦ Fin 3)).symm).trans
    ((Equiv.curry (Fin L) (Fin 2) (Fin 3)).symm.trans
      (Equiv.arrowCongr
        (((finCongr (by ring : 2 * L = L * 2)).trans finProdFinEquiv.symm).symm)
        (Equiv.refl (Fin 3))))

@[simp] theorem cellSiteEquiv_symm_apply (L : ℕ) (τ : Fin (2 * L) → Fin 3) (m : Fin L) :
    (cellSiteEquiv L).symm τ m = (τ (siteA L m), τ (siteB L m)) := by
  simp [cellSiteEquiv, siteA, siteB, piFinTwoEquiv, Equiv.curry, finProdFinEquiv, finCongr]
  congr 1
  exact Fin.ext (by simp only [Fin.val_mk]; omega)


/-- Reindexing along an equivalence preserves the trace. -/
theorem trace_reindex {m n : Type*} [Fintype m] [Fintype n] (e : m ≃ n) (M : Matrix m m ℂ) :
    trace (reindex e e M) = trace M := by
  simp only [trace, reindex_apply, diag_apply, submatrix_apply]
  exact Equiv.sum_comp e.symm fun i ↦ M i i

/-- Reindexing along an equivalence preserves the Gibbs trace. -/
theorem gibbsTrace_reindex {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (e : m ≃ n) (H : Matrix m m ℂ) (β : ℝ) :
    gibbsTrace (reindex e e H) β = gibbsTrace H β := by
  have hcont : Continuous (reindexAlgEquiv (R := ℂ) (A := ℂ) e : Matrix m m ℂ → Matrix n n ℂ) :=
    (reindexAlgEquiv (R := ℂ) (A := ℂ) e).toLinearMap.continuous_of_finiteDimensional
  have h2 := map_exp (reindexAlgEquiv (R := ℂ) (A := ℂ) e) hcont ((-(β : ℂ)) • H)
  simp only [coe_reindexAlgEquiv, map_smul] at h2
  rw [gibbsTrace, gibbsTrace]
  exact (congrArg (fun M : Matrix n n ℂ ↦ (trace M).re) h2.symm).trans
    (congrArg Complex.re (trace_reindex e (exp ((-(β : ℂ)) • H))))

/-- Away from the two sites of cell `j`, agreement of cell configurations is agreement of site
configurations (both sides in pair form). -/
theorem forall_cell_iff_forall_site (L : ℕ) (j : Fin L) (σ τ : Fin (2 * L) → Fin 3) :
    (∀ i : Fin L, i ≠ j → (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i))) ↔
      (∀ s : Fin (2 * L), s ≠ siteA L j → s ≠ siteB L j → σ s = τ s) := by
  constructor
  · intro h s hsA hsB
    have hc : cellOfSite L s ≠ j := by
      rcases site_eq_siteA_or_siteB L s with hs | hs
      · intro e; exact hsA (by rw [hs, e])
      · intro e; exact hsB (by rw [hs, e])
    have hcs := h (cellOfSite L s) hc
    rcases site_eq_siteA_or_siteB L s with hs | hs
    · rw [hs]; exact congrArg Prod.fst hcs
    · rw [hs]; exact congrArg Prod.snd hcs
  · intro h i hi
    have h1 : σ (siteA L i) = τ (siteA L i) :=
      h (siteA L i) (fun e ↦ hi (siteA_inj L e)) (fun e ↦ siteA_ne_siteB L i j e)
    have h2 : σ (siteB L i) = τ (siteB L i) :=
      h (siteB L i) (fun e ↦ siteB_ne_siteA L i j e) (fun e ↦ hi (siteB_inj L e))
    exact Prod.ext h1 h2

/-- The site condition (agreement away from the `b` spin of cell `j` and the `a` spin of cell `k`)
implies the cell condition (agreement away from cells `j, k`), in pair form. -/
theorem forall_site_implies_forall_cell2 (L : ℕ) (j k : Fin L) (σ τ : Fin (2 * L) → Fin 3)
    (h : ∀ s : Fin (2 * L), s ≠ siteB L j → s ≠ siteA L k → σ s = τ s)
    (i : Fin L) (hij : i ≠ j) (hik : i ≠ k) :
    (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i)) := by
  have h1 : σ (siteA L i) = τ (siteA L i) :=
    h (siteA L i) (fun e ↦ siteA_ne_siteB L i j e) (fun e ↦ hik (siteA_inj L e))
  have h2 : σ (siteB L i) = τ (siteB L i) :=
    h (siteB L i) (fun e ↦ hij (siteB_inj L e)) (fun e ↦ siteB_ne_siteA L i k e)
  exact Prod.ext h1 h2

/-- Cell agreement away from cells `j, k` gives site agreement at any site whose cell is not
`j` or `k`. -/
theorem cell2_agree_at_site (L : ℕ) (j k : Fin L) (σ τ : Fin (2 * L) → Fin 3)
    (hcell : ∀ i : Fin L, i ≠ j → i ≠ k →
      (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i)))
    (s : Fin (2 * L)) (hsj : cellOfSite L s ≠ j) (hsk : cellOfSite L s ≠ k) : σ s = τ s := by
  have h := hcell (cellOfSite L s) hsj hsk
  rcases site_eq_siteA_or_siteB L s with he | he
  · rw [he]; exact congrArg Prod.fst h
  · rw [he]; exact congrArg Prod.snd h

/-- Reindexing a site operator along the cell/site equivalence turns it into the bond operator on
the two sites of the cell. -/
theorem siteOp_reindex (L : ℕ) (j : Fin L) (A : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ) :
    (siteOp j A).reindex (cellSiteEquiv L) (cellSiteEquiv L) =
      bondOp (siteA L j) (siteB L j) A := by
  ext σ τ
  have hj : j ∈ (Finset.univ : Finset (Fin L)) := Finset.mem_univ j
  simp only [reindex_apply, submatrix_apply, siteOp, piTensor_apply, bondOp, of_apply,
    cellSiteEquiv_symm_apply]
  rw [← Finset.mul_prod_erase _ _ hj]
  rw [Finset.prod_congr rfl
    (g := fun i ↦ (1 : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ)
      (σ (siteA L i), σ (siteB L i)) (τ (siteA L i), τ (siteB L i))) (fun i hi ↦ by
      simp only [Finset.mem_erase, Finset.mem_univ, and_true] at hi
      simp [hi])]
  have hrest : (∏ i ∈ Finset.univ.erase j, (1 : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℂ)
      ((σ (siteA L i), σ (siteB L i))) ((τ (siteA L i), τ (siteB L i)))) =
      if (∀ i : Fin L, i ≠ j → (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i)))
        then 1 else 0 := by
    rw [Finset.prod_congr rfl
      (g := fun i ↦ if (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i))
        then (1 : ℂ) else 0) (fun i hi ↦ by simp [one_apply])]
    rw [Finset.prod_ite_zero, Finset.prod_const_one]
    by_cases h : ∀ i ∈ Finset.univ.erase j,
        (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i))
    · rw [if_pos h, if_pos fun i hi ↦ h i (Finset.mem_erase.2 ⟨hi, Finset.mem_univ i⟩)]
    · rw [if_neg h, if_neg (fun h' ↦ h fun i hi ↦ h' i (Finset.mem_erase.1 hi).1)]
  rw [hrest]
  have hjj : (if j = j then A else 1) (σ (siteA L j), σ (siteB L j))
      (τ (siteA L j), τ (siteB L j)) =
      A (σ (siteA L j), σ (siteB L j)) (τ (siteA L j), τ (siteB L j)) := by simp
  rw [hjj]
  by_cases h : ∀ s : Fin (2 * L), s ≠ siteA L j → s ≠ siteB L j → σ s = τ s
  · have hc := (forall_cell_iff_forall_site L j σ τ).2 h
    simp only [if_pos h, if_pos hc, mul_one]
  · have hc : ¬ (∀ i : Fin L, i ≠ j →
        (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i))) :=
      mt (forall_cell_iff_forall_site L j σ τ).1 h
    simp only [if_neg h, if_neg hc, mul_zero]

/-- The bond factors of the dimerized model: `C_α = √(1+δ) (1 ⊗ G_α)` and its reflection
`D_α = √(1+δ) (G_α ⊗ 1)`. -/
theorem Cc_eq (δ : ℝ) (α : Fin 3) : (model δ).Cc α = (Real.sqrt (1 + δ) : ℂ) • (1 ⊗ₖ Gc α) := by
  ext x y
  by_cases h : x.1 = y.1
  · simp [RPRing.Cc, model, Gc, kroneckerMap_apply, one_apply, map_apply, Matrix.smul_apply, h,
      Complex.ofReal_mul]
  · simp [RPRing.Cc, model, Gc, kroneckerMap_apply, one_apply, map_apply, Matrix.smul_apply, h]

/-- Mapping a real Kronecker product to `ℂ` commutes with the Kronecker product. -/
theorem map_kronecker_ofReal (A B : Matrix (Fin 3) (Fin 3) ℝ) :
    (A ⊗ₖ B).map Complex.ofReal = (A.map Complex.ofReal) ⊗ₖ (B.map Complex.ofReal) := by
  ext ⟨a, c⟩ ⟨b, d⟩
  simp [kroneckerMap_apply, map_apply, Complex.ofReal_mul]

/-- Mapping a real scalar multiple to `ℂ` commutes with the scalar multiple. -/
theorem map_smul_ofReal (c : ℝ) (M : Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℝ) :
    (c • M).map Complex.ofReal = (c : ℂ) • M.map Complex.ofReal := by
  ext ⟨a, c'⟩ ⟨b, d⟩
  simp [map_apply, Matrix.smul_apply, Complex.ofReal_mul]

theorem map_sum_ofReal {ι : Type*} (f : ι → Matrix (Fin 3 × Fin 3) (Fin 3 × Fin 3) ℝ)
    (s : Finset ι) :
    (∑ α ∈ s, f α).map Complex.ofReal = ∑ α ∈ s, (f α).map Complex.ofReal := by
  ext p q
  simp [map_apply, Matrix.sum_apply, Complex.ofReal_sum]

theorem Dc_eq (δ : ℝ) (α : Fin 3) : (model δ).Dc α = (Real.sqrt (1 + δ) : ℂ) • (Gc α ⊗ₖ 1) := by
  have hreal : (model δ).refl * (model δ).C α * (model δ).refl
      = Real.sqrt (1 + δ) • (G α ⊗ₖ (1 : Matrix (Fin 3) (Fin 3) ℝ)) := by
    rw [model, Matrix.mul_smul, Matrix.smul_mul, cellSwap_mul_kronecker_mul_cellSwap]
  have hDc : (model δ).Dc α =
      ((model δ).refl * (model δ).C α * (model δ).refl).map Complex.ofReal := by
    ext p q
    simp [RPRing.Dc, RPRing.Rc, RPRing.Cc, Matrix.map_mul, mul_apply, map_apply]
  rw [hDc, hreal, map_smul_ofReal, map_kronecker_ofReal]
  simp [Matrix.map_one, Gc]

/-- The entry of the cell bond matrix depends only on the `b` spin of the left cell and the `a`
spin of the right cell, up to Kronecker deltas in the other two components. -/
theorem bondMatrix_entry (δ : ℝ) (hmδ : -1 ≤ δ) (a b c d a' b' c' d' : Fin 3) :
    (model δ).bondMatrix ((a, b), (c, d)) ((a', b'), (c', d')) =
      -((1 + δ : ℝ) : ℂ) * (if a = a' then 1 else 0) * (if d = d' then 1 else 0) *
        ∑ α : Fin 3, Gc α b b' * Gc α c c' := by
  have hsq : (Real.sqrt (1 + δ) : ℂ) * Real.sqrt (1 + δ) = (1 + δ : ℝ) := by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (by nlinarith : (0 : ℝ) ≤ 1 + δ)]
  simp only [RPRing.bondMatrix, Matrix.sum_apply, Matrix.smul_apply, Matrix.neg_apply, Cc_eq δ,
    Dc_eq δ, kroneckerMap_apply, smul_kronecker, kronecker_smul, one_apply, hsq, smul_smul]
  ring_nf
  by_cases ha : a = a'
  · by_cases hd : d = d'
    · simp only [ha, hd, if_true, ite_mul, one_mul, zero_mul, Finset.mem_univ, ↓reduceIte]
      rw [← Finset.sum_mul]
      ring
    · simp [ha, hd]
  · simp [ha]

/-- Reindexing the inter-cell bond of the dimerized model gives the strong site bond from the `b`
spin of a cell to the `a` spin of the next cell. -/
theorem bond_reindex (δ : ℝ) (hmδ : -1 ≤ δ) (L : ℕ) (hL : 2 ≤ L) (j : Fin L) :
    ((model δ).bond L j).reindex (cellSiteEquiv L) (cellSiteEquiv L) =
      bondOp (siteB L j) (siteA L (finRotate L j)) (siteBondM (1 + δ)) := by
  have hjr : j ≠ finRotate L j := ne_finRotate hL j
  ext σ τ
  by_cases hsite : ∀ s : Fin (2 * L), s ≠ siteB L j → s ≠ siteA L (finRotate L j) → σ s = τ s
  · have hcell : ∀ i : Fin L, i ≠ j → i ≠ finRotate L j →
        (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i)) :=
      fun i hij hik ↦ forall_site_implies_forall_cell2 L j (finRotate L j) σ τ hsite i hij hik
    simp only [RPRing.bond, reindex_apply, submatrix_apply, bondOp, of_apply, if_pos hcell,
      if_pos hsite, cellSiteEquiv_symm_apply, siteBondM, Finset.sum_apply, Matrix.smul_apply,
      Matrix.neg_apply]
    rw [bondMatrix_entry δ hmδ (σ (siteA L j)) (σ (siteB L j)) (σ (siteA L (finRotate L j)))
      (σ (siteB L (finRotate L j))) (τ (siteA L j)) (τ (siteB L j)) (τ (siteA L (finRotate L j)))
      (τ (siteB L (finRotate L j)))]
    have haa : σ (siteA L j) = τ (siteA L j) :=
      hsite (siteA L j) (siteA_ne_siteB L j j) (fun e ↦ hjr (siteA_inj L e))
    have hdd : σ (siteB L (finRotate L j)) = τ (siteB L (finRotate L j)) :=
      hsite (siteB L (finRotate L j)) (fun e ↦ hjr (siteB_inj L e.symm))
        (siteB_ne_siteA L (finRotate L j) (finRotate L j))
    simp only [haa, hdd, if_true, one_mul, Matrix.sum_apply, kroneckerMap_apply]
    push_cast
    ring
  · have hnot : ¬ (∀ s : Fin (2 * L), s ≠ siteB L j → s ≠ siteA L (finRotate L j) →
        σ s = τ s) := hsite
    push_neg at hnot
    obtain ⟨s, hsB, hsA, hs⟩ := hnot
    simp only [RPRing.bond, reindex_apply, submatrix_apply, bondOp, of_apply, if_neg hsite,
      cellSiteEquiv_symm_apply, siteBondM, Finset.sum_apply, Matrix.smul_apply, Matrix.neg_apply]
    by_cases hcell : ∀ i : Fin L, i ≠ j → i ≠ finRotate L j →
        (σ (siteA L i), σ (siteB L i)) = (τ (siteA L i), τ (siteB L i))
    · rw [if_pos hcell, bondMatrix_entry δ hmδ (σ (siteA L j)) (σ (siteB L j))
        (σ (siteA L (finRotate L j))) (σ (siteB L (finRotate L j))) (τ (siteA L j))
        (τ (siteB L j)) (τ (siteA L (finRotate L j))) (τ (siteB L (finRotate L j)))]
      have hsd : s = siteA L j ∨ s = siteB L (finRotate L j) := by
        by_contra hnot2
        push_neg at hnot2
        obtain ⟨hnot1, hnot2'⟩ := hnot2
        have hcj : cellOfSite L s ≠ j := by
          rcases site_eq_siteA_or_siteB L s with he | he
          · intro e; exact hnot1 (by rw [he, e])
          · intro e; exact hsB (by rw [he, e])
        have hck : cellOfSite L s ≠ finRotate L j := by
          rcases site_eq_siteA_or_siteB L s with he | he
          · intro e; exact hsA (by rw [he, e])
          · intro e; exact hnot2' (by rw [he, e])
        exact hs (cell2_agree_at_site L j (finRotate L j) σ τ hcell s hcj hck)
      rcases hsd with hs' | hs'
      · rw [hs'] at hs; simp only [hs, if_false, zero_mul, mul_zero, Matrix.neg_apply,
          Finset.sum_const_zero]
      · rw [hs'] at hs; simp only [hs, if_false, zero_mul, mul_zero, Matrix.neg_apply,
          Finset.sum_const_zero]
    · rw [if_neg hcell]

/-- The on-site term of the dimerized model is the weak site bond matrix. -/
theorem oc_eq_siteBondM (δ : ℝ) : (model δ).oc = siteBondM (1 - δ) := by
  have ho : (model δ).o = -(1 - δ) • ∑ α : Fin 3, G α ⊗ₖ G α := rfl
  rw [RPRing.oc, ho, siteBondM, map_smul_ofReal, map_sum_ofReal]
  simp only [map_kronecker_ofReal, Gc, Complex.ofReal_neg]

/-- The site Hamiltonian is the cell model's Hamiltonian, reindexed cell/site. -/
theorem siteHam_reindex (δ : ℝ) (hmδ : -1 ≤ δ) (L : ℕ) (hL : 2 ≤ L) :
    siteHam δ L = ((model δ).hamiltonian L).reindex (cellSiteEquiv L) (cellSiteEquiv L) := by
  classical
  let e := cellSiteEquiv L
  have radd : ∀ (X Y : Matrix (Fin L → (Fin 3 × Fin 3)) (Fin L → (Fin 3 × Fin 3)) ℂ),
      (X + Y).reindex e e = X.reindex e e + Y.reindex e e :=
    fun X Y ↦ map_add (reindexAlgEquiv (R := ℂ) (A := ℂ) e) X Y
  have rsum : ∀ (s : Finset (Fin L)) (f : Fin L → Matrix (Fin L → (Fin 3 × Fin 3))
      (Fin L → (Fin 3 × Fin 3)) ℂ),
      (∑ j ∈ s, f j).reindex e e = ∑ j ∈ s, (f j).reindex e e :=
    fun s f ↦ map_sum (reindexAlgEquiv (R := ℂ) (A := ℂ) e) f s
  rw [show ((model δ).hamiltonian L) = (∑ j : Fin L, (model δ).bond L j) +
      ∑ j : Fin L, (model δ).site L j from rfl, radd, rsum, rsum]
  have hbond : ∀ j : Fin L, ((model δ).bond L j).reindex e e =
      bondOp (siteB L j) (siteA L (finRotate L j)) (siteBondM (1 + δ)) := fun j ↦
    bond_reindex δ hmδ L hL j
  have hsite : ∀ j : Fin L, ((model δ).site L j).reindex e e =
      bondOp (siteA L j) (siteB L j) (siteBondM (1 - δ)) := by
    intro j
    rw [RPRing.site, oc_eq_siteBondM δ, siteOp_reindex]
  rw [Finset.sum_congr rfl fun j _ ↦ hbond j, Finset.sum_congr rfl fun j _ ↦ hsite j]
  rw [siteHam]
  have hweak : ∑ j : Fin L, bondOp (siteA L j) (siteB L j) (siteBondM (1 - δ)) =
      (List.ofFn fun m ↦ bondOp (evenSite L m) (finRotate (2 * L) (evenSite L m))
        (siteBondM (1 - δ))).sum := by
    rw [List.sum_ofFn]
    refine Finset.sum_congr rfl fun m _ ↦ ?_
    have he : evenSite L m = siteA L m := rfl
    rw [he, finRotate_siteA]
  have hstrong : ∑ j : Fin L, bondOp (siteB L j) (siteA L (finRotate L j)) (siteBondM (1 + δ)) =
      (List.ofFn fun m ↦ bondOp (oddSite L m) (finRotate (2 * L) (oddSite L m))
        (siteBondM (1 + δ))).sum := by
    rw [List.sum_ofFn]
    refine Finset.sum_congr rfl fun m _ ↦ ?_
    have he : oddSite L m = siteB L m := rfl
    rw [he, finRotate_siteB]
  rw [add_comm, hstrong, hweak]

/-- **The dimerized spin-1 model has the all-lengths moment inequalities** (its cell transfer
operator is positive semidefinite, `X_β = K K* ⪰ 0`), so the two-purity criterion applies at odd
seed lengths. -/
theorem hasSpatialMomentBoundsAll (δ : ℝ) (hmδ : -1 ≤ δ) (hδ : δ ≤ 1) :
    HasSpatialMomentBoundsAll (model δ).spectrum 2 := by
  intro β hβ
  have hsite := siteHam_spatialMomentBoundsAll δ hmδ hδ β hβ
  have htr : ∀ L, 2 ≤ L →
      gibbsTrace (siteHam δ L) β = partitionFunction (model δ).spectrum L β := by
    intro L hL
    rw [siteHam_reindex δ hmδ L hL, gibbsTrace_reindex, ← (model δ).partitionFunction_eq_gibbsTrace]
    exact RPRing.partitionFunction_eq_sum (model δ) L β
  refine ⟨fun n hn ↦ ?_, fun n hn ↦ ?_, fun p q hp hpq ↦ ?_⟩
  · have h := hsite.sq_le n (by omega)
    rw [htr (2 * n) (by omega), htr n hn] at h
    exact h
  · have h := hsite.squaring n (by omega)
    rw [htr (2 * n) (by omega), htr (2 * (2 * n)) (by omega), htr n hn] at h
    exact h
  · have h := hsite.le_rpow p q (by omega) hpq
    rw [htr q (by omega), htr p (by omega)] at h
    exact h

end CertifiedGaps.SpinOneDimerized
