/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.SingleSiteLimit

/-!
# Open chains with two end fields in the single-site transfer matrix

The single-site transfer matrix `X` of `Transfer/SingleSite.lean` gives every ring trace as
`Tr X^L`. This file gives the open chains with arbitrary end fields as pairings of the same `X`
with two boundary vectors, which is the paper's Proposition tr:boundary for the Trotterized
kernel (OpenAI, *A boundary-field gap for the spin-one Heisenberg chain*, 2026).

The open chain on the `k + 1` sites `0, …, k` has the bonds `0, …, k + 1`: bond `j` joins the
sites `j - 1` and `j`, and the two *caps* `0` and `k + 1` touch only the sites `0` and `k`.
A cap carries a weight vector `u : ι → ℂ` and stands for the end term `-∑_α u_α C_α`. The time
step with one label per bond is

  `V(u) = (∑_e (∏ⱼ w_u(e j) κⱼ(e j)) ⨂ᵢ M(e i, e (i + 1))) · ⨂ᵢ (1 - u o)`,

with `κ₀(α) = u_L(α)`, `κ_{k+1}(α) = u_R(α)` and `κⱼ = 1` on the bonds of the chain. Its
derivative at `0` is `-H` for the open Hamiltonian

  `H = ∑ᵢ hᵢ,ᵢ₊₁ + ∑ᵢ oᵢ - ∑_α u_L(α) (C_α)₀ - ∑_α u_R(α) (C_α)ₖ`

(`hasDerivAt_openStep`), so `V(β/M)^M → e^{-βH}` (`tendsto_openStep_pow`). Exchanging the
space and time products as for the ring gives

  `Tr V(τ)^m = ⟨v_L, X^{k+1} v_R⟩`,  `v_u(γ) = ∏ₜ ω_u(γ t)`,  `ω_u(none) = 1`, `ω_u(α) = √τ u_α`

(`trace_openStep_pow`), with the bilinear pairing; the half weights of the two caps go into
the boundary vectors. The interior bonds are those of the ring, so `X` is the same matrix.

## Main results

* `trace_openStep_pow`: the transfer form of the open chain.
* `hasDerivAt_openStep`, `tendsto_openStep_pow`: the Gibbs limit of the time step.
* `tendsto_openPairing`: the pairings `⟨v_L, X_M^{k+1} v_R⟩` of the transfer matrices of `M`
  steps of width `β / M` converge to the open-chain partition function.
-/

open Matrix Filter Topology NormedSpace
open scoped Kronecker

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d]

/-! ### Boundary vectors -/

/-- The weight of one time label in a boundary vector: `1` for `none`, `√τ u_α` for `α`. -/
noncomputable def capWeight (τ : ℝ) (u : ι → ℂ) (a : Option ι) : ℂ :=
  a.elim 1 fun α ↦ (√τ : ℂ) * u α

/-- The boundary vector `v_u(γ) = ∏ₜ ω_u(γ t)` of a cap with weights `u`. -/
noncomputable def capVec {m : ℕ} (τ : ℝ) (u : ι → ℂ) (γ : Fin m → Option ι) : ℂ :=
  ∏ t, capWeight τ u (γ t)

/-- The conjugate of a boundary vector is the boundary vector of the conjugate weights. -/
theorem star_capVec {m : ℕ} (τ : ℝ) (u : ι → ℂ) :
    star (capVec (m := m) τ u) = capVec τ fun α ↦ star (u α) := by
  funext γ
  simp only [Pi.star_apply, capVec, star_prod]
  refine Finset.prod_congr rfl fun t _ ↦ ?_
  rcases γ t with _ | α
  · simp [capWeight]
  · simp [capWeight, Complex.conj_ofReal]

/-- The coefficient `κⱼ(a)` of a label on the bond `j` of the chain on `k + 1` sites. -/
def capCoeff {k : ℕ} (uL uR : ι → ℂ) (j : Fin (k + 2)) : Option ι → ℂ
  | none => 1
  | some α => if (j : ℕ) = 0 then uL α else if (j : ℕ) = k + 1 then uR α else 1

/-! ### The time step -/

section Step

variable [Fintype ι] (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

/-- The bond part of the open time step on `k + 1` sites. -/
noncomputable def openBondStep (k : ℕ) (uL uR : ι → ℂ) (u : ℂ) :
    Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ :=
  ∑ e : Fin (k + 2) → Option ι,
    (∏ j, labelWeight u (e j) * capCoeff uL uR j (e j)) •
      piTensor fun i : Fin (k + 1) ↦ (symFactor C (e i.castSucc) (e i.succ)).map (↑)

/-- The open time step `V(u)`: the bond part times the on-site part. -/
noncomputable def openStep (k : ℕ) (uL uR : ι → ℂ) (u : ℂ) :
    Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ :=
  openBondStep C k uL uR u * siteStep o (k + 1) u

/-- The end term `-∑_α u_α C_α` of a cap. -/
noncomputable def capTerm (u : ι → ℂ) : Matrix d d ℂ := -∑ α, u α • bondC C α

/-- The open chain on `k + 1` sites: the plain bonds, the on-site terms and the two end terms. -/
noncomputable def openHamiltonian (k : ℕ) (uL uR : ι → ℂ) :
    Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ :=
  ∑ i : Fin k, bondOp i.castSucc i.succ (signedBond C fun _ ↦ 1) +
    siteOp 0 (capTerm C uL) + siteOp (Fin.last k) (capTerm C uR) + siteSum (o.map (↑)) (k + 1)

/-- The tensor product `⨂ᵢ M(e i, e (i + 1)) (1 - u o)` of a label configuration. -/
noncomputable def openTerm {k : ℕ} (u : ℂ) (e : Fin (k + 2) → Option ι) :
    Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ :=
  piTensor fun i : Fin (k + 1) ↦ siteFactorC C o u (e i.castSucc) (e i.succ)

/-- The open time step as one sum over label configurations. -/
theorem openStep_eq_sum (k : ℕ) (uL uR : ι → ℂ) (u : ℂ) :
    openStep C o k uL uR u = ∑ e : Fin (k + 2) → Option ι,
      (∏ j, labelWeight u (e j) * capCoeff uL uR j (e j)) • openTerm C o u e := by
  simp only [openStep, openBondStep, siteStep, Finset.sum_mul, smul_mul_assoc, piTensor_mul,
    openTerm, siteFactorC]

end Step

/-! ### The transfer form -/

section TransferForm

variable [Fintype ι] (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- A label weight at real `τ ≥ 0` is the square of the half weight. -/
theorem labelWeight_ofReal {τ : ℝ} (hτ : 0 ≤ τ) (a : Option ι) :
    labelWeight (τ : ℂ) a = (((a.elim 1 fun _ ↦ √τ : ℝ) ^ 2 : ℝ) : ℂ) := by
  rcases a with _ | α
  · simp [labelWeight]
  · simp [labelWeight, Real.sq_sqrt hτ]

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The weight of a cap label is the half weight times the boundary weight. -/
theorem labelWeight_mul_capCoeff_zero {k : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ)
    (a : Option ι) :
    labelWeight (τ : ℂ) a * capCoeff uL uR (0 : Fin (k + 2)) a =
      ((a.elim 1 fun _ ↦ √τ : ℝ) : ℂ) * capWeight τ uL a := by
  rcases a with _ | α
  · simp [labelWeight, capCoeff, capWeight]
  · have h : (τ : ℂ) = (√τ : ℂ) * √τ := by rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτ]
    simp only [labelWeight, capCoeff, capWeight, Option.elim, Fin.val_zero, ite_true, h]
    ring

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The weight of a label on the right cap. -/
theorem labelWeight_mul_capCoeff_last {k : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ)
    (a : Option ι) :
    labelWeight (τ : ℂ) a * capCoeff uL uR (Fin.last (k + 1)) a =
      ((a.elim 1 fun _ ↦ √τ : ℝ) : ℂ) * capWeight τ uR a := by
  rcases a with _ | α
  · simp [labelWeight, capCoeff, capWeight]
  · have h : (τ : ℂ) = (√τ : ℂ) * √τ := by rw [← Complex.ofReal_mul, Real.mul_self_sqrt hτ]
    simp only [labelWeight, capCoeff, capWeight, Option.elim, Fin.val_last, h,
      show k + 1 ≠ 0 by omega, ite_false, ite_true]
    ring

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- On an interior bond the coefficient is `1`. -/
theorem capCoeff_interior {k : ℕ} (uL uR : ι → ℂ) {j : Fin (k + 2)} (h0 : (j : ℕ) ≠ 0)
    (hl : (j : ℕ) ≠ k + 1) (a : Option ι) : capCoeff uL uR j a = 1 := by
  rcases a with _ | α <;> simp [capCoeff, h0, hl]

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The labels of a sequence on the left cap: its half weight times its boundary vector. -/
theorem prod_capCoeff_zero {k m : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ)
    (γ : Fin m → Option ι) :
    ∏ t, labelWeight (τ : ℂ) (γ t) * capCoeff uL uR (0 : Fin (k + 2)) (γ t) =
      ((halfWeight τ γ : ℝ) : ℂ) * capVec τ uL γ := by
  simp only [labelWeight_mul_capCoeff_zero hτ, Finset.prod_mul_distrib, halfWeight,
    Complex.ofReal_prod, capVec]

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The labels of a sequence on the right cap. -/
theorem prod_capCoeff_last {k m : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ)
    (γ : Fin m → Option ι) :
    ∏ t, labelWeight (τ : ℂ) (γ t) * capCoeff uL uR (Fin.last (k + 1)) (γ t) =
      ((halfWeight τ γ : ℝ) : ℂ) * capVec τ uR γ := by
  simp only [labelWeight_mul_capCoeff_last hτ, Finset.prod_mul_distrib, halfWeight,
    Complex.ofReal_prod, capVec]

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The labels of a sequence on an interior bond: its squared half weight. -/
theorem prod_capCoeff_interior {k m : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ) (l : Fin k)
    (γ : Fin m → Option ι) :
    ∏ t, labelWeight (τ : ℂ) (γ t) * capCoeff uL uR l.castSucc.succ (γ t) =
      ((halfWeight τ γ : ℝ) : ℂ) ^ 2 := by
  have h0 : ((l.castSucc.succ : Fin (k + 2)) : ℕ) ≠ 0 := by simp
  have hl : ((l.castSucc.succ : Fin (k + 2)) : ℕ) ≠ k + 1 := by simp; omega
  simp only [capCoeff_interior uL uR h0 hl, mul_one, labelWeight_ofReal hτ, halfWeight,
    ← Finset.prod_pow, Complex.ofReal_prod, Complex.ofReal_pow]

omit [Fintype ι] in
/-- **The per-path identity** behind `trace_openStep_pow`: the time-ordered weights of a
space-time label configuration `q j t` of the open chain equal the transfer-matrix weights of
the open path `q 0, …, q (k + 1)` between the two boundary vectors. -/
theorem openPath_eq {k m : ℕ} {τ : ℝ} (hτ : 0 ≤ τ) (uL uR : ι → ℂ)
    (q : Fin (k + 2) → Fin m → Option ι) :
    (∏ t, ∏ j, labelWeight (τ : ℂ) (q j t) * capCoeff uL uR j (q j t)) *
        ∏ i : Fin (k + 1), (siteTrace C o τ (q i.castSucc) (q i.succ) : ℂ) =
      capVec τ uL (q 0) * (∏ i : Fin (k + 1),
        ((transfer C o m τ (q i.castSucc) (q i.succ) : ℝ) : ℂ)) * capVec τ uR (q (Fin.last _)) := by
  set H : Fin (k + 2) → ℂ := fun j ↦ ((halfWeight τ (q j) : ℝ) : ℂ)
  rw [Finset.prod_comm, Fin.prod_univ_succ, Fin.prod_univ_castSucc, prod_capCoeff_zero hτ,
    Fin.succ_last, prod_capCoeff_last hτ]
  simp only [prod_capCoeff_interior hτ]
  simp only [transfer, of_apply, Complex.ofReal_mul, Finset.prod_mul_distrib]
  rw [Fin.prod_univ_succ (fun i : Fin (k + 1) ↦ ((halfWeight τ (q i.castSucc) : ℝ) : ℂ)),
    Fin.prod_univ_castSucc (fun i : Fin (k + 1) ↦ ((halfWeight τ (q i.succ) : ℝ) : ℂ))]
  simp only [Fin.castSucc_zero, Fin.succ_last, ← Fin.succ_castSucc, Finset.prod_pow]
  ring

/-- **Transfer form of the open chain.** On the chain of `k + 1` sites,
`Tr V(τ)^m = ⟨v_L, X^{k+1} v_R⟩` (bilinear pairing) with the transfer matrix `X` of `m` time
steps and the boundary vectors of the two caps. -/
theorem trace_openStep_pow [DecidableEq ι] (k m : ℕ) (uL uR : ι → ℂ) {τ : ℝ} (hτ : 0 ≤ τ) :
    trace (openStep C o k uL uR τ ^ m) =
      capVec τ uL ⬝ᵥ ((transfer C o m τ ^ (k + 1)).map (↑) *ᵥ capVec τ uR) := by
  have hpow : (transfer C o m τ ^ (k + 1)).map (Complex.ofReal) =
      (List.ofFn fun _ : Fin (k + 1) ↦ (transfer C o m τ).map Complex.ofReal).prod := by
    rw [List.ofFn_const, List.prod_replicate]
    exact map_pow Complex.ofRealHom.mapMatrix _ _
  rw [hpow, dotProduct_list_prod_mulVec, ← List.prod_replicate, ← List.ofFn_const,
    openStep_eq_sum, list_prod_ofFn_sum, trace_sum]
  refine Fintype.sum_equiv (Equiv.piComm fun (_ : Fin m) (_ : Fin (k + 2)) ↦ Option ι) _ _
    fun E ↦ ?_
  rw [list_prod_ofFn_smul, trace_smul]
  simp only [openTerm, list_prod_ofFn_piTensor, trace_piTensor, trace_prod_siteFactorC,
    smul_eq_mul, map_apply]
  exact openPath_eq C o hτ uL uR (Equiv.piComm _ E)

end TransferForm

/-! ### First order -/

section FirstOrder

omit [Fintype d] [DecidableEq d] in
/-- The first-order coefficient of a configuration in the bond `j`, for per-bond coefficients
`cf` with `cf j none = 1`: nonzero exactly for the single-label configurations on `j`. -/
theorem firstOrder_eq_coeff [DecidableEq ι] [Fintype ι] {L : ℕ} (cf : Fin L → Option ι → ℂ)
    (hcf : ∀ j, cf j none = 1) (j : Fin L) (e : Fin L → Option ι) :
    (∏ k ∈ Finset.univ.erase j, labelWeight 0 (e k) * cf k (e k)) *
        (labelWeightDeriv (e j) * cf j (e j)) =
      ∑ α, if e = singleLabel j α then cf j (some α) else 0 := by
  by_cases h : ∀ k, k ≠ j → e k = none
  · have hprod : ∏ k ∈ Finset.univ.erase j, labelWeight 0 (e k) * cf k (e k) = 1 :=
      Finset.prod_eq_one fun k hk ↦ by
        rw [h k (Finset.ne_of_mem_erase hk), hcf]; simp [labelWeight]
    rw [hprod, one_mul]
    have hsingle : ∀ α, e = singleLabel j α ↔ e j = some α := fun α ↦ by
      constructor
      · rintro rfl; simp [singleLabel]
      · intro hj; ext1 k
        by_cases hk : k = j
        · subst hk; simp [singleLabel, hj]
        · simp [singleLabel, hk, h k hk]
    simp only [hsingle]
    rcases e j with _ | α₀
    · simp [labelWeightDeriv]
    · simp [labelWeightDeriv]
  · push Not at h
    obtain ⟨k, hkj, hk⟩ := h
    obtain ⟨β, hβ⟩ := Option.ne_none_iff_exists'.1 hk
    rw [Finset.prod_eq_zero (Finset.mem_erase.2 ⟨hkj, Finset.mem_univ k⟩)
      (by rw [hβ]; simp [labelWeight]), zero_mul]
    refine (Finset.sum_eq_zero fun α _ ↦ ?_).symm
    rw [ite_eq_right]
    rintro rfl
    simp [singleLabel, hkj] at hβ

variable (C : ι → Matrix d d ℝ)

omit [Fintype d] in
/-- A single label is `none` off its bond. -/
theorem singleLabel_of_ne {L : ℕ} {j k : Fin L} (h : k ≠ j) (α : ι) :
    singleLabel j α k = none := by
  simp [singleLabel, Function.update_of_ne h]

omit [Fintype d] in
/-- A single label on its bond. -/
theorem singleLabel_self {L : ℕ} (j : Fin L) (α : ι) : singleLabel j α j = some α := by
  simp [singleLabel]

/-- A single label on the left cap gives `C_α` on the site `0`. -/
theorem openTensor_single_zero {k : ℕ} (α : ι) :
    (piTensor fun i : Fin (k + 1) ↦ (symFactor C (singleLabel (0 : Fin (k + 2)) α i.castSucc)
      (singleLabel (0 : Fin (k + 2)) α i.succ)).map Complex.ofReal) = siteOp 0 (bondC C α) := by
  rw [siteOp]
  congr 1
  ext1 i
  rw [singleLabel_of_ne (Fin.succ_ne_zero i)]
  by_cases hi : i = 0
  · subst hi
    rw [Fin.castSucc_zero, singleLabel_self, symFactor_some_none]
    simp only [↓reduceIte]; rfl
  · rw [singleLabel_of_ne (fun h ↦ hi (Fin.castSucc_injective _ (h.trans Fin.castSucc_zero.symm))),
      symFactor_none_none, Matrix.map_one _ Complex.ofReal_zero Complex.ofReal_one]
    simp only [hi, ↓reduceIte]

/-- A single label on the right cap gives `C_α` on the site `k`. -/
theorem openTensor_single_last {k : ℕ} (α : ι) :
    (piTensor fun i : Fin (k + 1) ↦ (symFactor C (singleLabel (Fin.last (k + 1)) α i.castSucc)
      (singleLabel (Fin.last (k + 1)) α i.succ)).map Complex.ofReal) =
        siteOp (Fin.last k) (bondC C α) := by
  rw [siteOp]
  congr 1
  ext1 i
  rw [singleLabel_of_ne (Fin.castSucc_ne_last i)]
  by_cases hi : i = Fin.last k
  · subst hi
    rw [Fin.succ_last, singleLabel_self, symFactor_none_some]
    simp only [↓reduceIte]; rfl
  · rw [singleLabel_of_ne (fun h ↦ hi (Fin.succ_injective _ (h.trans (Fin.succ_last k).symm))),
      symFactor_none_none, Matrix.map_one _ Complex.ofReal_zero Complex.ofReal_one]
    simp only [hi, ↓reduceIte]

/-- A single label on the interior bond `l + 1` gives `C_α ⊗ C_α` on the sites `l, l + 1`. -/
theorem openTensor_single_interior {k : ℕ} (l : Fin k) (α : ι) :
    (piTensor fun i : Fin (k + 1) ↦ (symFactor C (singleLabel l.castSucc.succ α i.castSucc)
      (singleLabel l.castSucc.succ α i.succ)).map Complex.ofReal) =
        bondOp l.castSucc l.succ (bondC C α ⊗ₖ bondC C α) := by
  have hne : (l.castSucc : Fin (k + 1)) ≠ l.succ := by
    intro h; have := congrArg Fin.val h; simp at this
  rw [bondOp_kronecker hne]
  congr 1
  ext1 i
  by_cases h1 : i = l.castSucc
  · subst h1
    rw [singleLabel_of_ne (by intro h; have := congrArg Fin.val h; simp at this),
      singleLabel_self, symFactor_none_some]
    simp only [↓reduceIte]; rfl
  · by_cases h2 : i = l.succ
    · subst h2
      rw [← Fin.succ_castSucc, singleLabel_self,
        singleLabel_of_ne (by intro h; have := congrArg Fin.val h; simp at this),
        symFactor_some_none]
      simp only [h1, ↓reduceIte]; rfl
    · rw [singleLabel_of_ne (j := l.castSucc.succ) (k := i.castSucc)
          (fun h ↦ h2 (Fin.castSucc_injective _ (h.trans (Fin.succ_castSucc l)))),
        singleLabel_of_ne (fun h ↦ h1 (Fin.succ_injective _ h)), symFactor_none_none,
        Matrix.map_one _ Complex.ofReal_zero Complex.ofReal_one]
      simp only [h1, h2, ↓reduceIte]

end FirstOrder


/-! ### The Gibbs limit -/

section Limit

omit [Fintype d] in
/-- `A ↦ Aⱼ` is linear: a combination of one-site operators. -/
theorem siteOp_sum_smul {κ : Type*} [Fintype κ] {L : ℕ} (j : Fin L) (c : κ → ℂ)
    (A : κ → Matrix d d ℂ) : siteOp j (∑ x, c x • A x) = ∑ x, c x • siteOp j (A x) := by
  simp only [siteOp]
  exact (sum_smul_piTensor_single j c A).symm

omit [Fintype d] in
/-- `(-A)ⱼ = -Aⱼ`. -/
theorem siteOp_neg {L : ℕ} (j : Fin L) (A : Matrix d d ℂ) : siteOp j (-A) = -siteOp j A := by
  have := siteOp_sum_smul (κ := Unit) j (fun _ ↦ (-1 : ℂ)) (fun _ ↦ A)
  simpa using this

variable [Fintype ι] (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

omit [Fintype d] in
/-- A cap term on one site, label by label. -/
theorem siteOp_capTerm {L : ℕ} (j : Fin L) (u : ι → ℂ) :
    siteOp j (capTerm C u) = -∑ α, u α • siteOp j (bondC C α) := by
  rw [capTerm, siteOp_neg, siteOp_sum_smul]

/-- The first-order part of the open bond step is minus the bonds and the two end terms. -/
theorem openFirstOrder_sum (k : ℕ) (uL uR : ι → ℂ) :
    ∑ e : Fin (k + 2) → Option ι, (∑ j, (∏ l ∈ Finset.univ.erase j,
        labelWeight 0 (e l) * capCoeff uL uR l (e l)) *
          (labelWeightDeriv (e j) * capCoeff uL uR j (e j))) •
        (piTensor fun i : Fin (k + 1) ↦
          (symFactor C (e i.castSucc) (e i.succ)).map Complex.ofReal) =
      -(∑ i : Fin k, bondOp i.castSucc i.succ (signedBond C fun _ ↦ 1) +
        siteOp 0 (capTerm C uL) + siteOp (Fin.last k) (capTerm C uR)) := by
  classical
  simp only [firstOrder_eq_coeff (capCoeff uL uR) (fun _ ↦ rfl), Finset.sum_smul, ite_smul,
    zero_smul]
  rw [Finset.sum_comm]
  have hj : ∀ j : Fin (k + 2), (∑ e : Fin (k + 2) → Option ι, ∑ α,
      if e = singleLabel j α then capCoeff uL uR j (some α) •
        (piTensor fun i : Fin (k + 1) ↦
          (symFactor C (e i.castSucc) (e i.succ)).map Complex.ofReal) else 0) =
      ∑ α, capCoeff uL uR j (some α) • (piTensor fun i : Fin (k + 1) ↦
        (symFactor C (singleLabel j α i.castSucc) (singleLabel j α i.succ)).map
          Complex.ofReal) := fun j ↦ by
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun α _ ↦ by rw [Finset.sum_ite_eq']; simp
  simp only [hj]
  rw [Fin.sum_univ_succ, Fin.sum_univ_castSucc]
  simp only [openTensor_single_zero, openTensor_single_interior, Fin.succ_last,
    openTensor_single_last, siteOp_capTerm, signedBond, bondOp_neg, bondOp_sum, bondOp_smul]
  have h0 : ∀ α, capCoeff uL uR (0 : Fin (k + 2)) (some α) = uL α := fun α ↦ by simp [capCoeff]
  have hl : ∀ α, capCoeff uL uR (Fin.last (k + 1)) (some α) = uR α := fun α ↦ by
    simp [capCoeff]
  have hi : ∀ (l : Fin k) α, capCoeff uL uR l.castSucc.succ (some α) = 1 := fun l α ↦
    capCoeff_interior uL uR (by simp) (by simp; omega) _
  simp only [h0, hl, hi, one_smul, Complex.ofReal_one, Finset.sum_neg_distrib, neg_neg,
    neg_add_rev]
  abel

/-- The open bond step is `1` at `u = 0`. -/
theorem openBondStep_zero (k : ℕ) (uL uR : ι → ℂ) : openBondStep C k uL uR 0 = 1 := by
  rw [openBondStep, Finset.sum_eq_single (fun _ ↦ none)]
  · simp [labelWeight, capCoeff, Matrix.map_one, piTensor_one]
  · intro e _ he
    obtain ⟨j, hj⟩ := Function.ne_iff.1 he
    obtain ⟨β, hβ⟩ := Option.ne_none_iff_exists'.1 hj
    rw [Finset.prod_eq_zero (Finset.mem_univ j) (by rw [hβ]; simp [labelWeight]), zero_smul]
  · simp

section Derivative

attribute [local instance] Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- **Generator of the open time step.** `V'(0) = -H` for the open chain with end fields. -/
theorem hasDerivAt_openStep (k : ℕ) (uL uR : ι → ℂ) :
    HasDerivAt (fun u : ℂ ↦ openStep C o k uL uR u) (-openHamiltonian C o k uL uR) 0 := by
  classical
  have hb := hasDerivAt_sum_prod_smul (𝔸 := Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ)
    (fun e j u ↦ labelWeight u (e j) * capCoeff uL uR j (e j))
    (fun e j ↦ labelWeightDeriv (e j) * capCoeff uL uR j (e j))
    (fun e j ↦ hasDerivAt_labelWeight (e j) _ 0)
    (fun e : Fin (k + 2) → Option ι ↦ piTensor fun i : Fin (k + 1) ↦
      (symFactor C (e i.castSucc) (e i.succ)).map Complex.ofReal)
  rw [openFirstOrder_sum C k uL uR] at hb
  have hs : HasDerivAt (fun u : ℂ ↦ siteStep o (k + 1) u)
      (List.ofFn fun j : Fin (k + 1) ↦ -siteOp j (o.map Complex.ofReal)).sum 0 := by
    simp only [siteStep_eq_prod]
    exact hasDerivAt_list_prod_one_add_smul _
  have := hb.mul hs
  rw [siteStep_zero] at this
  have hb0 : (∑ e : Fin (k + 2) → Option ι,
      (∏ j, labelWeight 0 (e j) * capCoeff uL uR j (e j)) • piTensor fun i : Fin (k + 1) ↦
        (symFactor C (e i.castSucc) (e i.succ)).map Complex.ofReal) = 1 :=
    openBondStep_zero C k uL uR
  rw [hb0, mul_one, one_mul, List.sum_ofFn, Finset.sum_neg_distrib, ← neg_add] at this
  have e : openHamiltonian C o k uL uR = ∑ i : Fin k, bondOp i.castSucc i.succ
      (signedBond C fun _ ↦ 1) + siteOp 0 (capTerm C uL) + siteOp (Fin.last k) (capTerm C uR) +
        ∑ j : Fin (k + 1), siteOp j (o.map Complex.ofReal) := rfl
  rw [e]
  exact this

end Derivative

/-- **Gibbs limit of the open time step.** `V(β/M)^M → exp (-β H)`. -/
theorem tendsto_openStep_pow (k : ℕ) (uL uR : ι → ℂ) (β : ℝ) :
    Tendsto (fun M : ℕ ↦ openStep C o k uL uR ((β / M : ℝ) : ℂ) ^ M) atTop
      (𝓝 (exp ((-(β : ℂ)) • openHamiltonian C o k uL uR))) := by
  rcases isEmpty_or_nonempty (Fin (k + 1) → d) with hn | hn
  · exact tendsto_const_nhds.congr fun _ ↦ Subsingleton.elim _ _
  let _ : NormedRing (Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ) := Matrix.linftyOpNormedRing
  let _ : NormedAlgebra ℂ (Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ) :=
    Matrix.linftyOpNormedAlgebra
  let _ : NormedAlgebra ℚ (Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ) :=
    Matrix.linftyOpNormedAlgebra
  have h2 := (hasDerivAt_openStep C o k uL uR).scomp_of_eq (0 : ℂ)
    ((hasDerivAt_id (0 : ℂ)).mul_const (β : ℂ)) (by simp)
  have := tendsto_pow_of_hasDerivAt (𝔸 := Matrix (Fin (k + 1) → d) (Fin (k + 1) → d) ℂ)
    (g := fun u : ℂ ↦ openStep C o k uL uR (u * β))
    (A := (-(β : ℂ)) • openHamiltonian C o k uL uR)
    (by simp [openStep, openBondStep_zero, siteStep_zero])
    (by simpa [Function.comp_def, smul_neg, neg_smul] using h2)
  refine this.congr fun M ↦ ?_
  push_cast
  rw [div_eq_inv_mul]

/-- **The open-chain pairings converge to the Gibbs trace.** For `β ≥ 0`,
`Re ⟨v_L, X_M^{k+1} v_R⟩ → Re Tr e^{-βH}` for the transfer matrix `X_M` of `M` steps of width
`β / M`. -/
theorem tendsto_openPairing [DecidableEq ι] (k : ℕ) (uL uR : ι → ℂ) {β : ℝ} (hβ : 0 ≤ β) :
    Tendsto (fun M : ℕ ↦ (capVec (β / M) uL ⬝ᵥ
        ((transfer C o M (β / M) ^ (k + 1)).map (↑) *ᵥ capVec (β / M) uR)).re) atTop
      (𝓝 (gibbsTrace (openHamiltonian C o k uL uR) β)) := by
  have h := (Complex.continuous_re.comp (continuous_id.matrix_trace)).continuousAt.tendsto.comp
    (tendsto_openStep_pow C o k uL uR β)
  refine h.congr fun M ↦ ?_
  simp only [Function.comp_apply, id]
  rw [trace_openStep_pow C o k M uL uR (div_nonneg hβ (Nat.cast_nonneg M))]

end Limit

end CertifiedGaps.SingleSite
