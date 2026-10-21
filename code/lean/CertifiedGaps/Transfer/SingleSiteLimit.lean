/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Model
import CertifiedGaps.Transfer.ProductLimit
import CertifiedGaps.Transfer.SingleSite

/-!
# The single-site time step converges to the Gibbs operator

The time step `V(u) = stepOp C o L r π u` of `Transfer/SingleSite.lean` is a polynomial in `u`
with `V(0) = 1`. Its bond part has first-order term from the configurations with exactly one
nontrivial label: a label `α` on the bond `j - 1 → j` puts `C_α` on the left end of that bond and
`C_{πⱼ α}` on its right end, so the bond part contributes `∑ⱼ ∑_α sⱼ(α) C_α ⊗ C_{πⱼ α}`. The
on-site part `⨂ⱼ (1 - u o)` contributes `-∑ⱼ oⱼ`. Hence `V'(0) = -H_L(r, π)` with the ring
Hamiltonian

  `H_L(r, π) = ∑ⱼ (-∑_α sⱼ₊₁(α) C_α ⊗ C_{πⱼ₊₁ α})ⱼ,ⱼ₊₁ + ∑ⱼ oⱼ`,

whose closing bond `L - 1 → 0` is `-∑_α r_α C_α ⊗ C_{π α}`: the generic
`twistedRingHamiltonian (signedBond C 1) (twistBond C r π) L + siteSum o L`. The Chernoff lemma
`CertifiedGaps.tendsto_pow_of_hasDerivAt` gives `V(β/M)^M → exp (-β H_L(r, π))` on every ring of
`L ≥ 2` sites, odd lengths included. Together with `trace_stepOp_pow`, every twisted partition
function `Z_L(r, π, β)` is a limit of `∑_γ D_r(γ) (X^L) (π γ, γ)` with one real symmetric matrix
`X` for all `L`.

The sign twists (`π = id`) have their own names, `hamiltonian C o L r` and
`partitionFunction C o L r β`, and the trace form `Tr (diag(D_r) X^L)`.
A reflection-positive ring `RPRing` without an internal cell reflection is the case `r = 1`.

## Main results

* `hasDerivAt_stepOp`: `V'(0) = -H_L(r, π)`.
* `tendsto_trace_stepOp_pow`: `Re Tr V(β/M)^M → Z_L(r, π, β)`.
* `eventually_abs_twisted_transfer_sub_le`: for any finite set of lengths `L ≥ 2`, eventually in
  `M`, `|∑_γ D_r(γ) (X_M^L) (π γ, γ) - Z_L(r, π, β)| ≤ δ`, with `X_M` the transfer matrix of `M`
  steps of width `β / M`.
* `eventually_abs_trace_transfer_sub_le`: the same for a sign twist,
  `|Tr (diag(D_r) X_M^L) - Z_L(r, β)| ≤ δ`.
* `RPRing.hamiltonian_eq_singleSite`: an `RPRing` with `refl = 1` has `H_L = H_L(1)`.
* `RPRing.exists_transfer_approx`: such a ring has approximate transfer spectra at all lengths:
  real numbers `xᵢ` with `|∑ᵢ xᵢ^L - Z_L(β)| ≤ δ` for finitely many lengths `L ≥ 2`.
-/

open Matrix Filter Topology NormedSpace
open scoped Kronecker

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι]
  (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

/-! ### The ring Hamiltonian with a twisted closing bond -/

/-- The complexified bond factor `C_α`. -/
noncomputable def bondC (α : ι) : Matrix d d ℂ := (C α).map (↑)

/-- The bond `-∑_α r_α C_α ⊗ C_α` with the signs `r`; the plain bond is `r = 1`. -/
noncomputable def signedBond (r : ι → ℝ) : Matrix (d × d) (d × d) ℂ :=
  -∑ α, (r α : ℂ) • (bondC C α ⊗ₖ bondC C α)

/-- The bond `-∑_α r_α C_α ⊗ C_{π α}` with the signs `r` and the right end relabeled by `π`. -/
noncomputable def twistBond (r : ι → ℝ) (π : ι → ι) : Matrix (d × d) (d × d) ℂ :=
  -∑ α, (r α : ℂ) • (bondC C α ⊗ₖ bondC C (π α))

omit [Fintype d] [DecidableEq d] in
/-- Without a relabeling the twisted bond is the signed bond. -/
theorem twistBond_id (r : ι → ℝ) : twistBond C r id = signedBond C r := rfl

/-- The ring Hamiltonian `H_L(r, π) = ∑ⱼ hⱼ,ⱼ₊₁ + ∑ⱼ oⱼ` whose closing bond `L - 1 → 0` is
`-∑_α r_α C_α ⊗ C_{π α}`: the generic twisted ring of the plain bond and the twisted bond. -/
noncomputable def twistHamiltonian (L : ℕ) (r : ι → ℝ) (π : ι → ι) :
    Matrix (Fin L → d) (Fin L → d) ℂ :=
  twistedRingHamiltonian (signedBond C fun _ ↦ 1) (twistBond C r π) L + siteSum (o.map (↑)) L

/-- The twisted partition function `Z_L(r, π, β) = Re Tr exp (-β H_L(r, π))`. -/
noncomputable def twistPartitionFunction (L : ℕ) (r : ι → ℝ) (π : ι → ι) (β : ℝ) : ℝ :=
  gibbsTrace (twistHamiltonian C o L r π) β

/-- The ring Hamiltonian `H_L(r) = ∑ⱼ hⱼ,ⱼ₊₁ + ∑ⱼ oⱼ` whose closing bond `L - 1 → 0` is signed
by `r`: the generic twisted ring of the plain bond and the signed bond. -/
noncomputable def hamiltonian (L : ℕ) (r : ι → ℝ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  twistedRingHamiltonian (signedBond C fun _ ↦ 1) (signedBond C r) L + siteSum (o.map (↑)) L

/-- The partition function `Z_L(r, β) = Re Tr exp (-β H_L(r))`. -/
noncomputable def partitionFunction (L : ℕ) (r : ι → ℝ) (β : ℝ) : ℝ :=
  gibbsTrace (hamiltonian C o L r) β

omit [Fintype d] in
/-- A sign twist is a twist without a relabeling. -/
theorem twistHamiltonian_id (L : ℕ) (r : ι → ℝ) :
    twistHamiltonian C o L r id = hamiltonian C o L r := rfl

/-- A sign twist is a twist without a relabeling. -/
theorem twistPartitionFunction_id (L : ℕ) (r : ι → ℝ) (β : ℝ) :
    twistPartitionFunction C o L r id β = partitionFunction C o L r β := rfl

/-- The bond whose right end is the site `k`, in the form the label expansion of the time step
produces it: `-∑_α s_k(α) C_α ⊗ C_{π_k α}`, with the signs `r` and the relabeling `π` on the
closing bond (`k = 0`). -/
noncomputable def bondTerm (r : ι → ℝ) (π : ι → ι) {L : ℕ} (k : Fin L) :
    Matrix (d × d) (d × d) ℂ :=
  -∑ α, (bondSign r k (some α) : ℂ) • (bondC C α ⊗ₖ bondC C (labelAt π k α))

omit [Fintype d] in
/-- Summed over the bonds, the label form `bondTerm` is the bond part of `H_L(r, π)`. -/
theorem sum_bondOp_bondTerm (L : ℕ) (r : ι → ℝ) (π : ι → ι) :
    ∑ j : Fin L, bondOp j (finRotate L j) (bondTerm C r π (finRotate L j)) =
      twistedRingHamiltonian (signedBond C fun _ ↦ 1) (twistBond C r π) L := by
  rw [twistedRingHamiltonian]
  refine Finset.sum_congr rfl fun j _ ↦ congrArg _ ?_
  have hj : ((finRotate L j : ℕ) = 0) ↔ (j : ℕ) + 1 = L := by
    rw [val_finRotate]; split_ifs with h <;> simp [h]
  simp only [bondTerm, signedBond, twistBond, bondSign, labelAt, labelSign, Option.elim, hj]
  split_ifs <;> rfl

/-! ### First order in the bond part -/

/-- The configuration with the single label `α` on the bond `j - 1 → j`. -/
def singleLabel {L : ℕ} (j : Fin L) (α : ι) : Fin L → Option ι :=
  Function.update (fun _ ↦ none) j (some α)

omit [Fintype ι] in
/-- A single label `α` on the bond `j - 1 → j` gives `C_α ⊗ C_{πⱼ α}` on that bond. -/
theorem piTensor_singleLabel {L : ℕ} (hL : 2 ≤ L) (π : ι → ι) (j : Fin L) (α : ι) :
    (piTensor fun k ↦ (symFactor C (bondLabel π k (singleLabel j α k))
        (singleLabel j α (finRotate L k))).map Complex.ofReal) =
      bondOp ((finRotate L).symm j) j (bondC C α ⊗ₖ bondC C (labelAt π j α)) := by
  have hne : (finRotate L).symm j ≠ j := by
    have := ne_finRotate hL ((finRotate L).symm j)
    rwa [Equiv.apply_symm_apply] at this
  have hrot : finRotate L j ≠ j := (ne_finRotate hL j).symm
  have h1 : (symFactor C none none).map Complex.ofReal = 1 := by simp [Matrix.map_one]
  rw [bondOp_kronecker hne]
  congr 1
  ext1 k
  by_cases hk1 : k = (finRotate L).symm j
  · subst hk1
    simp only [singleLabel, Equiv.apply_symm_apply, Function.update_apply, hne, ite_true,
      ite_false, bondLabel, Option.map_none, symFactor_none_some, bondC]
  · by_cases hk2 : k = j
    · subst hk2
      simp only [singleLabel, Function.update_apply, hrot, hk1, ite_true, ite_false, bondLabel,
        Option.map_some, symFactor_some_none, bondC]
    · have hk3 : finRotate L k ≠ j := fun h ↦ hk1 (by rw [← h, Equiv.symm_apply_apply])
      simp only [singleLabel, Function.update_apply, hk2, hk3, ite_false, bondLabel,
        Option.map_none, h1, hk1]

/-- The derivative at `u = 0` of the weight `w_u(a)` of a label. -/
def labelWeightDeriv (a : Option ι) : ℂ := a.elim 0 fun _ ↦ 1

omit [Fintype ι] in
/-- The label weight has derivative `labelWeightDeriv`. -/
theorem hasDerivAt_labelWeight (a : Option ι) (c : ℂ) (u : ℂ) :
    HasDerivAt (fun u ↦ labelWeight u a * c) (labelWeightDeriv a * c) u := by
  rcases a with _ | α
  · simpa [labelWeight, labelWeightDeriv] using hasDerivAt_const u c
  · simpa [labelWeight, labelWeightDeriv] using (hasDerivAt_id u).mul_const c

/-- The first-order coefficient of a configuration, in the bond `j`: nonzero exactly for the
single-label configurations on `j`. -/
theorem firstOrder_eq [DecidableEq ι] {L : ℕ} (r : ι → ℝ) (j : Fin L)
    (e : Fin L → Option ι) :
    (∏ k ∈ Finset.univ.erase j, labelWeight 0 (e k) * (bondSign r k (e k) : ℂ)) *
        (labelWeightDeriv (e j) * (bondSign r j (e j) : ℂ)) =
      ∑ α, if e = singleLabel j α then (bondSign r j (some α) : ℂ) else 0 := by
  by_cases h : ∀ k, k ≠ j → e k = none
  · have hprod : ∏ k ∈ Finset.univ.erase j, labelWeight 0 (e k) * (bondSign r k (e k) : ℂ) = 1 :=
      Finset.prod_eq_one fun k hk ↦ by
        rw [h k (Finset.ne_of_mem_erase hk)]; simp [labelWeight, bondSign, labelSign]
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
    rw [Finset.prod_eq_zero (Finset.mem_univ_val k |> fun _ ↦
      Finset.mem_erase.2 ⟨hkj, Finset.mem_univ k⟩) (by rw [hβ]; simp [labelWeight]), zero_mul]
    refine (Finset.sum_eq_zero fun α _ ↦ ?_).symm
    rw [ite_eq_right]
    rintro rfl
    simp [singleLabel, hkj] at hβ


/-- The first-order part of the bond step is minus the bond part of `H_L(r, π)`. -/
theorem firstOrder_sum {L : ℕ} (hL : 2 ≤ L) (r : ι → ℝ) (π : ι → ι) :
    ∑ e : Fin L → Option ι, (∑ j, (∏ k ∈ Finset.univ.erase j,
        labelWeight 0 (e k) * (bondSign r k (e k) : ℂ)) *
          (labelWeightDeriv (e j) * (bondSign r j (e j) : ℂ))) •
        (piTensor fun k ↦
          (symFactor C (bondLabel π k (e k)) (e (finRotate L k))).map Complex.ofReal) =
      -∑ j : Fin L, bondOp j (finRotate L j) (bondTerm C r π (finRotate L j)) := by
  classical
  simp only [firstOrder_eq, Finset.sum_smul, ite_smul, zero_smul]
  rw [Finset.sum_comm]
  rw [← (Equiv.sum_comp (finRotate L) _), ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [Finset.sum_comm, bondTerm, bondOp_neg, neg_neg, bondOp_sum]
  refine Finset.sum_congr rfl fun α _ ↦ ?_
  rw [Finset.sum_ite_eq', ite_eq_left (Finset.mem_univ _), bondOp_smul,
    piTensor_singleLabel C hL, Equiv.symm_apply_apply]

/-- The bond step is `1` at `u = 0`. -/
theorem bondStep_zero (L : ℕ) (r : ι → ℝ) (π : ι → ι) : bondStep C L r π 0 = 1 := by
  rw [bondStep, Finset.sum_eq_single (fun _ ↦ none)]
  · simp [stepCoeff, labelWeight, bondSign, labelSign, bondLabel, Matrix.map_one, piTensor_one]
  · intro e _ he
    obtain ⟨k, hk⟩ := Function.ne_iff.1 he
    obtain ⟨β, hβ⟩ := Option.ne_none_iff_exists'.1 hk
    rw [stepCoeff, Finset.prod_eq_zero (Finset.mem_univ k) (by rw [hβ]; simp [labelWeight]),
      zero_smul]
  · simp

/-! ### The on-site part -/

omit [Fintype d] [Fintype ι] in
/-- One factor `1 - u o` at the site `j` is `1 - u oⱼ` on the ring. -/
theorem siteOp_one_sub (L : ℕ) (j : Fin L) (u : ℂ) (A : Matrix d d ℂ) :
    siteOp j (1 - u • A) = 1 + u • (-siteOp j A) := by
  have h := sum_smul_piTensor_single j ![(1 : ℂ), -u] ![(1 : Matrix d d ℂ), A]
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, one_smul,
    ite_self, piTensor_one, neg_smul] at h
  rw [siteOp, sub_eq_add_neg, ← h, smul_neg]
  rfl

omit [Fintype ι] in
/-- The on-site part is the ordered product of its one-site factors. -/
theorem siteStep_eq_prod (L : ℕ) (u : ℂ) :
    siteStep o L u = ((List.ofFn fun j : Fin L ↦ -siteOp j (o.map (↑))).map
      fun h ↦ 1 + u • h).prod := by
  have h : (fun h ↦ 1 + u • h) ∘ (fun j : Fin L ↦ -siteOp j (o.map Complex.ofReal)) =
      fun j ↦ siteOp j (1 - u • o.map Complex.ofReal) :=
    funext fun j ↦ (siteOp_one_sub L j u _).symm
  rw [List.map_ofFn, h]
  simp only [siteOp, list_prod_ofFn_piTensor, siteStep]
  congr 1
  funext i
  exact (list_prod_ofFn_ite (fun k ↦ i = k) (fun _ ↦ 1 - u • o.map Complex.ofReal) i
    fun k ↦ eq_comm).symm

omit [Fintype d] [Fintype ι] in
/-- The on-site part is `1` at `u = 0`. -/
theorem siteStep_zero (L : ℕ) : siteStep o L 0 = 1 := by
  simp [siteStep, piTensor_one]

/-! ### The derivative and the Gibbs limit -/

section Derivative

attribute [local instance] Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra

/-- **Generator of the time step.** `V'(0) = -H_L(r, π)` on every ring of `L ≥ 2` sites. -/
theorem hasDerivAt_stepOp {L : ℕ} (hL : 2 ≤ L) (r : ι → ℝ) (π : ι → ι) :
    HasDerivAt (fun u : ℂ ↦ stepOp C o L r π u) (-twistHamiltonian C o L r π) 0 := by
  have hb := hasDerivAt_sum_prod_smul (𝔸 := Matrix (Fin L → d) (Fin L → d) ℂ)
    (fun e j u ↦ labelWeight u (e j) * (bondSign r j (e j) : ℂ))
    (fun e j ↦ labelWeightDeriv (e j) * (bondSign r j (e j) : ℂ))
    (fun e j ↦ hasDerivAt_labelWeight (e j) _ 0)
    (fun e : Fin L → Option ι ↦ piTensor fun k ↦
      (symFactor C (bondLabel π k (e k)) (e (finRotate L k))).map Complex.ofReal)
  rw [firstOrder_sum C hL r π, sum_bondOp_bondTerm] at hb
  have hs := hasDerivAt_list_prod_one_add_smul
    (List.ofFn fun j : Fin L ↦ -siteOp j (o.map Complex.ofReal))
  simp only [← siteStep_eq_prod] at hs
  have := hb.mul hs
  rw [siteStep_zero] at this
  have hb0 : (∑ e : Fin L → Option ι, (∏ j, labelWeight 0 (e j) * (bondSign r j (e j) : ℂ)) •
      piTensor fun k ↦
        (symFactor C (bondLabel π k (e k)) (e (finRotate L k))).map Complex.ofReal) = 1 :=
    bondStep_zero C L r π
  rw [hb0, mul_one, one_mul, List.sum_ofFn, Finset.sum_neg_distrib, ← neg_add] at this
  exact this

end Derivative

/-- **Gibbs limit of the time step.** `V(β/M)^M → exp (-β H_L(r, π))` on every ring of `L ≥ 2`
sites. -/
theorem tendsto_stepOp_pow {L : ℕ} (hL : 2 ≤ L) (r : ι → ℝ) (π : ι → ι) (β : ℝ) :
    Tendsto (fun M : ℕ ↦ stepOp C o L r π ((β / M : ℝ) : ℂ) ^ M) atTop
      (𝓝 (exp ((-(β : ℂ)) • twistHamiltonian C o L r π))) := by
  rcases isEmpty_or_nonempty (Fin L → d) with hn | hn
  · exact tendsto_const_nhds.congr fun _ ↦ Subsingleton.elim _ _
  let _ : NormedRing (Matrix (Fin L → d) (Fin L → d) ℂ) := Matrix.linftyOpNormedRing
  let _ : NormedAlgebra ℂ (Matrix (Fin L → d) (Fin L → d) ℂ) := Matrix.linftyOpNormedAlgebra
  let _ : NormedAlgebra ℚ (Matrix (Fin L → d) (Fin L → d) ℂ) := Matrix.linftyOpNormedAlgebra
  have h2 := (hasDerivAt_stepOp C o hL r π).scomp_of_eq (0 : ℂ)
    ((hasDerivAt_id (0 : ℂ)).mul_const (β : ℂ)) (by simp)
  have := tendsto_pow_of_hasDerivAt (𝔸 := Matrix (Fin L → d) (Fin L → d) ℂ)
    (g := fun u : ℂ ↦ stepOp C o L r π (u * β)) (A := (-(β : ℂ)) • twistHamiltonian C o L r π)
    (by simp [stepOp, bondStep_zero, siteStep_zero])
    (by simpa [Function.comp_def, smul_neg, neg_smul] using h2)
  refine this.congr fun M ↦ ?_
  push_cast
  rw [div_eq_inv_mul]

/-- **Gibbs limit of the trace.** `Re Tr V(β/M)^M → Z_L(r, π, β)` on every ring of `L ≥ 2`
sites. -/
theorem tendsto_trace_stepOp_pow {L : ℕ} (hL : 2 ≤ L) (r : ι → ℝ) (π : ι → ι) (β : ℝ) :
    Tendsto (fun M : ℕ ↦ (trace (stepOp C o L r π ((β / M : ℝ) : ℂ) ^ M)).re) atTop
      (𝓝 (twistPartitionFunction C o L r π β)) :=
  (Complex.continuous_re.comp (continuous_id.matrix_trace)).continuousAt.tendsto.comp
    (tendsto_stepOp_pow C o hL r π β)

/-- **Transfer form of the twisted partition functions.** For every finite set of ring lengths
`L ≥ 2` and every `δ > 0`, eventually in the number `M` of time steps,
`|∑_γ D_r(γ) (X^L) (π γ, γ) - Z_L(r, π, β)| ≤ δ` for the transfer matrix `X` of `M` steps of
width `β / M`. -/
theorem eventually_abs_twisted_transfer_sub_le [DecidableEq ι] (r : ι → ℝ) (π : ι → ι) {β : ℝ}
    (hβ : 0 ≤ β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ M : ℕ in atTop, ∀ n ∈ ns,
      |∑ γ, labelChar r γ * (transfer C o M (β / M) ^ n) (fun t ↦ (γ t).map π) γ -
        twistPartitionFunction C o n r π β| ≤ δ := by
  rw [Filter.eventually_all_finset]
  intro n hn
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by have := hns n hn; omega⟩
  have hlim := tendsto_trace_stepOp_pow C o (hns _ hn) r π β
  refine ((Metric.tendsto_nhds.1 hlim δ hδ).mono fun M hM ↦ ?_)
  rw [Real.dist_eq] at hM
  have hτ : 0 ≤ β / M := div_nonneg hβ (Nat.cast_nonneg M)
  rw [trace_stepOp_pow C o k M r π hτ, Complex.ofReal_re] at hM
  exact hM.le

/-- **Transfer form of the partition functions.** For every finite set of ring lengths
`L ≥ 2` and every `δ > 0`, eventually in the number `M` of time steps,
`|Tr (diag(D_r) X^L) - Z_L(r, β)| ≤ δ` for the transfer matrix `X` of `M` steps of width
`β / M`. -/
theorem eventually_abs_trace_transfer_sub_le [DecidableEq ι] (r : ι → ℝ) {β : ℝ} (hβ : 0 ≤ β)
    (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ M : ℕ in atTop, ∀ n ∈ ns,
      |trace (diagonal (labelChar r) * transfer C o M (β / M) ^ n) -
        partitionFunction C o n r β| ≤ δ := by
  refine (eventually_abs_twisted_transfer_sub_le C o r id hβ ns hns hδ).mono fun M hM n hn ↦ ?_
  have h := hM n hn
  simp only [Option.map_id, id_eq, twistPartitionFunction_id] at h
  simpa only [trace, diag_apply, diagonal_mul] using h

omit [Fintype d] [DecidableEq d] [Fintype ι] in
/-- The trivial sign vector has the trivial character. -/
theorem labelChar_one {m : ℕ} (γ : Fin m → Option ι) : labelChar (fun _ ↦ (1 : ℝ)) γ = 1 :=
  Finset.prod_eq_one fun t _ ↦ by rcases γ t with _ | α <;> rfl

end CertifiedGaps.SingleSite

/-! ### Reflection-positive rings without an internal reflection -/

namespace CertifiedGaps.RPRing

open SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- An `RPRing` without an internal reflection is the single-site ring with the plain sign
`r = 1`. -/
theorem hamiltonian_eq_singleSite (h : R.refl = 1) (L : ℕ) :
    R.hamiltonian L = SingleSite.hamiltonian R.C R.o L (fun _ ↦ 1) := by
  have hb : R.bondMatrix = signedBond R.C fun _ ↦ 1 := by
    simp only [bondMatrix, signedBond, Dc_eq_of_refl R h, Complex.ofReal_one, one_smul]
    rfl
  rw [hamiltonian_eq_ringHamiltonian, SingleSite.hamiltonian, twistedRingHamiltonian_self, hb]
  rfl

/-- **All-length transfer spectra of a reflection-positive ring.** For an `RPRing` without an
internal reflection, every `β ≥ 0`, every finite set of ring lengths `L ≥ 2` (odd lengths
included) and every `δ > 0`, there are real numbers `xᵢ` (the eigenvalues of one real symmetric
transfer matrix) with `|∑ᵢ xᵢ^L - Z_L(β)| ≤ δ` for every `L` in the set. -/
theorem exists_transfer_approx (h : R.refl = 1) {β : ℝ} (hβ : 0 ≤ β) (ns : Finset ℕ)
    (hns : ∀ n ∈ ns, 2 ≤ n) {δ : ℝ} (hδ : 0 < δ) :
    ∃ (D : Type) (_ : Fintype D) (x : D → ℝ), ∀ n ∈ ns,
      |∑ i, x i ^ n - R.partitionFunction n β| ≤ δ := by
  classical
  obtain ⟨M, hM⟩ :=
    (eventually_abs_trace_transfer_sub_le R.C R.o (fun _ ↦ 1) hβ ns hns hδ).exists
  have hX := transfer_isHermitian R.C R.o M (β / M)
  refine ⟨_, inferInstance, hX.eigenvalues, fun n hn ↦ ?_⟩
  have h1 := hM n hn
  rw [show labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from funext labelChar_one, diagonal_one,
    Matrix.one_mul, trace_pow_eq_sum_eigenvalues_real hX] at h1
  rwa [partitionFunction, hamiltonian_eq_singleSite R h]

end CertifiedGaps.RPRing
