/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.ChainCheck
import CertifiedGaps.Ring.Transfer
import CertifiedGaps.Transfer.MomentForm
import CertifiedGaps.Transfer.ProductLimit

/-!
# The two-purity criterion for reflection-positive rings

For every reflection-positive ring (`RPRing`) the partition functions `Re Tr e^{-βH_L}` satisfy
the spatial moment inequalities above length `2` (`RPRing.spatialMomentBounds`): the Trotter
approximants factor through a positive semidefinite transfer matrix
(`RPRing.exists_posSemidef_transfer`), so they have the moment form, and they converge to the
partition function (`tendsto_list_prod_one_add_smul_pow`). The inequalities are closed, so they
pass to the limit.

Combined with a rational chain certificate this gives `RPRing.gapGT_of_chainCertificate`: the
gap bound `log r / β₀` on every even ring `L ≥ n₀` from the two seed inequalities alone.

## Main results

* `RPRing.spatialMomentBounds`: the moment inequalities at every `β > 0`.
* `RPRing.gapGT_of_chainCertificate`: the uniform gap from seeds and a valid certificate.
-/

open Filter Matrix
open scoped ComplexOrder

namespace CertifiedGaps

namespace RPRing

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- The Trotter approximants have the moment form, for every number of slices. -/
theorem hasMomentForm_trotterTrace (M : ℕ) {β : ℝ} (hβ : 0 < β) :
    HasMomentForm (R.trotterTrace M β) 2 := by
  obtain ⟨D, _, _, T, hT, htr⟩ := R.exists_posSemidef_transfer M (β / M) (by positivity)
  exact hasMomentForm_of_trace_pow hT fun N hN ↦ by rw [trotterTrace, htr N hN]

omit [Fintype d] [DecidableEq d] in
/-- The even and odd bonds together are all bonds. -/
theorem sum_evenBonds_append_oddBonds {M : Type*} [AddCommMonoid M] (L : ℕ) (f : Fin L → M) :
    ((evenBonds L ++ oddBonds L).map f).sum = ∑ j, f j := by
  have hperm : (evenBonds L ++ oddBonds L).Perm (List.finRange L) := by
    unfold evenBonds oddBonds
    simpa using List.filter_append_perm (fun j : Fin L ↦ decide (Even j.val)) (List.finRange L)
  rw [(hperm.map f).sum_eq, Fin.sum_univ_def]

/-- The Trotter approximants converge to the partition function. -/
theorem tendsto_trotterTrace (β : ℝ) (L : ℕ) :
    Tendsto (fun M : ℕ ↦ R.trotterTrace M β L) atTop (nhds (R.partitionFunction L β)) := by
  let hs : List (Matrix (Fin L → d) (Fin L → d) ℂ) :=
    (evenBonds L).map (fun j ↦ -R.bond L j) ++ List.ofFn (fun i ↦ -R.site L i) ++
      (oddBonds L).map (fun j ↦ -R.bond L j)
  have hslice : ∀ τ : ℝ, R.trotterSlice L τ = (hs.map fun h ↦ 1 + (τ : ℂ) • h).prod := by
    intro τ
    simp only [hs, trotterSlice, List.map_append, List.prod_append, List.map_map, List.map_ofFn,
      Function.comp_def, smul_neg, ← sub_eq_add_neg]
  have hsum : hs.sum = -R.hamiltonian L := by
    have h1 := sum_evenBonds_append_oddBonds L fun j ↦ -R.bond L j
    simp only [List.map_append, List.sum_append] at h1
    simp only [hs, List.sum_append, List.sum_ofFn]
    rw [show ∀ a b c : Matrix (Fin L → d) (Fin L → d) ℂ, a + b + c = (a + c) + b from
      fun a b c ↦ by abel, h1, Finset.sum_neg_distrib, Finset.sum_neg_distrib, hamiltonian,
      neg_add]
  have h := tendsto_list_prod_one_add_smul_pow hs β
  rw [hsum, smul_neg, ← neg_smul] at h
  have hc : Continuous fun A : Matrix (Fin L → d) (Fin L → d) ℂ ↦ (trace A).re :=
    Complex.continuous_re.comp (continuous_id.matrix_trace)
  have := (hc.tendsto _).comp h
  simpa [Function.comp_def, trotterTrace, hslice, partitionFunction] using this

/-- **The transfer-operator input.** For every `β > 0`, the traces `Re Tr e^{-βH_L}` satisfy the
spatial moment inequalities above length `2`. -/
theorem spatialMomentBounds (β : ℝ) (hβ : 0 < β) :
    SpatialMomentBounds (fun L ↦ R.partitionFunction L β) 2 :=
  SpatialMomentBounds.of_tendsto_hasMomentForm (Zs := fun M ↦ R.trotterTrace M β) two_pos
    (Eventually.of_forall fun M ↦ R.hasMomentForm_trotterTrace M hβ)
    fun L _ _ ↦ R.tendsto_trotterTrace β L

/-- The moment inequalities in the form of the abstract criterion, for the spectra `R.spectrum`. -/
theorem hasSpatialMomentBounds :
    HasSpatialMomentBounds (ι := fun L ↦ Fin L → d) R.spectrum 2 := fun β hβ ↦ by
  have h := R.spatialMomentBounds β hβ
  simp only [RPRing.partitionFunction_eq_sum] at h
  exact h

/-- The spatial purity `P_s(n, β) = Z_{2n}(β) / Z_n(β)²` of the ring. -/
noncomputable def spatialPurity (n : ℕ) (β : ℝ) : ℝ :=
  R.partitionFunction (2 * n) β / R.partitionFunction n β ^ 2

/-- The thermal purity `P_τ(L, β) = Z_L(2β) / Z_L(β)²` of the ring. -/
noncomputable def thermalPurity (L : ℕ) (β : ℝ) : ℝ :=
  R.partitionFunction L (2 * β) / R.partitionFunction L β ^ 2

/-- **The two-purity criterion for a reflection-positive ring.** Let `n₀ ≥ 2` be even, `β₀ > 0`,
and `c` a valid rational chain certificate with seed defects `s = c.sAt 0`, `t = c.tAt 0` and
rate `r = c.r`. If `P_s(n₀, β₀) ≥ 1 - s` and `P_τ(2n₀, β₀) ≥ 1 - t`, then every even ring of
`L ≥ n₀` cells has a nondegenerate ground state and gap `Δ_L > log r / β₀`. The defects and the
rate are passed as real numbers, so that the seeds can be stated with decimal literals. -/
theorem gapGT_of_chainCertificate [Nonempty d] {n₀ : ℕ} (hn₀ : 2 ≤ n₀) (hn₀e : Even n₀) {β₀ : ℝ}
    (hβ₀ : 0 < β₀) {c : ChainCertificate} (hc : c.Valid) {s t r : ℝ}
    (hcert : (c.sAt 0 : ℝ) = s ∧ (c.tAt 0 : ℝ) = t ∧ (c.r : ℝ) = r)
    (hs : 1 - s ≤ R.spatialPurity n₀ β₀) (ht : 1 - t ≤ R.thermalPurity (2 * n₀) β₀) {L : ℕ}
    (hL : n₀ ≤ L) (hLe : Even L) : GapGT (R.spectrum L) (Real.log r / β₀) := by
  obtain ⟨rfl, rfl, rfl⟩ := hcert
  refine CertifiedGaps.gapGT_of_chainCertificate (ι := fun L ↦ Fin L → d) R.spectrum
    (by omega) hn₀e hβ₀ (fun β hβ ↦ (R.hasSpatialMomentBounds β hβ).mono hn₀) c hc ?_ ?_ hL hLe
  · simpa only [CertifiedGaps.spatialPurity, CertifiedGaps.partitionFunction,
      spatialPurity, partitionFunction_eq_sum] using hs
  · simpa only [CertifiedGaps.thermalPurity, CertifiedGaps.partitionFunction,
      thermalPurity, partitionFunction_eq_sum] using ht

end RPRing

end CertifiedGaps
