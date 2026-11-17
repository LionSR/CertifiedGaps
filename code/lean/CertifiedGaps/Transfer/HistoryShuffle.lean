/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.TimeOrdering
import Mathlib.MeasureTheory.Integral.Prod

/-!
# Chronological interleavings of ordered histories

Finite families of precedence constraints select exactly the chronological chambers whose
permutations preserve those constraints. Specializing to the internal order of two bond
histories gives a genuine shuffle integral over the product of their ordered-time simplexes.
All identities are valid for Banach-valued integrands: no commutativity is required.
-/

open MeasureTheory Set

namespace CertifiedGaps.TimeOrdering

/-- The time cube with the strict precedence constraints specified by `r`. -/
def orderedTimeDomain (b : ℝ) {n : ℕ} (r : Fin n → Fin n → Prop) :
    Set (Fin n → ℝ) :=
  timeCube b n ∩ {t | ∀ i j, r i j → t i < t j}

/-- A chronological permutation preserves `r` when each required predecessor has an earlier
position. The inverse appears because a chamber permutation lists original event indices. -/
def PreservesTimeOrder {n : ℕ} (r : Fin n → Fin n → Prop)
    (σ : Equiv.Perm (Fin n)) : Prop :=
  ∀ i j, r i j → σ.symm i < σ.symm j

/-- An admissible chronological interleaving for the prescribed event order. -/
abbrev TimeShuffle {n : ℕ} (r : Fin n → Fin n → Prop) :=
  {σ : Equiv.Perm (Fin n) // PreservesTimeOrder r σ}
noncomputable instance {n : ℕ} (r : Fin n → Fin n → Prop) : Fintype (TimeShuffle r) :=
  Fintype.ofFinite _

/-- Finitely many strict precedence conditions are measurable. -/
theorem measurableSet_orderedTimeDomain (b : ℝ) {n : ℕ}
    (r : Fin n → Fin n → Prop) : MeasurableSet (orderedTimeDomain b r) := by
  refine (measurableSet_timeCube b n).inter ?_
  simp only [ofPred_forall]
  exact MeasurableSet.iInter fun i => MeasurableSet.iInter fun j =>
    MeasurableSet.iInter fun _ => measurableSet_lt (measurable_pi_apply i)
      (measurable_pi_apply j)

/-- On a collision-free chronological chamber, precedence of times is equivalent to precedence
of positions in its permutation. -/
theorem mem_orderedTimeDomain_iff {b : ℝ} {n : ℕ} {r : Fin n → Fin n → Prop}
    {σ : Equiv.Perm (Fin n)} {t : Fin n → ℝ} (ht : t ∈ timeChamber b n σ) :
    t ∈ orderedTimeDomain b r ↔ PreservesTimeOrder r σ := by
  constructor
  · intro h i j hij
    apply ht.2.lt_iff_lt.mp
    simpa only [Function.comp_apply, Equiv.apply_symm_apply] using h.2 i j hij
  · intro h
    refine ⟨ht.1, fun i j hij => ?_⟩
    simpa only [Function.comp_apply, Equiv.apply_symm_apply] using ht.2 (h i j hij)

/-- Every permitted chamber is contained in the constrained time domain. -/
theorem timeChamber_subset_orderedTimeDomain {b : ℝ} {n : ℕ}
    {r : Fin n → Fin n → Prop} (σ : TimeShuffle r) :
    timeChamber b n σ.val ⊆ orderedTimeDomain b r :=
  fun _ ht => (mem_orderedTimeDomain_iff ht).mpr σ.property

/-- Admissible interleavings partition the constrained domain outside the already-proved null
set of coinciding times. The precedence relation need not be a total order. -/
theorem iUnion_timeShuffle_ae_eq (b : ℝ) {n : ℕ} (r : Fin n → Fin n → Prop) :
    (⋃ σ : TimeShuffle r, timeChamber b n σ.val) =ᵐ[volume] orderedTimeDomain b r := by
  filter_upwards [ae_injective_times n] with t ht
  apply propext
  constructor
  · intro h
    obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp h
    exact timeChamber_subset_orderedTimeDomain σ hσ
  · intro h
    obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp (mem_iUnion_timeChamber h.1 ht)
    exact Set.mem_iUnion.mpr ⟨⟨σ, (mem_orderedTimeDomain_iff hσ).mp h⟩, hσ⟩

section Integral

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]

omit [CompleteSpace E] in
/-- A finite constrained-time integral is the sum of its admissible chronological chambers. -/
theorem integral_orderedTimeDomain_eq_sum_chambers
    {b : ℝ} {n : ℕ} {r : Fin n → Fin n → Prop} {f : (Fin n → ℝ) → E}
    (hf : IntegrableOn f (orderedTimeDomain b r)) :
    ∫ t in orderedTimeDomain b r, f t = ∑ σ : TimeShuffle r,
      ∫ t in timeChamber b n σ.val, f t := by
  classical
  rw [← setIntegral_congr_set (iUnion_timeShuffle_ae_eq b r)]
  refine integral_iUnion_fintype (fun σ => measurableSet_timeChamber b n σ.val) ?_
    (fun σ => hf.mono_set (timeChamber_subset_orderedTimeDomain σ))
  intro σ τ hστ
  exact pairwise_disjoint_timeChamber b n (fun h => hστ (Subtype.ext h))

/-- Genuine finite-degree reorganization into chronological coordinates, with only those
permutations that preserve the prescribed internal orders. -/
theorem integral_orderedTimeDomain_eq_sum_timeSimplex
    {b : ℝ} {n : ℕ} {r : Fin n → Fin n → Prop} {f : (Fin n → ℝ) → E}
    (hf : IntegrableOn f (orderedTimeDomain b r)) :
    ∫ t in orderedTimeDomain b r, f t = ∑ σ : TimeShuffle r,
      ∫ t in timeSimplex b n, f (t ∘ σ.val.symm) := by
  classical
  rw [integral_orderedTimeDomain_eq_sum_chambers hf]
  exact Finset.sum_congr rfl fun σ _ => integral_timeChamber b n σ.val f

end Integral

/-- Split a concatenated event tuple into its two original bond histories. -/
def splitHistoryTimes (n m : ℕ) :
    (Fin (n + m) → ℝ) ≃ᵐ ((Fin n → ℝ) × (Fin m → ℝ)) :=
  (MeasurableEquiv.piCongrLeft (fun _ : Fin n ⊕ Fin m => ℝ) finSumFinEquiv.symm).trans
    (MeasurableEquiv.sumPiEquivProdPi (fun _ : Fin n ⊕ Fin m => ℝ))

@[simp]
theorem splitHistoryTimes_apply (n m : ℕ) (t : Fin (n + m) → ℝ) :
    splitHistoryTimes n m t = (t ∘ Fin.castAdd m, t ∘ Fin.natAdd n) := by
  ext i <;> simp [splitHistoryTimes, MeasurableEquiv.coe_piCongrLeft,
    Equiv.piCongrLeft_apply_eq_cast, MeasurableEquiv.coe_sumPiEquivProdPi]

/-- Concatenation/splitting is an actual product-Lebesgue change of variables. -/
theorem measurePreserving_splitHistoryTimes (n m : ℕ) :
    MeasurePreserving (splitHistoryTimes n m) volume volume :=
  (volume_measurePreserving_sumPiEquivProdPi (fun _ : Fin n ⊕ Fin m => ℝ)).comp
    (volume_measurePreserving_piCongrLeft (fun _ : Fin n ⊕ Fin m => ℝ) finSumFinEquiv.symm)

/-- The predecessor relation for two histories with respectively `n` and `m` ordered events. -/
def twoHistoryOrder (n m : ℕ) (i j : Fin (n + m)) : Prop :=
  (∃ a b : Fin n, a < b ∧ Fin.castAdd m a = i ∧ Fin.castAdd m b = j) ∨
  (∃ a b : Fin m, a < b ∧ Fin.natAdd n a = i ∧ Fin.natAdd n b = j)

/-- A two-history shuffle preserves both original strict event orders. -/
theorem preserves_twoHistoryOrder_iff {n m : ℕ} (σ : Equiv.Perm (Fin (n + m))) :
    PreservesTimeOrder (twoHistoryOrder n m) σ ↔
      StrictMono (σ.symm ∘ Fin.castAdd m) ∧ StrictMono (σ.symm ∘ Fin.natAdd n) := by
  constructor
  · intro h
    exact ⟨fun a b hab => h _ _ (Or.inl ⟨a, b, hab, rfl, rfl⟩),
      fun a b hab => h _ _ (Or.inr ⟨a, b, hab, rfl, rfl⟩)⟩
  · rintro ⟨h₁, h₂⟩ i j (⟨a, b, hab, rfl, rfl⟩ | ⟨a, b, hab, rfl, rfl⟩)
    · exact h₁ hab
    · exact h₂ hab

/-- The flattened constrained domain is exactly the product of the two ordered histories,
not merely a domain of the same mass. -/
theorem preimage_prod_timeSimplex (b : ℝ) (n m : ℕ) :
    splitHistoryTimes n m ⁻¹' (timeSimplex b n ×ˢ timeSimplex b m) =
      orderedTimeDomain b (twoHistoryOrder n m) := by
  ext t
  simp only [Set.mem_preimage, splitHistoryTimes_apply, Set.mem_prod]
  change ((t ∘ Fin.castAdd m ∈ timeCube b n ∧ StrictMono (t ∘ Fin.castAdd m)) ∧
      (t ∘ Fin.natAdd n ∈ timeCube b m ∧ StrictMono (t ∘ Fin.natAdd n))) ↔ _
  constructor
  · rintro ⟨⟨hc₁, ho₁⟩, ⟨hc₂, ho₂⟩⟩
    refine ⟨?_, ?_⟩
    · intro i _
      obtain ⟨s, rfl⟩ := finSumFinEquiv.surjective i
      cases s with
      | inl a => exact hc₁ a (Set.mem_univ _)
      | inr a => exact hc₂ a (Set.mem_univ _)
    · rintro i j (⟨a, c, hac, rfl, rfl⟩ | ⟨a, c, hac, rfl, rfl⟩)
      · exact ho₁ hac
      · exact ho₂ hac
  · rintro ⟨hc, ho⟩
    exact ⟨⟨fun i _ => hc _ (Set.mem_univ _),
        fun a c hac => ho _ _ (Or.inl ⟨a, c, hac, rfl, rfl⟩)⟩,
      ⟨fun i _ => hc _ (Set.mem_univ _),
        fun a c hac => ho _ _ (Or.inr ⟨a, c, hac, rfl, rfl⟩)⟩⟩

section TwoHistoryIntegral

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]

/-- The product of two ordered bond-history simplexes is the sum over their chronological
interleavings. Each term uses the same globally ordered simplex, with the original histories
recovered by the inverse shuffle. This permits matrix- or operator-valued integrands without
commuting any events. -/
theorem integral_prod_timeSimplex_eq_sum_shuffles
    {b : ℝ} {n m : ℕ} {f : ((Fin n → ℝ) × (Fin m → ℝ)) → E}
    (hf : IntegrableOn f (timeSimplex b n ×ˢ timeSimplex b m)) :
    ∫ p in timeSimplex b n ×ˢ timeSimplex b m, f p =
      ∑ σ : TimeShuffle (twoHistoryOrder n m),
        ∫ t in timeSimplex b (n + m),
          f (t ∘ σ.val.symm ∘ Fin.castAdd m, t ∘ σ.val.symm ∘ Fin.natAdd n) := by
  have hmp := measurePreserving_splitHistoryTimes n m
  have hfi := (hmp.integrableOn_comp_preimage
    (splitHistoryTimes n m).measurableEmbedding).mpr hf
  rw [preimage_prod_timeSimplex] at hfi
  have h := hmp.setIntegral_preimage_emb (splitHistoryTimes n m).measurableEmbedding f
    (timeSimplex b n ×ˢ timeSimplex b m)
  rw [preimage_prod_timeSimplex] at h
  rw [← h]
  simpa only [Function.comp_apply, splitHistoryTimes_apply, Function.comp_assoc] using
    integral_orderedTimeDomain_eq_sum_timeSimplex hfi

/-- Iterated-integral form of the two-history shuffle identity. -/
theorem integral_integral_timeSimplex_eq_sum_shuffles
    {b : ℝ} {n m : ℕ} {f : ((Fin n → ℝ) × (Fin m → ℝ)) → E}
    (hf : IntegrableOn f (timeSimplex b n ×ˢ timeSimplex b m)) :
    (∫ x in timeSimplex b n, ∫ y in timeSimplex b m, f (x, y)) =
      ∑ σ : TimeShuffle (twoHistoryOrder n m),
        ∫ t in timeSimplex b (n + m),
          f (t ∘ σ.val.symm ∘ Fin.castAdd m, t ∘ σ.val.symm ∘ Fin.natAdd n) := by
  rw [← setIntegral_prod f hf]
  exact integral_prod_timeSimplex_eq_sum_shuffles hf

end TwoHistoryIntegral

section FiniteFamily

variable {ι : Type*} [Fintype ι] (counts : ι → ℕ)

/-- Flattening a finite family of time tuples preserves product Lebesgue measure. -/
theorem measurePreserving_uncurryHistoryTimes :
    MeasurePreserving (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm
      volume volume where
  measurable := (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm.measurable
  map_eq := by
    refine (Measure.pi_eq fun s _ => ?_).symm
    rw [MeasurableEquiv.map_apply]
    have hpre : (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ)).symm ⁻¹'
        (Set.univ.pi s) = Set.univ.pi (fun i => Set.univ.pi (fun j => s ⟨i, j⟩)) := by
      ext x
      simp [Set.mem_pi, MeasurableEquiv.coe_piCurry_symm, Sigma.forall, Sigma.uncurry]
    rw [hpre, volume_pi_pi]
    simp only [volume_pi_pi]
    rw [Fintype.prod_sigma]

/-- Enumerate the event slots of a finite family, including empty individual histories. -/
noncomputable def historyEventEquiv :
    (Σ i, Fin (counts i)) ≃ Fin (∑ i, counts i) :=
  Fintype.equivFinOfCardEq (by simp)

/-- Split globally indexed event slots into a finite family of original histories. -/
noncomputable def splitFamilyHistoryTimes :
    (Fin (∑ i, counts i) → ℝ) ≃ᵐ ((i : ι) → Fin (counts i) → ℝ) :=
  (MeasurableEquiv.piCongrLeft (fun _ : Σ i, Fin (counts i) => ℝ)
    (historyEventEquiv counts).symm).trans
      (MeasurableEquiv.piCurry (fun i (_ : Fin (counts i)) => ℝ))

@[simp]
theorem splitFamilyHistoryTimes_apply (t : Fin (∑ i, counts i) → ℝ)
    (i : ι) (j : Fin (counts i)) :
    splitFamilyHistoryTimes counts t i j = t (historyEventEquiv counts ⟨i, j⟩) := by
  simp [splitFamilyHistoryTimes, MeasurableEquiv.coe_piCongrLeft,
    Equiv.piCongrLeft_apply_eq_cast, Sigma.curry]

/-- The finite-family coordinate change uses actual product Lebesgue measure. -/
theorem measurePreserving_splitFamilyHistoryTimes :
    MeasurePreserving (splitFamilyHistoryTimes counts) volume volume :=
  (measurePreserving_uncurryHistoryTimes counts).symm.comp
    (volume_measurePreserving_piCongrLeft (fun _ : Σ i, Fin (counts i) => ℝ)
      (historyEventEquiv counts).symm)

/-- Only events belonging to the same history impose precedence constraints. -/
def familyHistoryOrder (p q : Fin (∑ i, counts i)) : Prop :=
  ∃ i, ∃ j k : Fin (counts i), j < k ∧
    historyEventEquiv counts ⟨i, j⟩ = p ∧ historyEventEquiv counts ⟨i, k⟩ = q

/-- The admissible shuffles are exactly the permutations preserving each history's order. -/
theorem preserves_familyHistoryOrder_iff (σ : Equiv.Perm (Fin (∑ i, counts i))) :
    PreservesTimeOrder (familyHistoryOrder counts) σ ↔
      ∀ i, StrictMono (fun j : Fin (counts i) =>
        σ.symm (historyEventEquiv counts ⟨i, j⟩)) := by
  constructor
  · intro h i j k hjk
    exact h _ _ ⟨i, j, k, hjk, rfl, rfl⟩
  · rintro h p q ⟨i, j, k, hjk, rfl, rfl⟩
    exact h i hjk

/-- The product of the bond-history simplexes is precisely the finite-family constrained domain. -/
theorem preimage_pi_timeSimplex (b : ℝ) :
    splitFamilyHistoryTimes counts ⁻¹'
        (Set.univ.pi fun i => timeSimplex b (counts i)) =
      orderedTimeDomain b (familyHistoryOrder counts) := by
  ext t
  simp only [Set.mem_preimage, Set.mem_pi, Set.mem_univ, forall_true_left]
  constructor
  · intro h
    refine ⟨?_, ?_⟩
    · intro p _
      obtain ⟨⟨i, j⟩, rfl⟩ := (historyEventEquiv counts).surjective p
      simpa only [splitFamilyHistoryTimes_apply] using (h i).1 j (Set.mem_univ _)
    · rintro p q ⟨i, j, k, hjk, rfl, rfl⟩
      simpa only [splitFamilyHistoryTimes_apply] using (h i).2 hjk
  · rintro ⟨hc, ho⟩ i
    refine ⟨?_, ?_⟩
    · intro j _
      simpa only [splitFamilyHistoryTimes_apply] using hc _ (Set.mem_univ _)
    · intro j k hjk
      simpa only [splitFamilyHistoryTimes_apply] using
        ho _ _ ⟨i, j, k, hjk, rfl, rfl⟩

/-- Fixed finite event counts on any finite family of bonds reorganize into exactly the global
chronological interleavings preserving every bond's internal order. This is a finite-degree
identity, separate from any subsequent infinite summation over the event counts. -/
theorem integral_pi_timeSimplex_eq_sum_shuffles
    {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [CompleteSpace E]
    {b : ℝ} {f : ((i : ι) → Fin (counts i) → ℝ) → E}
    (hf : IntegrableOn f (Set.univ.pi fun i => timeSimplex b (counts i))) :
    ∫ p in Set.univ.pi (fun i => timeSimplex b (counts i)), f p =
      ∑ σ : TimeShuffle (familyHistoryOrder counts),
        ∫ t in timeSimplex b (∑ i, counts i),
          f (fun i j => t (σ.val.symm (historyEventEquiv counts ⟨i, j⟩))) := by
  have hmp := measurePreserving_splitFamilyHistoryTimes counts
  have hfi := (hmp.integrableOn_comp_preimage
    (splitFamilyHistoryTimes counts).measurableEmbedding).mpr hf
  rw [preimage_pi_timeSimplex] at hfi
  have h := hmp.setIntegral_preimage_emb (splitFamilyHistoryTimes counts).measurableEmbedding f
    (Set.univ.pi fun i => timeSimplex b (counts i))
  rw [preimage_pi_timeSimplex] at h
  rw [← h]
  convert integral_orderedTimeDomain_eq_sum_timeSimplex hfi using 1
  · rfl
  · apply Finset.sum_congr rfl
    intro σ _
    congr 1
    funext t
    congr 1
    funext i j
    simp only [splitFamilyHistoryTimes_apply, Function.comp_apply]

end FiniteFamily

end CertifiedGaps.TimeOrdering
