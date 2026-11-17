import Mathlib.Data.Fin.Tuple.Sort
import Mathlib.Data.Fintype.Perm
import Mathlib.MeasureTheory.Integral.Bochner.Set
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar

/-!
# Chronological chambers

Permutation chambers partition the finite-dimensional time cube outside the null set of
colliding times. This is the finite-degree analytic ingredient in reorganizing continuous
Poisson event expansions. The chronological convention agrees with the history-space
construction in the pinned spatial-transfer proof:
the September 24, 2026 periodic spin-one Haldane-gap preprint, `transfer.tex`,
commit `adc7f1241b42e322a6451854ab7e4b4c146bf78a` in `openai/math`.
-/

open MeasureTheory Set
open scoped ENNReal

namespace CertifiedGaps.TimeOrdering

/-- The open time cube, including its one-point zero-dimensional case. -/
def timeCube (b : ℝ) (n : ℕ) : Set (Fin n → ℝ) :=
  Set.univ.pi fun _ => Set.Ioo 0 b

/-- The increasing-time simplex inside the time cube. -/
def timeSimplex (b : ℝ) (n : ℕ) : Set (Fin n → ℝ) :=
  timeCube b n ∩ {t | StrictMono t}

/-- The chamber in which the coordinates indexed by `σ` occur in increasing order. -/
def timeChamber (b : ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) : Set (Fin n → ℝ) :=
  timeCube b n ∩ {t | StrictMono (t ∘ σ)}

/-- Time cubes are measurable. -/
theorem measurableSet_timeCube (b : ℝ) (n : ℕ) : MeasurableSet (timeCube b n) :=
  MeasurableSet.univ_pi fun _ => measurableSet_Ioo

/-- Strictly increasing tuples form a measurable set. -/
theorem measurableSet_strictMono (n : ℕ) :
    MeasurableSet {t : Fin n → ℝ | StrictMono t} := by
  simp only [StrictMono, ofPred_forall]
  exact MeasurableSet.iInter fun i => MeasurableSet.iInter fun j =>
    MeasurableSet.iInter fun _ => measurableSet_lt (measurable_pi_apply i)
      (measurable_pi_apply j)

/-- The time simplex is measurable. -/
theorem measurableSet_timeSimplex (b : ℝ) (n : ℕ) : MeasurableSet (timeSimplex b n) :=
  (measurableSet_timeCube b n).inter (measurableSet_strictMono n)

/-- Every chronological chamber is measurable. -/
theorem measurableSet_timeChamber (b : ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) :
    MeasurableSet (timeChamber b n σ) := by
  exact (measurableSet_timeCube b n).inter
    ((measurableSet_strictMono n).preimage (measurable_pi_iff.mpr fun i =>
      measurable_pi_apply (σ i)))

/-- A collision of two distinct coordinates has Lebesgue measure zero. -/
theorem volume_timeCollision {n : ℕ} {i j : Fin n} (hij : i ≠ j) :
    volume {t : Fin n → ℝ | t i = t j} = 0 := by
  let S : Submodule ℝ (Fin n → ℝ) :=
    { carrier := {t | t i = t j}
      zero_mem' := rfl
      add_mem' := by intro x y hx hy; change x i + y i = x j + y j; rw [hx, hy]
      smul_mem' := by intro c x hx; change c * x i = c * x j; rw [hx] }
  apply Measure.addHaar_submodule volume S
  intro htop
  have hm : Pi.single i (1 : ℝ) ∈ S := by rw [htop]; trivial
  change (Pi.single i 1 : Fin n → ℝ) i = (Pi.single i 1 : Fin n → ℝ) j at hm
  simp [hij.symm] at hm

/-- Almost every tuple has distinct times. -/
theorem ae_injective_times (n : ℕ) :
    ∀ᵐ t : Fin n → ℝ, Function.Injective t := by
  have h : ∀ᵐ t : Fin n → ℝ, ∀ i j, t i = t j → i = j := by
    apply ae_all_iff.mpr; intro i
    apply ae_all_iff.mpr; intro j
    by_cases hij : i = j
    · exact ae_of_all _ fun _ _ => hij
    · have hn : ∀ᵐ t : Fin n → ℝ, t i ≠ t j :=
        compl_mem_ae_iff.mpr (volume_timeCollision hij)
      exact hn.mono fun _ h he => (h he).elim
  exact h

/-- Distinct chronological chambers do not intersect. -/
theorem pairwise_disjoint_timeChamber (b : ℝ) (n : ℕ) :
    Pairwise (fun σ τ => Disjoint (timeChamber b n σ) (timeChamber b n τ)) := by
  intro σ τ hστ
  apply Set.disjoint_left.mpr
  intro t htσ htτ
  have heq := Tuple.unique_monotone htσ.2.monotone htτ.2.monotone
  have hinj : Function.Injective t :=
    htσ.2.injective.of_comp_right σ.surjective
  apply hστ
  apply Equiv.ext; intro i
  exact hinj (congrFun heq i)

/-- Every collision-free tuple in the time cube belongs to a chronological chamber. -/
theorem mem_iUnion_timeChamber {b : ℝ} {n : ℕ} {t : Fin n → ℝ}
    (ht : t ∈ timeCube b n) (hinj : Function.Injective t) :
    t ∈ ⋃ σ, timeChamber b n σ := by
  refine Set.mem_iUnion.mpr ⟨Tuple.sort t, ht, ?_⟩
  exact (Tuple.monotone_sort t).strictMono_of_injective (hinj.comp (Tuple.sort t).injective)

/-- The chronological chambers cover the time cube up to null collisions. -/
theorem iUnion_timeChamber_ae_eq (b : ℝ) (n : ℕ) :
    (⋃ σ, timeChamber b n σ) =ᵐ[volume] timeCube b n := by
  filter_upwards [ae_injective_times n] with t ht
  apply propext
  constructor
  · intro h; obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp h; exact hσ.1
  · exact fun h => mem_iUnion_timeChamber h ht

/-- Coordinate permutation with the convention `t ↦ t ∘ σ`. -/
def permuteTimes {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    (Fin n → ℝ) ≃ᵐ (Fin n → ℝ) :=
  MeasurableEquiv.piCongrLeft (fun _ => ℝ) σ.symm

@[simp]
theorem permuteTimes_apply {n : ℕ} (σ : Equiv.Perm (Fin n)) (t : Fin n → ℝ) :
    permuteTimes σ t = t ∘ σ := by
  ext i
  simp [permuteTimes, MeasurableEquiv.coe_piCongrLeft, Equiv.piCongrLeft_apply]

/-- Permuting event coordinates preserves product Lebesgue measure. -/
theorem measurePreserving_permuteTimes {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    MeasurePreserving (permuteTimes σ) volume volume :=
  volume_measurePreserving_piCongrLeft (fun _ => ℝ) σ.symm

@[simp]
theorem comp_perm_mem_timeCube {n : ℕ} (σ : Equiv.Perm (Fin n))
    (b : ℝ) (t : Fin n → ℝ) : t ∘ σ ∈ timeCube b n ↔ t ∈ timeCube b n := by
  simp only [timeCube, Set.mem_pi, Set.mem_univ, forall_true_left,
    Function.comp_apply]
  exact (σ.surjective.forall (p := fun i => t i ∈ Ioo 0 b)).symm

/-- The increasing simplex pulls back to the corresponding chronological chamber. -/
theorem preimage_timeSimplex (b : ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) :
    permuteTimes σ ⁻¹' timeSimplex b n = timeChamber b n σ := by
  ext t
  simp only [Set.mem_preimage, timeSimplex, timeChamber, Set.mem_inter_iff,
    Set.mem_ofPred_eq, permuteTimes_apply, comp_perm_mem_timeCube]

/-- All chronological chambers have the same volume. -/
theorem volume_timeChamber (b : ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) :
    volume (timeChamber b n σ) = volume (timeSimplex b n) := by
  rw [← preimage_timeSimplex]
  exact (measurePreserving_permuteTimes σ).measure_preimage
    (measurableSet_timeSimplex b n).nullMeasurableSet

/-- The finite time cube has volume `b ^ n` for nonnegative `b`. -/
theorem volume_timeCube (b : ℝ) (n : ℕ) :
    volume (timeCube b n) = ENNReal.ofReal b ^ n := by
  simp [timeCube, volume_pi_pi, Real.volume_Ioo]

/-- In particular the finite time cube always has finite volume. -/
theorem volume_timeCube_ne_top (b : ℝ) (n : ℕ) : volume (timeCube b n) ≠ ⊤ := by
  rw [volume_timeCube]
  finiteness

/-- The increasing simplex has finite volume, including degree zero. -/
theorem volume_timeSimplex_ne_top (b : ℝ) (n : ℕ) : volume (timeSimplex b n) ≠ ⊤ :=
  ne_top_of_le_ne_top (volume_timeCube_ne_top b n) (measure_mono Set.inter_subset_left)

/-- There are `n!` chronological chambers, so the cube volume is `n!` times simplex volume. -/
theorem factorial_mul_volume_timeSimplex (b : ℝ) (n : ℕ) :
    (n.factorial : ℝ≥0∞) * volume (timeSimplex b n) = ENNReal.ofReal b ^ n := by
  calc
    (n.factorial : ℝ≥0∞) * volume (timeSimplex b n) =
        ∑' σ : Equiv.Perm (Fin n), volume (timeChamber b n σ) := by
      simp [volume_timeChamber, Fintype.card_perm, nsmul_eq_mul]
    _ = volume (⋃ σ, timeChamber b n σ) :=
      (measure_iUnion (pairwise_disjoint_timeChamber b n) (measurableSet_timeChamber b n)).symm
    _ = volume (timeCube b n) := measure_congr (iUnion_timeChamber_ae_eq b n)
    _ = ENNReal.ofReal b ^ n := volume_timeCube b n

/-- The factorial normalization of the ordered-time Lebesgue simplex. -/
theorem volume_timeSimplex (b : ℝ) (n : ℕ) :
    volume (timeSimplex b n) = ENNReal.ofReal b ^ n / n.factorial := by
  apply (ENNReal.eq_div_iff (by simp [Nat.factorial_ne_zero]) (by simp)).mpr
  simpa [mul_comm] using factorial_mul_volume_timeSimplex b n

/-- The real-valued simplex mass used in the ordered exponential series. -/
theorem volumeReal_timeSimplex {b : ℝ} (hb : 0 ≤ b) (n : ℕ) :
    volume.real (timeSimplex b n) = b ^ n / n.factorial := by
  simp [measureReal_def, volume_timeSimplex, ENNReal.toReal_div, ENNReal.toReal_pow,
    ENNReal.toReal_ofReal hb]

section Integral

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- An integrable time-cube integrand is the finite sum of its chronological chamber integrals. -/
theorem integral_timeCube_eq_sum_chambers [CompleteSpace E]
    {b : ℝ} {n : ℕ} {f : (Fin n → ℝ) → E}
    (hf : IntegrableOn f (timeCube b n)) :
    ∫ t in timeCube b n, f t = ∑ σ : Equiv.Perm (Fin n),
      ∫ t in timeChamber b n σ, f t := by
  rw [← setIntegral_congr_set (iUnion_timeChamber_ae_eq b n)]
  exact integral_iUnion_fintype (measurableSet_timeChamber b n)
    (pairwise_disjoint_timeChamber b n) (fun _ => hf.mono_set Set.inter_subset_left)

/-- Transport a chamber integral to increasing chronological coordinates. -/
theorem integral_timeChamber [CompleteSpace E] (b : ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n))
    (f : (Fin n → ℝ) → E) :
    ∫ t in timeChamber b n σ, f t =
      ∫ t in timeSimplex b n, f (t ∘ σ.symm) := by
  have h := (measurePreserving_permuteTimes σ).setIntegral_preimage_emb
    (permuteTimes σ).measurableEmbedding (fun t => f (t ∘ σ.symm)) (timeSimplex b n)
  rw [preimage_timeSimplex] at h
  simpa [Function.comp_def] using h

/-- Time ordering is an honest finite decomposition of the cube integral. Collisions have
already been proved null, rather than assumed away. No symmetry of the integrand is required. -/
theorem integral_timeCube_eq_sum_timeSimplex [CompleteSpace E]
    {b : ℝ} {n : ℕ} {f : (Fin n → ℝ) → E}
    (hf : IntegrableOn f (timeCube b n)) :
    ∫ t in timeCube b n, f t = ∑ σ : Equiv.Perm (Fin n),
      ∫ t in timeSimplex b n, f (t ∘ σ.symm) := by
  rw [integral_timeCube_eq_sum_chambers hf]
  exact Finset.sum_congr rfl fun σ _ => integral_timeChamber b n σ f

end Integral

end CertifiedGaps.TimeOrdering
