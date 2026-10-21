import CertifiedGaps.Transfer.HistoryShuffle


/-!
# Chronological chambers for atomless time measures

The chamber decomposition uses product symmetry and null collisions, not translation
invariance. The time domain can be any measurable subset of the real line, including all of it.
-/

open MeasureTheory Set
open scoped ENNReal

namespace CertifiedGaps.AtomlessTimeOrdering

/-- The finite time cube over an arbitrary time domain. -/
def timeCube (s : Set ℝ) (n : ℕ) : Set (Fin n → ℝ) :=
  Set.univ.pi fun _ => s

/-- Strictly chronological times in a specified domain. -/
def timeSimplex (s : Set ℝ) (n : ℕ) : Set (Fin n → ℝ) :=
  timeCube s n ∩ {t | StrictMono t}

/-- The chronological chamber indexed by a coordinate permutation. -/
def timeChamber (s : Set ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) : Set (Fin n → ℝ) :=
  timeCube s n ∩ {t | StrictMono (t ∘ σ)}

/-- Independent coordinates with the same event-time measure. -/
noncomputable def timeMeasure (ν : Measure ℝ) (n : ℕ) : Measure (Fin n → ℝ) :=
  Measure.pi fun _ => ν

instance (ν : Measure ℝ) [SigmaFinite ν] (n : ℕ) : SigmaFinite (timeMeasure ν n) :=
  inferInstanceAs (SigmaFinite (Measure.pi fun _ : Fin n => ν))

instance (ν : Measure ℝ) [IsFiniteMeasure ν] (n : ℕ) : IsFiniteMeasure (timeMeasure ν n) :=
  inferInstanceAs (IsFiniteMeasure (Measure.pi fun _ : Fin n => ν))

theorem measurableSet_timeCube {s : Set ℝ} (hs : MeasurableSet s) (n : ℕ) :
    MeasurableSet (timeCube s n) := MeasurableSet.univ_pi fun _ => hs

theorem measurableSet_timeSimplex {s : Set ℝ} (hs : MeasurableSet s) (n : ℕ) :
    MeasurableSet (timeSimplex s n) :=
  (measurableSet_timeCube hs n).inter (TimeOrdering.measurableSet_strictMono n)

theorem measurableSet_timeChamber {s : Set ℝ} (hs : MeasurableSet s) (n : ℕ)
    (σ : Equiv.Perm (Fin n)) : MeasurableSet (timeChamber s n σ) :=
  (measurableSet_timeCube hs n).inter
    ((TimeOrdering.measurableSet_strictMono n).preimage (measurable_pi_iff.mpr fun i =>
      measurable_pi_apply (σ i)))

/-- Distinct product coordinates collide with probability zero for any atomless time law. -/
theorem measure_timeCollision (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    {n : ℕ} {i j : Fin n} (hij : i ≠ j) :
    timeMeasure ν n {t | t i = t j} = 0 := by
  cases n with
  | zero => exact Fin.elim0 i
  | succ n =>
    obtain ⟨k, rfl⟩ := Fin.exists_succAbove_eq hij.symm
    have hs : MeasurableSet {p : ℝ × (Fin n → ℝ) | p.1 = p.2 k} :=
      measurableSet_eq_fun measurable_fst ((measurable_pi_apply k).comp measurable_snd)
    have hnull : (ν.prod (Measure.pi fun _ : Fin n => ν))
        {p : ℝ × (Fin n → ℝ) | p.1 = p.2 k} = 0 := by
      apply Measure.measure_prod_null_of_ae_null hs
      filter_upwards [] with x
      change (Measure.pi fun _ : Fin n => ν) {t : Fin n → ℝ | x = t k} = 0
      have hset : {t : Fin n → ℝ | x = t k} =
          (fun t : Fin n → ℝ => t k) ⁻¹' ({x} : Set ℝ) := by
        ext t
        simp [eq_comm]
      rw [hset]
      exact Measure.pi_eval_preimage_null (fun _ : Fin n => ν) (measure_singleton x)
    have hmp := measurePreserving_piFinSuccAbove (fun _ : Fin (n + 1) => ν) i
    have he := hmp.measure_preimage hs.nullMeasurableSet
    simpa [timeMeasure, MeasurableEquiv.piFinSuccAbove, Fin.insertNthEquiv, Fin.removeNth] using
      he.trans hnull

/-- Almost every finite tuple has distinct coordinates. -/
theorem ae_injective_times (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν] (n : ℕ) :
    ∀ᵐ t ∂timeMeasure ν n, Function.Injective t := by
  apply ae_all_iff.mpr
  intro i
  apply ae_all_iff.mpr
  intro j
  by_cases hij : i = j
  · exact ae_of_all _ fun _ _ => hij
  · have hn : ∀ᵐ t ∂timeMeasure ν n, t i ≠ t j :=
      compl_mem_ae_iff.mpr (measure_timeCollision ν hij)
    exact hn.mono fun _ h he => (h he).elim

/-- Chronological chambers are disjoint before taking any measure. -/
theorem pairwise_disjoint_timeChamber (s : Set ℝ) (n : ℕ) :
    Pairwise (fun σ τ => Disjoint (timeChamber s n σ) (timeChamber s n τ)) := by
  intro σ τ hστ
  apply Set.disjoint_left.mpr
  intro t htσ htτ
  have heq := Tuple.unique_monotone htσ.2.monotone htτ.2.monotone
  have hinj : Function.Injective t := htσ.2.injective.of_comp_right σ.surjective
  apply hστ
  apply Equiv.ext
  intro i
  exact hinj (congrFun heq i)

theorem mem_iUnion_timeChamber {s : Set ℝ} {n : ℕ} {t : Fin n → ℝ}
    (ht : t ∈ timeCube s n) (hinj : Function.Injective t) :
    t ∈ ⋃ σ, timeChamber s n σ := by
  refine Set.mem_iUnion.mpr ⟨Tuple.sort t, ht, ?_⟩
  exact (Tuple.monotone_sort t).strictMono_of_injective (hinj.comp (Tuple.sort t).injective)

/-- Chronological chambers partition the product time domain outside a proved null set. -/
theorem iUnion_timeChamber_ae_eq (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    (s : Set ℝ) (n : ℕ) :
    (⋃ σ, timeChamber s n σ) =ᵐ[timeMeasure ν n] timeCube s n := by
  filter_upwards [ae_injective_times ν n] with t ht
  apply propext
  constructor
  · intro h
    obtain ⟨σ, hσ⟩ := Set.mem_iUnion.mp h
    exact hσ.1
  · exact fun h => mem_iUnion_timeChamber h ht

/-- Coordinate permutations preserve arbitrary equal-coordinate product measures. -/
theorem measurePreserving_permuteTimes (ν : Measure ℝ) [SigmaFinite ν]
    {n : ℕ} (σ : Equiv.Perm (Fin n)) :
    MeasurePreserving (TimeOrdering.permuteTimes σ) (timeMeasure ν n) (timeMeasure ν n) :=
  MeasureTheory.measurePreserving_piCongrLeft (fun _ => ν) σ.symm

@[simp]
theorem comp_perm_mem_timeCube {n : ℕ} (σ : Equiv.Perm (Fin n))
    (s : Set ℝ) (t : Fin n → ℝ) : t ∘ σ ∈ timeCube s n ↔ t ∈ timeCube s n := by
  simp only [timeCube, Set.mem_pi, Set.mem_univ, forall_true_left, Function.comp_apply]
  exact (σ.surjective.forall (p := fun i => t i ∈ s)).symm

theorem preimage_timeSimplex (s : Set ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) :
    TimeOrdering.permuteTimes σ ⁻¹' timeSimplex s n = timeChamber s n σ := by
  ext t
  simp only [Set.mem_preimage, timeSimplex, timeChamber, Set.mem_inter_iff,
    Set.mem_ofPred_eq, TimeOrdering.permuteTimes_apply, comp_perm_mem_timeCube]

theorem measure_timeChamber (ν : Measure ℝ) [SigmaFinite ν]
    {s : Set ℝ} (hs : MeasurableSet s) (n : ℕ) (σ : Equiv.Perm (Fin n)) :
    timeMeasure ν n (timeChamber s n σ) = timeMeasure ν n (timeSimplex s n) := by
  rw [← preimage_timeSimplex]
  exact (measurePreserving_permuteTimes ν σ).measure_preimage
    (measurableSet_timeSimplex hs n).nullMeasurableSet

theorem measure_timeCube (ν : Measure ℝ) [SigmaFinite ν] (s : Set ℝ) (n : ℕ) :
    timeMeasure ν n (timeCube s n) = ν s ^ n := by
  simp [timeMeasure, timeCube, Measure.pi_pi]

/-- The factorial law holds for every atomless sigma-finite time measure. -/
theorem measure_timeSimplex (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    {s : Set ℝ} (hs : MeasurableSet s) (n : ℕ) :
    timeMeasure ν n (timeSimplex s n) = ν s ^ n / n.factorial := by
  apply (ENNReal.eq_div_iff (by simp [Nat.factorial_ne_zero]) (by simp)).mpr
  calc
    n.factorial * timeMeasure ν n (timeSimplex s n) =
        ∑' σ : Equiv.Perm (Fin n), timeMeasure ν n (timeChamber s n σ) := by
      simp [measure_timeChamber ν hs, Fintype.card_perm, nsmul_eq_mul]
    _ = timeMeasure ν n (⋃ σ, timeChamber s n σ) :=
      (measure_iUnion (pairwise_disjoint_timeChamber s n)
        (measurableSet_timeChamber hs n)).symm
    _ = timeMeasure ν n (timeCube s n) := measure_congr (iUnion_timeChamber_ae_eq ν s n)
    _ = ν s ^ n := measure_timeCube ν s n

section Integral

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- Transport a chamber integral without changing the integrand's event order. -/
theorem integral_timeChamber (ν : Measure ℝ) [SigmaFinite ν]
    (s : Set ℝ) (n : ℕ) (σ : Equiv.Perm (Fin n)) (f : (Fin n → ℝ) → E) :
    ∫ t in timeChamber s n σ, f t ∂timeMeasure ν n =
      ∫ t in timeSimplex s n, f (t ∘ σ.symm) ∂timeMeasure ν n := by
  have h := (measurePreserving_permuteTimes ν σ).setIntegral_preimage_emb
    (TimeOrdering.permuteTimes σ).measurableEmbedding (fun t => f (t ∘ σ.symm))
      (timeSimplex s n)
  rw [preimage_timeSimplex] at h
  simpa [Function.comp_def] using h

/-- Finite chronological decomposition for arbitrary Banach-valued integrands. -/
theorem integral_timeCube_eq_sum_timeSimplex (ν : Measure ℝ) [SigmaFinite ν] [NullSingletonClass ν]
    {s : Set ℝ} (hs : MeasurableSet s) {n : ℕ} {f : (Fin n → ℝ) → E}
    (hf : IntegrableOn f (timeCube s n) (timeMeasure ν n)) :
    ∫ t in timeCube s n, f t ∂timeMeasure ν n = ∑ σ : Equiv.Perm (Fin n),
      ∫ t in timeSimplex s n, f (t ∘ σ.symm) ∂timeMeasure ν n := by
  rw [← setIntegral_congr_set (iUnion_timeChamber_ae_eq ν s n)]
  rw [integral_iUnion_fintype (measurableSet_timeChamber hs n)
    (pairwise_disjoint_timeChamber s n) (fun _ => hf.mono_set Set.inter_subset_left)]
  exact Finset.sum_congr rfl fun σ _ => integral_timeChamber ν s n σ f

end Integral

end CertifiedGaps.AtomlessTimeOrdering
