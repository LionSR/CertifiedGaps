/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Open.Gap
import CertifiedGaps.Heisenberg.Gap
import CertifiedGaps.Models.SpinOneLargeD

/-!
# The boundary-field gap of the open spin-one Heisenberg chain

This file states and proves, from four finite inputs, Theorem 2 of the thermal-MPO note (the
boundary-field gap of OpenAI's *A boundary-field gap for the spin-one Heisenberg chain*, 2026):
for every odd `n ≥ 1921` the open chain

  `H_n^{same} = ∑ᵢ Sᵢ · Sᵢ₊₁ - h S₀ᶻ - h S_{n-1}ᶻ`,  `h = 3/5`,

on `(ℂ³)^{⊗n}` (`endFieldHamiltonian h h n`) has a unique ground state and
`Δ_n > (log 10)/392 ≈ 0.005874`. The same inputs give the same gap for the opposite-field chain
at every even `n ≥ 1920`.

The chain is the open chain of the Heisenberg ring `SpinOneLargeD.model 0` with the boundary
weights `w_z = -i h` (`boundary h`), whose end term is `-h Sᶻ`. At odd length that open chain is
`H_n^{same}`, and at even length it is the opposite-field chain
`H_n^{opp} = ∑ᵢ Sᵢ · Sᵢ₊₁ - h S₀ᶻ + h S_{n-1}ᶻ` (`endFieldHamiltonian h (-h) n`;
`openHamiltonian_odd`, `openHamiltonian_even`),
which is the chain the boundary inputs are stated for.

## The inputs

`BoundaryFieldInputs h` collects the four inequalities of the note's eq. (C1), each a statement
about traces `Re Tr e^{-βH}` of Lean Hamiltonians at `β = 49`:
* `ring_spatial`, `ring_thermal`: `P_s(120, 49) ≥ 0.948`, `P_τ(240, 49) ≥ 0.948` for the ring;
* `boundary`: the boundary defect `e₁₄₂ ≤ 0.14` of the opposite-field chain;
* `open_thermal`: the opposite-field thermal purity `P_τ(240, 49) ≥ 0.79`.

The values come from Ref. [boundary] of the note (the finite initialization of OpenAI's
boundary-field preprint).
`spin_one_boundary_field_gap` uses them as the `sorry` lemmas `seed_*` below.

## Main results

* `openGapGT_of_boundaryFieldInputs_certificate`: the gap of the open chain at every length from
  the four inputs and any valid certificate for their seeds, without `sorry`. Its odd and even
  cases are `gapGT_of_boundaryFieldInputs_certificate` (equal fields) and
  `oppositeFieldGapGT_of_boundaryFieldInputs_certificate` (opposite fields), with the note's and
  the sharper certificate as `gapGT_of_boundaryFieldInputs(_sharp)` and
  `oppositeFieldGapGT_of_boundaryFieldInputs(_sharp)`.
* `spin_one_boundary_field_gap(_sharp)`: the equal-field gaps `(log 10)/392` and
  `log(10⁶/2241)/392` at odd `n ≥ 1921`, and `spin_one_opposite_field_gap(_sharp)` the same gaps
  for the opposite-field chain at even `n ≥ 1920`, with the four inputs as open `sorry` lemmas.
-/

open Matrix NormedSpace Real

namespace CertifiedGaps.SpinOneHeisenberg

/-! ### The open chains -/

/-- The chain `∑ᵢ Sᵢ · Sᵢ₊₁ - h₀ S₀ᶻ - h₁ S_{n-1}ᶻ` with end fields `h₀`, `h₁` on `n` sites
(`0` for `n = 0`). The equal-field chain is `h₁ = h₀`, the opposite-field chain `h₁ = -h₀`. -/
def endFieldHamiltonian (h₀ h₁ : ℝ) : (n : ℕ) → Matrix (Fin n → Fin 3) (Fin n → Fin 3) ℂ
  | 0 => 0
  | n + 1 => ∑ j : Fin n, bondOp j.castSucc j.succ SpinOne.exchange +
      siteOp 0 (-(h₀ : ℂ) • SpinOne.S 2) + siteOp (Fin.last n) (-(h₁ : ℂ) • SpinOne.S 2)

/-- A real multiple of `Sᶻ` is Hermitian. -/
theorem smul_spinOne_isHermitian (h : ℝ) : ((h : ℂ) • SpinOne.S 2).IsHermitian := by
  change ((h : ℂ) • SpinOne.S 2)ᴴ = (h : ℂ) • SpinOne.S 2
  rw [conjTranspose_smul, (SpinOne.S_isHermitian 2).eq, Complex.star_def, Complex.conj_ofReal]

/-- The chain with end fields is Hermitian. -/
theorem endFieldHamiltonian_isHermitian (h₀ h₁ : ℝ) (n : ℕ) :
    (endFieldHamiltonian h₀ h₁ n).IsHermitian := by
  cases n with
  | zero => exact isHermitian_zero
  | succ n =>
    have hz : ∀ h : ℝ, (-(h : ℂ) • SpinOne.S 2).IsHermitian := fun h ↦ by
      rw [← Complex.ofReal_neg]; exact smul_spinOne_isHermitian _
    exact ((isSelfAdjoint_sum _ fun j _ ↦ bondOp_isHermitian _ _ SpinOne.exchange_isHermitian).add
      (siteOp_isHermitian _ (hz h₀))).add (siteOp_isHermitian _ (hz h₁))

/-- The spectrum of the equal-field chain with multiplicity. -/
noncomputable def equalFieldSpectrum (h : ℝ) (n : ℕ) : (Fin n → Fin 3) → ℝ :=
  (endFieldHamiltonian_isHermitian h h n).eigenvalues

/-- The spectrum of the opposite-field chain with multiplicity. -/
noncomputable def oppositeFieldSpectrum (h : ℝ) (n : ℕ) : (Fin n → Fin 3) → ℝ :=
  (endFieldHamiltonian_isHermitian h (-h) n).eigenvalues

/-- The partition function `Re Tr e^{-β H_n^{opp}}` of the opposite-field chain. -/
noncomputable def oppositeFieldPartitionFunction (h : ℝ) (n : ℕ) (β : ℝ) : ℝ :=
  (trace (exp ((-(β : ℂ)) • endFieldHamiltonian h (-h) n))).re

/-- The thermal purity `P_τ(n, β)` of the opposite-field chain. -/
noncomputable def oppositeFieldThermalPurity (h : ℝ) (n : ℕ) (β : ℝ) : ℝ :=
  oppositeFieldPartitionFunction h n (2 * β) / oppositeFieldPartitionFunction h n β ^ 2

/-! ### The chains as open chains of the Heisenberg ring -/

/-- The boundary weights `w_z = -i h` of the Heisenberg ring `SpinOneLargeD.model 0`: the end term
is `-h Sᶻ` (`boundary_mat`). -/
noncomputable def boundary (h : ℝ) : (SpinOneLargeD.model 0).Boundary where
  w α := if α = 2 then -(Complex.I * h) else 0
  w_symm α := .inr ⟨SpinOne.G_transpose α, by
    split_ifs <;> simp [Complex.conj_ofReal]⟩

/-- The end term is `-h Sᶻ`. -/
theorem boundary_mat (h : ℝ) : (boundary h).mat = -(h : ℂ) • SpinOne.S 2 := by
  ext a c
  simp only [RPRing.Boundary.mat, boundary, Fin.sum_univ_three, RPRing.Cc, SpinOneLargeD.model,
    SpinOne.G, SpinOne.S, Matrix.neg_apply, Matrix.add_apply, Matrix.smul_apply, map_apply,
    of_apply, smul_eq_mul]
  simp only [Fin.reduceEq, ite_false, ite_true, zero_mul, zero_add]
  push_cast
  ring_nf

/-- The transposed end term is `h Sᶻ`. -/
theorem boundary_mat_transpose (h : ℝ) : (boundary h).matᵀ = (h : ℂ) • SpinOne.S 2 := by
  rw [boundary_mat]
  ext a c
  simp only [transpose_apply, Matrix.smul_apply, SpinOne.S, of_apply, smul_eq_mul,
    SpinOne.leviCivita_swap 2 c a]
  push_cast
  ring

/-- The ring has no on-site term at `D = 0`. -/
theorem oc_zero : (SpinOneLargeD.model 0).oc = 0 := by
  rw [SpinOneLargeD.oc_eq]; simp

/-- At odd length the open chain of the Heisenberg ring is the equal-field chain. -/
theorem openHamiltonian_odd (h : ℝ) {n : ℕ} (hn : Odd n) :
    (SpinOneLargeD.model 0).openHamiltonian (boundary h) n = endFieldHamiltonian h h n := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by obtain ⟨k, hk⟩ := hn; omega⟩
  have hne : ¬ Even (m + 1) := Nat.not_even_iff_odd.2 hn
  simp only [RPRing.openHamiltonian, endFieldHamiltonian, RPRing.Boundary.endMat, hne,
    ite_false, SpinOneLargeD.bondMatrix_eq, oc_zero, siteOp_zero, Finset.sum_const_zero,
    add_zero, boundary_mat]

/-- At even length the open chain of the Heisenberg ring is the opposite-field chain. -/
theorem openHamiltonian_even (h : ℝ) {n : ℕ} (hn : Even n) (hn0 : n ≠ 0) :
    (SpinOneLargeD.model 0).openHamiltonian (boundary h) n = endFieldHamiltonian h (-h) n := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  simp only [RPRing.openHamiltonian, endFieldHamiltonian, RPRing.Boundary.endMat, hn,
    ite_true, SpinOneLargeD.bondMatrix_eq, oc_zero, siteOp_zero, Finset.sum_const_zero,
    add_zero]
  rw [boundary_mat_transpose, boundary_mat, Complex.ofReal_neg, neg_neg]

/-- The spectrum of the odd open chain of the ring is that of the equal-field chain. -/
theorem openSpectrum_odd (h : ℝ) {n : ℕ} (hn : Odd n) :
    (SpinOneLargeD.model 0).openSpectrum (boundary h) n = equalFieldSpectrum h n :=
  IsHermitian.eigenvalues_eq_of_eq _ _ (openHamiltonian_odd h hn)

/-- The spectrum of the even open chain of the ring is that of the opposite-field chain. -/
theorem openSpectrum_even (h : ℝ) {n : ℕ} (hn : Even n) (hn0 : n ≠ 0) :
    (SpinOneLargeD.model 0).openSpectrum (boundary h) n = oppositeFieldSpectrum h n :=
  IsHermitian.eigenvalues_eq_of_eq _ _ (openHamiltonian_even h hn hn0)

/-- At even length the partition functions of the ring's open chain and the opposite-field
chain agree. -/
theorem openPartitionFunction_even (h : ℝ) (N : ℕ) (β : ℝ) :
    (SpinOneLargeD.model 0).openPartitionFunction (boundary h) (2 * N) β =
      oppositeFieldPartitionFunction h (2 * N) β := by
  rcases Nat.eq_zero_or_pos N with rfl | hN
  · simp [RPRing.openPartitionFunction, oppositeFieldPartitionFunction, RPRing.openHamiltonian,
      endFieldHamiltonian]
  · rw [RPRing.openPartitionFunction, openHamiltonian_even h (even_two_mul N) (by omega)]
    rfl

/-- The boundary defect only sees the even-length traces. -/
theorem boundaryDefect_congr {Z G G' : ℕ → ℝ} (hG : ∀ N, G (2 * N) = G' (2 * N)) (K : ℕ) :
    boundaryDefect Z G (2 * K) = boundaryDefect Z G' (2 * K) := by
  simp only [boundaryDefect, leadingOverlap, hG]

/-! ### The theorem -/

/-- The four inputs of the boundary-field route at `β₀ = 49`, the note's eq. (C1). -/
structure BoundaryFieldInputs (h : ℝ) : Prop where
  /-- The ring spatial purity `P_s(120, 49) ≥ 0.948`. -/
  ring_spatial : 948 / 1000 ≤ spatialPurity 120 49
  /-- The ring thermal purity `P_τ(240, 49) ≥ 0.948`. -/
  ring_thermal : 948 / 1000 ≤ thermalPurity 240 49
  /-- The boundary defect `e₁₄₂ ≤ 0.14` of the opposite-field chain at `β = 49`. -/
  boundary : boundaryDefect (fun L ↦ partitionFunction L 49)
    (fun n ↦ oppositeFieldPartitionFunction h n 49) 142 ≤ 14 / 100
  /-- The opposite-field thermal purity `P_τ(240, 49) ≥ 0.79`. -/
  open_thermal : 79 / 100 ≤ oppositeFieldThermalPurity h 240 49

/-- **The open-chain gap from the inputs and a certificate.** Under `BoundaryFieldInputs h`, a
valid certificate `c` for the seeds of the inputs (`N₀ = 60`, `K₀ = 71`, `u₀ ≥ 0.052`,
`E₀ ≥ 0.14`, `η₀ ≥ 0.21`) gives `Δ_n > log(1/(10x)) / (49 · 2^m)` for the open chain of the
Heisenberg ring with end term `-h Sᶻ` at every length `n ≥ 240 · 2^m`: the equal-field chain at
odd `n`, the opposite-field chain at even `n`. -/
theorem openGapGT_of_boundaryFieldInputs_certificate {h : ℝ} (hin : BoundaryFieldInputs h)
    (c : OpenCertificate) (hc : c.Valid) (hN : c.N₀ = 60) (hK : c.K₀ = 71)
    (hu : 13 / 250 ≤ c.u₀) (hE : 14 / 100 ≤ c.EAt 0) (hη : 21 / 100 ≤ c.ηAt 0) {n : ℕ}
    (hn : 240 * 2 ^ c.m ≤ n) :
    GapGT ((SpinOneLargeD.model 0).openSpectrum (boundary h) n)
      (Real.log (1 / (10 * (c.x : ℝ))) / (49 * 2 ^ c.m)) := by
  have hZ : ∀ L β, (SpinOneLargeD.model 0).partitionFunction L β = partitionFunction L β :=
    SpinOneLargeD.partitionFunction_zero
  have hu' : (13 / 250 : ℝ) ≤ c.u₀ := by simpa using (Rat.cast_le (K := ℝ)).2 hu
  have hE' : (14 / 100 : ℝ) ≤ c.EAt 0 := by simpa using (Rat.cast_le (K := ℝ)).2 hE
  have hη' : (21 / 100 : ℝ) ≤ c.ηAt 0 := by simpa using (Rat.cast_le (K := ℝ)).2 hη
  refine (SpinOneLargeD.model 0).openGapGT_of_certificate rfl (boundary h) c hc (β₀ := 49)
    (by norm_num) ?_ ?_ ?_ ?_ (by rw [hN]; omega) <;> (try rw [hN]) <;> (try rw [hK])
  · have := hin.ring_spatial
    simp only [RPRing.spatialPurity, hZ, spatialPurity] at this ⊢
    linarith
  · have := hin.ring_thermal
    simp only [RPRing.thermalPurity, hZ, thermalPurity] at this ⊢
    linarith
  · simp only [hZ]
    rw [boundaryDefect_congr (G' := fun n ↦ oppositeFieldPartitionFunction h n 49)
      (openPartitionFunction_even h · 49)]
    linarith [hin.boundary]
  · have := hin.open_thermal
    simp only [RPRing.openThermalPurity, openPartitionFunction_even]
    simp only [oppositeFieldThermalPurity] at this
    linarith

/-- **The boundary-field gap from its inputs and a certificate**, at odd lengths: the
equal-field chain. -/
theorem gapGT_of_boundaryFieldInputs_certificate {h : ℝ} (hin : BoundaryFieldInputs h)
    (c : OpenCertificate) (hc : c.Valid) (hN : c.N₀ = 60) (hK : c.K₀ = 71)
    (hu : 13 / 250 ≤ c.u₀) (hE : 14 / 100 ≤ c.EAt 0) (hη : 21 / 100 ≤ c.ηAt 0) {n : ℕ}
    (hn : 240 * 2 ^ c.m ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum h n) (Real.log (1 / (10 * (c.x : ℝ))) / (49 * 2 ^ c.m)) := by
  rw [← openSpectrum_odd h hno]
  exact openGapGT_of_boundaryFieldInputs_certificate hin c hc hN hK hu hE hη hn

/-- **The opposite-field gap from the inputs and a certificate**, at even lengths. -/
theorem oppositeFieldGapGT_of_boundaryFieldInputs_certificate {h : ℝ}
    (hin : BoundaryFieldInputs h) (c : OpenCertificate) (hc : c.Valid) (hN : c.N₀ = 60)
    (hK : c.K₀ = 71) (hu : 13 / 250 ≤ c.u₀) (hE : 14 / 100 ≤ c.EAt 0) (hη : 21 / 100 ≤ c.ηAt 0)
    {n : ℕ} (hn : 240 * 2 ^ c.m ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum h n) (Real.log (1 / (10 * (c.x : ℝ))) / (49 * 2 ^ c.m)) := by
  rw [← openSpectrum_even h hne (by have := Nat.one_le_two_pow (n := c.m); omega)]
  exact openGapGT_of_boundaryFieldInputs_certificate hin c hc hN hK hu hE hη hn

/-- **The boundary-field gap from its inputs.** Under `BoundaryFieldInputs h`, every odd
equal-field chain of `n ≥ 1921` sites has a unique ground state and gap `Δ_n > (log 10)/392`. -/
theorem gapGT_of_boundaryFieldInputs {h : ℝ} (hin : BoundaryFieldInputs h) {n : ℕ}
    (hn : 1921 ≤ n) (hno : Odd n) : GapGT (equalFieldSpectrum h n) (Real.log 10 / 392) := by
  have h := gapGT_of_boundaryFieldInputs_certificate hin boundaryFieldCertificate
    boundaryFieldCertificate_valid rfl rfl (by decide +kernel) (by decide +kernel)
    (by decide +kernel) (n := n)
    (by rw [show boundaryFieldCertificate.m = 3 from rfl]; omega) hno
  rw [show boundaryFieldCertificate.m = 3 from rfl,
    show (boundaryFieldCertificate.x : ℝ) = 1 / 100 by norm_num [boundaryFieldCertificate]] at h
  convert h using 2 <;> norm_num

/-- **A sharper boundary-field gap from the same inputs.** Under `BoundaryFieldInputs h`, the
certificate `sharpBoundaryFieldCertificate` gives every odd equal-field chain of `n ≥ 1921` sites
a unique ground state and gap `Δ_n > log(10⁶/2241)/392 ≈ 0.01556`. -/
theorem gapGT_of_boundaryFieldInputs_sharp {h : ℝ} (hin : BoundaryFieldInputs h) {n : ℕ}
    (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum h n) (Real.log (1000000 / 2241) / 392) := by
  have h := gapGT_of_boundaryFieldInputs_certificate hin sharpBoundaryFieldCertificate
    sharpBoundaryFieldCertificate_valid rfl rfl (by decide +kernel) (by decide +kernel)
    (by decide +kernel) (n := n)
    (by rw [show sharpBoundaryFieldCertificate.m = 3 from rfl]; omega) hno
  rw [show sharpBoundaryFieldCertificate.m = 3 from rfl,
    show (sharpBoundaryFieldCertificate.x : ℝ) = 2241 / 10000000 by
      norm_num [sharpBoundaryFieldCertificate]] at h
  convert h using 2 <;> norm_num

/-- **The opposite-field gap from the inputs.** Under `BoundaryFieldInputs h`, every even
opposite-field chain of `n ≥ 1920` sites has a unique ground state and gap
`Δ_n > (log 10)/392`. -/
theorem oppositeFieldGapGT_of_boundaryFieldInputs {h : ℝ} (hin : BoundaryFieldInputs h)
    {n : ℕ} (hn : 1920 ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum h n) (Real.log 10 / 392) := by
  have h := oppositeFieldGapGT_of_boundaryFieldInputs_certificate hin boundaryFieldCertificate
    boundaryFieldCertificate_valid rfl rfl (by decide +kernel) (by decide +kernel)
    (by decide +kernel) (n := n)
    (by rw [show boundaryFieldCertificate.m = 3 from rfl]; omega) hne
  rw [show boundaryFieldCertificate.m = 3 from rfl,
    show (boundaryFieldCertificate.x : ℝ) = 1 / 100 by norm_num [boundaryFieldCertificate]] at h
  convert h using 2 <;> norm_num

/-- **A sharper opposite-field gap from the same inputs**: `Δ_n > log(10⁶/2241)/392` for every
even `n ≥ 1920`. -/
theorem oppositeFieldGapGT_of_boundaryFieldInputs_sharp {h : ℝ} (hin : BoundaryFieldInputs h)
    {n : ℕ} (hn : 1920 ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum h n) (Real.log (1000000 / 2241) / 392) := by
  have h := oppositeFieldGapGT_of_boundaryFieldInputs_certificate hin
    sharpBoundaryFieldCertificate sharpBoundaryFieldCertificate_valid rfl rfl (by decide +kernel)
    (by decide +kernel) (by decide +kernel) (n := n)
    (by rw [show sharpBoundaryFieldCertificate.m = 3 from rfl]; omega) hne
  rw [show sharpBoundaryFieldCertificate.m = 3 from rfl,
    show (sharpBoundaryFieldCertificate.x : ℝ) = 2241 / 10000000 by
      norm_num [sharpBoundaryFieldCertificate]] at h
  convert h using 2 <;> norm_num

/-- **Open input 1.** The ring spatial purity `P_s(120, 49) ≥ 0.948`. -/
theorem seed_ring_spatial : 948 / 1000 ≤ spatialPurity 120 49 := by
  sorry

/-- **Open input 2.** The ring thermal purity `P_τ(240, 49) ≥ 0.948`. -/
theorem seed_ring_thermal : 948 / 1000 ≤ thermalPurity 240 49 := by
  sorry

/-- **Open input 3.** The boundary defect `e₁₄₂ ≤ 0.14` of the opposite-field chain with
`h = 3/5` at `β = 49`. -/
theorem seed_boundary : boundaryDefect (fun L ↦ partitionFunction L 49)
    (fun n ↦ oppositeFieldPartitionFunction (3 / 5) n 49) 142 ≤ 14 / 100 := by
  sorry

/-- **Open input 4.** The opposite-field thermal purity `P_τ(240, 49) ≥ 0.79` with `h = 3/5`. -/
theorem seed_open_thermal : 79 / 100 ≤ oppositeFieldThermalPurity (3 / 5) 240 49 := by
  sorry

/-- **Theorem 2** of the thermal-MPO note (the boundary-field gap). For every odd `n ≥ 1921`
the open spin-one Heisenberg chain with equal end fields `h = 3/5`,
`∑ᵢ Sᵢ · Sᵢ₊₁ - (3/5)(S₀ᶻ + S_{n-1}ᶻ)`, has a unique ground state and every other level lies
more than `(log 10)/392 ≈ 0.005874` above it. -/
theorem spin_one_boundary_field_gap {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum (3 / 5) n) (Real.log 10 / 392) :=
  gapGT_of_boundaryFieldInputs ⟨seed_ring_spatial, seed_ring_thermal, seed_boundary,
    seed_open_thermal⟩ hn hno

/-- **The sharper boundary-field gap.** From the same four inputs as
`spin_one_boundary_field_gap`, for every odd `n ≥ 1921` the open spin-one Heisenberg chain with
equal end fields `h = 3/5` has a unique ground state and every other level lies more than
`log(10⁶/2241)/392 ≈ 0.01556` above it. -/
theorem spin_one_boundary_field_gap_sharp {n : ℕ} (hn : 1921 ≤ n) (hno : Odd n) :
    GapGT (equalFieldSpectrum (3 / 5) n) (Real.log (1000000 / 2241) / 392) :=
  gapGT_of_boundaryFieldInputs_sharp ⟨seed_ring_spatial, seed_ring_thermal, seed_boundary,
    seed_open_thermal⟩ hn hno

/-- **The opposite-field gap.** From the same four inputs, for every even `n ≥ 1920` the open
spin-one Heisenberg chain with opposite end fields,
`∑ᵢ Sᵢ · Sᵢ₊₁ - (3/5)(S₀ᶻ - S_{n-1}ᶻ)`, has a unique ground state and every other level lies
more than `(log 10)/392` above it. -/
theorem spin_one_opposite_field_gap {n : ℕ} (hn : 1920 ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum (3 / 5) n) (Real.log 10 / 392) :=
  oppositeFieldGapGT_of_boundaryFieldInputs ⟨seed_ring_spatial, seed_ring_thermal, seed_boundary,
    seed_open_thermal⟩ hn hne

/-- **The sharper opposite-field gap**: `Δ_n > log(10⁶/2241)/392 ≈ 0.01556` for every even
`n ≥ 1920`, from the same four inputs. -/
theorem spin_one_opposite_field_gap_sharp {n : ℕ} (hn : 1920 ≤ n) (hne : Even n) :
    GapGT (oppositeFieldSpectrum (3 / 5) n) (Real.log (1000000 / 2241) / 392) :=
  oppositeFieldGapGT_of_boundaryFieldInputs_sharp ⟨seed_ring_spatial, seed_ring_thermal,
    seed_boundary, seed_open_thermal⟩ hn hne

end CertifiedGaps.SpinOneHeisenberg
