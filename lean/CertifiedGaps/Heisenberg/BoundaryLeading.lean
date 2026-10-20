/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Heisenberg.BoundaryInitTable
import CertifiedGaps.Heisenberg.BoundaryTraces
import CertifiedGaps.Heisenberg.SectorTransferC
import CertifiedGaps.Local.SpinOneRotation
import CertifiedGaps.Open.LeadingLimit
import CertifiedGaps.Transfer.SingleSiteLeading

/-!
# The leading boundary weight at `b₀` from the fifteen thermal traces

This file proves `BoundaryInit.ThermalCore` of `Heisenberg/BoundaryInitTable.lean` from the
seed-36 trace table, the table of the fifteen short-chain thermal traces and the bound
`Z̃₇₂(b₀) > 1` of the periodic trial state (`BoundaryInit.thermalCore_of_tables`). It follows
the reduction of the thermal appendix of OpenAI's *A boundary-field gap for the spin-one
Heisenberg chain* (2026), eqs. (appth:average) and (appth:polynomial), with the leading-line
invariance of Proposition per:concentration, in the single-site transfer matrix `X_M` of `M` time
steps (`Transfer/SingleSite.lean`) and then in the limit `M → ∞`.

## The argument

Write `β = b₀ = 49/4`, `s = e^{-βa}` for the shift and `X̃ = s X_M`.

1. **Pairings** (`Transfer/OpenSingleSite.lean`). With the boundary vector `v` of the end field
   `+h Sᶻ` (`h = 3/5`), whose conjugate is the vector of `-h Sᶻ`, the pairing
   `⟨v, X̃ⁿ w⟩` of `v` with the vector `w` of an end field `c S^α` tends to `sⁿ` times the
   partition function of the chain with end fields `-h Sᶻ` and `c S^α`. This gives the even
   opposite-field traces (`w = v`), the odd equal-field traces (`w = D_x v`, the `π` rotation
   about `x`) and, with the average `Π` over the rotation group `A₄`, the twelve chains with the
   right field along `±x, ±y, ±z`. Global rotations about `z` identify the four transverse
   chains with `Z_{n,0}`, so `⟨v, X̃ⁿ Π v⟩ → I(n) = (Z_{n,1} + Z_{n,-1} + 4 Z_{n,0})/6`.
2. **The leading line** (`Transfer/SingleSiteLeading.lean`). For large `M` the sector
   eigenvalues of `X̃` obey the seed-36 table, whose filters and masses give one eigenvalue
   `|λ₀| > ℓ = .999` of the trivial `D₂` block, the other eigenvalues of that block below
   `(m_J - ℓ¹²)^{1/12} < .873126` and those of the block `N` below the cap `.87125`. So `X̃` has an
   isolated leading line with `r = .873126`, fixed by the rotation group.
3. **The witness** (`Open/TwoBoundary.lean`). The approximants then satisfy the boundary moment
   inequalities with `μ = |λ₀| > ℓ`, `ρ ≤ r / |λ₀| < r / ℓ = .874`, and
   `c P_*(λ₀)² ≤ ∑ₙ aₙ ⟨v, X̃ⁿ Π v⟩`, `P_*² = ∑ₙ aₙ xⁿ`.
4. **The limit** (`Open/LeadingLimit.lean`). These closed facts pass to a witness of the exact
   shifted traces. Its sign is `+1` (or `c = 0`), and `∑ aₙ I(n) < .003745` is
   `BoundaryTraces.moment_polynomial_lt_of_table`.

## Main results

* `pairing_tendsto`, `openM_tendsto`, `ringM_tendsto`: the pairings and traces of `X̃_M`
  converge to the shifted open-chain and ring partition functions.
* `avgM_tendsto`: the averaged pairings converge to `I(n)`.
* `eventually_leadingLine`: the leading line of `X̃_M` for large `M`, from the trace table.
* `witness_of_tables`: the witness with `λ₀ ≥ .999`, `ρ₀ ≤ .874`, `c₀ P_*(λ₀)² ≤ .003745`.
* `BoundaryInit.thermalCore_of_tables`: `ThermalCore` from the two tables and
  `Z̃₇₂(b₀) > 1`; `BoundaryInit.thermalCore_of_tables_mps` with the periodic trial state.

## Differences from the paper

* The paper's transfer operator `X_b` acts on continuous-time histories; here it is replaced by
  the single-site transfer matrices `X_M` of `M` time steps, whose pairings converge to the same
  partition functions. Every spectral statement (leading line, its invariance, the moment
  form) is proved for `X_M` and passed to the traces by limits.
* The sector bounds are those of the Lean seed-36 certificate (`seed36`: caps `1.0013641597`,
  `.87125`, masses `1.0657405`, `.2783105`, `ℓ = .999`) instead of the paper's
  (`U = 1.00139`, `V = .872`, `m_J = 1.0657205`, `m_N = .2783155`). They give the residual ratio
  `ρ₀ ≤ .873126/.999 < .874` directly from the eigenvalue caps; the paper's `D(k)` is not used.
* The witness is obtained at `λ₀ ≥ ℓ = .999`; `ThermalCore` asks only `λ₀ ≥ 9/10`.
* The sign of the leading eigenvalue is proved in the limit from the positivity of the odd
  open-chain traces (`OpenMomentWitness.eq_zero_or_sign_eq_one`) instead of the ring traces.
* Leading-line invariance: the paper lets `C` act on the simple real leading line by `+1` and
  the three `π` rotations by a common sign equal to its square. Here the table places the leading
  eigenvalue in the trivial `D₂` block, so `D₂` fixes the line outright; the cyclic relabeling
  acts by a real cube root of `1`.
* The witness of `ThermalCore` is now built from the single-site transfer matrices, where the
  odd open-chain traces are the same-field pairings `⟨v, X̃ⁿ D_x v⟩` (the paper's eq. (tr:open)
  with `U_P`), not from the checkerboard transfer matrix of `Open/Transfer.lean`. Both
  approximate the same partition functions, so the statement is unchanged.
-/

open Matrix Filter Topology Finset
open scoped Kronecker

namespace CertifiedGaps.SpinOneHeisenberg

namespace BoundaryLeading

open SingleSite SpinOne

/-! ### End fields as cap weights -/

/-- The cap weights `i c` on the axis `α`, for the end term `c S^α`. -/
noncomputable def fieldW (α : Fin 3) (c : ℝ) : Fin 3 → ℂ :=
  fun β ↦ if β = α then Complex.I * c else 0

/-- The left end field `-(3/5) Sᶻ`. -/
noncomputable def wL : Fin 3 → ℂ := fieldW 2 (-(3 / 5))

/-- The right end field `+(3/5) Sᶻ` of the opposite-field chain. -/
noncomputable def wR : Fin 3 → ℂ := fieldW 2 (3 / 5)

/-- The end term of the weights `fieldW α c` is `c S^α`. -/
theorem capTerm_fieldW (α : Fin 3) (c : ℝ) : capTerm G (fieldW α c) = (c : ℂ) • S α := by
  ext a b
  simp only [capTerm, fieldW, bondC, Matrix.neg_apply, Matrix.sum_apply, Matrix.smul_apply,
    map_apply, smul_eq_mul, ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ, ite_true,
    S, G, of_apply]
  push_cast
  ring

/-- The open chain with end terms `-(3/5) Sᶻ` and `c S^α` on `k + 1` sites. -/
noncomputable def fieldH (k : ℕ) (α : Fin 3) (c : ℝ) :
    Matrix (Fin (k + 1) → Fin 3) (Fin (k + 1) → Fin 3) ℂ :=
  ∑ j : Fin k, bondOp j.castSucc j.succ exchange + siteOp 0 (-((3 / 5 : ℝ) : ℂ) • S 2) +
    siteOp (Fin.last k) ((c : ℂ) • S α)

/-- The single-site open chain of the spin-one generators is `fieldH`. -/
theorem openHamiltonian_fieldW (k : ℕ) (α : Fin 3) (c : ℝ) :
    openHamiltonian G 0 k wL (fieldW α c) = fieldH k α c := by
  rw [openHamiltonian, signedBond_G_one, wL, capTerm_fieldW, capTerm_fieldW,
    Matrix.map_zero _ Complex.ofReal_zero, siteSum_zero, add_zero, fieldH]
  push_cast; rfl

/-- `fieldH` along `z` with `c = ±3/5` and along `x` with `c = 3/5` are the paper's `H_{n,g}`. -/
theorem fieldH_eq_hamiltonian (k : ℕ) :
    fieldH k 2 (3 / 5) = BoundaryTraces.hamiltonian 1 (k + 1) ∧
      fieldH k 2 (-(3 / 5)) = BoundaryTraces.hamiltonian (-1) (k + 1) ∧
        fieldH k 0 (3 / 5) = BoundaryTraces.hamiltonian 0 (k + 1) := by
  refine ⟨?_, ?_, ?_⟩ <;>
    simp only [fieldH, BoundaryTraces.hamiltonian, BoundaryTraces.rightField] <;> push_cast <;>
    simp

/-! ### Rotations about `z` -/

/-- Conjugating a one-site operator by the same unitary on every site. -/
private theorem conj_siteOp_const' {L : ℕ} (V : Matrix (Fin 3) (Fin 3) ℂ) (hV : Vᴴ * V = 1)
    (j : Fin L) (A : Matrix (Fin 3) (Fin 3) ℂ) :
    (piTensor fun _ : Fin L ↦ V)ᴴ * siteOp j A * piTensor (fun _ : Fin L ↦ V) =
      siteOp j (Vᴴ * A * V) := by
  have h := conj_siteOp (fun _ : Fin L ↦ Vᴴ) (fun _ ↦ by rw [conjTranspose_conjTranspose, hV]) j A
  simp only [conjTranspose_conjTranspose] at h
  rw [conjTranspose_piTensor, h]

/-- A rotation `Q` fixing `Sᶻ` and mapping `S^α` to `σ S^{α'}` maps `fieldH α c` to
`fieldH α' (σ c)`, so the two chains have the same Gibbs trace. -/
theorem gibbsTrace_fieldH_rot {Q : Matrix (Fin 3) (Fin 3) ℝ} (hQ : PreservesCross Q)
    (hQo : Qᵀ * Q = 1) (k : ℕ) {α α' : Fin 3} {σ : ℝ}
    (hz : (Q.map (↑))ᴴ * S 2 * Q.map (↑) = S 2)
    (hα : (Q.map (↑))ᴴ * S α * Q.map (↑) = (σ : ℂ) • S α') (c : ℝ) (β : ℝ) :
    gibbsTrace (fieldH k α c) β = gibbsTrace (fieldH k α' (σ * c)) β := by
  set V : Matrix (Fin 3) (Fin 3) ℂ := Q.map (↑)
  have hV : Vᴴ * V = 1 := map_conjTranspose_mul_map hQo
  have hex : (V ⊗ₖ V)ᴴ * exchange * (V ⊗ₖ V) = exchange := by
    rw [← rotExchange_one, conj_rotExchange hQ, Matrix.mul_one, hQo]
  rw [← gibbsTrace_conj (piTensor fun _ : Fin (k + 1) ↦ V) _
    (piTensor_const_conjTranspose_mul_self V hV)]
  congr 1
  simp only [fieldH, Matrix.mul_add, Matrix.add_mul, Finset.mul_sum, Finset.sum_mul]
  congr 1
  · congr 1
    · refine Finset.sum_congr rfl fun j _ ↦ ?_
      rw [conj_bondOp_const V hV (by simp [Fin.ext_iff, Fin.val_succ]), hex]
    · rw [conj_siteOp_const' V hV, Matrix.mul_smul, Matrix.smul_mul, hz]
  · rw [conj_siteOp_const' V hV, Matrix.mul_smul, Matrix.smul_mul, hα, smul_smul]
    push_cast; ring_nf

/-- The `π` rotation about `z`. -/
def rotPiZ : Matrix (Fin 3) (Fin 3) ℝ := !![-1, 0, 0; 0, -1, 0; 0, 0, 1]

/-- The rotation by `π/2` about `z`, with rows `Q_{α ·}`: `S^x ↦ -S^y`, `S^y ↦ S^x`. -/
def rotHalfPiZ : Matrix (Fin 3) (Fin 3) ℝ := !![0, -1, 0; 1, 0, 0; 0, 0, 1]

/-- The `π` rotation about `z` is proper: it preserves the cross product. -/
theorem rotPiZ_preservesCross : PreservesCross rotPiZ := by
  intro α b c
  fin_cases α <;> fin_cases b <;> fin_cases c <;>
    simp [rotPiZ, Fin.sum_univ_three, leviCivita]

/-- The `π/2` rotation about `z` is proper: it preserves the cross product. -/
theorem rotHalfPiZ_preservesCross : PreservesCross rotHalfPiZ := by
  intro α b c
  fin_cases α <;> fin_cases b <;> fin_cases c <;>
    simp [rotHalfPiZ, Fin.sum_univ_three, leviCivita]

/-- The `π` rotation about `z` is orthogonal. -/
theorem rotPiZ_orthogonal : rotPiZᵀ * rotPiZ = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [rotPiZ, mul_apply, Fin.sum_univ_three]

/-- The `π/2` rotation about `z` is orthogonal. -/
theorem rotHalfPiZ_orthogonal : rotHalfPiZᵀ * rotHalfPiZ = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [rotHalfPiZ, mul_apply, Fin.sum_univ_three]

/-- **The transverse chains agree.** The right fields `c S^x`, `-c S^x`, `c S^y` and `-c S^y`
give the same partition function: rotations about `z` fix the bonds and the left field. -/
theorem gibbsTrace_fieldH_transverse (k : ℕ) (c β : ℝ) :
    gibbsTrace (fieldH k 0 (-c)) β = gibbsTrace (fieldH k 0 c) β ∧
      gibbsTrace (fieldH k 1 c) β = gibbsTrace (fieldH k 0 c) β ∧
        gibbsTrace (fieldH k 1 (-c)) β = gibbsTrace (fieldH k 0 c) β := by
  have hS : ∀ (Q : Matrix (Fin 3) (Fin 3) ℝ), PreservesCross Q → ∀ α,
      (Q.map (↑))ᴴ * S α * Q.map (↑) = ∑ β, (Q α β : ℂ) • S β :=
    fun Q hQ α ↦ conj_S_of_preservesCross hQ α
  have z1 : (rotPiZ.map (↑))ᴴ * S 2 * rotPiZ.map (↑) = S 2 := by
    rw [hS _ rotPiZ_preservesCross]; simp [rotPiZ, Fin.sum_univ_three]
  have z2 : (rotHalfPiZ.map (↑))ᴴ * S 2 * rotHalfPiZ.map (↑) = S 2 := by
    rw [hS _ rotHalfPiZ_preservesCross]; simp [rotHalfPiZ, Fin.sum_univ_three]
  have x1 : (rotPiZ.map (↑))ᴴ * S 0 * rotPiZ.map (↑) = ((-1 : ℝ) : ℂ) • S 0 := by
    rw [hS _ rotPiZ_preservesCross]; simp [rotPiZ, Fin.sum_univ_three]
  have y1 : (rotHalfPiZ.map (↑))ᴴ * S 1 * rotHalfPiZ.map (↑) = ((1 : ℝ) : ℂ) • S 0 := by
    rw [hS _ rotHalfPiZ_preservesCross]; simp [rotHalfPiZ, Fin.sum_univ_three]
  have e1 : ∀ c, gibbsTrace (fieldH k 0 (-c)) β = gibbsTrace (fieldH k 0 c) β := fun c ↦ by
    rw [gibbsTrace_fieldH_rot rotPiZ_preservesCross rotPiZ_orthogonal k z1 x1, neg_one_mul,
      neg_neg]
  have e2 : ∀ c, gibbsTrace (fieldH k 1 c) β = gibbsTrace (fieldH k 0 c) β := fun c ↦ by
    rw [gibbsTrace_fieldH_rot rotHalfPiZ_preservesCross rotHalfPiZ_orthogonal k z2 y1, one_mul]
  exact ⟨e1 c, e2 c, (e2 (-c)).trans (e1 c)⟩

/-! ### The approximating pairings -/

/-- The shift factor `s = e^{-b₀ a}` per site. -/
noncomputable abbrev sf : ℝ := BoundaryInit.shiftFactor (49 / 4)

/-- The time step `b₀ / M`. -/
noncomputable abbrev τM (M : ℕ) : ℝ := (49 / 4 : ℝ) / M

/-- The shifted transfer matrix `X̃_M = s X_M` of `M` time steps. -/
noncomputable def tX (M : ℕ) : Matrix (Fin M → Option (Fin 3)) (Fin M → Option (Fin 3)) ℂ :=
  (sf • transfer G 0 M (τM M)).map (↑)

/-- The boundary vector of the right end field `+(3/5) Sᶻ`. -/
noncomputable def vM (M : ℕ) : (Fin M → Option (Fin 3)) → ℂ := capVec (τM M) wR

/-- The pairing `Re ⟨v, X̃ⁿ w⟩` with the boundary vector of the end field `c S^α`. -/
noncomputable def pairingM (M n : ℕ) (α : Fin 3) (c : ℝ) : ℝ :=
  (star (vM M) ⬝ᵥ (tX M ^ n *ᵥ capVec (τM M) (fieldW α c))).re

/-- `v` is conjugate to the vector of the left end field. -/
theorem star_vM (M : ℕ) : star (vM M) = capVec (τM M) wL := by
  rw [vM, SingleSite.star_capVec]
  congr 1
  funext α
  simp only [wR, wL, fieldW]
  split_ifs <;> simp

/-- `X̃ⁿ = sⁿ Xⁿ`. -/
theorem tX_pow (M n : ℕ) :
    tX M ^ n = ((sf ^ n : ℝ) : ℂ) • (transfer G 0 M (τM M) ^ n).map (↑) := by
  rw [tX, show ((sf • transfer G 0 M (τM M)).map Complex.ofReal) =
    Complex.ofRealHom.mapMatrix (sf • transfer G 0 M (τM M)) from rfl, ← map_pow, smul_pow]
  ext i j
  simp

/-- The pairing is `sⁿ` times the unshifted bilinear pairing. -/
theorem pairingM_eq (M n : ℕ) (α : Fin 3) (c : ℝ) :
    pairingM M n α c = sf ^ n * (capVec (τM M) wL ⬝ᵥ
      ((transfer G 0 M (τM M) ^ n).map (↑) *ᵥ capVec (τM M) (fieldW α c))).re := by
  rw [pairingM, star_vM, tX_pow, smul_mulVec, dotProduct_smul, smul_eq_mul,
    Complex.re_ofReal_mul]

/-- **The pairings converge** to the shifted partition function of the chain with end fields
`-(3/5) Sᶻ` and `c S^α`. -/
theorem pairing_tendsto (k : ℕ) (α : Fin 3) (c : ℝ) :
    Tendsto (fun M ↦ pairingM M (k + 1) α c) atTop
      (𝓝 (sf ^ (k + 1) * gibbsTrace (fieldH k α c) (49 / 4))) := by
  simp only [pairingM_eq]
  have h := tendsto_openPairing G 0 k wL (fieldW α c) (β := 49 / 4) (by norm_num)
  rw [openHamiltonian_fieldW] at h
  exact h.const_mul _

/-- The ring approximants `Re Tr X̃^L`. -/
noncomputable def ringM (M L : ℕ) : ℝ := (trace (tX M ^ L)).re

/-- The trace of a complexified real matrix. -/
theorem trace_map_ofReal {D : Type*} [Fintype D] (A : Matrix D D ℝ) :
    trace (A.map (↑) : Matrix D D ℂ) = ((trace A : ℝ) : ℂ) := by
  simp [trace, Complex.ofReal_sum]

/-- **The ring approximants converge** to the shifted ring traces. -/
theorem ringM_tendsto {L : ℕ} (hL : 2 ≤ L) :
    Tendsto (fun M ↦ ringM M L) atTop (𝓝 (BoundaryInit.shiftedRingZ (49 / 4) L)) := by
  obtain ⟨n, rfl⟩ : ∃ n, L = n + 1 := ⟨L - 1, by omega⟩
  have h := tendsto_trace_stepOp_pow G 0 hL (fun _ ↦ 1) id (49 / 4)
  have e : ∀ M : ℕ, (trace (stepOp G 0 (n + 1) (fun _ ↦ 1) id (((49 / 4 : ℝ) / M : ℝ) : ℂ) ^
      M)).re = trace (transfer G 0 M (τM M) ^ (n + 1)) := fun M ↦ by
    rw [trace_stepOp_pow_id G 0 n M _ (by positivity), Complex.ofReal_re,
      show labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from funext labelChar_one, diagonal_one,
      Matrix.one_mul]
  simp only [e] at h
  have hZ : twistPartitionFunction G 0 (n + 1) (fun _ ↦ 1) id (49 / 4) =
      partitionFunction (n + 1) (49 / 4) := by
    rw [twistPartitionFunction_id, SingleSite.partitionFunction, hamiltonian_plain]; rfl
  rw [hZ] at h
  have := h.const_mul (sf ^ (n + 1))
  refine this.congr fun M ↦ ?_
  rw [ringM, tX_pow, trace_smul, trace_map_ofReal, smul_eq_mul, Complex.re_ofReal_mul,
    Complex.ofReal_re]

/-- A sign symmetry flips the end field on its axis by its sign. -/
theorem charMat_mulVec_fieldW {m : ℕ} (r : Fin 3 → ℝ) (τ : ℝ) (α : Fin 3) (c : ℝ) :
    charMat (m := m) r *ᵥ capVec τ (fieldW α c) = capVec τ (fieldW α (r α * c)) := by
  rw [charMat_mulVec_capVec]
  congr 1
  funext β
  simp only [fieldW]
  split_ifs with h
  · subst h; push_cast; ring
  · simp

/-- The open approximants: the opposite-field pairing at even lengths, the same-field pairing
`⟨v, X̃ⁿ D_x v⟩` at odd lengths. -/
noncomputable def openM (M n : ℕ) : ℝ :=
  if Even n then pairingM M n 2 (3 / 5) else pairingM M n 2 (-(3 / 5))

/-- **The open approximants converge** to the shifted open-chain traces of
`Heisenberg/BoundaryInit.lean` (opposite fields at even, equal fields at odd lengths). -/
theorem openM_tendsto {n : ℕ} (hn : 1 ≤ n) :
    Tendsto (fun M ↦ openM M n) atTop (𝓝 (BoundaryInit.shiftedOpenZ (49 / 4) n)) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  obtain ⟨h1, h2, -⟩ := fieldH_eq_hamiltonian k
  by_cases he : Even (k + 1)
  · simp only [openM, he, ite_true]
    convert pairing_tendsto k 2 (3 / 5) using 2
    have hopp : (SpinOneLargeD.model 0).openPartitionFunction (boundary (3 / 5)) (k + 1)
        (49 / 4) = oppositeFieldPartitionFunction (3 / 5) (k + 1) (49 / 4) := by
      obtain ⟨N, hN⟩ := he
      rw [hN, ← two_mul]
      exact openPartitionFunction_even _ N _
    rw [h1, BoundaryInit.shiftedOpenZ, BoundaryInit.openZ, hopp, BoundaryTraces.hamiltonian_one]
    rfl
  · simp only [openM, he, ite_false]
    convert pairing_tendsto k 2 (-(3 / 5)) using 2
    rw [h2, BoundaryInit.shiftedOpenZ, BoundaryInit.openZ, RPRing.openPartitionFunction,
      openHamiltonian_odd _ (Nat.not_even_iff_odd.1 he), ← BoundaryTraces.hamiltonian_neg_one]
    rfl

/-! ### The averaged pairing -/

/-- The averaged approximants `Re ⟨v, X̃ⁿ Π v⟩`. -/
noncomputable def avgM (M n : ℕ) : ℝ :=
  (star (vM M) ⬝ᵥ (tX M ^ n *ᵥ (avgProj spinOneA4 M *ᵥ vM M))).re

/-- The cyclic relabeling moves the end field from `z` to `y` and from `y` to `x`. -/
theorem fieldW_comp_finRotate (c : ℝ) :
    fieldW 2 c ∘ finRotate 3 = fieldW 1 c ∧ fieldW 1 c ∘ finRotate 3 = fieldW 0 c := by
  constructor <;> (funext β; fin_cases β <;> simp [fieldW])

/-- **The average of the right boundary vector over the rotation group**: the twelve rotated
end fields, `±x`, `±y`, `±z`, each twice (written with the signs of the `π` rotations). -/
theorem avgProj_mulVec_vM (M : ℕ) :
    avgProj spinOneA4 M *ᵥ vM M = ((1 / 12 : ℝ) : ℂ) • (
      (capVec (τM M) (fieldW 2 (3 / 5)) + capVec (τM M) (fieldW 2 (rotX 2 * (3 / 5))) +
        capVec (τM M) (fieldW 2 (rotY 2 * (3 / 5))) +
        capVec (τM M) (fieldW 2 (rotX 2 * rotY 2 * (3 / 5)))) +
      (capVec (τM M) (fieldW 1 (3 / 5)) + capVec (τM M) (fieldW 1 (rotX 1 * (3 / 5))) +
        capVec (τM M) (fieldW 1 (rotY 1 * (3 / 5))) +
        capVec (τM M) (fieldW 1 (rotX 1 * rotY 1 * (3 / 5)))) +
      (capVec (τM M) (fieldW 0 (3 / 5)) + capVec (τM M) (fieldW 0 (rotX 0 * (3 / 5))) +
        capVec (τM M) (fieldW 0 (rotY 0 * (3 / 5))) +
        capVec (τM M) (fieldW 0 (rotX 0 * rotY 0 * (3 / 5))))) := by
  obtain ⟨hz, hy⟩ := fieldW_comp_finRotate (3 / 5)
  have hP : ∀ u, cycleMat (m := M) spinOneA4 *ᵥ capVec (τM M) u =
      capVec (τM M) (u ∘ finRotate 3) := fun u ↦ relabelMat_mulVec_capVec _ _ _
  rw [avgProj, ← mulVec_mulVec, projC_mulVec, vM, wR, hP, hz, hP, hy, mulVec_smul,
    mulVec_add, mulVec_add, projD_mulVec, projD_mulVec, projD_mulVec]
  have hx : spinOneA4.x.sign = rotX := rfl
  have hy' : spinOneA4.y.sign = rotY := rfl
  simp only [hx, hy', charMat_mulVec_fieldW]
  simp only [mul_assoc]
  push_cast
  module

/-- The averaged approximant is the average of the twelve rotated pairings. -/
theorem avgM_eq (M n : ℕ) : avgM M n = 1 / 12 * (
    (pairingM M n 2 (3 / 5) + pairingM M n 2 (rotX 2 * (3 / 5)) +
      pairingM M n 2 (rotY 2 * (3 / 5)) + pairingM M n 2 (rotX 2 * rotY 2 * (3 / 5))) +
    (pairingM M n 1 (3 / 5) + pairingM M n 1 (rotX 1 * (3 / 5)) +
      pairingM M n 1 (rotY 1 * (3 / 5)) + pairingM M n 1 (rotX 1 * rotY 1 * (3 / 5))) +
    (pairingM M n 0 (3 / 5) + pairingM M n 0 (rotX 0 * (3 / 5)) +
      pairingM M n 0 (rotY 0 * (3 / 5)) + pairingM M n 0 (rotX 0 * rotY 0 * (3 / 5)))) := by
  rw [avgM, avgProj_mulVec_vM, mulVec_smul, dotProduct_smul, smul_eq_mul,
    Complex.re_ofReal_mul]
  simp only [mulVec_add, dotProduct_add, Complex.add_re]
  rfl

/-- **The averaged pairings converge to `I(n)`** (eq. (appth:average)): the twelve rotated
chains are the chains `Z_{n,1}`, `Z_{n,-1}` twice each and eight transverse chains equal to
`Z_{n,0}`. -/
theorem avgM_tendsto {n : ℕ} (hn : 1 ≤ n) :
    Tendsto (fun M ↦ avgM M n) atTop (𝓝 (BoundaryTraces.moment n)) := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  simp only [avgM_eq]
  have T := fun α c ↦ pairing_tendsto k α c
  have h := ((((((T 2 (3 / 5)).add (T 2 (rotX 2 * (3 / 5)))).add (T 2 (rotY 2 * (3 / 5)))).add
    (T 2 (rotX 2 * rotY 2 * (3 / 5)))).add ((((T 1 (3 / 5)).add (T 1 (rotX 1 * (3 / 5)))).add
      (T 1 (rotY 1 * (3 / 5)))).add (T 1 (rotX 1 * rotY 1 * (3 / 5))))).add
    ((((T 0 (3 / 5)).add (T 0 (rotX 0 * (3 / 5)))).add (T 0 (rotY 0 * (3 / 5)))).add
      (T 0 (rotX 0 * rotY 0 * (3 / 5))))).const_mul (1 / 12)
  convert h using 2
  obtain ⟨h1, h2, h0⟩ := fieldH_eq_hamiltonian k
  obtain ⟨e1, e2, e3⟩ := gibbsTrace_fieldH_transverse k (3 / 5) (49 / 4)
  have r : rotX 2 * (3 / 5 : ℝ) = -(3 / 5) ∧ rotY 2 * (3 / 5 : ℝ) = -(3 / 5) ∧
      rotX 2 * rotY 2 * (3 / 5 : ℝ) = 3 / 5 ∧ rotX 1 * (3 / 5 : ℝ) = -(3 / 5) ∧
      rotY 1 * (3 / 5 : ℝ) = 3 / 5 ∧ rotX 1 * rotY 1 * (3 / 5 : ℝ) = -(3 / 5) ∧
      rotX 0 * (3 / 5 : ℝ) = 3 / 5 ∧ rotY 0 * (3 / 5 : ℝ) = -(3 / 5) ∧
      rotX 0 * rotY 0 * (3 / 5 : ℝ) = -(3 / 5) := by
    norm_num [rotX, rotY]
  obtain ⟨r1, r2, r3, r4, r5, r6, r7, r8, r9⟩ := r
  rw [r1, r2, r3, r4, r5, r6, r7, r8, r9, e1, e2, e3, h1, h2, h0, BoundaryTraces.moment,
    BoundaryTraces.trace_eq, BoundaryTraces.trace_eq, BoundaryTraces.trace_eq]
  ring

/-! ### The sector eigenvalues of `X̃` and the trace table -/

section Sectors

open D2Symmetry

/-- The real transfer matrix `X_M`. -/
noncomputable abbrev TM (M : ℕ) : Matrix (Fin M → Option (Fin 3)) (Fin M → Option (Fin 3)) ℝ :=
  transfer G 0 M (τM M)

/-- The eigenvalues of the trivial `D₂` block of `X_M`. -/
noncomputable abbrev eJ (M : ℕ) : (Fin M → Option (Fin 3)) → ℝ :=
  (spinOneA4.block_isHermitian M (τM M) sectorJ).eigenvalues

/-- The eigenvalues of the block `N` of `X_M`. -/
noncomputable abbrev eN (M : ℕ) : (Fin M → Option (Fin 3)) → ℝ :=
  (spinOneA4.block_isHermitian M (τM M) sectorN).eigenvalues

/-- The sector data of `X̃_M`: the scaled eigenvalues of the two blocks. -/
noncomputable def xData (M : ℕ) :
    (Fin M → Option (Fin 3)) ⊕ (Fin M → Option (Fin 3)) → ℝ :=
  Sum.elim (fun j ↦ sf * eJ M j) (fun j ↦ sf * eN M j)

/-- The two sectors `J` and `N` of the sector data. -/
noncomputable def secData (M : ℕ) :
    Fin (seed36.scheme.m + 1) → Finset ((Fin M → Option (Fin 3)) ⊕ (Fin M → Option (Fin 3))) :=
  ![Finset.univ.map Function.Embedding.inl, Finset.univ.map Function.Embedding.inr]

/-- The sector `J` of the sector data. -/
theorem secData_zero (M : ℕ) : secData M 0 = Finset.univ.map Function.Embedding.inl := rfl

/-- The sector `N` of the sector data. -/
theorem secData_one (M : ℕ) : secData M 1 = Finset.univ.map Function.Embedding.inr := rfl

/-- The power sums of the two sectors are the scaled traces of the blocks. -/
theorem sum_xData_pow (M n : ℕ) :
    ∑ i ∈ secData M 0, xData M i ^ n = sf ^ n * trace (spinOneA4.block M (τM M) sectorJ ^ n) ∧
      ∑ i ∈ secData M 1, xData M i ^ n =
        sf ^ n * trace (spinOneA4.block M (τM M) sectorN ^ n) := by
  rw [secData_zero, secData_one, Finset.sum_map, Finset.sum_map,
    trace_pow_eq_sum_eigenvalues_real (spinOneA4.block_isHermitian M (τM M) sectorJ),
    trace_pow_eq_sum_eigenvalues_real (spinOneA4.block_isHermitian M (τM M) sectorN)]
  simp only [Function.Embedding.inl_apply, Function.Embedding.inr_apply, xData, Sum.elim_inl,
    Sum.elim_inr, mul_pow, ← Finset.mul_sum, and_self]

/-- **The sector traces of `X̃_M`** are the shifted plain and twisted traces of `X_M`. -/
theorem trace_secData (M k : ℕ) :
    seed36.scheme.trace (secData M) (xData M) 0 (k + 1) = sf ^ (k + 1) * trace (TM M ^ (k + 1)) ∧
      seed36.scheme.trace (secData M) (xData M) 1 (k + 1) =
        sf ^ (k + 1) * trace (diagonal (labelChar rotZ) * TM M ^ (k + 1)) := by
  obtain ⟨h0, h1⟩ := sum_xData_pow M (k + 1)
  have hs := spinOneD2.trace_transfer_pow_sectors M k (τM M)
  have ht := spinOneD2.trace_twist_transfer_pow_sectors M k (τM M)
  rw [spinOneD2_sz] at ht
  constructor
  · refine (SectorScheme.d2_trace_zero (secData M) (xData M) (k + 1)).trans ?_
    change ∑ i ∈ secData M 0, xData M i ^ (k + 1) + 3 * ∑ i ∈ secData M 1, xData M i ^ (k + 1) = _
    rw [h0, h1, TM, hs]
    change _ = sf ^ (k + 1) * (trace (spinOneA4.block M (τM M) sectorJ ^ (k + 1)) +
      3 * trace (spinOneA4.block M (τM M) sectorN ^ (k + 1)))
    ring
  · refine (SectorScheme.d2_trace_one (secData M) (xData M) (k + 1)).trans ?_
    change ∑ i ∈ secData M 0, xData M i ^ (k + 1) - ∑ i ∈ secData M 1, xData M i ^ (k + 1) = _
    rw [h0, h1, TM, ht]
    change _ = sf ^ (k + 1) * (trace (spinOneA4.block M (τM M) sectorJ ^ (k + 1)) -
      trace (spinOneA4.block M (τM M) sectorN ^ (k + 1)))
    ring

/-- The certificate's shift factor is `s`. -/
theorem seed36_shift : Real.exp (-((seed36.beta : ℝ) * (seed36.kappa : ℝ))) = sf := by
  simp only [sf, BoundaryInit.shiftFactor]
  congr 1
  norm_num [seed36, kappa]

/-- The shift factor is positive. -/
theorem sf_pos : 0 < sf := BoundaryInit.shiftFactor_pos _

/-- The shift factor is at most `1`. -/
theorem sf_le_one : sf ≤ 1 := by
  simp only [sf, BoundaryInit.shiftFactor]
  exact Real.exp_le_one_iff.2 (by norm_num)

/-- The plain and twisted traces of `X_M` converge to the ring partition functions. -/
theorem eventually_traces {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ M in atTop, ∀ n ∈ seedLengths,
      |trace (TM M ^ n) - partitionFunction n (49 / 4)| ≤ δ ∧
        |trace (diagonal (labelChar rotZ) * TM M ^ n) - twistedPartitionFunction n (49 / 4)| ≤
          δ := by
  have hns : ∀ n ∈ seedLengths, 2 ≤ n :=
    fun _ hn ↦ TableCertificate.two_le_of_mem_lengths seed36_valid hn
  filter_upwards [eventually_abs_trace_transfer_sub_le G 0 (fun _ ↦ 1) (β := 49 / 4)
    (by norm_num) seedLengths hns hδ, eventually_abs_trace_transfer_sub_le G 0 rotZ
    (β := 49 / 4) (by norm_num) seedLengths hns hδ] with M h1 h2 n hn
  refine ⟨?_, ?_⟩
  · have := h1 n hn
    rwa [show labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from funext labelChar_one, diagonal_one,
      Matrix.one_mul, SingleSite.partitionFunction, hamiltonian_plain] at this
  · have := h2 n hn
    rwa [SingleSite.partitionFunction, hamiltonian_twist] at this

/-- **The sector data of `X̃_M` obey the seed-36 table** for large `M`, and their trace at length
`72` exceeds `1`. -/
theorem eventually_holds (htab : TraceTable) (h72 : 1 < shiftedZ 72 (49 / 4)) :
    ∀ᶠ M in atTop, seed36.table.Holds (seed36.scheme.trace (secData M) (xData M)) ∧
      1 < seed36.scheme.trace (secData M) (xData M) 0 (2 * seed36.n₀) := by
  have hstrict := seed36_holdsStrict htab
  obtain ⟨μ, hμ, hmargin⟩ := seed36.table.exists_margin seed36.table.lengths hstrict
  set δ := min μ ((shiftedZ 72 (49 / 4) - 1) / 2)
  have hδ : 0 < δ := lt_min hμ (by linarith)
  filter_upwards [eventually_traces hδ] with M hM
  have hsh : ∀ (n : ℕ) (z : ℝ), shiftedTrace (seed36.beta : ℝ) (seed36.kappa : ℝ) n z =
      sf ^ n * z := fun n z ↦ by rw [shiftedTrace_eq_pow, seed36_shift]
  -- Every tabulated or seed length `n` is approximated within `sⁿ δ ≤ δ`.
  have happ : ∀ n ∈ seedLengths, ∀ g : Fin 2,
      |seed36.scheme.trace (secData M) (xData M) g n -
        shiftedTrace (seed36.beta : ℝ) (seed36.kappa : ℝ) n (seedTraces g n seed36.beta)| ≤ δ := by
    intro n hn g
    obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 :=
      ⟨n - 1, by have := TableCertificate.two_le_of_mem_lengths seed36_valid hn; omega⟩
    obtain ⟨e0, e1⟩ := trace_secData M k
    have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
    have hsn : sf ^ (k + 1) ≤ 1 := pow_le_one₀ sf_pos.le sf_le_one
    have hs0 : 0 ≤ sf ^ (k + 1) := pow_nonneg sf_pos.le _
    obtain ⟨a1, a2⟩ := hM (k + 1) hn
    fin_cases g
    · simp only [Fin.zero_eta, hsh, seedTraces, Matrix.cons_val_zero]
      rw [hb]
      change |seed36.scheme.trace (secData M) (xData M) 0 (k + 1) - _| ≤ δ
      rw [e0, ← mul_sub, abs_mul, abs_of_nonneg hs0]
      nlinarith [abs_nonneg (trace (TM M ^ (k + 1)) - partitionFunction (k + 1) (49 / 4))]
    · simp only [Fin.mk_one, hsh, seedTraces, Matrix.cons_val_one, Matrix.cons_val_zero]
      rw [hb]
      change |seed36.scheme.trace (secData M) (xData M) 1 (k + 1) - _| ≤ δ
      rw [e1, ← mul_sub, abs_mul, abs_of_nonneg hs0]
      nlinarith [abs_nonneg (trace (diagonal (labelChar rotZ) * TM M ^ (k + 1)) -
        twistedPartitionFunction (k + 1) (49 / 4))]
  refine ⟨fun n hn g b ↦ ?_, ?_⟩
  · have hn' : n ∈ seedLengths := Finset.mem_union_left _ (List.mem_toFinset.2 hn)
    have a := abs_le.1 (happ n hn' g)
    have hm := hmargin n hn g b
    have hδμ : δ ≤ μ := min_le_left _ _
    exact ⟨fun hb ↦ by linarith [hm.1 hb], fun hb ↦ by linarith [hm.2 hb]⟩
  · have h72' : (72 : ℕ) ∈ seedLengths := by decide
    have a := abs_le.1 (happ 72 h72' 0)
    have hδ2 : δ ≤ (shiftedZ 72 (49 / 4) - 1) / 2 := min_le_right _ _
    have e : shiftedTrace (seed36.beta : ℝ) (seed36.kappa : ℝ) 72 (seedTraces 0 72 seed36.beta) =
        shiftedZ 72 (49 / 4) := by
      have hb : ((seed36.beta : ℚ) : ℝ) = 49 / 4 := by norm_num [seed36]
      rw [hb]; rfl
    change 1 < seed36.scheme.trace (secData M) (xData M) (0 : Fin 2) 72
    linarith

/-- The bound `r = .999 × .874 = .873126` on the rest of the spectrum of `X̃_M`. -/
noncomputable def rBound : ℝ := 873126 / 1000000

/-- **The leading line of `X̃_M` for large `M`.** The eigenvalue caps and masses of the seed-36
table give an eigenvalue `|λ₀| > .999` of the trivial block, with every other eigenvalue of `X̃_M`
(in the trivial block by the mass, in the three nontrivial blocks by the cap `.87125`) at most
`r = .873126` in modulus. Its eigenvector lies in the trivial `D₂` character and is fixed by
the cyclic relabeling. -/
theorem eventually_leadingLine (htab : TraceTable) (h72 : 1 < shiftedZ 72 (49 / 4)) :
    ∀ᶠ M in atTop, ∃ (e₀ : ℝ) (ψ : (Fin M → Option (Fin 3)) → ℂ),
      IsLeadingLine (tX M) ψ e₀ rBound ∧ 999 / 1000 < |e₀| ∧
        (∀ γ, spinOneA4.key γ ≠ sectorJ → ψ γ = 0) ∧ cycleMat spinOneA4 M *ᵥ ψ = ψ := by
  filter_upwards [eventually_holds htab h72] with M hM
  obtain ⟨htable, h2⟩ := hM
  have hcap := TableCertificate.abs_lt_cap_of_table seed36_valid htable
  have hmass := TableCertificate.sum_pow_le_mass_of_table seed36_valid htable 0
  obtain ⟨i₀, hi₀, hℓ, -, -⟩ := TableCertificate.exists_window_of_table seed36_valid htable h2
  rw [secData_zero, Finset.mem_map] at hi₀
  obtain ⟨j₀, -, rfl⟩ := hi₀
  have hℓ' : (999 / 1000 : ℝ) < |sf * eJ M j₀| := by
    have : seed36.seed.ell = 999 / 1000 := rfl
    rw [this] at hℓ; push_cast at hℓ; exact hℓ
  have hr0 : (0 : ℝ) ≤ rBound := by norm_num [rBound]
  -- The block `N`: the cap `.87125 < r`.
  have hN : ∀ j, |sf * eN M j| ≤ rBound := fun j ↦ by
    have h := hcap 1 (4356234043 / 5000000000) rfl (Sum.inr j) (by
      rw [secData_one]; exact Finset.mem_map_of_mem _ (Finset.mem_univ j))
    simp only [xData, Sum.elim_inr] at h
    have : ((4356234043 / 5000000000 : ℚ) : ℝ) ≤ rBound := by norm_num [rBound]
    linarith
  -- The trivial block: the mass leaves at most `m_J - ℓ¹² < r¹²` to every other eigenvalue.
  have hJ : ∀ j, j ≠ j₀ → |sf * eJ M j| ≤ rBound := fun j hj ↦ by
    have hm : ∑ i ∈ secData M 0, xData M i ^ 12 ≤ 2131481 / 2000000 := by
      have := hmass
      have e : seed36.seed.mass 0 = 2131481 / 2000000 := rfl
      have eq : seed36.seed.q = 12 := rfl
      rw [e, eq] at this; push_cast at this; exact this
    rw [secData_zero, Finset.sum_map] at hm
    simp only [Function.Embedding.inl_apply, xData, Sum.elim_inl] at hm
    have he : Even 12 := ⟨6, rfl⟩
    have hpair : (sf * eJ M j) ^ 12 + (sf * eJ M j₀) ^ 12 ≤ ∑ i, (sf * eJ M i) ^ 12 :=
      calc _ = ∑ i ∈ {j, j₀}, (sf * eJ M i) ^ 12 := (Finset.sum_pair hj).symm
        _ ≤ _ := Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
          fun i _ _ ↦ he.pow_nonneg _
    have h0 : (999 / 1000 : ℝ) ^ 12 < (sf * eJ M j₀) ^ 12 := by
      rw [← he.pow_abs (sf * eJ M j₀)]
      exact pow_lt_pow_left₀ hℓ' (by norm_num) (by norm_num : (12 : ℕ) ≠ 0)
    have h1 : |sf * eJ M j| ^ 12 ≤ rBound ^ 12 := by
      rw [he.pow_abs (sf * eJ M j)]
      have : (2131481 / 2000000 : ℝ) - (999 / 1000) ^ 12 ≤ rBound ^ 12 := by norm_num [rBound]
      linarith
    exact (pow_le_pow_iff_left₀ (abs_nonneg _) hr0 (by norm_num)).1 h1
  have hr : rBound < |sf * eJ M j₀| := by
    have : rBound < 999 / 1000 := by norm_num [rBound]
    linarith
  obtain ⟨ψ, hL, hsupp, hP⟩ := isLeadingLine_of_blocks spinOneA4 M (τM M) sf_pos j₀ hr0 hJ hN hr
  exact ⟨sf * eJ M j₀, ψ, hL, hℓ', hsupp, hP⟩

end Sectors

/-! ### The witness and its limit -/

/-- The lengths `2, …, 6` of the leading-weight polynomial. -/
def polyLengths : Finset ℕ := {2, 3, 4, 5, 6}

/-- The coefficients of `P_*(x)² = x⁶ - .28x⁵ - .3404x⁴ + .0504x³ + .0324x²`. -/
noncomputable def polyCoef : ℕ → ℝ
  | 2 => 324 / 10000
  | 3 => 504 / 10000
  | 4 => -(3404 / 10000)
  | 5 => -(28 / 100)
  | 6 => 1
  | _ => 0

/-- The coefficients are those of `P_*²`. -/
theorem sum_polyCoef (x : ℝ) :
    ∑ n ∈ polyLengths, polyCoef n * x ^ n = (x ^ 3 - 7 / 50 * x ^ 2 - 9 / 50 * x) ^ 2 := by
  norm_num [polyLengths, polyCoef]
  ring

/-- The averaged traces enter through eq. (appth:polynomial). -/
theorem sum_polyCoef_moment :
    ∑ n ∈ polyLengths, polyCoef n * BoundaryTraces.moment n =
      BoundaryTraces.moment 6 - 28 / 100 * BoundaryTraces.moment 5 -
        3404 / 10000 * BoundaryTraces.moment 4 + 504 / 10000 * BoundaryTraces.moment 3 +
          324 / 10000 * BoundaryTraces.moment 2 := by
  norm_num [polyLengths, polyCoef]
  ring

/-- **The approximants satisfy the bundled inequalities** for every large `M`: the boundary
moment inequalities with `μ ≥ .999` and `ρ ≤ .874`, and the leading-weight bound under the
average `Π` over the rotation group. -/
theorem eventually_leadingWitness (htab : TraceTable) (h72 : 1 < shiftedZ 72 (49 / 4)) :
    ∀ᶠ M in atTop, ∃ μ c ρ s, LeadingWitness (ringM M) (openM M) (avgM M) polyLengths polyCoef
      (999 / 1000) (874 / 1000) μ c ρ s := by
  filter_upwards [eventually_leadingLine htab h72] with M hM
  obtain ⟨e₀, ψ, hL, hℓ, hsupp, hP⟩ := hM
  set T : Matrix (Fin M → Option (Fin 3)) (Fin M → Option (Fin 3)) ℝ := sf • TM M
  have htX : tX M = T.map (↑) := rfl
  -- The `π` rotation about `x`.
  set U := charMat (m := M) spinOneA4.x.sign
  have hsign := spinOneA4.x.sign_sq
  have hU : Uᴴ * U = 1 := by rw [charMat_conjTranspose, charMat_mul_self hsign]
  have hTx : ∀ γ γ', T γ γ' = (labelChar spinOneA4.x.sign γ * labelChar spinOneA4.x.sign γ') *
      T γ γ' := fun γ γ' ↦ by
    have h : TM M γ γ' = (labelChar spinOneA4.x.sign γ * labelChar spinOneA4.x.sign γ') *
        TM M γ γ' :=
      transfer_eq_char_mul spinOneA4.x.trace_map spinOneA4.x.map_C spinOneA4.x.map_o M
        (τM M) γ γ'
    simp only [T, Matrix.smul_apply, smul_eq_mul]
    conv_lhs => rw [h]
    ring
  have hUX : U * tX M = tX M * U := by rw [htX]; exact charMat_mul_comm hsign hTx
  have hUψ : U *ᵥ ψ = ψ := by
    rw [charMat_mulVec]
    funext γ
    by_cases hγ : spinOneA4.key γ = D2Symmetry.sectorJ
    · have : labelChar spinOneA4.x.sign γ = 1 := by
        simp only [D2Symmetry.key, D2Symmetry.sectorJ, Prod.mk.injEq, decide_eq_true_eq] at hγ
        exact hγ.1
      rw [this, Complex.ofReal_one, one_mul]
    · rw [hsupp γ hγ, mul_zero]
  -- The boundary moment inequalities.
  have hUv : U *ᵥ vM M = capVec (τM M) (fieldW 2 (-(3 / 5))) := by
    rw [vM, wR, charMat_mulVec_fieldW]
    congr 2
    norm_num [show spinOneA4.x.sign = rotX from rfl, rotX]
  have hr : rBound ≤ |e₀| := by
    have : rBound < 999 / 1000 := by norm_num [rBound]
    linarith
  obtain ⟨ρ, hw, hρ⟩ := hL.openMomentWitness hr (vM M) hU hUX hUψ (Z := ringM M)
    (G := openM M) (fun N _ ↦ rfl)
    (fun N _ ↦ by simp only [openM, even_two_mul, ite_true]; rfl)
    (fun N _ ↦ by
      have hodd : ¬ Even (2 * N + 1) := by simp
      simp only [openM, hodd, ite_false, hUv]; rfl)
  -- The leading-weight bound under `Π`.
  have hPX : avgProj spinOneA4 M * tX M = tX M * avgProj spinOneA4 M := by
    rw [htX]
    refine avgProj_mul_comm spinOneA4 M (fun γ γ' h ↦ ?_) (fun γ γ' ↦ ?_)
    · simp [T, spinOneA4.transfer_eq_zero_of_key_ne (τM M) h]
    · have h : TM M (relabelSeq spinOneA4.cycle.perm γ) (relabelSeq spinOneA4.cycle.perm γ') =
          TM M γ γ' := spinOneA4.cycle.transfer_map M (τM M) γ γ'
      simp only [T, Matrix.smul_apply, smul_eq_mul, h]
  have hpoly := hL.mul_poly_le (vM M) (avgProj_conjTranspose spinOneA4 M)
    (avgProj_mul_self spinOneA4 M) hPX (avgProj_mulVec_of_fixed spinOneA4 M hsupp hP)
    polyLengths polyCoef (fun n hn ↦ by
      simp only [polyLengths, Finset.mem_insert, Finset.mem_singleton] at hn; omega)
    (fun x ↦ by rw [sum_polyCoef]; positivity)
  have hsμ : (if 0 ≤ e₀ then (1 : ℝ) else -1) * |e₀| = e₀ := by
    split_ifs with h
    · rw [abs_of_nonneg h, one_mul]
    · rw [abs_of_neg (not_le.1 h)]; ring
  refine ⟨|e₀|, _, ρ, _, LeadingWitness.mk hw hℓ.le ?_ ?_⟩
  · have hρ0 := hw.rho_nonneg
    have : rBound = 999 / 1000 * (874 / 1000) := by norm_num [rBound]
    nlinarith
  · rw [hsμ]; exact hpoly

/-- The odd open-chain traces are positive. -/
theorem shiftedOpenZ_pos (n : ℕ) : 0 < BoundaryInit.shiftedOpenZ (49 / 4) n :=
  mul_pos (pow_pos (BoundaryInit.shiftFactor_pos _) _) (BoundaryInit.openZ_pos n _)

/-- The ring traces are positive. -/
theorem shiftedRingZ_pos (L : ℕ) : 0 < BoundaryInit.shiftedRingZ (49 / 4) L :=
  mul_pos (pow_pos (BoundaryInit.shiftFactor_pos _) _) (BoundaryInit.partitionFunction_pos' L _)

/-- **The leading boundary weight at `b₀`.** From the seed-36 trace table, the fifteen-trace
table and `Z̃₇₂(b₀) > 1`: a witness of the shifted boundary moment inequalities at `b₀` with
`λ₀ ≥ .999`, `ρ₀ ≤ .874` and `c₀ P_*(λ₀)² ≤ .003745`. -/
theorem witness_of_tables (htab : TraceTable) (hb : BoundaryTraces.BoundaryTraceTable)
    (h72 : 1 < shiftedZ 72 (49 / 4)) :
    ∃ μ c ρ s : ℝ, OpenMomentWitness (BoundaryInit.shiftedRingZ (49 / 4))
      (BoundaryInit.shiftedOpenZ (49 / 4)) μ c ρ s ∧ 999 / 1000 ≤ μ ∧ ρ ≤ 874 / 1000 ∧
        c * (μ ^ 3 - 7 / 50 * μ ^ 2 - 9 / 50 * μ) ^ 2 ≤ 3745 / 10 ^ 6 := by
  obtain ⟨μ, c, ρ, s, w⟩ := LeadingWitness.of_tendsto (eventually_leadingWitness htab h72)
    (fun n _ ↦ ringM_tendsto (by omega)) (fun n _ ↦ openM_tendsto (by omega))
    (fun n hn ↦ avgM_tendsto (by
      simp only [polyLengths, Finset.mem_insert, Finset.mem_singleton] at hn; omega))
    (shiftedRingZ_pos 2) (shiftedRingZ_pos 4)
  have hI : ∑ n ∈ polyLengths, polyCoef n * BoundaryTraces.moment n ≤ 3745 / 10 ^ 6 := by
    rw [sum_polyCoef_moment]
    exact (BoundaryTraces.moment_polynomial_lt_of_table hb).le
  have h := w.mul_poly_le (by norm_num) (by norm_num) (fun N _ ↦ shiftedOpenZ_pos _)
    (by norm_num) hI
  rw [sum_polyCoef] at h
  exact ⟨μ, c, ρ, s, w.toOpenMomentWitness, w.ell_le, w.rho_le, h⟩

end BoundaryLeading

namespace BoundaryInit

/-- **`ThermalCore` from the tables.** The core H1 data of the boundary initialization follow
from the seed-36 trace table, the table of the fifteen short-chain thermal traces and
`e^{-72 b₀ a} Z₇₂(b₀) > 1` (the periodic trial state). -/
theorem thermalCore_of_tables (htab : TraceTable) (hb : BoundaryTraces.BoundaryTraceTable)
    (h72 : 1 < shiftedZ 72 (49 / 4)) : ThermalCore := by
  obtain ⟨μ, c, ρ, s, w, hμ, hρ, hc⟩ := BoundaryLeading.witness_of_tables htab hb h72
  exact ⟨⟨μ, c, ρ, s, w, by linarith, hρ, hc⟩, (BoundaryTraces.sixth_of_table hb).le⟩

/-- **`ThermalCore` from the tables and a trial state**: the form of `thermalCore_of_tables`
with the periodic matrix-product trial state of `Seed36.lean`. -/
theorem thermalCore_of_tables_mps (htab : TraceTable) (hb : BoundaryTraces.BoundaryTraceTable)
    {D : Type*} [Fintype D] [DecidableEq D] (A : Fin 3 → Matrix D D ℤ)
    (hU : 500000 * trace (mpsTransfer A ^ 70 * mpsInsertion A) <
      -700741 * trace (mpsTransfer A ^ 72)) : ThermalCore :=
  thermalCore_of_tables htab hb (one_lt_shiftedZ_of_mps A hU (by norm_num))

end BoundaryInit

end CertifiedGaps.SpinOneHeisenberg
