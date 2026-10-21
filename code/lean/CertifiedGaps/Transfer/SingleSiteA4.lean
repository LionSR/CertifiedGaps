/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.CyclicSectors
import CertifiedGaps.Transfer.SingleSiteD2

/-!
# `A₄` sectors of the single-site transfer matrix

A `D2Symmetry` (`Transfer/SingleSiteD2.lean`) splits the transfer matrix `X` into the trivial
`D₂` block `Y_J` and three nontrivial blocks with equal traces of powers, `N`. If its cyclic
relabeling `π` has order three, the symmetry group is `A₄ = D₂ ⋊ ℤ₃`, and a third twisted trace
resolves `Y_J` further. The ring twisted by `π` (closing bond `-∑_α C_α ⊗ C_{π α}`) has, by
`trace_stepOp_pow`, the transfer form `∑_γ (X^L) (π γ, γ) = Tr (Π X^L)` with `Π` the permutation
matrix of the relabeling. Since `π` moves the three nontrivial `D₂` characters cyclically, only the
trivial block contributes (`sum_mul_pow_apply_eq_blockOf`):

  `Tr (Π X^L) = Tr (Π Y_J^L)`.

On `Y_J` the relabeling is a symmetry of order three, so `exists_cyclic_sectors` gives real
numbers `aᵢ` (sector `I`, fixed by `π`) and `bᵢ` (sector `O`, one of the two complex-conjugate
characters `ω`, `ω̄`) with `Tr Y_J^L = I_L + 2 O_L` and `Tr (Π Y_J^L) = I_L - O_L`. With the `D₂`
identities this is the scheme `SectorScheme.ion`:

  `Z_L = I_L + 2 O_L + 3 N_L`,  `Z_L(P) = I_L + 2 O_L - N_L`,  `Z_L(C) = I_L - O_L`,

where `P` is the sign twist `s_z` and `C` the relabeling twist `π`.

## Main results

* `SectorApprox.congr`: sector data transfer along traces equal at the tabulated lengths.
* `A4Symmetry`: a `D2Symmetry` whose cyclic relabeling has order three.
* `A4Symmetry.sum_twist_transfer_pow_eq_block`: the `π`-twisted trace sees only the trivial block.
* `A4Symmetry.sectorApprox_ion`: the sector representation `SectorApprox SectorScheme.ion` of the
  plain, `s_z`-twisted and `π`-twisted partition functions.
* `RPRing.sectorApprox_ion`: the same for a reflection-positive ring without an internal
  reflection.
-/

open Matrix Filter Finset

namespace CertifiedGaps

/-- Sector data for some traces are sector data for any traces equal to them at the lengths `ns`. -/
theorem SectorApprox.congr {S : SectorScheme} {Z Z' : Fin (S.t + 1) → ℕ → ℝ} {ns : Finset ℕ}
    (h : SectorApprox S Z ns) (hZ : ∀ g, ∀ n ∈ ns, Z g n = Z' g n) : SectorApprox S Z' ns :=
  fun δ hδ ↦ by
    obtain ⟨ι, i, sec, x, hx⟩ := h δ hδ
    exact ⟨ι, i, sec, x, fun n hn g ↦ by rw [← hZ g n hn]; exact hx n hn g⟩

namespace SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d]

/-- An **`A₄` symmetry** of the bond factors `C` and the on-site term `o`: a `D₂` symmetry with
a cyclic relabeling `π` of order three (for spin chains, the cyclic permutation `x → y → z` of
the axes, a rotation by `2π/3`). -/
structure A4Symmetry (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ) extends D2Symmetry C o where
  cycle_cube : ∀ α, cycle.perm (cycle.perm (cycle.perm α)) = α

/-- The action of the cyclic relabeling fixes only the trivial `D₂` character. -/
theorem d2KeyPerm_eq_self_iff (c : Bool × Bool) : d2KeyPerm c = c ↔ c = D2Symmetry.sectorJ := by
  revert c; decide

namespace A4Symmetry

variable {C : ι → Matrix d d ℝ} {o : Matrix d d ℝ} (S : A4Symmetry C o)

/-- The cyclic relabeling of sequences relabels every time label. -/
theorem relabel_apply {m : ℕ} (γ : Fin m → Option ι) :
    S.relabel m γ = fun t ↦ (γ t).map S.cycle.perm := rfl

/-- The cyclic relabeling of sequences has order three. -/
theorem relabel_cube (m : ℕ) (γ : Fin m → Option ι) :
    S.relabel m (S.relabel m (S.relabel m γ)) = γ := by
  funext t
  simp only [relabel_apply, Option.map_map]
  rcases γ t with _ | α
  · rfl
  · exact congrArg some (S.cycle_cube α)

/-- The cyclic relabeling preserves the trivial block. -/
theorem block_relabel (m : ℕ) (τ : ℝ) (γ γ' : Fin m → Option ι) :
    S.block m τ D2Symmetry.sectorJ (S.relabel m γ) (S.relabel m γ') =
      S.block m τ D2Symmetry.sectorJ γ γ' := by
  have hJ : ∀ γ : Fin m → Option ι, S.key (S.relabel m γ) = D2Symmetry.sectorJ ↔
      S.key γ = D2Symmetry.sectorJ := fun γ ↦ by
    rw [S.key_relabel]
    constructor
    · intro h
      have := d2KeyPerm_injective (h.trans ((d2KeyPerm_eq_self_iff _).2 rfl).symm)
      exact this
    · intro h; rw [h]; rfl
  simp only [blockOf, of_apply, hJ]
  congr 1
  exact S.cycle.transfer_map m τ γ γ'

variable [Fintype ι] [DecidableEq ι]

/-- **The relabeling twist sees only the trivial block.** The `π`-twisted marked trace of the
transfer matrix is that of its trivial `D₂` block: `∑_γ (X^(n+1)) (π γ, γ) = Tr (Π Y_J^(n+1))`. -/
theorem sum_twist_transfer_pow_eq_block (m n : ℕ) (τ : ℝ) :
    ∑ γ, labelChar (fun _ ↦ (1 : ℝ)) γ *
        (transfer C o m τ ^ (n + 1)) (fun t ↦ (γ t).map S.cycle.perm) γ =
      ∑ γ, (S.block m τ D2Symmetry.sectorJ ^ (n + 1)) (S.relabel m γ) γ := by
  have hr : ∀ γ : Fin m → Option ι,
      (transfer C o m τ ^ (n + 1)) (fun t ↦ (γ t).map S.cycle.perm) γ =
        (transfer C o m τ ^ (n + 1)) (S.relabel m γ) γ := fun _ ↦ rfl
  simp only [hr]
  rw [sum_mul_pow_apply_eq_blockOf (T := transfer C o m τ)
    (fun _ _ h ↦ S.transfer_eq_zero_of_key_ne τ h) (g := S.relabel m) (κ' := d2KeyPerm)
    S.key_relabel d2KeyPerm_eq_self_iff]
  simp only [labelChar_one, one_mul]

omit [DecidableEq ι] in
/-- **`A₄` sector moment form.** For an `A₄` symmetry, every `β ≥ 0` and every finite set of
ring lengths `n ≥ 2` (odd lengths included), the plain partition functions, those twisted by the
signs `s_z` and those twisted by the cyclic relabeling `π` are, to any accuracy,
`I_n + 2 O_n + 3 N_n`, `I_n + 2 O_n - N_n` and `I_n - O_n`: power sums of real numbers in three
sectors. -/
theorem sectorApprox_ion {β : ℝ} (hβ : 0 ≤ β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    SectorApprox SectorScheme.ion
      ![fun n ↦ partitionFunction C o n (fun _ ↦ 1) β, fun n ↦ partitionFunction C o n S.sz β,
        fun n ↦ twistPartitionFunction C o n (fun _ ↦ 1) S.cycle.perm β] ns := by
  classical
  intro δ hδ
  obtain ⟨M, hM⟩ := (((eventually_abs_trace_transfer_sub_le C o (fun _ ↦ 1) hβ ns hns hδ).and
    (eventually_abs_trace_transfer_sub_le C o S.sz hβ ns hns hδ)).and
    (eventually_abs_twisted_transfer_sub_le C o (fun _ ↦ 1) S.cycle.perm hβ ns hns hδ)).exists
  set τ := β / M
  obtain ⟨a, b, hab⟩ := exists_cyclic_sectors (S.block_isHermitian M τ D2Symmetry.sectorJ)
    (S.relabel M) (S.relabel_cube M) (S.block_relabel M τ)
  have hN := S.block_isHermitian M τ D2Symmetry.sectorN
  -- One copy of the index set of the transfer matrix for each of the sectors `I`, `O`, `N`.
  let E := (Fin M → Option ι) ⊕ (Fin M → Option ι) ⊕ (Fin M → Option ι)
  let x : E → ℝ := Sum.elim a (Sum.elim b hN.eigenvalues)
  let sec : Fin (SectorScheme.ion.m + 1) → Finset E :=
    ![univ.map Function.Embedding.inl,
      (univ.map Function.Embedding.inl).map Function.Embedding.inr,
      (univ.map Function.Embedding.inr).map Function.Embedding.inr]
  refine ⟨E, inferInstance, sec, x, fun n hn g ↦ ?_⟩
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by have := hns n hn; omega⟩
  have sI : ∑ i ∈ sec 0, x i ^ (k + 1) = ∑ i, a i ^ (k + 1) := by
    change ∑ i ∈ univ.map Function.Embedding.inl, x i ^ (k + 1) = _
    rw [Finset.sum_map]; rfl
  have sO : ∑ i ∈ sec 1, x i ^ (k + 1) = ∑ i, b i ^ (k + 1) := by
    change ∑ i ∈ (univ.map Function.Embedding.inl).map Function.Embedding.inr, x i ^ (k + 1) = _
    rw [Finset.sum_map, Finset.sum_map]; rfl
  have sN : ∑ i ∈ sec 2, x i ^ (k + 1) =
      trace (S.block M τ D2Symmetry.sectorN ^ (k + 1)) := by
    rw [trace_pow_eq_sum_eigenvalues_real hN]
    change ∑ i ∈ (univ.map Function.Embedding.inr).map Function.Embedding.inr, x i ^ (k + 1) = _
    rw [Finset.sum_map, Finset.sum_map]; rfl
  obtain ⟨hJ, hC⟩ := hab (k + 1) (by omega)
  obtain ⟨⟨h1, h2⟩, h3⟩ := hM
  have e1 := h1 (k + 1) hn
  have e2 := h2 (k + 1) hn
  have e3 := h3 (k + 1) hn
  rw [show labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from funext labelChar_one, diagonal_one,
    Matrix.one_mul, S.trace_transfer_pow_sectors, hJ] at e1
  rw [S.trace_twist_transfer_pow_sectors, hJ] at e2
  rw [S.sum_twist_transfer_pow_eq_block, hC] at e3
  fin_cases g
  · change |SectorScheme.ion.trace sec x 0 (k + 1) - partitionFunction C o (k + 1) (fun _ ↦ 1) β|
      ≤ δ
    rw [SectorScheme.ion_trace_zero, sI, sO, sN]
    simpa [add_assoc] using e1
  · change |SectorScheme.ion.trace sec x 1 (k + 1) - partitionFunction C o (k + 1) S.sz β| ≤ δ
    rw [SectorScheme.ion_trace_one, sI, sO, sN]
    simpa [add_assoc] using e2
  · change |SectorScheme.ion.trace sec x 2 (k + 1) -
      twistPartitionFunction C o (k + 1) (fun _ ↦ 1) S.cycle.perm β| ≤ δ
    rw [SectorScheme.ion_trace_two, sI, sO]
    simpa using e3

end A4Symmetry

end SingleSite

/-! ### Reflection-positive rings without an internal reflection -/

namespace RPRing

open SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- **`A₄` sector moment form of a reflection-positive ring.** For an `RPRing` without an
internal reflection whose bond factors and on-site term carry an `A₄` symmetry, its partition
functions and those of the rings twisted by the signs `s_z` and by the cyclic relabeling `π`
satisfy `SectorApprox SectorScheme.ion`. -/
theorem sectorApprox_ion (h : R.refl = 1) (S : A4Symmetry R.C R.o) {β : ℝ} (hβ : 0 ≤ β)
    (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    SectorApprox SectorScheme.ion
      ![fun n ↦ R.partitionFunction n β, fun n ↦ SingleSite.partitionFunction R.C R.o n S.sz β,
        fun n ↦ twistPartitionFunction R.C R.o n (fun _ ↦ 1) S.cycle.perm β] ns := by
  have hZ : (fun n ↦ R.partitionFunction n β) =
      fun n ↦ SingleSite.partitionFunction R.C R.o n (fun _ ↦ 1) β := by
    ext n
    rw [partitionFunction, hamiltonian_eq_singleSite R h]
    rfl
  rw [hZ]
  exact S.sectorApprox_ion hβ ns hns

end RPRing

end CertifiedGaps
