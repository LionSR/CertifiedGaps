/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Seed.Certificate
import CertifiedGaps.Transfer.SingleSiteLimit
import CertifiedGaps.Transfer.SingleSiteSectors

/-!
# `D₂` sectors of the single-site transfer matrix

The seed certificates of `Seed/Certificate.lean` take as analytic input `D2Approx Z ZP ns`: the
untwisted and twisted ring traces are, to any accuracy, the power sums `J_n + 3 N_n` and
`J_n - N_n` of real numbers in two sectors. This file derives that input for every single-site
ring (`Transfer/SingleSite.lean`) whose bond factors carry a `D₂ = {1, D_x, D_y, D_z}` symmetry
with a cyclic relabeling, as the spin-`S` Heisenberg chains in a Cartesian basis do.

A model supplies a `D2Symmetry C o`:

* two sign symmetries `x`, `y` (`SignSymmetry`): `Φ_x (C_α) = s_x(α) C_α`, and so on. They make
  the transfer matrix `X` block diagonal by the character
  `key γ = (D_x(γ) = 1, D_y(γ) = 1)`, with `D_z = D_x D_y` the character of `s_z = s_x s_y`;
* a label symmetry `cycle` (`LabelSymmetry`) with `s_x ∘ π = s_z` and `s_y ∘ π = s_x`. It
  preserves `X` and permutes the three nontrivial characters cyclically, so their blocks have
  equal traces of powers.

With `J` the trivial character and `N` the one with `D_z = 1`, `D_x = D_y = -1`:

  `Tr X^L = Tr Y_J^L + 3 Tr Y_N^L`,  `Tr (D_z X^L) = Tr Y_J^L - Tr Y_N^L`,

and letting the number of time steps grow (`eventually_abs_trace_transfer_sub_le`) gives
`D2Approx` for the plain ring and the ring with closing-bond signs `s_z`.

## Main results

* `D2Symmetry`: the data a model supplies.
* `D2Symmetry.trace_transfer_pow_sectors`, `D2Symmetry.trace_twist_transfer_pow_sectors`: the two
  sector identities above.
* `D2Symmetry.d2Approx`: the sector moment form of the plain and `s_z`-twisted partition
  functions.
* `RPRing.d2Approx`: the same for a reflection-positive ring without an internal reflection.
-/

open Matrix Filter

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d]

/-- A **`D₂` symmetry with a cyclic relabeling** of the bond factors `C` and the on-site term
`o`: two commuting sign symmetries `x`, `y` (the third is `z = x y`) and a label symmetry whose
permutation `π` turns the `x` signs into the `z` signs and the `y` signs into the `x` signs, so
that it permutes the three nontrivial characters cyclically. -/
structure D2Symmetry (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ) where
  /-- The first sign symmetry (the `π` rotation about `x` for spin chains). -/
  x : SignSymmetry C o
  /-- The second sign symmetry (the `π` rotation about `y` for spin chains). -/
  y : SignSymmetry C o
  /-- The label symmetry cycling the three nontrivial characters (`x → y → z → x`). -/
  cycle : LabelSymmetry C o
  x_cycle : ∀ α, x.sign (cycle.perm α) = x.sign α * y.sign α
  y_cycle : ∀ α, y.sign (cycle.perm α) = x.sign α

/-- The action of the cyclic relabeling on characters `(D_x = 1, D_y = 1)`:
`(D_x, D_y) ↦ (D_z, D_x)`. -/
def d2KeyPerm (c : Bool × Bool) : Bool × Bool := (c.1 == c.2, c.1)

/-- The action on characters is a bijection. -/
theorem d2KeyPerm_injective : Function.Injective d2KeyPerm := by decide

namespace D2Symmetry

variable {C : ι → Matrix d d ℝ} {o : Matrix d d ℝ} (S : D2Symmetry C o)

/-- The signs `s_z = s_x s_y` of the third symmetry (the twist of the twisted ring). -/
def sz (α : ι) : ℝ := S.x.sign α * S.y.sign α

/-- `D_z = D_x D_y`. -/
theorem labelChar_sz {m : ℕ} (γ : Fin m → Option ι) :
    labelChar S.sz γ = labelChar S.x.sign γ * labelChar S.y.sign γ :=
  labelChar_mul _ _ γ

/-- The `D₂` character of a time-label sequence: whether `D_x = 1` and whether `D_y = 1`. -/
noncomputable def key {m : ℕ} (γ : Fin m → Option ι) : Bool × Bool :=
  (decide (labelChar S.x.sign γ = 1), decide (labelChar S.y.sign γ = 1))

/-- The cyclic relabeling of time-label sequences. -/
def relabel (m : ℕ) : (Fin m → Option ι) ≃ (Fin m → Option ι) :=
  Equiv.arrowCongr (Equiv.refl _) (Equiv.optionCongr S.cycle.perm)

/-- The character of a relabeled sequence. -/
theorem key_relabel {m : ℕ} (γ : Fin m → Option ι) :
    S.key (S.relabel m γ) = d2KeyPerm (S.key γ) := by
  have h : ∀ r, labelChar r (S.relabel m γ) = labelChar (r ∘ S.cycle.perm) γ :=
    fun r ↦ labelChar_map_labels r S.cycle.perm γ
  have hx : S.x.sign ∘ S.cycle.perm = S.sz := funext S.x_cycle
  have hy : S.y.sign ∘ S.cycle.perm = S.x.sign := funext S.y_cycle
  simp only [key, d2KeyPerm, h, hx, hy, labelChar_sz]
  rcases labelChar_eq_one_or_neg_one S.x.sign_sq γ with h1 | h1 <;>
    rcases labelChar_eq_one_or_neg_one S.y.sign_sq γ with h2 | h2 <;> simp [h1, h2]

/-- The transfer matrix is block diagonal by the `D₂` character. -/
theorem transfer_eq_zero_of_key_ne {m : ℕ} (τ : ℝ) {γ γ' : Fin m → Option ι}
    (h : S.key γ ≠ S.key γ') : transfer C o m τ γ γ' = 0 := by
  simp only [key, ne_eq, Prod.mk.injEq, not_and_or, decide_eq_decide] at h
  rcases h with h | h
  · exact S.x.transfer_eq_zero m τ fun e ↦ h (by rw [e])
  · exact S.y.transfer_eq_zero m τ fun e ↦ h (by rw [e])

/-- The block of the transfer matrix on one `D₂` character, zero elsewhere. -/
noncomputable abbrev block (m : ℕ) (τ : ℝ) (c : Bool × Bool) :
    Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ :=
  blockOf (transfer C o m τ) S.key c

/-- Each block is real symmetric. -/
theorem block_isHermitian (m : ℕ) (τ : ℝ) (c : Bool × Bool) : (S.block m τ c).IsHermitian :=
  blockOf_isHermitian (transfer_isHermitian _ _ m τ) _ c

variable [Fintype ι] [DecidableEq ι]

/-- The cyclic relabeling preserves the traces of powers of the blocks. -/
theorem trace_block_keyPerm_pow (m : ℕ) (τ : ℝ) (c : Bool × Bool) (L : ℕ) :
    trace (S.block m τ (d2KeyPerm c) ^ L) = trace (S.block m τ c ^ L) :=
  trace_blockOf_pow_of_reindex _ _ (S.relabel m) d2KeyPerm_injective
    (S.cycle.transfer_map m τ) S.key_relabel c L

/-- The trivial character. -/
def sectorJ : Bool × Bool := (true, true)

/-- The nontrivial character with `D_z = 1`, `D_x = D_y = -1`. -/
def sectorN : Bool × Bool := (false, false)

/-- `Tr (diag(φ ∘ key) X^(n+1))` in blocks, with the three nontrivial blocks identified. -/
theorem trace_diagonal_mul_transfer_pow (m n : ℕ) (τ : ℝ) (φ : Bool × Bool → ℝ) :
    trace (diagonal (fun γ ↦ φ (S.key γ)) * transfer C o m τ ^ (n + 1)) =
      φ sectorJ * trace (S.block m τ sectorJ ^ (n + 1)) +
        (φ (false, true) + φ (true, false) + φ sectorN) *
          trace (S.block m τ sectorN ^ (n + 1)) := by
  rw [trace_diagonal_mul_pow_eq_sum_blockOf (transfer C o m τ) S.key
    (fun _ _ h ↦ S.transfer_eq_zero_of_key_ne τ h) φ n]
  have h1 := S.trace_block_keyPerm_pow m τ (false, false) (n + 1)
  have h2 := S.trace_block_keyPerm_pow m τ (true, false) (n + 1)
  simp only [d2KeyPerm] at h1 h2
  norm_num at h1 h2
  simp only [Fintype.sum_prod_type, Fintype.sum_bool, sectorJ, sectorN, h2, h1]
  ring

/-- **Plain ring in sectors.** `Tr X^(n+1) = Tr Y_J^(n+1) + 3 Tr Y_N^(n+1)`. -/
theorem trace_transfer_pow_sectors (m n : ℕ) (τ : ℝ) :
    trace (transfer C o m τ ^ (n + 1)) =
      trace (S.block m τ sectorJ ^ (n + 1)) + 3 * trace (S.block m τ sectorN ^ (n + 1)) := by
  have h := S.trace_diagonal_mul_transfer_pow m n τ fun _ ↦ 1
  rw [diagonal_one, Matrix.one_mul] at h
  rw [h]
  ring

/-- **Twisted ring in sectors.** `Tr (D_z X^(n+1)) = Tr Y_J^(n+1) - Tr Y_N^(n+1)`. -/
theorem trace_twist_transfer_pow_sectors (m n : ℕ) (τ : ℝ) :
    trace (diagonal (labelChar S.sz) * transfer C o m τ ^ (n + 1)) =
      trace (S.block m τ sectorJ ^ (n + 1)) - trace (S.block m τ sectorN ^ (n + 1)) := by
  let φ : Bool × Bool → ℝ := fun c ↦ (if c.1 then 1 else -1) * (if c.2 then 1 else -1)
  have hφ : labelChar S.sz = fun γ : Fin m → Option ι ↦ φ (S.key γ) := by
    ext γ
    simp only [φ, key, labelChar_sz]
    rcases labelChar_eq_one_or_neg_one S.x.sign_sq γ with h1 | h1 <;>
      rcases labelChar_eq_one_or_neg_one S.y.sign_sq γ with h2 | h2 <;> norm_num [h1, h2]
  rw [hφ, S.trace_diagonal_mul_transfer_pow]
  norm_num [φ, sectorJ, sectorN]
  ring

omit [DecidableEq ι] in
/-- **Sector moment form.** For a `D₂` symmetry with a cyclic relabeling, every `β ≥ 0` and
every finite set of ring lengths `n ≥ 2` (odd lengths included), the plain partition functions
`Z_n(β)` and the partition functions `Z_n(β, s_z)` of the ring twisted by `s_z` are, to any
accuracy, `J_n + 3 N_n` and `J_n - N_n`: power sums of the eigenvalues of the two sector blocks
of the transfer matrix of finitely many time steps. -/
theorem d2Approx {β : ℝ} (hβ : 0 ≤ β) (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    D2Approx (fun n ↦ partitionFunction C o n (fun _ ↦ 1) β)
      (fun n ↦ partitionFunction C o n S.sz β) ns := by
  classical
  intro δ hδ
  obtain ⟨M, hM⟩ := ((eventually_abs_trace_transfer_sub_le C o (fun _ ↦ 1) hβ ns hns hδ).and
    (eventually_abs_trace_transfer_sub_le C o S.sz hβ ns hns hδ)).exists
  have hJ := S.block_isHermitian M (β / M) sectorJ
  have hN := S.block_isHermitian M (β / M) sectorN
  refine ⟨(Fin M → Option ι) ⊕ (Fin M → Option ι), inferInstance,
    Finset.univ.map Function.Embedding.inl, Finset.univ.map Function.Embedding.inr,
    Sum.elim hJ.eigenvalues hN.eigenvalues, fun n hn ↦ ?_⟩
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by have := hns n hn; omega⟩
  have sJ : ∑ i ∈ Finset.univ.map Function.Embedding.inl,
      Sum.elim hJ.eigenvalues hN.eigenvalues i ^ (k + 1) =
        trace (S.block M (β / M) sectorJ ^ (k + 1)) := by
    rw [Finset.sum_map, trace_pow_eq_sum_eigenvalues_real hJ]; rfl
  have sN : ∑ i ∈ Finset.univ.map Function.Embedding.inr,
      Sum.elim hJ.eigenvalues hN.eigenvalues i ^ (k + 1) =
        trace (S.block M (β / M) sectorN ^ (k + 1)) := by
    rw [Finset.sum_map, trace_pow_eq_sum_eigenvalues_real hN]; rfl
  have h1 := hM.1 (k + 1) hn
  have h2 := hM.2 (k + 1) hn
  rw [show labelChar (fun _ ↦ (1 : ℝ)) = fun _ ↦ 1 from funext labelChar_one, diagonal_one,
    Matrix.one_mul, S.trace_transfer_pow_sectors] at h1
  rw [S.trace_twist_transfer_pow_sectors] at h2
  rw [sJ, sN]
  exact ⟨h1, h2⟩

end D2Symmetry

end CertifiedGaps.SingleSite

/-! ### Reflection-positive rings without an internal reflection -/

namespace CertifiedGaps.RPRing

open SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d] [Fintype ι] (R : RPRing d ι)

/-- **Sector moment form of a reflection-positive ring.** For an `RPRing` without an internal
reflection whose bond factors and on-site term carry a `D₂` symmetry with a cyclic relabeling,
its partition functions and those of the ring twisted by `s_z` satisfy `D2Approx`. -/
theorem d2Approx (h : R.refl = 1) (S : D2Symmetry R.C R.o) {β : ℝ} (hβ : 0 ≤ β)
    (ns : Finset ℕ) (hns : ∀ n ∈ ns, 2 ≤ n) :
    D2Approx (fun n ↦ R.partitionFunction n β)
      (fun n ↦ SingleSite.partitionFunction R.C R.o n S.sz β) ns := by
  classical
  have hZ : (fun n ↦ R.partitionFunction n β) =
      fun n ↦ SingleSite.partitionFunction R.C R.o n (fun _ ↦ 1) β := by
    ext n
    rw [partitionFunction, hamiltonian_eq_singleSite R h]
    rfl
  rw [hZ]
  exact S.d2Approx hβ ns hns

end CertifiedGaps.RPRing
