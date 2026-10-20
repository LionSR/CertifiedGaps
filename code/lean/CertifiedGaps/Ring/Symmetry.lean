/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Ring.Hamiltonian

/-!
# Signed-permutation symmetries of ring Hamiltonians

Thermal traces computed from representative columns (`Trace/Frobenius.lean`) need symmetries
that map a basis vector `e_τ` to a phase times another basis vector `e_{e τ}`: site permutations
(rotation and reflection of the ring), on-site relabelings (such as `m ↦ -m` for spin one) and
diagonal sign gauges. This file states them entrywise, for matrices over any commutative
semiring, so that an integer Hamiltonian can be checked on its integer entries and then cast to
`ℂ`.

* `Intertwines e u H H'`: `H' (e σ) (e τ) u_τ = u_σ H σ τ`; it composes (`Intertwines.trans`),
  adds up (`Intertwines.sum`) and maps along ring homomorphisms (`Intertwines.map`).
* `signedPerm e u`, the matrix of `e_τ ↦ u_τ e_{e τ}`: it intertwines `H` with `H'` as
  matrices (`signedPerm_mul`), commutes with `H` when `H' = H` (`commute_signedPerm`) and is
  unitary over `ℂ` for unit phases (`signedPerm_conjTranspose_mul_self`).
* On a ring of `L` sites: `relabel π` (`σ ↦ σ ∘ π`) moves `bondOp j k B` to
  `bondOp (π j) (π k) B` (`bondOp_relabel`); `onsite f` (`σ ↦ f ∘ σ`) conjugates the bond
  (`bondOp_onsite`); `bondOp j k B = bondOp k j (B with its two sites swapped)` (`bondOp_swap`);
  and a sign `p (σ s)` at one site of the bond turns `B` into `B'` when
  `B' (a, b) (c, e) p_e = p_b B (a, b) (c, e)` (`intertwines_bondOp_gauge_right`, `_left`).
* `seamHamiltonian B Bc L m`: the ring whose bond `m → m + 1` carries `Bc`; the twisted ring
  `twistedRingHamiltonian B Bc L` is the case `m = L - 1` (`twistedRingHamiltonian_eq_seam`).
  The rotation moves the seam by one bond (`intertwines_seam_rotate`), the reflection maps the
  seam `m` to `rev (m + 1)` (`intertwines_seam_reverse`), and a sign gauge on the site `s` moves
  the seam between the two bonds at `s` (`intertwines_seam_gauge`). Composing the rotation with
  the gauge gives the signed rotation symmetry of a twisted ring (`intertwines_twisted_rotate`),
  as in `thermal72.tex` (the unitary `G_k`).
-/

open Matrix

namespace CertifiedGaps

/-! ### Intertwining by signed permutations -/

section Intertwines

variable {ι R : Type*}

/-- `(e, u)` **intertwines** `H` with `H'`: `H' (e σ) (e τ) u_τ = u_σ H σ τ` for all labels.
For a bijection `e` this says `P H = H' P` for the signed permutation `P e_τ = u_τ e_{e τ}`
(`signedPerm_mul`); for `H' = H` and unit `u`, `P` is a unitary symmetry of `H`. -/
def Intertwines [Mul R] (e : ι → ι) (u : ι → R) (H H' : Matrix ι ι R) : Prop :=
  ∀ σ τ, H' (e σ) (e τ) * u τ = u σ * H σ τ

/-- Intertwiners compose: first `(e₁, u₁)`, then `(e₂, u₂)`. -/
theorem Intertwines.trans [CommSemigroup R] {e₁ e₂ : ι → ι} {u₁ u₂ : ι → R}
    {H₁ H₂ H₃ : Matrix ι ι R} (h₁ : Intertwines e₁ u₁ H₁ H₂) (h₂ : Intertwines e₂ u₂ H₂ H₃) :
    Intertwines (e₂ ∘ e₁) (fun σ ↦ u₂ (e₁ σ) * u₁ σ) H₁ H₃ := fun σ τ ↦ by
  simp only [Function.comp_apply]
  rw [← mul_assoc, h₂, mul_assoc, h₁, mul_assoc]

/-- A sign gauge (`e = id`) with `u_σ² = 1` intertwines in both directions. -/
theorem Intertwines.symm_of_id [CommMonoid R] {u : ι → R} (hu : ∀ σ, u σ * u σ = 1)
    {H H' : Matrix ι ι R} (h : Intertwines id u H H') : Intertwines id u H' H := fun σ τ ↦ by
  have h' := h σ τ
  simp only [id] at h' ⊢
  calc H σ τ * u τ = u σ * u σ * H σ τ * u τ := by rw [hu, one_mul]
    _ = u σ * (u σ * H σ τ) * u τ := by rw [mul_assoc (u σ)]
    _ = u σ * (H' σ τ * u τ) * u τ := by rw [h']
    _ = u σ * H' σ τ * (u τ * u τ) := by simp only [mul_assoc]
    _ = u σ * H' σ τ := by rw [hu, mul_one]

/-- Intertwiners add up. -/
theorem Intertwines.sum [CommSemiring R] {κ : Type*} (s : Finset κ) {e : ι → ι} {u : ι → R}
    {H H' : κ → Matrix ι ι R} (h : ∀ k ∈ s, Intertwines e u (H k) (H' k)) :
    Intertwines e u (∑ k ∈ s, H k) (∑ k ∈ s, H' k) := fun σ τ ↦ by
  simp only [Matrix.sum_apply, Finset.sum_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun k hk ↦ h k hk σ τ

/-- Intertwiners map along ring homomorphisms, for example from `ℤ` to `ℂ`. -/
theorem Intertwines.map {S : Type*} [Semiring R] [Semiring S] (f : R →+* S) {e : ι → ι}
    {u : ι → R} {H H' : Matrix ι ι R} (h : Intertwines e u H H') :
    Intertwines e (fun σ ↦ f (u σ)) (H.map f) (H'.map f) := fun σ τ ↦ by
  simp only [map_apply, ← map_mul, h σ τ]

/-- The signed permutation matrix `e_τ ↦ u_τ e_{e τ}`. -/
def signedPerm [Zero R] [DecidableEq ι] (e : ι ≃ ι) (u : ι → R) : Matrix ι ι R :=
  of fun σ τ ↦ if σ = e τ then u τ else 0

variable [Fintype ι] [DecidableEq ι]

/-- `P e_r = u_r e_{e r}`. -/
theorem signedPerm_mulVec_single [Semiring R] (e : ι ≃ ι) (u : ι → R) (r : ι) :
    signedPerm e u *ᵥ Pi.single r 1 = u r • Pi.single (e r) 1 := by
  ext σ
  simp [signedPerm, Pi.single_apply]

/-- **Intertwining as a matrix identity**: `P H = H' P`. -/
theorem signedPerm_mul [CommSemiring R] {e : ι ≃ ι} {u : ι → R} {H H' : Matrix ι ι R}
    (h : Intertwines e u H H') : signedPerm e u * H = H' * signedPerm e u := by
  ext σ τ
  have hσ : ∀ ρ, (σ = e ρ) ↔ (e.symm σ = ρ) := fun ρ ↦ e.symm_apply_eq.symm
  simp only [mul_apply, signedPerm, of_apply, ite_mul, zero_mul, mul_ite, mul_zero, hσ,
    Finset.sum_ite_eq, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  have := h (e.symm σ) τ
  rw [Equiv.apply_symm_apply] at this
  rw [this]

/-- A signed permutation that intertwines `H` with itself commutes with `H`. -/
theorem commute_signedPerm [CommSemiring R] {e : ι ≃ ι} {u : ι → R} {H : Matrix ι ι R}
    (h : Intertwines e u H H) : Commute (signedPerm e u) H :=
  signedPerm_mul h

/-- With unit phases the signed permutation is unitary. -/
theorem signedPerm_conjTranspose_mul_self {e : ι ≃ ι} {u : ι → ℂ} (hu : ∀ σ, ‖u σ‖ = 1) :
    (signedPerm e u)ᴴ * signedPerm e u = 1 := by
  ext σ τ
  simp only [mul_apply, conjTranspose_apply, signedPerm, of_apply]
  rw [Finset.sum_eq_single (e σ) (fun ρ _ hρ ↦ by simp [hρ]) (by simp)]
  simp only [ite_true, e.injective.eq_iff, one_apply]
  split_ifs with h
  · subst h
    rw [RCLike.star_def, RCLike.conj_mul, hu]; simp
  · simp

end Intertwines

/-! ### Site permutations, on-site relabelings and sign gauges of a bond -/

section Bond

variable {d : Type} [DecidableEq d] {L : ℕ} {R : Type*}

/-- Relabeling the sites by `π`: the configuration `σ ↦ σ ∘ π`. -/
def relabel (π : Equiv.Perm (Fin L)) : Equiv.Perm (Fin L → d) where
  toFun σ := σ ∘ π
  invFun σ := σ ∘ π.symm
  left_inv σ := by ext i; simp
  right_inv σ := by ext i; simp

omit [DecidableEq d] in
/-- The relabeled configuration at site `i` is the letter at `π i`. -/
@[simp] theorem relabel_apply (π : Equiv.Perm (Fin L)) (σ : Fin L → d) (i : Fin L) :
    relabel π σ i = σ (π i) := rfl

/-- Relabeling the letters on every site by `f`: the configuration `σ ↦ f ∘ σ`. -/
def onsite (f : Equiv.Perm d) : Equiv.Perm (Fin L → d) where
  toFun σ := f ∘ σ
  invFun σ := f.symm ∘ σ
  left_inv σ := by ext i; simp
  right_inv σ := by ext i; simp

omit [DecidableEq d] in
/-- The relabeled configuration at site `i` is `f` of the letter at `i`. -/
@[simp] theorem onsite_apply (f : Equiv.Perm d) (σ : Fin L → d) (i : Fin L) :
    onsite f σ i = f (σ i) := rfl

/-- **A site permutation moves the bond**: `B` on `j, k` read through `σ ∘ π` is `B` on
`π j, π k`. -/
theorem bondOp_relabel [Zero R] (π : Equiv.Perm (Fin L)) (j k : Fin L)
    (B : Matrix (d × d) (d × d) R) (σ τ : Fin L → d) :
    bondOp j k B (relabel π σ) (relabel π τ) = bondOp (π j) (π k) B σ τ := by
  have hc : (∀ i, i ≠ j → i ≠ k → σ (π i) = τ (π i)) ↔
      (∀ i, i ≠ π j → i ≠ π k → σ i = τ i) := by
    refine ⟨fun h i hj hk ↦ ?_, fun h i hj hk ↦ h (π i) (π.injective.ne hj) (π.injective.ne hk)⟩
    have := h (π.symm i) (fun e ↦ hj (by rw [← e, Equiv.apply_symm_apply]))
      (fun e ↦ hk (by rw [← e, Equiv.apply_symm_apply]))
    rwa [Equiv.apply_symm_apply] at this
  simp only [bondOp, of_apply, relabel_apply]
  exact if_congr hc rfl rfl

/-- **An on-site relabeling conjugates the bond**: `B` read through `f ∘ σ` is the matrix
`B (f a, f b) (f c, f e)`. -/
theorem bondOp_onsite [Zero R] (f : Equiv.Perm d) (j k : Fin L) (B : Matrix (d × d) (d × d) R)
    (σ τ : Fin L → d) :
    bondOp j k B (onsite f σ) (onsite f τ) =
      bondOp j k (B.submatrix (Prod.map f f) (Prod.map f f)) σ τ := by
  simp only [bondOp, of_apply, onsite_apply, f.injective.eq_iff, submatrix_apply, Prod.map]

/-- The two sites of a bond can be listed in either order, with the factors of `B` swapped. -/
theorem bondOp_swap [Zero R] (j k : Fin L) (B : Matrix (d × d) (d × d) R) :
    bondOp j k B = bondOp k j (B.submatrix Prod.swap Prod.swap) := by
  ext σ τ
  simp only [bondOp, of_apply, submatrix_apply, Prod.swap_prod_mk]
  exact if_congr ⟨fun h i hk hj ↦ h i hj hk, fun h i hj hk ↦ h i hk hj⟩ rfl rfl

/-- A site permutation intertwines the moved bond with the original one. -/
theorem intertwines_bondOp_relabel [CommSemiring R] (π : Equiv.Perm (Fin L)) (j k : Fin L)
    (B : Matrix (d × d) (d × d) R) :
    Intertwines (relabel π) (fun _ ↦ 1) (bondOp (π j) (π k) B) (bondOp j k B) := fun σ τ ↦ by
  rw [bondOp_relabel, mul_one, one_mul]

/-- A sign `p (σ s)` at a site `s` off the bond commutes with the bond. -/
theorem intertwines_bondOp_gauge_of_ne [CommSemiring R] (p : d → R) {s j k : Fin L}
    (hj : s ≠ j) (hk : s ≠ k) (B : Matrix (d × d) (d × d) R) :
    Intertwines id (fun σ ↦ p (σ s)) (bondOp j k B) (bondOp j k B) := fun σ τ ↦ by
  simp only [id, bondOp, of_apply]
  split_ifs with h
  · rw [h s hj hk, mul_comm]
  · simp

/-- A sign `p (σ k)` on the second site of the bond turns `B` into `B'` when
`B' (a, b) (c, e) p_e = p_b B (a, b) (c, e)`, that is `B' = (1 ⊗ P) B (1 ⊗ P)⁻¹` with
`P = diag p`. -/
theorem intertwines_bondOp_gauge_right [CommSemiring R] (p : d → R) {j k : Fin L}
    {B B' : Matrix (d × d) (d × d) R}
    (h : ∀ a b c e, B' (a, b) (c, e) * p e = p b * B (a, b) (c, e)) :
    Intertwines id (fun σ ↦ p (σ k)) (bondOp j k B) (bondOp j k B') := fun σ τ ↦ by
  simp only [id, bondOp, of_apply]
  split_ifs
  · exact h _ _ _ _
  · simp

/-- A sign `p (σ j)` on the first site of the bond turns `B` into `B'` when
`B' (a, b) (c, e) p_c = p_a B (a, b) (c, e)`. -/
theorem intertwines_bondOp_gauge_left [CommSemiring R] (p : d → R) {j k : Fin L}
    {B B' : Matrix (d × d) (d × d) R}
    (h : ∀ a b c e, B' (a, b) (c, e) * p c = p a * B (a, b) (c, e)) :
    Intertwines id (fun σ ↦ p (σ j)) (bondOp j k B) (bondOp j k B') := fun σ τ ↦ by
  simp only [id, bondOp, of_apply]
  split_ifs
  · exact h _ _ _ _
  · simp

end Bond

/-! ### Rings with a seam -/

section Seam

variable {d : Type} [DecidableEq d] {L : ℕ} {R : Type*}

/-- The ring `∑ⱼ Bⱼ,ⱼ₊₁` whose bond `m → m + 1` carries `Bc` in place of `B` (the seam). -/
def seamHamiltonian [AddCommMonoid R] (B Bc : Matrix (d × d) (d × d) R) (L : ℕ) (m : Fin L) :
    Matrix (Fin L → d) (Fin L → d) R :=
  ∑ j : Fin L, bondOp j (finRotate L j) (if j = m then Bc else B)

/-- The twisted ring is the ring with its seam on the closing bond `L - 1 → 0`. -/
theorem twistedRingHamiltonian_eq_seam [AddCommMonoid R] (B Bc : Matrix (d × d) (d × d) R)
    (m : Fin L) (hm : m.val + 1 = L) :
    twistedRingHamiltonian B Bc L = seamHamiltonian B Bc L m := by
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  have : (j.val + 1 = L) ↔ j = m := by rw [Fin.ext_iff]; omega
  simp only [this]

/-- **Rotation moves the seam by one bond**: the ring read through `σ ∘ finRotate` has its
seam one bond further on. -/
theorem intertwines_seam_rotate [CommSemiring R] (B Bc : Matrix (d × d) (d × d) R) (m : Fin L) :
    Intertwines (relabel (finRotate L)) (fun _ ↦ 1) (seamHamiltonian B Bc L (finRotate L m))
      (seamHamiltonian B Bc L m) := fun σ τ ↦ by
  simp only [seamHamiltonian, Matrix.sum_apply, bondOp_relabel, mul_one, one_mul]
  conv_rhs => rw [← Equiv.sum_comp (finRotate L)]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [(finRotate L).injective.eq_iff]

/-- `(j + 1)` reflected, plus one, is `j` reflected: the reflection reverses each bond. -/
theorem finRotate_rev_finRotate (j : Fin L) :
    finRotate L (Fin.rev (finRotate L j)) = Fin.rev j := by
  rcases L with _ | L
  · exact j.elim0
  apply Fin.ext
  rw [val_finRotate, Fin.val_rev, val_finRotate, Fin.val_rev]
  have := j.isLt
  split_ifs <;> omega

/-- **Reflection reverses the ring**: the ring read through `σ ∘ rev` has its seam on the bond
`rev (m + 1)`, provided both bond operators are symmetric under swapping their two sites. -/
theorem intertwines_seam_reverse [CommSemiring R] {B Bc : Matrix (d × d) (d × d) R}
    (hB : B.submatrix Prod.swap Prod.swap = B) (hBc : Bc.submatrix Prod.swap Prod.swap = Bc)
    (m : Fin L) :
    Intertwines (relabel Fin.revPerm) (fun _ ↦ 1)
      (seamHamiltonian B Bc L (Fin.rev (finRotate L m))) (seamHamiltonian B Bc L m) :=
  fun σ τ ↦ by
  simp only [seamHamiltonian, Matrix.sum_apply, bondOp_relabel, mul_one, one_mul,
    Fin.revPerm_apply]
  -- The bond `j` becomes the bond `k = rev (j + 1)`, traversed backwards.
  let e : Equiv.Perm (Fin L) := (finRotate L).trans Fin.revPerm
  conv_rhs => rw [← Equiv.sum_comp e]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  have he : e j = Fin.rev (finRotate L j) := rfl
  rw [he, finRotate_rev_finRotate, bondOp_swap (Fin.rev j)]
  have hj : (Fin.rev (finRotate L j) = Fin.rev (finRotate L m)) ↔ j = m := by
    rw [Fin.rev_inj, (finRotate L).injective.eq_iff]
  rw [if_congr hj rfl rfl]
  split_ifs <;> simp only [hB, hBc]

/-- **An on-site relabeling** that preserves both bond operators commutes with the ring. -/
theorem intertwines_seam_onsite [CommSemiring R] (f : Equiv.Perm d)
    {B Bc : Matrix (d × d) (d × d) R} (hB : B.submatrix (Prod.map f f) (Prod.map f f) = B)
    (hBc : Bc.submatrix (Prod.map f f) (Prod.map f f) = Bc) (m : Fin L) :
    Intertwines (onsite f) (fun _ ↦ 1) (seamHamiltonian B Bc L m) (seamHamiltonian B Bc L m) :=
  fun σ τ ↦ by
  simp only [seamHamiltonian, Matrix.sum_apply, bondOp_onsite, mul_one, one_mul]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  split_ifs <;> simp only [hB, hBc]

/-- **A sign gauge moves the seam across a site.** The sign `p (σ s)` turns the ring with its
seam on the bond `s → s + 1` into the ring with its seam on the bond `s - 1 → s`, provided
`Bc (a, b) (c, e) p_e = p_b B (a, b) (c, e)` and `B (a, b) (c, e) p_c = p_a Bc (a, b) (c, e)`
(for the `π` twist of spin one, `p = (-1)^m`). -/
theorem intertwines_seam_gauge [CommSemiring R] (hL : 2 ≤ L) (p : d → R)
    {B Bc : Matrix (d × d) (d × d) R}
    (h₁ : ∀ a b c e, Bc (a, b) (c, e) * p e = p b * B (a, b) (c, e))
    (h₂ : ∀ a b c e, B (a, b) (c, e) * p c = p a * Bc (a, b) (c, e)) (s : Fin L) :
    Intertwines id (fun σ ↦ p (σ s)) (seamHamiltonian B Bc L s)
      (seamHamiltonian B Bc L ((finRotate L).symm s)) := by
  refine Intertwines.sum _ fun j _ ↦ ?_
  have hprev : (j = (finRotate L).symm s) ↔ finRotate L j = s := Equiv.eq_symm_apply _
  rw [if_congr hprev rfl rfl]
  by_cases hjs : j = s
  · subst hjs
    have hne : finRotate L j ≠ j := (ne_finRotate hL j).symm
    simp only [hne, ite_false, ite_true]
    exact intertwines_bondOp_gauge_left p h₂
  · by_cases hks : finRotate L j = s
    · subst hks
      simp only [hjs, ite_false, ite_true]
      exact intertwines_bondOp_gauge_right p h₁
    · simp only [hjs, hks, ite_false]
      exact intertwines_bondOp_gauge_of_ne p (Ne.symm hjs) (Ne.symm hks) B

/-- **The signed rotation of a twisted ring.** Rotating the twisted ring by one site and then
applying the sign `p (σ 0)` maps it to itself. This is the paper's signed shift `G_k`
(`thermal72.tex`, §"Representative columns"), one step at a time. -/
theorem intertwines_twisted_rotate [CommSemiring R] (hL : 2 ≤ L) (p : d → R)
    (hp : ∀ a, p a * p a = 1) {B Bc : Matrix (d × d) (d × d) R}
    (h₁ : ∀ a b c e, Bc (a, b) (c, e) * p e = p b * B (a, b) (c, e))
    (h₂ : ∀ a b c e, B (a, b) (c, e) * p c = p a * Bc (a, b) (c, e)) :
    Intertwines (relabel (finRotate L))
      (fun σ ↦ 1 * p (σ ⟨0, by omega⟩)) (twistedRingHamiltonian B Bc L)
      (twistedRingHamiltonian B Bc L) := by
  set last : Fin L := ⟨L - 1, by omega⟩
  set zero : Fin L := ⟨0, by omega⟩
  have hlast : finRotate L last = zero := by
    apply Fin.ext; rw [val_finRotate]; simp [last, zero]; omega
  have hzero : (finRotate L).symm zero = last := by rw [Equiv.symm_apply_eq, hlast]
  rw [twistedRingHamiltonian_eq_seam B Bc last (by simp [last]; omega)]
  -- Gauge the seam from `L - 1 → 0` to `0 → 1`, then rotate it back.
  have hg := (intertwines_seam_gauge hL p h₁ h₂ zero).symm_of_id (fun σ ↦ hp _)
  rw [hzero] at hg
  have hr := intertwines_seam_rotate B Bc last
  rw [hlast] at hr
  exact hg.trans hr

/-- The untwisted ring (`Bc = B`) is invariant under the rotation, with no sign. -/
theorem intertwines_ring_rotate [CommSemiring R] (hL : 0 < L) (B : Matrix (d × d) (d × d) R) :
    Intertwines (relabel (finRotate L)) (fun _ ↦ 1) (twistedRingHamiltonian B B L)
      (twistedRingHamiltonian B B L) := by
  have hs : ∀ m m' : Fin L, seamHamiltonian B B L m = seamHamiltonian B B L m' := fun _ _ ↦ by
    simp only [seamHamiltonian, ite_self]
  have h := intertwines_seam_rotate B B (⟨L - 1, by omega⟩ : Fin L)
  rwa [hs (finRotate L _) ⟨L - 1, by omega⟩,
    ← twistedRingHamiltonian_eq_seam B B _ (by simp; omega)] at h

/-- The twisted ring is invariant under the reflection `j ↦ L - 1 - j`, which maps the closing
bond to itself, when both bond operators are symmetric under swapping their sites. -/
theorem intertwines_twisted_reverse [CommSemiring R] (hL : 0 < L)
    {B Bc : Matrix (d × d) (d × d) R} (hB : B.submatrix Prod.swap Prod.swap = B)
    (hBc : Bc.submatrix Prod.swap Prod.swap = Bc) :
    Intertwines (relabel Fin.revPerm) (fun _ ↦ 1) (twistedRingHamiltonian B Bc L)
      (twistedRingHamiltonian B Bc L) := by
  set last : Fin L := ⟨L - 1, by omega⟩
  have h0 : finRotate L last = ⟨0, hL⟩ := by
    apply Fin.ext; rw [val_finRotate]; simp [last]; omega
  have hlast : Fin.rev (finRotate L last) = last := by
    rw [h0]; apply Fin.ext; simp [last]
  have h := intertwines_seam_reverse hB hBc last
  rwa [hlast, ← twistedRingHamiltonian_eq_seam B Bc last (by simp [last]; omega)] at h

/-- An on-site relabeling that preserves both bond operators is a symmetry of the twisted
ring. -/
theorem intertwines_twisted_onsite [CommSemiring R] (f : Equiv.Perm d)
    {B Bc : Matrix (d × d) (d × d) R} (hB : B.submatrix (Prod.map f f) (Prod.map f f) = B)
    (hBc : Bc.submatrix (Prod.map f f) (Prod.map f f) = Bc) :
    Intertwines (onsite f) (fun _ ↦ 1) (twistedRingHamiltonian B Bc L)
      (twistedRingHamiltonian B Bc L) := fun σ τ ↦ by
  simp only [twistedRingHamiltonian, Matrix.sum_apply, bondOp_onsite, mul_one, one_mul]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  split_ifs <;> simp only [hB, hBc]

end Seam

end CertifiedGaps
