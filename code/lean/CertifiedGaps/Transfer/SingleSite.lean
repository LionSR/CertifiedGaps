/-
Copyright (c) 2026 thermal-mpo-gap contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: thermal-mpo-gap contributors
-/
import CertifiedGaps.Transfer.MatrixLemmas

/-!
# A single-site time step and its transfer matrix, at every ring length

The checkerboard factorization of `CertifiedGaps.Ring.Transfer` needs a ring of even length.
Spectral filters also use odd lengths, so here the imaginary-time step is built from one label
per bond and a site factor that does not depend on the order of the two bonds meeting at a site.

The data are a finite site space `d`, real bond factors `C : ι → Matrix d d ℝ` (the bond is
`h = -∑_α C_α ⊗ C_α`), a real on-site term `o : Matrix d d ℝ`, and a twist of the closing bond
`L - 1 → 0`: a sign vector `r : ι → ℝ` and a relabeling `π : ι → ι` of its right end, so that the
closing bond is `-∑_α r_α C_α ⊗ C_{π α}` (`r = 1`, `π = id` is the plain ring). A sign twist
(`π = id`) is a ring twisted by a symmetry that is diagonal on the bond factors; a relabeling
twist is one that permutes them, such as the cyclic permutation of the spin axes. No symmetry of
`C` or `o` is needed in this file.

A bond label `a : Option ι` carries `C_a`, with `C_none = 1`, and a site whose two bonds carry
`a` and `b` gets the symmetrized factor `M(a, b) = (C_a C_b + C_b C_a) / 2`. The label `e j` sits
on the bond `j - 1 → j`, so the site `j` carries `M(e j, e (j + 1))`, except that the site `0`
sees the closing-bond label relabeled, `M(π (e 0), e 1)`. The time step with weight `u` per
nontrivial label is the bond part times the on-site part,

  `V(u) = (∑_e (∏ⱼ w_u(e j) sⱼ(e j)) ⨂ⱼ M(πⱼ (e j), e (j + 1))) · ⨂ⱼ (1 - u o)`,

where `w_u(none) = 1`, `w_u(some α) = u`, and the sign `sⱼ` is `r α` on the closing bond (label
`e 0`) and `1` elsewhere, and `πⱼ` is `π` for `j = 0` and the identity elsewhere.

Expanding `V(τ)^m` and exchanging the space and time products gives a transfer matrix on
time-label sequences `γ : Fin m → Option ι`,

  `X(γ, γ') = h(γ) h(γ') Tr ∏ₜ M(γ t, γ' t) (1 - τ o)`,  `h(γ) = ∏ₜ √τ^[γ t ≠ none]`,

which is real and symmetric (because `M` is symmetric in its two labels) and does not depend on
the ring length.

## Main results

* `transfer_isHermitian`: `X` is real symmetric.
* `trace_stepOp_pow`: `Tr V(τ)^m = ∑_γ D_r(γ) (X^(n+1)) (π γ, γ)` on a ring of `n + 1` sites,
  where `D_r(γ) = ∏ₜ r(γ t)` and `π γ` relabels every time label of `γ`.
* `trace_stepOp_pow_id`: for a sign twist, `Tr V(τ)^m = Tr (diag(D_r) X^(n+1))`.
-/

open Matrix

namespace CertifiedGaps.SingleSite

variable {d ι : Type} [Fintype d] [DecidableEq d]

/-! ### Labels and site factors -/

section Factors

variable (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

/-- The bond factor carried by a label, with `C_none = 1`. -/
def labelGen : Option ι → Matrix d d ℝ
  | none => 1
  | some α => C α

/-- The symmetrized product `M(a, b) = (C_a C_b + C_b C_a) / 2` of two labels. -/
noncomputable def symFactor (a b : Option ι) : Matrix d d ℝ :=
  (1 / 2 : ℝ) • (labelGen C a * labelGen C b + labelGen C b * labelGen C a)

/-- The site factor `M(a, b) (1 - τ o)` of one time step. -/
noncomputable def siteFactor (τ : ℝ) (a b : Option ι) : Matrix d d ℝ :=
  symFactor C a b * (1 - τ • o)

/-- The complex site factor `M(a, b) (1 - u o)`. -/
noncomputable def siteFactorC (u : ℂ) (a b : Option ι) : Matrix d d ℂ :=
  (symFactor C a b).map (↑) * (1 - u • o.map (↑))

/-- The symmetrized product is symmetric in its two labels; this makes the transfer matrix
symmetric. -/
theorem symFactor_comm (a b : Option ι) : symFactor C a b = symFactor C b a := by
  simp only [symFactor, add_comm]

/-- Two trivial labels give the identity. -/
@[simp] theorem symFactor_none_none : symFactor C none none = 1 := by
  simp only [symFactor, labelGen, mul_one, ← two_smul ℝ (1 : Matrix d d ℝ), smul_smul]
  norm_num

/-- A single nontrivial label gives its bond factor. -/
@[simp] theorem symFactor_none_some (α : ι) : symFactor C none (some α) = C α := by
  simp only [symFactor, labelGen, mul_one, one_mul, ← two_smul ℝ (C α), smul_smul]
  norm_num

/-- A single nontrivial label gives its bond factor. -/
@[simp] theorem symFactor_some_none (α : ι) : symFactor C (some α) none = C α := by
  rw [symFactor_comm, symFactor_none_some]

/-- At a real time step the complex site factor is the complexified real one. -/
theorem siteFactorC_ofReal (τ : ℝ) (a b : Option ι) :
    siteFactorC C o τ a b = (siteFactor C o τ a b).map (↑) := by
  have h : ∀ M : Matrix d d ℝ, M.map Complex.ofReal = Complex.ofRealHom.mapMatrix M := fun _ ↦ rfl
  rw [siteFactorC, siteFactor, h (_ * _), map_mul, map_sub, map_one, ← h, ← h]
  congr 2
  ext i j
  simp

/-- The trace of an ordered product of complex site factors at a real time step is real. -/
theorem trace_prod_siteFactorC {m : ℕ} (τ : ℝ) (a b : Fin m → Option ι) :
    trace (List.ofFn fun t ↦ siteFactorC C o τ (a t) (b t)).prod =
      ((trace (List.ofFn fun t ↦ siteFactor C o τ (a t) (b t)).prod : ℝ) : ℂ) := by
  simp only [siteFactorC_ofReal, trace_list_prod_ofFn_map_ofReal]

end Factors

/-! ### The time step -/

/-- The complex weight `u` of a nontrivial label (`1` for `none`). -/
def labelWeight (u : ℂ) (a : Option ι) : ℂ := a.elim 1 fun _ ↦ u

/-- The sign `r α` of a nontrivial label (`1` for `none`). -/
def labelSign (r : ι → ℝ) (a : Option ι) : ℝ := a.elim 1 r

/-- The sign of the label `a` on the bond `j - 1 → j`: `r` on the closing bond `L - 1 → 0`
(`j = 0`), trivial elsewhere. -/
def bondSign (r : ι → ℝ) {L : ℕ} (j : Fin L) (a : Option ι) : ℝ :=
  if (j : ℕ) = 0 then labelSign r a else 1

/-- The coefficient `∏ⱼ w_u(e j) sⱼ(e j)` of a label configuration. -/
def stepCoeff {L : ℕ} (r : ι → ℝ) (u : ℂ) (e : Fin L → Option ι) : ℂ :=
  ∏ j, labelWeight u (e j) * (bondSign r j (e j) : ℂ)

/-- The relabeling `πⱼ` of the bond `j - 1 → j` at its right end, the site `j`: `π` on the
closing bond (`j = 0`), the identity elsewhere. -/
def labelAt (π : ι → ι) {L : ℕ} (j : Fin L) (α : ι) : ι := if (j : ℕ) = 0 then π α else α

/-- The label of the bond `j - 1 → j` as the site `j` sees it, relabeled by `πⱼ`. -/
def bondLabel (π : ι → ι) {L : ℕ} (j : Fin L) (a : Option ι) : Option ι := a.map (labelAt π j)

omit [Fintype d] [DecidableEq d] in
/-- The closing bond is relabeled by `π`. -/
theorem labelAt_zero (π : ι → ι) {n : ℕ} : labelAt π (0 : Fin (n + 1)) = π :=
  funext fun _ ↦ by simp [labelAt]

omit [Fintype d] [DecidableEq d] in
/-- The other bonds are not relabeled. -/
theorem labelAt_of_ne (π : ι → ι) {L : ℕ} {j : Fin L} (hj : (j : ℕ) ≠ 0) : labelAt π j = id :=
  funext fun _ ↦ by simp [labelAt, hj]

omit [Fintype d] [DecidableEq d] in
/-- Without a relabeling every site sees the bond label itself. -/
@[simp] theorem bondLabel_id {L : ℕ} (j : Fin L) (a : Option ι) : bondLabel id j a = a := by
  have h : labelAt (id : ι → ι) j = id := funext fun _ ↦ ite_self _
  simp [bondLabel, h]

omit [Fintype d] [DecidableEq d] in
/-- Read along a time-label sequence on each bond of a ring of `n + 1` sites, the relabeled
labels are the configuration with its closing-bond sequence (`q 0`) relabeled by `π`. -/
theorem bondLabel_seq (π : ι → ι) {n m : ℕ} (q : Fin (n + 1) → Fin m → Option ι)
    (k : Fin (n + 1)) :
    (fun t ↦ bondLabel π k (q k t)) = Function.update q 0 (fun t ↦ (q 0 t).map π) k := by
  by_cases hk : k = 0
  · subst hk
    simp [bondLabel, labelAt_zero]
  · simp [bondLabel, labelAt_of_ne π (fun h ↦ hk (Fin.ext h)), Function.update_of_ne hk]

variable [Fintype ι] [DecidableEq ι] (C : ι → Matrix d d ℝ) (o : Matrix d d ℝ)

/-- The bond part `∑_e (∏ⱼ w_u(e j) sⱼ(e j)) ⨂ⱼ M(πⱼ (e j), e (j + 1))` of the time step on a ring
of `L` sites, with the closing bond signed by `r` and relabeled by `π` at its right end. -/
noncomputable def bondStep (L : ℕ) (r : ι → ℝ) (π : ι → ι) (u : ℂ) :
    Matrix (Fin L → d) (Fin L → d) ℂ :=
  ∑ e : Fin L → Option ι,
    stepCoeff r u e • piTensor fun j ↦
      (symFactor C (bondLabel π j (e j)) (e (finRotate L j))).map (↑)

/-- The on-site part `⨂ⱼ (1 - u o)` of the time step. -/
noncomputable def siteStep (L : ℕ) (u : ℂ) : Matrix (Fin L → d) (Fin L → d) ℂ :=
  piTensor fun _ : Fin L ↦ 1 - u • o.map (↑)

/-- The time step `V(u)`: the bond part times the on-site part. -/
noncomputable def stepOp (L : ℕ) (r : ι → ℝ) (π : ι → ι) (u : ℂ) :
    Matrix (Fin L → d) (Fin L → d) ℂ :=
  bondStep C L r π u * siteStep o L u

/-- The tensor product `⨂ⱼ M(πⱼ (e j), e (j + 1)) (1 - u o)` of a label configuration. -/
noncomputable def stepTerm {L : ℕ} (π : ι → ι) (u : ℂ) (e : Fin L → Option ι) :
    Matrix (Fin L → d) (Fin L → d) ℂ :=
  piTensor fun j ↦ siteFactorC C o u (bondLabel π j (e j)) (e (finRotate L j))

omit [DecidableEq ι] in
/-- The time step as one sum over label configurations. -/
theorem stepOp_eq_sum (L : ℕ) (r : ι → ℝ) (π : ι → ι) (u : ℂ) :
    stepOp C o L r π u = ∑ e : Fin L → Option ι, stepCoeff r u e • stepTerm C o π u e := by
  simp only [stepOp, bondStep, siteStep, Finset.sum_mul, smul_mul_assoc, piTensor_mul,
    stepTerm, siteFactorC]

/-! ### The transfer matrix -/

/-- The half weight `h(γ) = ∏ₜ √τ^[γ t ≠ none]` of a time-label sequence. -/
noncomputable def halfWeight (τ : ℝ) {m : ℕ} (γ : Fin m → Option ι) : ℝ :=
  ∏ t, (γ t).elim 1 fun _ ↦ √τ

/-- The site trace `Tr ∏ₜ M(γ t, γ' t) (1 - τ o)` of two time-label sequences. -/
noncomputable def siteTrace (τ : ℝ) {m : ℕ} (γ γ' : Fin m → Option ι) : ℝ :=
  trace (List.ofFn fun t ↦ siteFactor C o τ (γ t) (γ' t)).prod

/-- The transfer matrix `X(γ, γ') = h(γ) h(γ') Tr ∏ₜ M(γ t, γ' t) (1 - τ o)` on time-label
sequences. -/
noncomputable def transfer (m : ℕ) (τ : ℝ) :
    Matrix (Fin m → Option ι) (Fin m → Option ι) ℝ :=
  of fun γ γ' ↦ halfWeight τ γ * halfWeight τ γ' * siteTrace C o τ γ γ'

/-- The character `D_r(γ) = ∏ₜ r(γ t)` of a time-label sequence. -/
def labelChar (r : ι → ℝ) {m : ℕ} (γ : Fin m → Option ι) : ℝ := ∏ t, labelSign r (γ t)

omit [Fintype ι] [DecidableEq ι] in
/-- The site trace is symmetric in its two sequences. -/
theorem siteTrace_comm (τ : ℝ) {m : ℕ} (γ γ' : Fin m → Option ι) :
    siteTrace C o τ γ γ' = siteTrace C o τ γ' γ := by
  simp only [siteTrace, siteFactor, symFactor_comm C (γ _)]

omit [Fintype ι] [DecidableEq ι] in
/-- The transfer matrix is symmetric. -/
theorem transfer_symm (m : ℕ) (τ : ℝ) (γ γ' : Fin m → Option ι) :
    transfer C o m τ γ γ' = transfer C o m τ γ' γ := by
  simp only [transfer, of_apply, siteTrace_comm C o τ γ, mul_comm (halfWeight τ γ)]

omit [Fintype ι] [DecidableEq ι] in
/-- The transfer matrix is real symmetric. -/
theorem transfer_isHermitian (m : ℕ) (τ : ℝ) : (transfer C o m τ).IsHermitian := by
  ext γ γ'
  simp only [conjTranspose_apply, star_trivial, transfer_symm C o m τ γ']

/-! ### The trace of the time step -/

/-- The real weight `τ` of a nontrivial label. -/
def labelWeightR (τ : ℝ) (a : Option ι) : ℝ := a.elim 1 fun _ ↦ τ

omit [Fintype d] [DecidableEq d] [Fintype ι] [DecidableEq ι] in
/-- At real `u = τ` the coefficient is real. -/
theorem stepCoeff_ofReal {L : ℕ} (r : ι → ℝ) (τ : ℝ) (e : Fin L → Option ι) :
    stepCoeff r (τ : ℂ) e = ((∏ j, labelWeightR τ (e j) * bondSign r j (e j) : ℝ) : ℂ) := by
  rw [stepCoeff, Complex.ofReal_prod]
  refine Finset.prod_congr rfl fun j _ ↦ ?_
  rcases e j with _ | α <;> simp [labelWeight, labelWeightR]

omit [Fintype d] [DecidableEq d] [Fintype ι] [DecidableEq ι] in
/-- The weight of a label is the square of its half weight. -/
theorem labelWeightR_eq_sq {τ : ℝ} (hτ : 0 ≤ τ) (a : Option ι) :
    labelWeightR τ a = (a.elim 1 fun _ ↦ √τ) ^ 2 := by
  rcases a with _ | α
  · simp [labelWeightR]
  · simp [labelWeightR, Real.sq_sqrt hτ]

omit [Fintype d] [DecidableEq d] [Fintype ι] [DecidableEq ι] in
/-- Relabeling the time labels of a sequence preserves its half weight. -/
theorem halfWeight_map (τ : ℝ) (π : ι → ι) {m : ℕ} (γ : Fin m → Option ι) :
    halfWeight τ (fun t ↦ (γ t).map π) = halfWeight τ γ :=
  Finset.prod_congr rfl fun t _ ↦ by
    change ((γ t).map π).elim _ _ = (γ t).elim _ _
    rcases γ t with _ | α <;> rfl

omit [Fintype ι] [DecidableEq ι] in
/-- The per-path identity behind `trace_stepOp_pow`: the time-ordered weights of a space-time
label configuration `q j t` equal the transfer-matrix weights of the closed path `q`, whose
first step leaves from the relabeled closing-bond sequence. -/
theorem stepPath_eq {n m : ℕ} (r : ι → ℝ) (π : ι → ι) {τ : ℝ} (hτ : 0 ≤ τ)
    (q : Fin (n + 1) → Fin m → Option ι) :
    (∏ t, ∏ j, labelWeightR τ (q j t) * bondSign r j (q j t)) *
        ∏ j, siteTrace C o τ (fun t ↦ bondLabel π j (q j t)) (q (finRotate (n + 1) j)) =
      labelChar r (q 0) * ∏ k, transfer C o m τ (fun t ↦ bondLabel π k (q k t))
        (q (finRotate (n + 1) k)) := by
  have hsign : ∏ j : Fin (n + 1), ∏ t, bondSign r j (q j t) = labelChar r (q 0) := by
    rw [Fin.prod_univ_succ]
    simp [bondSign, labelChar]
  have hw : ∀ j, ∏ t, labelWeightR τ (q j t) = halfWeight τ (q j) ^ 2 := fun j ↦ by
    simp only [labelWeightR_eq_sq hτ, halfWeight, Finset.prod_pow]
  have hrot : ∏ k, halfWeight τ (q (finRotate (n + 1) k)) = ∏ k, halfWeight τ (q k) :=
    Equiv.prod_comp (finRotate (n + 1)) fun k ↦ halfWeight τ (q k)
  have hh : ∀ k, halfWeight τ (fun t ↦ bondLabel π k (q k t)) = halfWeight τ (q k) :=
    fun k ↦ halfWeight_map τ _ (q k)
  rw [Finset.prod_comm, ← hsign]
  simp only [Finset.prod_mul_distrib, hw, transfer, of_apply, hh, hrot, Finset.prod_pow]
  ring

/-- **Transfer form of the time step.** On a ring of `n + 1` sites,
`Tr V(τ)^m = ∑_γ D_r(γ) (X^(n+1)) (π γ, γ)` with the transfer matrix `X` of `m` time steps,
where `π γ` relabels every time label of `γ`: the relabeling twist is the permutation matrix
`Π (γ, γ') = [γ' = π γ]` in `Tr (diag(D_r) Π X^(n+1))`. -/
theorem trace_stepOp_pow (n m : ℕ) (r : ι → ℝ) (π : ι → ι) {τ : ℝ} (hτ : 0 ≤ τ) :
    trace (stepOp C o (n + 1) r π τ ^ m) =
      ((∑ γ, labelChar r γ * (transfer C o m τ ^ (n + 1)) (fun t ↦ (γ t).map π) γ : ℝ) : ℂ) := by
  rw [← List.prod_replicate, ← List.ofFn_const, stepOp_eq_sum, list_prod_ofFn_sum, trace_sum,
    sum_mul_pow_apply_eq_sum_path, Complex.ofReal_sum]
  refine Fintype.sum_equiv (Equiv.piComm fun (_ : Fin m) (_ : Fin (n + 1)) ↦ Option ι) _ _
    fun E ↦ ?_
  rw [list_prod_ofFn_smul, trace_smul]
  simp only [stepTerm, list_prod_ofFn_piTensor, trace_piTensor, trace_prod_siteFactorC,
    stepCoeff_ofReal, ← Complex.ofReal_prod, smul_eq_mul, ← Complex.ofReal_mul]
  congr 1
  rw [← funext (bondLabel_seq π (Equiv.piComm _ E))]
  exact stepPath_eq C o r π hτ (Equiv.piComm _ E)

/-- **Transfer form of a sign-twisted time step.** On a ring of `n + 1` sites,
`Tr V(τ)^m = Tr (diag(D_r) X^(n+1))` when the closing bond carries only the signs `r`. -/
theorem trace_stepOp_pow_id (n m : ℕ) (r : ι → ℝ) {τ : ℝ} (hτ : 0 ≤ τ) :
    trace (stepOp C o (n + 1) r id τ ^ m) =
      ((trace (diagonal (labelChar r) * transfer C o m τ ^ (n + 1)) : ℝ) : ℂ) := by
  rw [trace_stepOp_pow C o n m r id hτ]
  congr 1
  simp only [Option.map_id, id_eq, trace, diag_apply, diagonal_mul]

end CertifiedGaps.SingleSite
