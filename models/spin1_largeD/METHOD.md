# Method: model, enclosures, transfer squaring and Weyl transport

Everything below is conditional on the repository's two-purity (spacetime bootstrap) theorem.
Stage 1 (`certify.jl`, `seed.jl`) produces complete verified spectra at L = 6, 12 and the ED
seed n₀ = 6. Stage 2 (`src/TransferCertificate.jl`, `transfer_certify.jl`) reuses those
spectra to certify n₀ = 12, …, 96 and nearby D intervals without building a larger ring.

## 1. Exact model and reflection positivity

In the Cartesian basis take the integer matrices

```
Gx = [0 0 0; 0 0 1; 0 -1 0]
Gy = [0 0 -1; 0 0 0; 1 0 0]
Gz = [0 1 0; -1 0 0; 0 0 0]
```

and set Sₐ = −iGₐ. `rp_identity` (in `model.jl`) verifies in integer arithmetic that
[Sx, Sy] = iSz, Σₐ Sₐ² = 2I, that the independently built bond h = P − |Ω⟩⟨Ω| (P the swap,
Ω = |xx⟩ + |yy⟩ + |zz⟩) equals Σₐ Sₐ ⊗ Sₐ, and that

    −h = Gx ⊗ Gx + Gy ⊗ Gy + Gz ⊗ Gz.

Every coefficient is 1 and every Gₐ is real, so Gₐ = conj(Gₐ). The on-site term
D(Sᶻ)² = D·diag(1,1,0) is real symmetric for every rational D. This is the local RP
condition; it does not by itself establish the seed hypothesis.

In the normalised Sᶻ basis m = −1, 0, 1, each allowed move (mⱼ, mⱼ₊₁) → (mⱼ+1, mⱼ₊₁−1) has
amplitude 1 and the diagonal is Σ mⱼmⱼ₊₁ + D Σ mⱼ². `test/reviewer.jl` checks this basis
change exactly, including the L = 2 convention of two periodic bonds.

## 2. Symmetry accounting

Product states are split by total Sᶻ = q; q and −q have equal spectra (m → −m). In q = 0
the blocks use translation together with spin inversion (character e^{2πikt/L} pᶠ, p = ±1);
for q ≠ 0, translation alone. An orbit contributes only if the character is trivial on its
stabiliser. With normalised orbit vectors, a transition from representative a to TᵘFᶠb has
matrix element

    h · sqrt(size(a)/size(b)) · exp(−2πiku/L) · pᶠ.

Complex conjugation pairs k with L − k, so each computed block has charge multiplicity 1 or 2
and momentum multiplicity 1 or 2. Every run asserts Σ multiplicity × dimension = 3^L; at
L = 12 this is 531441 states in 92 blocks, the largest of dimension 5798. `test/small.jl`
compares all eigenvalues with the full product-space Hamiltonian through L = 6, and
`test/symmetry.jl` checks the covariances and orbit counts through L = 12.

## 3. Verified eigenvalue enclosures (stage 1)

Every rational transition, phase and square root in a block is also accumulated in
outward-rounded interval arithmetic. The Float64 Hermitian matrix A is built from the upper
triangle; the entrywise deviations |Re| + |Im| summed along rows give c ≥ ‖H_block − A‖₂.
Physical norms give ‖H_block‖₂ ≤ L(2 + |D|) (bond eigenvalues −2, −1, 1), so
M = L(2 + |D|) + c bounds ‖A‖₂.

For the full dense eigensystem AQ ≈ QΛ of a block of dimension m, with u = 2⁻⁵³ and
γ = (8m + 32)u / (1 − (8m + 32)u), the standard binary64/BLAS dot-product model gives

    ‖AQ − QΛ‖₂ ≤ r = ‖fl(AQ − QΛ)‖_F + γ(√m M ‖Q‖_F + ‖Q‖_F ‖Λ‖₂) + 10⁻²⁸⁰,
    ‖Q*Q − I‖₂ ≤ e = ‖fl(Q*Q − I)‖_F + γ(‖Q‖_F² + √m) + 10⁻²⁸⁰.

A floating sum of N squared magnitudes is enlarged by γ_{4N+16} before an outward square
root; the tiny additive terms dominate subnormal effects at these sizes. The code asserts
e < 1 and runs with one BLAS thread. It does not trust LAPACK's eigenvalues.

**Polar/Weyl theorem.** If e < 1, Q is invertible, its polar factor U is unitary and
‖Q − U‖₂ ≤ d = e/(1 + √(1 − e)). Hence ‖AU − UΛ‖₂ ≤ r + (M + ‖Λ‖₂)d, which bounds
‖A − UΛU*‖₂. By Weyl's inequality (min-max), each sorted eigenvalue of H_block lies within

    η = c + r + (M + ‖Λ‖₂) e/(1 + √(1 − e))

of the corresponding sorted diagonal entry of Λ, with multiplicity and in degenerate
eigenspaces. Each block's centres and common radius η are stored; the merged, sorted list is
the complete spectrum with every multiplicity.

## 4. Partition functions and the ED seed (stage 1)

For a complete spectrum, Σ exp(−β[λ − η, λ + η]) is summed with outward rounding after an
exact binary energy shift; log and unshift give an interval containing log Z. Then

    εₛ = upper(1 − exp(log Z_{2n}(β) − 2 log Z_n(β))),
    εₜ = upper(1 − exp(log Z_{2n}(2β) − 2 log Z_{2n}(β))),

whose binary upper endpoints are converted exactly to `Rational{BigInt}` and passed with the
rational β to the unmodified `certified_gap_bound`. Only even n₀ are eligible.

## 5. Spatial squaring and the missing thermal ring (stage 2)

For even n, σₙ = X_βⁿ / Zₙ(β) is positive with trace one. With eigenvalues pᵢ, P = Σpᵢ² and
P′ = Σpᵢ⁴/P² are the spatial purities at n and 2n, and

    1 − P′ = 2 Σ_{i<j} pᵢ²pⱼ² / P² ≤ (1 − P)² / (2P²),

because (1 − P)² = (2Σ_{i<j} pᵢpⱼ)² ≥ 4Σ_{i<j} pᵢ²pⱼ². A defect bound s therefore lifts to
f(s) = s²/[2(1 − s)²], capped at 1 (`square_defect`). The four finite inputs
log Z₆(β), log Z₁₂(β), log Z₆(2β), log Z₁₂(2β) bound s = 1 − Pₛ(6, β), s₂ = 1 − Pₛ(6, 2β) and
t = 1 − P_τ(12, β). The exact identity

    P_τ(4n, β) = Pₛ(2n, 2β) P_τ(2n, β)² / Pₛ(2n, β)²

with Pₛ ≤ 1 gives the update

    s′ = f(s),   s₂′ = f(s₂),   t′ = 1 − (1 − s₂′)(1 − t)²,

valid seeds at 2n and the **same** β. Iterating gives n₀ = 12, 24, 48, … from L ≤ 12
(`seed_ladder`). This is distinct from the core's own space-and-temperature recurrence, which
then runs unchanged. Every acceptance goes through interval endpoints converted exactly to
rationals; a floating-point search value never accepts a seed.

## 6. Parameter transport (stage 2)

Write H_L(D) = H_L(D*) + δV with V = Σ(Sᶻ)², 0 ≤ V ≤ L, so ‖δ(V − L/2)‖ ≤ L|δ|/2. By Weyl,
every eigenvalue of H_L(D) − δL/2 lies within L|δ|/2 of the corresponding one at D*, so log Z
of the shifted Hamiltonian lies within βL|δ|/2 of the reference enclosure; for an interval
[D_lo, D_hi] use max |D − D*|. The shift cancels exactly in both purities. No extensive
perturbation is applied to the infinite-volume gap: the transported finite seed is passed to
the bootstrap at each D of the interval separately (`transported_logz`).

## 7. Partition tails and finite gaps (stage 2)

Stage 2 sums the lowest K = 8192 sorted levels with interval exponentials. The remaining
N − K true levels are ≥ E_{K+1} − max_{i>K} rᵢ, so their shifted Boltzmann sum lies in
[0, (N − K) exp(−β(E_{K+1} − r_max − shift))]. This tail comes from a **complete** enclosure,
unlike a Krylov cutoff, and is never dropped in an accepting inequality.

A finite-ring gap uses the global maximum radius for both of the two lowest sorted centres,
since blocks may cross after sorting; transport adds L|δ|. Dividing a certified bound B by this
positive Δ₁₂ interval gives the plotted ratio interval. Δ₁₂ is a finite-volume quantity.

## 8. Even lengths and uniqueness

With x = L/n ∈ [1, 2], log-convexity of the even transfer moments and monotonicity of ℓᵖ norms
give log P_τ(L, β) ≥ log P_τ(2n, β) + 2 log Pₛ(n, β), so every even L in [n, 2n] has thermal
defect at most 1 − (1 − t)(1 − s)²; the dyadic bands from n₀ cover all even L ≥ n₀. Thermal
purity above 1/2 forces the largest Gibbs weight above 1/2, hence a unique ground state, and
the core bounds the gap from the ratio of the next weight. At D = 4 the n₀ = 12 conclusion is
extended to L = 6, 8, 10 by direct verified gaps (`bridge_certificate`), each of which exceeds
the uniform bound.
