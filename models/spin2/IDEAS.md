# Spin-2 Heisenberg ring: routes to a gap certificate

Model: `H_L = Σ_j S_j·S_{j+1}` on an even ring, `S = 2`, `S² = 6`. Numerical reference values
(Todo–Kato, not certified): bulk gap `Δ ≈ 0.08917`, correlation length `ξ ≈ 49.49`, energy per
bond `e_∞ ≈ −4.761249`, hence a magnon velocity `v = Δξ ≈ 4.41`. The spin-1 chain has
`Δ ≈ 0.4105`, `ξ ≈ 6.0`: spin 2 is 4.6 times softer and 8.2 times longer-ranged.

This memo records the candidate routes to the two seed inequalities of the doubling criterion,
why each does or does not work for this model, how its cost scales with the seed `(n₀, β₀)`, and
what was done. The repository already contains a long certified-MPO programme for this model
(see `README.md`); the memo is the missing scale analysis that explains its outcome and ranks
the alternatives. Nothing here certifies a spin-2 gap.

## 1. What the criterion demands for spin 2

Write `Z_L(β) = e^{−βE₀(L)} (1 + R_L(β))`, `R_L(β) = Σ_{m≥1} e^{−β(E_m(L) − E₀(L))}`, and
`δ_n = E₀(2n) − 2E₀(n) > 0` (small rings have the lower energy per site). Then exactly

    Pₛ(n, β) = e^{−βδ_n} (1 + R_{2n}(β)) / (1 + R_n(β))²,
    P_τ(L, β) = (1 + R_L(2β)) / (1 + R_L(β))².

So the two seeds are controlled by two numbers of very different character:

* the **ring-closure energy** `δ_n ∼ e^{−n/ξ}`, which forces `β δ_n ≲ 0.2` (spatial defect);
* the **excitation sum** `R_{2n}(β) ∼ 3·(number of magnon modes within 1/β of the gap)·e^{−βΔ}`,
  which forces `βΔ ≳ 6–8` (thermal defect).

The certified D = 295 SU(2) trial energies in `results/trial-295-symmetric-certificate.json`
(upper bounds on `E₀(N)`, but their *differences* are the best available estimate of `δ_n`) give

| n | δ_n from trial energies | β allowed by β δ_n ≤ 0.2 | β needed for 1 − P_τ(2n, β) ≤ 0.1 |
|---|---|---|---|
| 128 | 0.01064 | ≤ 19 | ≳ 48 |
| 256 | 4.69·10⁻⁴ | ≤ 430 | ≳ 54 |
| 512 | 5.0·10⁻⁷ | ≤ 4·10⁵ | ≳ 61 |

(the last column is the magnon-gas estimate of `experiments/seed_scale.jl`). Pushing the model
defects through the exact recurrence, the window is empty for every `n ≤ 160` and opens at
`n = 176–192` (depending on how `δ_n` is extrapolated beyond the trial lengths), with
`β ∈ [48, 70–110]`. The float map in `experiments/seed-scale.json` therefore puts the smallest
admissible seed near **(n₀, β₀) ≈ (192–256, 50–70)**. The gap bound that such a seed would
certify is `0.008–0.03` (10–35 % of the bulk gap); the best point on the grid is
`(256, 64–68)` with a bound near `0.03`.

Equivalently (the one-trace proposition already in `scripts/one_trace_target.jl`): a single
certified **upper** bound `Ẑ_256(64) ≤ 1.08` on the shifted trace with the certified shift
`a = 4.7612463` would give `Δ_L > 0.0030` for all even `L ≥ 512`. Everything below is about how
to obtain an upper bound of that kind. Upper bounds on `Z` are the hard direction: every
variational principle (Peierls–Bogoliubov, Gibbs, Golden–Thompson with a trial state) bounds
`Z` from below. An upper bound needs control of the *whole* spectrum of a 200-site ring.

Unequal defects do not change the scale. The exact recurrence (`certified_gap_bound`, outward
rounded) admits at most `εₛ ≈ 0.2176` when `εₜ → 0` and at most `εₜ ≈ 0.2006` when `εₛ → 0`
(table in `experiments/seed-scale.json`, exact rationals); relaxing `βδ_n ≤ 0.13` to `≤ 0.24`
moves `n₀` by about `ξ ln 1.8 ≈ 25` sites. The alternative seed pair of survey Theorem 2.24(v),
`P_τ(n₀, β₀)` and `Pₛ(n₀, 2β₀)`, needs the same upper bound on `Z_{n₀}` at the doubled `β`.

## 2. Why the spin-1 recipe does not transfer

The spin-1 certificate (`models/spin1`) never computed a 2304-site ring. It used exact thermal
traces of rings of **4–12 sites** at `β₀ = 49/4`, read as moments `Z_k(β₀) = Σ_i λ_i(β₀)^k` of
the self-adjoint spatial transfer operator `X_β`, and a variational trial at length 72 or 120.
Polynomial filters bound the top eigenvalue from above, the trial bounds it from below, and the
capped-mass inequality `Σ_{i≥1} λ_i^{jk} ≤ (Σ_{i≥1} λ_i^k)^j` pushes the residual mass to the
seed lengths. The whole argument rests on the residual `k`-th moment being below one:

    Σ_{i≥1} (λ_i/λ₀)^k = Ẑ_k(β₀)/λ̂₀^k − 1 ≈ exp[β₀ k (e_∞ − e₀(k))] (1 + R_k(β₀)) − 1 < 1.

For spin 1 this is `12.25 · 12 · 4.4·10⁻³ ≈ 0.65` (`e_∞ − e₀(12) ≈ 4.4·10⁻³`): it works, barely.
For spin 2 at `β₀ = 64` and the largest ring that exact diagonalisation can reach, `k = 16`
(`e₀(16) = −4.7735` by DMRG, `e_∞ − e₀(16) ≈ 0.012`), the exponent is `64 · 16 · 0.012 ≈ 12.5`.
The residual mass is `e^{12.5}`, not `0.65`. Using `k(e_∞ − e₀(k)) ≈ δ_k/2` and the trial
`δ_n`, the residual mass drops below one only for `k ≳ 100–128`. Temperature lifting
(`Ẑ_k(rβ₀) ≤ Ẑ_k(β₀)^r`) does not rescue smaller `k`: it is lossless only when the `k`-ring is
already in its ground state at `β₀`, which is the same `β₀ ≳ 5/Δ_k` requirement, and the
ring-closure exponent scales with `r` anyway. So the moment route needs the same object as every
other route: a certified thermal trace of a ring of at least a hundred spin-2 sites at
`β ≈ 50–64`. Exact diagonalisation (`5^100`) is out by ninety orders of magnitude.

## 3. Routes considered

Each entry: idea, what it would need for spin 2, verdict.

**A. Closed forms.** No free-fermion, commuting or frustration-free point. Ruled out.

**B. High-temperature cluster expansion.** Converges for `β ≲ 1/(2e·‖h‖) ≈ 0.03` with
`‖h‖ = 6`; we need `β ≈ 60`. Out by three orders of magnitude; the Trotter/polynomial seeds of
the MPO programme already play this role at `β ≈ 10⁻⁶`. Ruled out.

**C. Perturbation from a solvable point.** The spin-2 AKLT model (projector `P₃ + P₄` on the
bond) is frustration free with a short `ξ`, but `h_Heis + 6 = P₁ + 3P₂ + 6P₃ + 10P₄`: the
difference is not small in any norm, and no interpolation keeps a certified gap. Large-`D`,
dimerised or ladder limits leave the Heisenberg point. Ruled out for *this* model (it is the
right route for anisotropic spin-2 models with `ξ ≲ 5`).

**D. Transfer-operator gap by positivity.** In the Cartesian basis `X_β` has signed kernel
entries; in the Marshall basis the Trotter gates are entrywise nonnegative, but the column
transfer matrix of a Trotterised torus has dimension `25^m` with `m ≳ 10⁵` slices needed for a
3 % Trotter loss at `β = 64`, and its Birkhoff projective diameter is astronomically large.
Ruled out.

**E. Small-ring moments and filters (the spin-1 route).** Needs rings of `k ≳ 100` sites
(section 2). Ruled out at exact-diagonalisation sizes; at `k ≈ 100–128` it is just another form of
the certified-trace problem H.

**F. Variational sandwich.** `Z_{2n} ≥ e^{−βE_trial(2n)}` is certified (trial energies at
128–1024 sites). The other side needs a **lower** bound on `E₀(n)` within about `10⁻³` in total
at `n ≈ 256` (`4·10⁻⁶` per site) plus a bound on `R_n(β)`. Anderson-type and local
sum-of-squares bounds have a boundary error `O(1/w)` per site for a window of `w` sites
(`≈ 0.6/w`), so `w ≈ 10⁵`. Ruled out.

**G. Verified eigenvalue enclosures (Temple, Kato, Lehmann–Maehly).** All give lower bounds
on `E₀` only in terms of a lower bound on the energy *above* the trial subspace, i.e. a
gap: circular for the ring whose gap is the goal. Ruled out as a seed source (useful later
for finite rings once a gap is known).

**H. Certified thermal MPO (the repository's main programme).** Exact polynomial seed at
`β/2^K`, `K` certified squarings of an `N = 256` MPO, signed backward trace correction.
Status in `README.md`: two squarings completed at `β ≈ 10⁻⁶` with HS error `1.7·10⁻¹⁰` at bond
dimension 16–32 (`K = 28`, degree-2 seed). The normalised relative HS error at least doubles
per squaring (the `2b/κ ≥ 2` factor is intrinsic to norm-only error control), so the residual
of step `j` is charged `2^{K−j}` times at the end. Two consequences:

* The seed degree matters more than anything else: a degree-12 polynomial seed needs only
  `K = 16` squarings for the same `β = 64` (README table), so the total amplification is
  `2^16 ≈ 6.6·10⁴` instead of `2^28 ≈ 2.7·10⁸`. `experiments/seed-scale.json` lists the
  step-resolved budget: with a final slack of 0.02, the residual allowed at step `j` is
  `0.02/(16·2^{16−j})`, from `4·10⁻⁸` at `β = 0.002` to `1.3·10⁻³` at `β = 64`.
* The expensive steps are the last ones (`β ≥ 1`), where the operator entanglement of
  `e^{−βH/2}` is largest and then saturates at the square of the ground-state Schmidt
  spectrum; there the budget is loose (`10⁻⁵`–`10⁻³`). What decides feasibility is therefore the
  bond dimension at which the central-cut operator Schmidt spectrum of `e^{−βH/2}` has
  discarded weight `10⁻¹⁰`–`10⁻⁶` for `β = 1…64`. `experiments/operator_entanglement.jl`
  (open chain, numerical) and `experiments/ground_schmidt.jl` (DMRG) measure exactly this;
  see section 6.

The mathematics is complete; whether the computation is feasible is a question of bond
dimension, answered in section 6. **Not ruled out.** The one methodological change that would
relax the budget further: an error functional weighted by the thermal state (trace norm against
`e^{−βH}`) instead of HS or operator norms, so that compression errors in high-energy
components are not charged at full weight.

**I. SU(2) sector resolution.** The leading spatial defect is the triplet magnon,
`(λ₁/λ₀)^n = e^{−n/ξ}`. The singlet sector of `X_β` has the two-magnon gap, `ξ_eff ≈ ξ/2`, so
a criterion phrased for the twist-averaged torus `∫dg Z_n(β; g)` would halve `n₀`. The
thermal purity of that ensemble is the purity of the direct integral `⊕_g e^{−βH^g}`, which is
positive and trace class, so the squaring lemma survives; but the gap conclusion then
concerns a Haar-average over twists and needs a Lipschitz bound in the twist to return to
the untwisted ring, costing `O(‖h‖θ)` in the gap. Speculative, at best a factor two in `n₀`,
and the certified-trace bottleneck remains. Not pursued.

**J. Quantum Monte Carlo.** The even ring is bipartite, so the purities at `(256, 64)` can be
sampled without a sign problem (replica/swap estimators). A statistical error bar is not a
certificate, and the criterion rejects it. Worth running only to confirm the float scale map.

**K. Alternative seed pair (survey Theorem 2.24(v)).** Same bottleneck; see section 1.

## 4. Choice and outcome

No route avoids a certified upper bound on the partition function of a ring of at least about
two hundred spin-2 sites at `β ≳ 50`, accurate to a few percent. Route H is the only one with
complete mathematics and an implementation; its remaining cost is set by bond dimension, not
by ideas. This session therefore did not attempt a certificate. It delivered

1. this memo;
2. `experiments/seed_scale.jl`: the float seed-scale map (magnon-gas excitation sums with
   `Δ = 0.08917`, `v = 4.41`; `δ_n` from the certified trial energies, two extrapolations in
   `n`), the exact admissible region of the recurrence for unequal defects, and the error
   budget of the squaring route; output `experiments/seed-scale.json`, every float clearly
   labelled as an estimate;
3. `test/seed_scale.jl`: fast checks of the model functions and of the admissibility table;
4. `experiments/ground_schmidt.jl` and `experiments/operator_entanglement.jl`: the numerical
   bond-dimension measurements of section 6, which turn route H's cost into numbers.

**Negative result.** At exact-diagonalisation sizes (`n ≤ 8`, `2n ≤ 16`) no admissible seed
exists for the spin-2 Heisenberg ring: the spatial window `β ≤ 0.22/δ_n` and the thermal window
`β ≥ 6/Δ_{2n}` are disjoint (`β ≤ 0.15` against `β ≥ 11` at `n = 4`, `β ≤ 0.33` against
`β ≥ 19` at `n = 8`), and the model windows stay empty for every `n ≤ 160`. The estimated seed
scale is `(n₀, β₀) ≈ (192–256, 50–70)`.

## 5. What would change the verdict

* A certified `Ẑ_n(β) ≤ U` for any `(n, β, U)` in the feasible region of `seed-scale.json`,
  with the matching trial at `2hn` (available for `2hn ≤ 1024`).
* A compression error certificate measured against the thermal weight rather than in HS or
  operator norm (route H), bringing the required `D` back to the MPS scale `D₀ ≈ 300`.
* Not a longer trial: `h = 4` in the one-trace scheme (a trial at `N = 2048`) would allow
  `β = 32` on the 256-ring, but there the excitation factor `1 + R_256(32) ≈ 1.28` already
  exceeds what that scheme tolerates. The thermal input must be taken at `β ≳ 50`.

## 6. Measured bond dimensions for the squaring route (numerical, October 12)

Two uncertified ITensors runs on open chains (`experiments/ground_schmidt.jl`,
`experiments/operator_entanglement.jl`; outputs `ground-schmidt-n64.json` and
`operator-entanglement-n32-d200.json`).

**Large β.** For `β ≫ 1/Δ` the thermal amplitude is the ground-state projector and its
operator Schmidt values at a cut are the products `λ_a λ_b` of the state's Schmidt values.
DMRG on 64 open sites (8 sweeps to D = 400, energy not fully converged) gives centre-bond
discarded weights `7.4·10⁻⁸` at D = 100, `5.5·10⁻⁹` at D = 200, `4.7·10⁻¹⁰` at D = 300 and
entanglement entropy 1.22. Keeping `D` state values on each side gives an operator rank `D²` with
discarded weight `≈ 2·tail(D)`. Against the step budget of route H (`K = 16`, slack 0.02):

| step | β | allowed relative HS residual | state D | operator D² |
|---|---|---|---|---|
| 16 | 64 | 1.3·10⁻³ | 64 | 4 100 |
| 15 | 32 | 6.3·10⁻⁴ | 81 | 6 600 |
| 14 | 16 | 3.1·10⁻⁴ | 112 | 12 500 |
| 13 | 8 | 1.6·10⁻⁴ | 162 | 26 000 |
| 12 | 4 | 7.8·10⁻⁵ | 227 | 52 000 |
| 11 | 2 | 3.9·10⁻⁵ | 283 | 80 000 |
| 10 | 1 | 2.0·10⁻⁵ | 330 | 109 000 |

The product form overestimates the operator rank at `β ≲ 8`, where the amplitude is still a
thermal mixture; the thermal run below measures the true spectrum there. These ranks are for a
dense bond; U(1) charge blocks (already implemented) and SU(2) multiplets reduce the arithmetic
by one to two orders of magnitude.

**Thermal run.** Open chain of 32 sites, `e^{−βH/2}` by second-order Trotter, `maxdim = 200`,
central-cut operator Schmidt spectrum at `β = 2^k`; discarded weight `w(D) = Σ_{i>D} s_i²`, so the
relative HS residual of a rank-`D` compression is `√w(D)`:

| β | χ reached | w(32) | w(64) | w(128) | budget (√w) |
|---|---|---|---|---|---|
| 1/8 | 42 | 8·10⁻¹¹ | 0 | 0 | 2.4·10⁻⁶ |
| 1/4 | 67 | 3.7·10⁻⁸ | 4.9·10⁻¹³ | 0 | 4.9·10⁻⁶ |
| 1/2 | 155 | 5.8·10⁻⁶ | 1.4·10⁻⁸ | 5.2·10⁻¹¹ | 9.8·10⁻⁶ |
| 1 | 200 (cap) | 2.2·10⁻⁴ | 3.8·10⁻⁶ | 1.5·10⁻⁷ | 2.0·10⁻⁵ |
| 2 | 200 (cap) | 2.4·10⁻³ | 1.7·10⁻⁴ | 1.3·10⁻⁵ | 3.9·10⁻⁵ |
| 4 | 200 (cap) | 1.1·10⁻² | 1.8·10⁻³ | 1.5·10⁻⁴ | 7.8·10⁻⁵ |

Rows at the cap understate the tails. Up to `β = 1/2` the budget is met by `D ≈ 128`
(`√w(128) ≈ 7·10⁻⁶` against `10⁻⁵`). From `β = 1` on it is not: the budget weight at `β = 1` is
`4·10⁻¹⁰` and at `β = 2` it is `1.5·10⁻⁹`, while `w(128)` is `10⁻⁷` and `10⁻⁵`. Each doubling of
`D` lowers `w` by a factor of about 13 in the measured range, so meeting the budget at
`β = 1–2` needs `D` of a few thousand even on this optimistic extrapolation, and the
product-form estimate from the ground state gives `10⁴–10⁵` for `β = 2–8`. The operator
entanglement entropy keeps rising (0.62, 1.22, 1.78, 2.21, 2.55 at `β = 1/4, 1/2, 1, 2, 4`); on
the open chain the free spin-1 edge modes add to the ground-state limit of about 2.4 (twice the
state entropy), so the ring value will be somewhat lower.

**Conclusion for route H.** The certified-squaring route needs bond dimensions of
`2·10³–10⁵` on the six most expensive squarings (`β = 1` to `64`), against the `16–34` at
which the exact-witness machinery has run at `N = 256` (about 160 s per squaring at
`D = 16`). With the contraction cost growing like `D³`–`D⁴`, that is a factor `10⁵`–`10⁹` in
time: two to four orders of magnitude beyond a feasible campaign, before symmetry blocks (one
to two orders back). The route is not wrong, it is expensive. What would make it feasible is a
compression certificate that charges errors against the thermal weight rather than in HS norm,
which would let the late steps run at the state bond dimension `D ≈ 100–300` instead of its
square. That is the one open problem this memo identifies as worth attacking.
