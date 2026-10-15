# Route memo: spin-1 chain with single-ion anisotropy D

Model H = Σ Sⱼ·Sⱼ₊₁ + D Σ (Sᶻⱼ)², one-site cells, Cartesian RP basis. This memo was written in
two rounds, each before its computations: first for the ED baseline, then for the follow-up
that goes beyond ED. The outcomes are recorded at the end.

## Candidate routes

* **Closed forms.** The exchange terms commute neither with the single-ion term nor with
  adjacent bonds. No free-particle or frustration-free partition function is known for this
  D range, and the D = ∞ product limit alone does not give finite-D seeds.
* **Cluster / low-temperature expansion.** The product vacuum is a natural reference, but pair
  creation by exchange needs a polymer estimate uniform in volume and β. An extensive
  remainder cannot control log Z₂ₙ − 2 log Zₙ, and a convergent, numerically useful remainder
  near D ≈ 3 would be substantial new analysis. Not pursued; cost after such a derivation could
  be independent of n up to exponentially small winding polymers.
* **Symmetry-resolved dense ED (baseline).** Sᶻ and translation reduce the largest L = 12 block
  from 73789 to about 6000; spin inversion pairs ±Sᶻ and conjugation pairs ±k. Cost is cubic
  in block dimension and exponential in L; L = 14 blocks (~42000) are infeasible densely on this
  host, and n₀ = 7 would be odd anyway. Certification needs every level of every block, with
  verified residual/orthogonality bounds and a polar/Weyl enclosure; Krylov low spectra are
  only for screening.
* **Transfer moments (chosen beyond-ED route).** For even n, X_βⁿ/Zₙ is a positive density, and
  squaring improves its purity defect s to at most s²/[2(1 − s)²]. This bounds Pₛ at
  n = 12, 24, … from L = 6, 12 alone. The thermal input follows from the exact identity
  P_τ(2n, b) = Pₛ(n, 2b) P_τ(n, b)² / Pₛ(n, b)² and Pₛ ≤ 1, tracking spatial defects at both b
  and 2b. Cost: O(log(n/6)) interval operations after four stored log Z enclosures. Risk: the
  thermal defect grows with every lift. A Birkhoff contraction estimate is less direct because
  no strictly positive finite kernel with controlled Trotter error is available.
* **Perturbation in D (second chosen route).** V = Σ(Sᶻ)² satisfies 0 ≤ V ≤ L, so
  ‖V − L/2‖ ≤ L/2 and Weyl moves each verified finite spectrum to nearby D with radius
  L|δD|/2 and a cancelling energy shift. Its defect loss is O(nβ|δD|): good for certifying
  parameter intervals without another ED run, not a uniform-volume perturbation theorem.
* **Smaller seeds (survey Thm 2.24(v)).** Replace P_τ(2n₀) by P_τ(n₀) and Pₛ(n₀, 2β₀); the
  coupling identity above is the same mechanism, used inside the transfer route.
* **Variational sandwiches.** A trial MPS bounds E₀(2n) from above, but no matching lower bound
  on E₀(n) with exponentially small surface error is available; an Anderson-type bound with
  O(1) error cannot resolve δₙ near threshold. Left as a separate route.

## Choice

ED baseline: exact RP, a floating-point screen over D, β and n₀ ≤ 6, then complete verified
dense spectra at L = 6, 12 for the best reachable seeds. Follow-up: transfer squaring and
centred Weyl transport on those stored spectra, a verified small-ring bridge, and a refined
certified D interval, with every failure recorded rather than hidden.

## Outcome

* ED seeds n₀ = 6 certify D = 4 (0.694694192) and D = 3.206 (0.002213456); the screen finds
  no even seed with n₀ ≤ 6 at D ≤ 3.205.
* Transfer squaring succeeds at D = 4: n₀ = 12, β = 29/10, bound 0.899143399, from L ≤ 12 only;
  verified L = 8, 10 gaps extend it to every even L ≥ 6. Further lifts (n₀ = 24, 48, 96) still
  close but with smaller bounds, because thermal purity is lost at each step.
* At D = 3.206 no lift through n₀ = 192 closes on β = 1.5:0.025:6 (plus 381/125): thermal-purity
  loss prevents it. This negative applies to these estimates only.
* Centred Weyl transport certifies D ∈ [3.2059998, 3.206] with n₀ = 6, β = 381/125.
* Pushing below D ≈ 3.2 needs seeds of roughly n₀ ≈ 10–26 at β ≈ 4–30 (a crude extrapolation
  of δ(n) from n = 4, 6), i.e. a route that bounds spatial purity at large n without losing
  thermal purity, such as a polymer expansion or a verified MPS/cluster sandwich.
