# Transfer moments and centered Weyl transport

All statements here are conditional on the repository's RP spatial-moment theorem and
certified finite spectra. The new calculation never constructs a ring larger than 12.
Exact rational parameters are enclosed outward; the final gap uses the unmodified
`CertifiedGaps.certified_gap_bound`.

## Spatial squaring and the missing thermal ring

For even n, σₙ = Xᵦⁿ / Zₙ(β) is positive with trace one. If its eigenvalues are pᵢ,
P = Σpᵢ² and P′ = Σpᵢ⁴/P². Then

    1−P′ = 2 Σᵢ<ⱼ pᵢ²pⱼ²/P² ≤ (1−P)²/(2P²).

Thus a defect upper bound s gives f(s) = s²/[2(1−s)²] after spatial squaring.
The code caps this by the trivial upper bound one. The bound is monotone for 0 ≤ s < 1.

Four finite inputs, log Z₆(β), log Z₁₂(β), log Z₆(2β), log Z₁₂(2β), supply bounds on
s = 1−Pₛ(6,β), s₂ = 1−Pₛ(6,2β), t = 1−Pτ(12,β). The exact identity

    Pτ(4n,β) = Pₛ(2n,2β) Pτ(2n,β)² / Pₛ(2n,β)²

and Pₛ ≤ 1 imply the update

    s′ = min(1,f(s)),  s₂′ = min(1,f(s₂)),
    t′ = 1−(1−s₂′)(1−t)².

These are upper bounds for the seed at 2n, at the **same** β. Iteration gives seeds
n₀ = 12,24,48,… from the L = 6/12 inputs. It is distinct from the subsequent core
space-and-temperature doubling recurrence. Every acceptance uses outward-rounded interval
endpoints converted exactly to rationals; a floating-point search result cannot accept a seed.

At D = 4 the n₀ = 12, β = 29/10 seed gives gap ≥ 0.899143399258073 for even L ≥ 12.
Verified positive gaps at L = 6,8,10 all exceed that value, hence the stronger uniform
conclusion holds for every even L ≥ 6. The small finite bridge does not change the seed.
At D = 3.206 none of the five larger seeds closes on the tested β grid. Repeated length
squaring reduces the spatial defect but eventually loses too much thermal purity.
This is a failure of these estimates on this grid, not a no-go theorem for larger seeds.

## Parameter transport

Write H_L(D) = H_L(D*) + δV, where V = Σ(Sᶻ)² and 0 ≤ V ≤ L·I. Therefore

    ‖δ(V−L·I/2)‖ ≤ L|δ|/2.

Weyl's inequality bounds every eigenvalue of H_L(D)−δL·I/2 by the corresponding
H_L(D*) eigenvalue plus/minus L|δ|/2. Consequently log Z for this shifted Hamiltonian
lies within βL|δ|/2 of its reference enclosure. For an interval [Dlo,Dhi], use
max(|Dlo−D*|,|Dhi−D*|). The scalar energy shift cancels exactly in **both** purity ratios.
No extensive perturbation estimate is applied directly to the infinite-volume gap.
The transported finite seed is passed to the bootstrap separately at each D in the interval.

At D* = 3.206 this certifies [3.2059998,3.206], β = 381/125, n₀ = 6, with bound
0.001001626234429964. The next 10⁻⁷-grid point fails this particular transport search.
Its physical seed might still be admissible; the method makes no exclusion claim.

## Fast partition enclosures from a complete spectrum

For centers Eᵢ and verified radii rᵢ, sorted by centers, the first K = 8192 levels are
summed with interval exponentials (all levels if there are fewer). For the remaining N−K,

    Eᵢ(true) ≥ E[K+1] − maxᵢ>K rᵢ,

so their shifted Boltzmann sum lies in [0, (N−K) exp(−β(E[K+1]−rmax−shift))]. This is a
rigorous tail from **complete** spectral enclosures, unlike an assumed cutoff after a
Krylov computation. Each sector eigenvalue retains its multiplicity. The omitted positive
tail is never discarded in an accepting inequality.

For a finite gap the global maximum radius is used for both ordered eigenvalues. Different
blocks can cross when their centers are sorted; using only the radii attached to the first
two centers would require a separate ordering proof. Parameter transport adds L|δ| to the
gap radius. Dividing a certified constant B by the positive Δ₁₂ interval gives the plotted
ratio interval. This denominator is finite-volume and does not assert a thermodynamic gap.

## Reproducibility and checker limits

`certificate.json` and `results/finite_bridge.json` fingerprint proof sources and the pinned
project; the former also fingerprints the historical complete spectral files. The checker
reconstructs those files byte for byte, verifies the original RP/symmetry checks, recomputes
the four partition enclosures, transport, spatial lifts and core recurrence, and rebuilds
both bridge spectra. It rejects changed defect bounds, lengths, parameter intervals and
claims. The L = 12 eigensystems and residuals are the archived verified inputs, not recomputed
by the fast checker. Their enclosure proof and independent replay are in METHOD.md and audit/.
