Certifies the periodic spin-1/2 Heisenberg ladder using exact checkerboard reflection positivity, complete verified sector spectra, outward-rounded partition sums, and the shared certified recurrence. All changes stay in `models/heisenberg_ladder/`, with a standalone Julia project and pinned manifest.

| J⊥ | Certified gap ≥ | Seed (n₀, β₀), in rungs |
|---|---:|---|
| 4 | 1.102403091 | (4, 2) |
| 3.5 | 0.778648190 | (4, 9/4) |
| 3 | 0.454189334 | (4, 5/2) |
| 2.5 | 0.042433721 | (4, 3) |

All bounds hold for **N = 8, 12, 16, … spin sites**, with unique ground states. Exact defect bounds for every row are in the README and certificate.

The local transition on n₀∈{2,4}, β₀=k/4, k=1,…,64 lies between 326133/131072 ≈ 2.48819733 (all grid seeds rigorously rejected) and 163069/65536 ≈ 2.48823547 (Δ≥0.001531165). The latter seed is (4, 3, 4693199/25000000, 1317379/25000000). This is not a continuous-β optimum or physical critical coupling; monotonicity in J is not assumed. At J⊥=2 and 1.5, interval coverage and an infinite-β tail reject n₀=2,4 for **every β>0**.

The pre-computation memo weighs every AGENTS.md route. An analytic negative result proves that the centered independent-rung norm sandwich fails at arbitrary seed sizes for 0<J⊥≤4. Local perturbation does certify small coupling intervals around the stronger seeds. All six historical scans through 10 rungs/20 spins are replayed with freshly reconstructed exact sector multiplicities; the long ED runs are not repeated.

Validation: 57 fast tests, all 1024 grid decisions resolved rigorously, historical scan replay, and a separate fresh directed-MPFR/trace-moment audit (37.6 s). The full fresh checker passed in 308.34 s: RP, complete spectra, seed inequalities, recurrence, all-temperature rejection, coupling continuation, and source/environment hashes. No core changes.

**Certificate summary:** at J⊥=4, **Δ≥1.102403091 for N=8,12,16,… sites**, seed **(n₀,β₀,εₛ,εₜ)=(4,2,3508321/100000000,123121/5000000)**. This calculation took **46.1 s**; the full suite took **366.7 s**, excluding startup/install, with each run capped at 45 minutes. The retained ten-rung log-Z grids/finite-ring gaps, heuristic seed-scale extrapolations, and timings remain floating point. Fresh eigensolver outputs are covered by rigorous eigenvalue radii.
