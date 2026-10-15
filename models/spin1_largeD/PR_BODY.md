The existing large-D certificate stopped at an ED seed n₀ = 6. Transfer-operator squaring
now certifies n₀ = 12 using only the stored L = 6/12 spectra, without diagonalizing L = 24.
Fresh verified L = 8/10 spectra bridge the stronger D = 4 bound to every even L ≥ 6.
Centered Weyl transport also certifies D ∈ [3.2059998, 3.206].

All changes are confined to `models/spin1_largeD/`. The route memo was committed before
computation. The model uses its own pinned Julia project depending on `julia/`; no core
addition was needed. The checker verifies exact RP, partition enclosures, transfer and
parameter bounds, the core recurrence, source hashes, and the finite-ring bridge.

Validation: 3,381 tests, including altered-certificate rejection. Repeated the quarter-step
D grid and refined β near D = 3.205/3.206; retained every scan row and a reproducible ratio
curve. L = 14 block counting found a 28.0 GB matrix requirement on a 16 GiB host. The
original survey seed (4,6,3) now has the rigorous lower bound 0.543812926; the quoted
floating-point 0.543842 is not assumed.

Certified summary (bounds rounded down; lengths are sites):

| D | Bound and validity | Seed (n₀, β₀, εₛ, εₜ), exact |
|---|---|---|
| 4 | 0.899143399, every even L ≥ 6 | (12, 29/10, 8376656523711495/2305843009213693952, 283096258190373/9007199254740992), plus finite L = 6/8/10 bridge |
| [3.2059998, 3.206] | 0.001001626, every even L ≥ 6 | (6, 381/125, 863115144109137/4503599627370496, 418622220052377/9007199254740992) |

Runtime: about 9 s for the new scalar search, 7 s for the bridge, and 40 s for the checker;
the reused verified L = 12 inputs originally took 35–36 min each. Every run was below
45 min. Floating point remains in the screening scan, unsuccessful-grid diagnostics,
plot coordinates and timing. Certified markers also enclose the finite Δ₁₂ denominator;
it is not identified with the thermodynamic gap. D = 3.2059998 is the smallest certified
value in this search, not a proven global threshold. Spatial lifts at D = 3.206 fail on
the tested β grid; that limitation is documented without asserting gaplessness.
