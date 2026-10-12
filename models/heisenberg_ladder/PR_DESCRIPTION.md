This adds a verified two-purity certificate for the spin-1/2 two-leg Heisenberg ladder at J_parallel=1, J_perp=4. The outward-rounded lower bound is **1.10240309174412 for every even L >= 4 one-rung cells**: N=2L >= 8 spin sites, with N divisible by four.

The seed is (n0,beta)=(4,2), with rational upper defects `3508321//100000000` and `123121//5000000`. Exact rational RP identities, complete dense symmetry spectra, residual/orthogonality error bounds, outward partition-function enclosures, and the recurrence are implemented only under `certificates/heisenberg_ladder/`. Seeds are rounded upward on a 10^-8 grid to retain slack across fresh residual calculations; the historical tighter snapshot is preserved unchanged.

### Completed numerical search

All six couplings {4,3.5,3,2.5,2,1.5} are complete through 10 rungs (20 spins). Every 10-rung spectrum accounts for 1,048,576 states; the largest dense block has dimension 8,398. All individual dense runs stayed below approximately 45 minutes. No scan or supervisor remains running.

For valid even n0 in {2,4}, beta=0.25:0.25:16, numerical bounds are 1.10240444, 0.77864899, 0.45418980, 0.04243459 at the first four couplings. No even seed was found at 2 or 1.5. The observed finite-grid transition is bracketed by 2.4881973266601562 and 2.4882354736328125. These are **numerical, not certified** threshold results. The J=4 10-rung reference gap is 3.1375441605; the certificate is 35.136% of this finite-ring reference.

The report records failed-grid purity maxima, delta(n), and explicitly tentative extrapolations (n0,beta)=(6,5.25) at J=2 and (8,8.5) at J=1.5. Odd n0=3,5 are diagnostics only: the theorem requires even n0, and odd original ladder rings are frustrated.

One J=3.5 attempt was OOM-killed. The lower-memory retry checkpointed every sector. Its final dimension assertion caught a threaded aggregation race, which was corrected; all saved sectors were independently checked and reassembled without repeating eigensolves or accepting a partial spectrum.

### Validation

On Julia 1.12.4 against main `22b2759` (`CertifiedGaps`), with OPENBLAS_NUM_THREADS=1:

- `julia --project=julia certificates/heisenberg_ladder/test/runtests.jl`: passed in 22.70 s, including fresh verified 4/8-rung eigensolves, hashes, complete sector counts, exact RP, strict rational seed inequalities and bound recomputation.
- `julia --project=julia certificates/heisenberg_ladder/smoke.jl`: passed in 12.71 s; independently constructed full Hamiltonians and spectra through 5 rungs and exact second moments.
- `julia --project=julia certificates/heisenberg_ladder/audit.jl`: passed in 12.47 s; fresh complete spectra without saved radii, exact trace moments, and independently implemented directed 256-bit MPFR partition sums.
- `julia --project=julia certificates/heisenberg_ladder/check_scans.jl`: passed; all six compacted numerical grids, dimensions and best-seed calculations replay correctly.
- All 25 core sources match the pre-rename source after namespace replacement. The original validated ThermalMPOGap snapshot and source hashes are preserved byte-for-byte under `legacy/`; current proof outputs were regenerated with CertifiedGaps.

### Reproduction and limits

Regenerate with `julia --project=julia certificates/heisenberg_ladder/certify.jl`, then run the checker and audit above. Numerical scans use `scan_lowmem.jl '4//1'` and the other rational couplings. REPORT.md gives the exact RP decomposition, uniqueness argument, enclosure theorem, normalization, timings and provenance. Numerical JSONs retain complete sector counts, low levels and logZ grids; verified spectra retain all eigenvalue enclosures.

The fresh audit shares the symmetry builder and residual-bound implementation; it is not an independent eigensolver or formal verification of the entire pipeline. The repository's abstract theorem and ordinary IEEE/OpenBLAS arithmetic model remain assumptions. No rigorous thermodynamic reference gap, coupling threshold, monotonicity, or absence of seeds beyond the grid is claimed. A different arithmetic backend may require regenerated enclosures.

No shared-core changes or CI additions. This PR remains **draft**; no merge is requested.
