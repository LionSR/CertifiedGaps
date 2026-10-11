This adds a verified two-purity certificate for the spin-1/2 two-leg Heisenberg ladder at J_parallel=1, J_perp=4. The outward-rounded lower bound is **1.102403170532509 for every even L >= 4 one-rung cells**: N=2L >= 8 spin sites, with N divisible by four.

The seed is (n0,beta)=(4,2), with rational upper defects `17541603577//500000000000` and `3078024217//125000000000`. Exact rational RP identities, complete dense symmetry spectra, residual/orthogonality error bounds, outward partition-function enclosures, and the recurrence are implemented only under `certificates/heisenberg_ladder/`.

### Validation

On Julia 1.12.4 against current main (`22b2759`, `CertifiedGaps`):

- `julia --project=julia certificates/heisenberg_ladder/test/runtests.jl`: passed; 24.6 s including fresh verified 4/8-rung eigensolves, source hashes, sector completeness, stored-seed replay, exact RP, and bound recomputation.
- `julia --project=julia certificates/heisenberg_ladder/smoke.jl`: passed; independent full Hamiltonians and spectra through 5 rungs, plus exact second-moment normalization checks; 14.3 s.
- `julia --project=julia certificates/heisenberg_ladder/audit.jl`: passed; 12.8 s. Replays all 256 and 65,536 eigenvalues without reading saved radii; independently checks exact trace moments and the three log partition functions using directed 256-bit MPFR arithmetic.
- The `CertifiedGaps` migration was checked against all core source files: they match the old source after namespace replacement. No core file is changed by this PR. The original validated `ThermalMPOGap` snapshot and matching source hashes are preserved byte-for-byte under `legacy/`.

### Numerical scope and work still running

The valid even-seed grid is complete for J_perp in {4,3.5,3,2.5,2,1.5}, n0 in {2,4}, beta=0.25:0.25:16. Numerical bounds are 1.10240444, 0.77864899, 0.45418980, 0.04243459 for the first four values; no even seed was found at 2 or 1.5. The observed grid threshold is bracketed by 2.4881973267 and 2.4882354736. These results are **numerical, not certified**.

Full 10-rung dense scans are complete for J=4,3,2.5,2,1.5. Only the checkpointed retry for J=3.5 remains active. The earlier J=3.5 10-rung attempt was OOM-killed; no partial spectrum from that process is represented as complete. This PR stays draft and will be updated with the remaining results. The completed 10-rung J=4 reference gap is numerical: 3.1375441605.

Odd n0=3,5 are diagnostics only: the theorem requires even n0, and odd original ladder rings are frustrated. A 10-rung ring therefore does not provide a valid new n0=5 certificate.

### Audit limitations and reproduction

The fresh audit shares the symmetry builder and residual-bound implementation; it is not an independently implemented eigensolver or a formal verification of the full pipeline. The repository's abstract two-purity theorem and ordinary IEEE/OpenBLAS arithmetic model remain assumptions. No rigorous thermodynamic reference gap, coupling threshold, monotonicity result, or absence of seeds beyond the searched grid is claimed. Fresh replay is tested with Julia 1.12.4 and the checked-in manifest; a different BLAS implementation may require regenerated enclosures rather than bit-identical replay.

Set `OPENBLAS_NUM_THREADS=1`. Regenerate the certificate with `julia --project=julia certificates/heisenberg_ladder/certify.jl`, then run the checker and audit above. Numerical scans use `scan_lowmem.jl '4//1'` (and the other rational couplings). `REPORT.md` gives the theorem, exact decomposition, normalization, validity range, provenance, and full reproduction details. Numerical JSONs retain complete sector dimensions, low levels and logZ grids; the rigorous complete eigenvalue enclosures are retained in full.

No CI, merge, or shared-core changes are included.
