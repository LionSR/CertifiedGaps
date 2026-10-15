# Candidate routes (non-certifying numerical investigation)

1. Dense diagonalization in each Sz sector through n=8 gives an independent
   multiplicity check and complete thermal sums. Cost is cubic in sector size.
2. Sparse restarted Lanczos/Arnoldi through n=12 gives low spectra at cost roughly
   nnz(H) times Krylov iterations, with explicit residual checks. Repeated
   eigenvalues and missing states are risks: validate against dense spectra and
   request guard eigenvalues above every counting threshold. No Float64 count is
   a certified spectral anchor.
3. Chebyshev block filtering from individual computational-basis states tests the
   proposed small recipe. Oversampling and QR are necessary to avoid losing low
   states. Measure total vector Hv calls, not just polynomial degree. Compare
   floor-rounded integer recurrences with Float64, separating vector generation
   from numerical QR and from the eventual exact Gram checker.
4. A sharper tail bound using sector thresholds may save work; retain the proposed
   dimension-times-exponential bound as the baseline. Quantify the difference
   between absolute shifted-trace errors and relative errors, particularly odd n.
5. Packed ternary basis has cheap fixed shifts but 3^n lanes; compressed Sz basis
   has fewer lanes but irregular permutations. Do not assume both advantages at
   once. Report both layouts and treat kernel runtime as conditional on an actual
   packing implementation.

Chosen route: dense/sparse numerical spectra, followed by measured Chebyshev
filtering and exact integer recurrence probes. This is an experiment and cost
audit, not a certificate. All new files remain in this task directory.
