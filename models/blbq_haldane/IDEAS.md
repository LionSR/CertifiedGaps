# Routes considered

The requested deliverable is a finite-size scaling survey, not an assumed positive certificate.
One cell is one spin-1 site. The theorem requires even n0; odd n0 are diagnostics only.

* Free/commuting closed forms: unavailable for the generic interacting BLBQ chain.
* Rigorous cluster expansions: no applicable convergent remainder supplied near the Haldane
  regime; developing one is outside a finite-size survey.
* Full dense ED with magnetization alone: largest blocks grow as 3^L/sqrt(L), infeasible
  at L=14 on this machine.
* Translation, magnetization and spin reversal: translation reduces blocks by about L;
  conjugate momenta and opposite magnetizations share spectra. This gives an independently
  checkable complete spectrum at moderate sizes. Dense diagonalization remains cubic per
  block. Preflight counts and measured timings determine the largest affordable size.
* SU(2) highest weights would reduce blocks further, but constructing an orthonormal,
  translation-resolved highest-weight basis is an additional algorithm; simple subtraction
  of magnetization-sector spectra does not reduce the cost of their dense diagonalization.
* Truncated Lanczos spectra can miss degeneracies and cannot certify Z. They will not be
  used for partition functions in this survey.
* Extrapolation: fit the even-n ground energy mismatch and finite-size triplet gap. A dilute
  triplet partition ansatz can estimate scales, but these are numerical model-dependent
  forecasts, especially near theta=-pi/4.

Chosen route: exact symbolic RP coefficient identities; complete dense magnetization and
momentum ED with multiplicities and small-system crosschecks; numerical purity grid and
explicitly conditional extrapolation. Only if an eligible seed appears, implement and run
verified full eigenvalue enclosures. A resource-limited survey is reported as such.

## Follow-up: constructive SU(2) route

In magnetization M=S and momentum k, D=S+ has full row rank and kernel dimension
r=d(M=S,k)-d(M=S+1,k). The spectrum of D'*D is exactly
J(J+1)-S(S+1), J=S,...,L. Thus the polynomial product over J>S of
(I-D'*D/[J(J+1)-S(S+1)]) is the highest-weight projector. Project deterministic Gaussian
columns and QR-orthonormalize them; then densely diagonalize Q'*H*Q with multiplet
weight 2S+1. Sparse D multiplication avoids dense Casimir diagonalization or a large
sparse-QR factorization. This is a numerical reduction, not a rigorous enclosure.

The implementation orders polynomial roots to limit intermediate amplification, checks
raising residuals and orthogonality, and first compares complete weighted spectra with
unreduced small rings. Benchmarks determine whether projected bases and dense eigensolves
fit the 16 GiB memory and 45-minute per-run budgets before attempting L=14.

## Outcome

No certificate. The exact RP identity holds for −π/4 < θ ≤ 0, but on complete spectra through
L = 14 no even seed n₀ ≤ 6 comes near the admissible region at any of the six angles (best
min(Pₛ, P_τ) = 0.329 at θ = 0, n₀ = 6, against about 0.873 needed). The SU(2) route above made
L = 14 complete; L = 16 (n₀ = 8) exceeds the 16 GiB machine. Conditional extrapolations put the
first admissible seed at n₀ ≈ 20–48, so a certificate needs a non-ED route (README, "Result").
