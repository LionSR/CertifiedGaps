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
